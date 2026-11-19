#!/usr/bin/env python3
"""Offline construction of a GRCh38 HapMap3+ / GWAS panel with ancestry and build exports."""
import argparse
import base64
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
from collections import defaultdict
from datetime import datetime, timezone
import gzip
import getpass
import hashlib
from http.client import HTTPException
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from time import sleep
from urllib.request import Request, urlopen
import zipfile

import pandas as pd

ROOT = Path(os.environ.get('HVS_ROOT', Path(__file__).resolve().parent)).resolve()
DNA = re.compile('[ACGT]+')
RSID = re.compile(r'rs[0-9]+')
KEY = ['CHROM', 'POS', 'REF', 'ALT']
CANDIDATES = ['CHROM', 'POS', 'RSID', 'A0', 'A1', 'SOURCE', 'POS37_NATIVE']
REFERENCE_CHR1_LENGTHS = {'reference_grch38': 248956422, 'reference_grch37': 249250621}


def local(value):
    p = Path(os.path.expandvars(os.path.expanduser(value)))
    return p.resolve() if p.is_absolute() else (ROOT / p).resolve()


def chrom_number(value):
    value = str(value).strip().removeprefix('chr')
    if value.startswith('NC_'):
        value = value[3:].split('.')[0].lstrip('0')
    # R numeric columns may arrive as integer-valued floats (e.g. 22.0).
    if re.fullmatch(r'[0-9]+(?:\.0+)?', value):
        number = int(value.split('.')[0])
        return number if 1 <= number <= 22 else None
    return None


def identity(filename):
    p = Path(filename)
    s = p.stat()
    return {'path': str(p.resolve()), 'size': s.st_size, 'mtime_ns': s.st_mtime_ns}


def digest(filename):
    with Path(filename).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def run(*args):
    print('RUN', ' '.join(map(str, args)), flush=True)
    subprocess.run(list(map(str, args)), check=True)


def picard_command(cfg):
    """Resolve module environment variables without invoking a shell."""
    command = cfg['picard_command']
    if not isinstance(command, list) or not command or not all(
            isinstance(x, str) and x for x in command):
        raise ValueError('picard_command must be a nonempty argument list')
    command = [os.path.expandvars(os.path.expanduser(arg)) for arg in command]
    if any(re.search(r'\$(?:[A-Za-z_]|\{)', arg) for arg in command):
        raise ValueError('Unresolved variable in picard_command; load the picard module to set EBROOTPICARD')
    if '-jar' in command:
        index = command.index('-jar') + 1
        if index == len(command):
            raise ValueError('picard_command requires a JAR path after -jar')
        jar = local(command[index])
        if not jar.is_file() or not jar.stat().st_size:
            raise FileNotFoundError(f'Picard JAR not found or empty: {jar}; check the picard module/path')
        command[index] = str(jar)
    return command


def catalog_member(archive):
    """Find the full Catalog table without reading/decompressing its contents."""
    members = [m for m in archive.namelist() if not m.endswith('/')
               and not m.startswith('__MACOSX/') and Path(m).name != '.DS_Store']
    if len(members) != 1:
        raise ValueError('Use the single full Catalog TSV ZIP, not the publication-year split ZIP')
    member = members[0]
    if not archive.getinfo(member).file_size:
        raise ValueError('GWAS Catalog archive contains an empty table')
    return member


def catalog_snapshot(filename):
    """Identify the staged file offline; never assign today's release to an older download."""
    filename = Path(filename)
    if not filename.is_file() or not filename.stat().st_size:
        raise FileNotFoundError(f'GWAS Catalog missing or empty: {filename}. Run bash run.sh --download-sources or stage it locally.')
    names = [filename.name]
    metadata = {}
    if filename.suffix.lower() == '.zip':
        with zipfile.ZipFile(filename) as archive:
            member = catalog_member(archive)
        names.append(Path(member).name)
        metadata['archive_member'] = member
    # Official filenames use rYYYY-MM-DD for the release and eNNN for Ensembl.
    dates = {match for name in names for match in re.findall(r'(?:^|[_-])r(\d{4}-\d{2}-\d{2})(?=[_.-]|$)', name)}
    ensembl = {match for name in names for match in re.findall(r'(?:^|[_-])e(\d+)(?=[_.-]|$)', name)}
    if len(dates) > 1 or len(ensembl) > 1:
        raise ValueError(f'Conflicting GWAS Catalog release metadata in filename/archive: {filename}')
    if dates:
        date = next(iter(dates))
        datetime.strptime(date, '%Y-%m-%d')
        metadata.update(release=f'GWAS Catalog r{date}', identification='versioned_filename')
    else:
        # Generic filenames still have an exact reproducible identity. This is a
        # local snapshot identifier, not a claim about the upstream release date.
        checksum = digest(filename)
        metadata.update(release=f'local-sha256:{checksum}', sha256=checksum,
                        identification='content_sha256', upstream_release='unknown')
    if ensembl:
        metadata['ensembl_release'] = next(iter(ensembl))
        metadata['release'] += f" e{metadata['ensembl_release']}"
    return metadata


def resolve_catalog_release(cfg):
    release = cfg.get('gwas_catalog_release')
    if release is not None and not isinstance(release, str):
        raise ValueError('gwas_catalog_release must be a string or null (use auto for automatic identification)')
    if release is None or not release.strip() or release.strip().lower() == 'auto' or 'REPLACE_' in release:
        snapshot = catalog_snapshot(local(cfg['gwas_catalog']))
        cfg['gwas_catalog_release'] = snapshot['release']
        cfg['gwas_catalog_snapshot'] = snapshot
    else:
        cfg['gwas_catalog_release'] = release.strip()
    return cfg


def load_config(filename):
    cfg = json.loads(local(filename).read_text())
    cfg['picard_command'] = picard_command(cfg)
    if cfg['gwas_catalog_build'] != 'GRCh38' or cfg['gnomad_data_type'] != 'genomes':
        raise ValueError('Initial workflow requires GRCh38 Catalog mappings and gnomAD genomes')
    if cfg['gnomad_release'] != '4.1':
        raise ValueError('This workflow is pinned to gnomAD genomes v4.1; audit headers before changing release')
    if not 0 < cfg['maf_gt'] < 0.5 or not isinstance(cfg['mac_gt'], int) or cfg['mac_gt'] < 0:
        raise ValueError('Invalid MAF/MAC thresholds')
    if not 0 < cfg['gwas_p_threshold'] <= 5e-8:
        raise ValueError('Use a positive genome-wide-significance threshold <= 5e-8')
    if not cfg['ancestries'] or not cfg['require_pass']:
        raise ValueError('Configure ancestries and require PASS reference variants')
    for name, groups in cfg['ancestries'].items():
        if not re.fullmatch('[A-Z]+', name) or not groups or len(set(groups)) != len(groups):
            raise ValueError('Ancestry names must be uppercase; source groups must be distinct')
        if not set(groups) <= {'afr', 'amr', 'asj', 'eas', 'fin', 'mid', 'nfe', 'sas'}:
            raise ValueError(f'Unsupported gnomAD groups for {name}')
    cfg['worker_network'] = 'offline'
    cfg.setdefault('scratch_dir', '/scratch/szabad')
    cfg.setdefault('download_workers', 4)
    cfg['wall_time'] = wall_time(cfg.get('wall_time', '01:00:00'))
    if not isinstance(cfg['download_workers'], int) or not 1 <= cfg['download_workers'] <= 32:
        raise ValueError('download_workers must be between 1 and 32')
    return resolve_catalog_release(cfg)


def scalar(value):
    if isinstance(value, (tuple, list)):
        value = value[0] if len(value) == 1 else None
    return value


def frequency(counts, groups, maf_gt, mac_gt):
    """Pool disjoint group AC/AN, then apply strict minor-frequency/count thresholds."""
    ac = an = 0
    for group in groups:
        a, n = counts.get(f'AC_{group}'), counts.get(f'AN_{group}')
        if a is None or n is None:
            return {'AC': None, 'AN': None, 'MAF': None, 'MAC': None, 'STATUS': 'missing_counts'}
        a, n = scalar(a), scalar(n)
        if a is None or n is None or not math.isfinite(a) or not math.isfinite(n) or int(a) != a or int(n) != n or not 0 <= a <= n:
            return {'AC': None, 'AN': None, 'MAF': None, 'MAC': None, 'STATUS': 'invalid_counts'}
        ac += int(a)
        an += int(n)
    if an == 0:
        return {'AC': ac, 'AN': an, 'MAF': None, 'MAC': None, 'STATUS': 'no_called_alleles'}
    mac = min(ac, an - ac)
    maf = mac / an
    status = 'pass' if maf > maf_gt and mac > mac_gt else 'below_maf_or_mac'
    return {'AC': ac, 'AN': an, 'MAF': maf, 'MAC': mac, 'STATUS': status}


def parse_catalog(chunk, threshold):
    """Keep unambiguous single-variant curated associations; audit unresolved significant rows."""
    required = {'P-VALUE', 'CHR_ID', 'CHR_POS', 'SNP_ID_CURRENT'}
    if not required <= set(chunk):
        raise ValueError(f'GWAS Catalog missing columns: {required - set(chunk)}')
    p = pd.to_numeric(chunk['P-VALUE'], errors='coerce')
    chunk = chunk.loc[p.ge(0) & p.lt(threshold)].copy()
    # A multi-variant/haplotype/interaction row must not be expanded into independent associations.
    ids = chunk['SNP_ID_CURRENT'].str.strip().str.replace(r'^rs', '', regex=True)
    chromosome = chunk['CHR_ID'].map(chrom_number)
    pos = pd.to_numeric(chunk['CHR_POS'], errors='coerce')
    valid = ids.str.fullmatch(r'[0-9]+') & chromosome.notna() & pos.gt(0) & pos.mod(1).eq(0)
    if 'SNPS' in chunk:
        valid &= ~chunk.SNPS.str.contains(r';|,|\s[xX]\s', regex=True)
    accepted = pd.DataFrame({'CHROM': chromosome[valid].astype(int), 'POS': pos[valid].astype('int64'),
        'RSID': 'rs' + ids[valid], 'A0': '', 'A1': '', 'SOURCE': 'GWAS_CATALOG', 'POS37_NATIVE': ''})
    rejected = chunk.loc[~valid].copy()
    rejected['REASON'] = 'unresolved_or_non_autosomal_or_multi_variant_association'
    return accepted, rejected


def read_hm3(filename, *, with_rejections=False):
    if filename.suffix == '.rds':
        import pyreadr
        frame = pyreadr.read_r(str(filename))[None]
    else:
        frame = pd.read_csv(filename, sep='\t', dtype=str, keep_default_na=False)
    required = {'chr', 'pos', 'pos_hg38', 'rsid', 'a0', 'a1'}
    if not required <= set(frame):
        raise ValueError(f'HapMap3+ file missing columns: {required - set(frame)}')
    frame = frame.rename(columns={'chr': 'CHROM', 'pos_hg38': 'POS', 'rsid': 'RSID',
                                  'a0': 'A0', 'a1': 'A1', 'pos': 'POS37_NATIVE'})
    frame['SOURCE'] = 'HAPMAP3_PLUS'
    frame = frame[CANDIDATES].copy()
    chromosome = frame.CHROM.map(chrom_number)
    pos38 = pd.to_numeric(frame.POS, errors='coerce')
    pos37 = pd.to_numeric(frame.POS37_NATIVE, errors='coerce')
    reason = pd.Series('', index=frame.index)
    checks = [
        (chromosome.notna(), 'non_autosomal_chromosome'),
        (pos38.notna(), 'missing_or_invalid_grch38_coordinate'),
        (pos38.gt(0) & pos38.mod(1).eq(0), 'invalid_grch38_coordinate'),
        (pos37.gt(0) & pos37.mod(1).eq(0), 'invalid_native_grch37_coordinate'),
        (frame.RSID.str.fullmatch(RSID.pattern).fillna(False), 'invalid_or_missing_rsid'),
        ((frame.A0.str.fullmatch('[ACGT]') & frame.A1.str.fullmatch('[ACGT]')
          & frame.A0.ne(frame.A1)).fillna(False), 'invalid_snp_alleles'),
    ]
    for valid, label in checks:
        reason.loc[reason.eq('') & ~valid] = label
    rejected = frame.loc[reason.ne('')].copy()
    rejected['REASON'] = reason.loc[rejected.index]
    frame = frame.loc[reason.eq('')].copy()
    if frame.empty:
        raise ValueError('No valid autosomal HapMap3+ variants with GRCh38 coordinates and rsIDs')
    frame['CHROM'] = chromosome.loc[frame.index].astype(int)
    frame['POS'] = pos38.loc[frame.index].astype('int64')
    frame['POS37_NATIVE'] = pos37.loc[frame.index].astype('int64')
    return (frame, rejected) if with_rejections else frame


def candidates(cfg, directory):
    hm3, hm3_rejected = read_hm3(local(cfg['hapmap3_plus']), with_rejections=True)
    hm3_rejected.to_csv(directory / 'hapmap3_unresolved.tsv.gz', sep='\t', index=False)
    frames = [hm3]
    filename = local(cfg['gwas_catalog'])
    rejected_count = significant = 0
    def consume(handle):
        nonlocal rejected_count, significant
        first = True
        for chunk in pd.read_csv(handle, sep='\t', dtype=str, keep_default_na=False, chunksize=50000):
            p = pd.to_numeric(chunk['P-VALUE'], errors='coerce')
            # Preserve study/trait/accession evidence, including unresolved significant records.
            chunk.loc[p.ge(0) & p.lt(cfg['gwas_p_threshold'])].to_csv(
                directory / 'gwas_significant_associations.tsv.gz', sep='\t', index=False,
                mode='w' if first else 'a', header=first, compression='gzip')
            accepted, rejected = parse_catalog(chunk, cfg['gwas_p_threshold'])
            frames.append(accepted)
            significant += len(accepted) + len(rejected)
            rejected_count += len(rejected)
            rejected.to_csv(directory / 'gwas_unresolved.tsv.gz', sep='\t', index=False,
                            mode='w' if first else 'a', header=first, compression='gzip')
            first = False
    if filename.suffix.lower() == '.zip':
        with zipfile.ZipFile(filename) as archive:
            with archive.open(catalog_member(archive)) as handle:
                consume(handle)
    else:
        consume(filename)
    table = pd.concat(frames, ignore_index=True).drop_duplicates()
    for chrom in range(1, 23):
        table.loc[table.CHROM.eq(chrom)].to_csv(directory / f'chr_{chrom}.tsv.gz', sep='\t', index=False)
    return {'candidate_rows': len(table), 'hapmap3_accepted_rows': len(hm3),
            'hapmap3_rejected_rows': len(hm3_rejected),
            'hapmap3_rejection_reasons': hm3_rejected.REASON.value_counts().to_dict(),
            'significant_catalog_rows': significant,
            'unresolved_catalog_rows': rejected_count,
            'source_files': [dict(identity(local(cfg[k])), sha256=digest(local(cfg[k])))
                             for k in ('hapmap3_plus', 'gwas_catalog')]}


def variant_key(chrom, pos, ref, alt):
    return f'{chrom}:{int(pos)}:{ref}:{alt}'


def record_key(record):
    return (record.contig, record.pos, record.ref, record.alts[0])


def write_vcf(table, filename, fasta):
    import pysam
    header = pysam.VariantHeader()
    for name, length in zip(fasta.references, fasta.lengths):
        header.contigs.add(name, length=length)
    with pysam.VariantFile(str(filename), 'w', header=header) as handle:
        for row in table.itertuples(index=False):
            handle.write(handle.new_record(contig=row.CHROM, start=int(row.POS)-1,
                stop=int(row.POS)-1+len(row.REF), id=row.VARIANT_KEY, alleles=(row.REF, row.ALT), filter='PASS'))


def normalize(source, destination, reference, directory):
    normalized = directory / (destination.name + '.unsorted.vcf')
    run('bcftools', 'norm', '-f', reference, '-c', 'e', '-Ov', '-o', normalized, source)
    run('bcftools', 'sort', '-Oz', '-o', destination, normalized)
    run('bcftools', 'index', '--force', '--tbi', destination)
    normalized.unlink()


def subset_vcf(source, positions, destination, directory):
    import pysam
    with pysam.VariantFile(str(source)) as vcf:
        available = {}
        for contig in vcf.header.contigs:
            chrom = chrom_number(contig)
            if chrom is not None:
                if chrom in available:
                    raise ValueError(f'Multiple primary contigs for chromosome {chrom} in {source}')
                available[chrom] = contig
    rows = [(available[int(c)], int(p), int(p)) for c, p in positions.drop_duplicates().itertuples(index=False, name=None)]
    regions = directory / (destination.name + '.regions.tsv')
    pd.DataFrame(rows).to_csv(regions, sep='\t', header=False, index=False)
    run('bcftools', 'view', '--regions-overlap', '0', '-R', regions, '-Oz', '-o', destination, source)
    run('bcftools', 'index', '--force', '--tbi', destination)


def resolve_candidates(table, dbsnp, fasta):
    """HM3 supplies an allele pair; Catalog rows require an unambiguous dbSNP allele pair."""
    import pysam
    lookup = defaultdict(set)
    if dbsnp:
        with pysam.VariantFile(str(dbsnp)) as vcf:
            for record in vcf:
                if not record.alts or len(record.alts) != 1 or not DNA.fullmatch(record.ref) or not DNA.fullmatch(record.alts[0]):
                    continue
                for rsid in (record.id or '').split(';'):
                    lookup[(record.pos, rsid)].add((record.ref, record.alts[0]))
    resolved, rejected = [], []
    for row in table.itertuples(index=False):
        chrom = f'chr{int(row.CHROM)}'
        if row.SOURCE == 'HAPMAP3_PLUS':
            ref = fasta.fetch(chrom, int(row.POS)-1, int(row.POS)).upper()
            pair = (row.A0, row.A1)
            if ref not in pair:
                complement = str.maketrans('ACGT', 'TGCA')
                pair = tuple(a.translate(complement) for a in pair)
            alleles = {(ref, next(a for a in pair if a != ref))} if ref in pair else set()
        else:
            alleles = lookup.get((int(row.POS), row.RSID), set())
        if len(alleles) != 1:
            rejected.append({**row._asdict(), 'REASON': 'unresolved_or_ambiguous_alleles'})
            continue
        ref, alt = next(iter(alleles))
        if ref == alt or fasta.fetch(chrom, int(row.POS)-1, int(row.POS)-1+len(ref)).upper() != ref:
            rejected.append({**row._asdict(), 'REASON': 'grch38_reference_mismatch'})
            continue
        resolved.append({'CHROM': chrom, 'POS': int(row.POS), 'REF': ref, 'ALT': alt,
            'RSID': row.RSID, 'SOURCE': row.SOURCE, 'POS37_NATIVE': row.POS37_NATIVE,
            'VARIANT_KEY': variant_key(chrom, row.POS, ref, alt)})
    return (pd.DataFrame(resolved),
            pd.DataFrame(rejected, columns=[*CANDIDATES, 'REASON']))


def vcf_table(filename):
    import pysam
    rows = []
    with pysam.VariantFile(str(filename)) as handle:
        for r in handle:
            rows.append((r.id, r.contig, r.pos, r.ref, r.alts[0]))
    return pd.DataFrame(rows, columns=['VARIANT_KEY', *KEY])


def mapping_checks(original, lifted, returned):
    """Require one-to-one mappings, exact allele-aware normalized round trips, and native HM3 agreement."""
    forward = lifted.loc[~lifted.VARIANT_KEY.duplicated(keep=False) & ~lifted.duplicated(KEY, keep=False)]
    reverse = returned.loc[~returned.VARIANT_KEY.duplicated(keep=False) & ~returned.duplicated(KEY, keep=False)]
    roundtrip = original[['VARIANT_KEY', *KEY]].merge(reverse, on='VARIANT_KEY', suffixes=('_original', '_returned'))
    matched = roundtrip.VARIANT_KEY.loc[(roundtrip[[k+'_original' for k in KEY]].to_numpy()
                                       == roundtrip[[k+'_returned' for k in KEY]].to_numpy()).all(axis=1)]
    forward = forward.loc[forward.VARIANT_KEY.isin(matched)].copy()
    expected = original.set_index('VARIANT_KEY')
    native_text = expected.POS37_NATIVE.fillna('').astype(str).reindex(forward.VARIANT_KEY)
    native = pd.to_numeric(native_text, errors='coerce').to_numpy()
    # Conflicting native positions (e.g. "123;456") must fail, rather than become missing.
    native_ok = native_text.eq('').to_numpy() | (native == forward.POS.to_numpy())
    same_chr = forward.CHROM.to_numpy() == expected.CHROM.reindex(forward.VARIANT_KEY).to_numpy()
    forward = forward.loc[native_ok & same_chr]
    # Cross-chromosome mappings are conservatively excluded in this autosomal reference panel.
    return forward


def build37_export(selected):
    """Give each export build-specific primary coordinates, retaining the stable GRCh38 key."""
    portable = selected.loc[selected.BUILD37_STATUS.eq('roundtrip_pass')].copy()
    portable = portable.rename(columns={**{k: k+'38' for k in KEY},
        'AC': 'AC_GRCH38_ALT', 'AN': 'AN_GRCH38', 'MAF': 'MAF_GRCH38', 'MAC': 'MAC_GRCH38'})
    portable = portable.rename(columns={k+'37': k for k in KEY})
    portable['POS'] = portable.POS.astype('int64')
    portable['TARGET_VARIANT_KEY'] = [variant_key(*r) for r in portable[KEY].itertuples(index=False, name=None)]
    return portable


def liftover(cfg, source, destination, chain, reference, rejects):
    run(*picard_command(cfg), 'LiftoverVcf', f'I={source}', f'O={destination}',
        f'CHAIN={chain}', f'R={reference}', f'REJECT={rejects}',
        'RECOVER_SWAPPED_REF_ALT=true', 'CREATE_INDEX=false', 'WARN_ON_MISSING_CONTIG=true')


def process_chromosome(cfg, chrom, directory, output):
    import pysam
    candidates_file = output / 'candidates' / f'chr_{chrom}.tsv.gz'
    table = pd.read_csv(candidates_file, sep='\t', keep_default_na=False,
                        dtype={'RSID': str, 'A0': str, 'A1': str, 'POS37_NATIVE': str})
    if table.empty:
        raise ValueError(f'No candidate variants for chr{chrom}')
    ref38, ref37 = local(cfg['reference_grch38']), local(cfg['reference_grch37'])
    with tempfile.TemporaryDirectory(prefix=f'hvs-chr{chrom}-', dir=os.environ['SLURM_TMPDIR']) as tmp:
        tmp = Path(tmp)
        catalog = table.loc[table.SOURCE.eq('GWAS_CATALOG'), ['CHROM', 'POS']]
        dbsnp = None
        if not catalog.empty:
            dbsnp = tmp / 'dbsnp.vcf.gz'
            subset_vcf(local(cfg['dbsnp_vcf']), catalog, dbsnp, tmp)
        with pysam.FastaFile(str(ref38)) as fasta:
            resolved, rejected = resolve_candidates(table, dbsnp, fasta)
            rejected.to_csv(directory / 'unresolved_candidates.tsv.gz', sep='\t', index=False)
            if resolved.empty:
                raise ValueError('No resolved canonical candidates')
            # Normalize before deduplicating by biological allele identity.
            raw = tmp / 'raw.vcf'
            write_vcf(resolved.drop_duplicates('VARIANT_KEY'), raw, fasta)
        normalized = tmp / 'normalized.vcf.gz'
        normalize(raw, normalized, ref38, tmp)
        norm = vcf_table(normalized).rename(columns={'VARIANT_KEY': 'OLD_KEY'})
        annotations = resolved.rename(columns={'VARIANT_KEY': 'OLD_KEY'}).drop(columns=KEY)
        joined = norm.merge(annotations, on='OLD_KEY', validate='one_to_many')
        joined['VARIANT_KEY'] = [variant_key(*r) for r in joined[KEY].itertuples(index=False, name=None)]
        def combine(series):
            return ';'.join(sorted({str(x) for x in series if str(x)}))
        master = joined.groupby(['VARIANT_KEY', *KEY], as_index=False).agg(
            RSIDS=('RSID', combine), SOURCES=('SOURCE', combine), POS37_NATIVE=('POS37_NATIVE', combine))
        master = master.sort_values(KEY).reset_index(drop=True)
        canonical = tmp / 'canonical.vcf'
        with pysam.FastaFile(str(ref38)) as fasta:
            write_vcf(master, canonical, fasta)
        frequencies, duplicate = candidate_counts(cfg, chrom, master, tmp, output)
        # Build mapping is independent of ancestry filtering and preserves stable GRCh38 IDs.
        lifted, back = tmp / 'lifted37.vcf', tmp / 'returned38.vcf'
        liftover(cfg, canonical, lifted, local(cfg['chain_38_to_37']), ref37, directory / 'rejected37.vcf')
        norm37 = tmp / 'normalized37.vcf.gz'
        normalize(lifted, norm37, ref37, tmp)
        liftover(cfg, norm37, back, local(cfg['chain_37_to_38']), ref38, directory / 'rejected_roundtrip.vcf')
        normback = tmp / 'normalized_back.vcf.gz'
        normalize(back, normback, ref38, tmp)
        mapped = mapping_checks(master, vcf_table(norm37), vcf_table(normback))
        map37 = mapped.rename(columns={k: k+'37' for k in KEY})
        master = master.merge(map37, on='VARIANT_KEY', how='left', validate='one_to_one')
        master['POS37'] = master.POS37.astype('Int64')
        master['BUILD37_STATUS'] = master.CHROM37.notna().map({True: 'roundtrip_pass', False: 'unmapped_or_ambiguous_or_failed_validation'})
        master.to_csv(directory / 'master.GRCh38.tsv.gz', sep='\t', index=False)
        summary = {'gnomad_counts_table': str(output / 'gnomad_counts' / f'chr_{chrom}' / 'counts.tsv.gz'),
                   'candidates': len(table), 'resolved_source_rows': len(resolved), 'canonical_variants': len(master),
                   'candidate_rows_by_source': table.SOURCE.value_counts().to_dict(),
                   'unresolved_rows_by_source': rejected.SOURCE.value_counts().to_dict(),
                   'resolution_outcomes': rejected.REASON.value_counts().to_dict(),
                   'canonical_variants_by_source': {s: int(master.SOURCES.str.split(';').map(lambda x: s in x).sum())
                       for s in ('HAPMAP3_PLUS', 'GWAS_CATALOG')},
                   'portable_grch37': len(mapped), 'ancestries': {}}
        for ancestry, source_groups in cfg['ancestries'].items():
            stats = []
            for key in master.VARIANT_KEY:
                counts = frequencies.get(key)
                if counts is None:
                    result = {'STATUS': 'no_frequency_record'}
                elif key in duplicate:
                    result = {'STATUS': 'duplicate_frequency_record'}
                elif not counts['PASS']:
                    result = {'STATUS': 'reference_filter_failed'}
                else:
                    result = frequency(counts, source_groups, cfg['maf_gt'], cfg['mac_gt'])
                stats.append(result)
            audit = pd.concat([master.reset_index(drop=True),
                pd.DataFrame(stats, columns=['AC', 'AN', 'MAF', 'MAC', 'STATUS'])], axis=1)
            for column in ('AC', 'AN', 'MAC'):
                audit[column] = audit[column].astype('Int64')
            audit.to_csv(directory / f'{ancestry}.audit.tsv.gz', sep='\t', index=False)
            selected = audit.loc[audit.STATUS.eq('pass')]
            selected.to_csv(directory / f'{ancestry}.GRCh38.tsv.gz', sep='\t', index=False)
            portable = build37_export(selected)
            portable.to_csv(directory / f'{ancestry}.GRCh37.tsv.gz', sep='\t', index=False)
            summary['ancestries'][ancestry] = {'selected_grch38': len(selected), 'selected_grch37': len(portable),
                                             'outcomes': audit.STATUS.value_counts().to_dict()}
        return summary


def input_inventory(cfg, chroms):
    required = [local(cfg[k]) for k in ('hapmap3_plus', 'gwas_catalog', 'dbsnp_vcf',
                'reference_grch38', 'reference_grch37', 'chain_38_to_37', 'chain_37_to_38')]
    for key in ('reference_grch38', 'reference_grch37'):
        reference = local(cfg[key])
        required += [Path(str(reference)+'.fai'), reference.with_suffix('.dict')]
    vcfs = [local(cfg['dbsnp_vcf'])] + [p for c in chroms
        if nonempty(p := local(cfg['gnomad_vcf'].format(chrom=c)))]
    required += vcfs
    missing_indexes = []
    for vcf in vcfs:
        index = next((Path(str(vcf)+ext) for ext in ('.tbi', '.csi')
                      if Path(str(vcf)+ext).is_file() and Path(str(vcf)+ext).stat().st_size), None)
        if index is None:
            missing_indexes.append(str(vcf)+'.tbi/.csi')
        else:
            required.append(index)
    required = sorted(set(required))
    missing = [str(p) for p in required if not p.is_file() or not p.stat().st_size]
    missing += missing_indexes
    return required, missing


def check_inputs(cfg, chroms):
    required, missing = input_inventory(cfg, chroms)
    if missing:
        raise FileNotFoundError('References/indexes are not ready. Run bash run.sh to download and prepare them:\n' + '\n'.join(missing))
    command = picard_command(cfg)
    if '-jar' in command:
        required.append(Path(command[command.index('-jar') + 1]))
    for tool in ('bcftools', command[0]):
        if not shutil.which(tool):
            raise FileNotFoundError(f'Missing executable: {tool}')
    # Inspect small FASTA indexes, never the multi-GB FASTAs themselves on login.
    for key, chr1_length in REFERENCE_CHR1_LENGTHS.items():
        fai = pd.read_csv(str(local(cfg[key]))+'.fai', sep='\t', header=None, usecols=[0, 1])
        lengths = fai.set_index(0)[1].to_dict()
        if lengths.get('chr1') != chr1_length or not all(f'chr{c}' in lengths for c in chroms):
            raise ValueError(f'{key} must contain chr-prefixed primary contigs for its configured build')
    return {'config': cfg, 'inputs': [identity(p) for p in sorted(set(required))],
            'code_sha256': digest(Path(__file__)),
            'packages': {p: version(p) for p in ('pandas', 'numpy', 'pysam', 'pyreadr')},
            'bcftools_version': subprocess.check_output(['bcftools', '--version'], text=True).strip(),
            'picard_command': command}



def require_login_network():
    if os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Network access is forbidden in SLURM workers; the login-node controller must stage inputs first')


def managed_gnomad(cfg, chrom, output):
    namespace = hashlib.sha256(str(output.resolve()).encode()).hexdigest()[:16]
    return local(cfg.get('scratch_dir', '/scratch/szabad')) / 'harmonized-variant-set' / namespace / f'chr_{chrom}' / f'gnomad.chr{chrom}.vcf.bgz'


def count_cache_ready(output, chrom):
    cache = output / 'gnomad_counts' / f'chr_{chrom}'
    return nonempty(cache / 'counts.tsv.gz') and nonempty(cache / 'manifest.json')


def staged_chromosomes(cfg, chroms):
    output = local(cfg['output_dir'])
    return {chrom for chrom in chroms if count_cache_ready(output, chrom)
            or nonempty(local(cfg['gnomad_vcf'].format(chrom=chrom)))
            or (nonempty(managed_gnomad(cfg, chrom, output))
                and nonempty(Path(str(managed_gnomad(cfg, chrom, output)) + '.tbi'))
                and nonempty(Path(str(managed_gnomad(cfg, chrom, output)) + '.download.json')))}


def stage_gnomad(cfg, chrom):
    """Login-node I/O only. No scientific filtering or normalization occurs here."""
    require_login_network()
    output = local(cfg['output_dir'])
    cache = output / 'gnomad_counts'
    cache.mkdir(parents=True, exist_ok=True)
    with (cache / f'.chr_{chrom}.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if chrom in staged_chromosomes(cfg, [chrom]):
            return
        source = managed_gnomad(cfg, chrom, output)
        url = cfg.get('reference_urls', {}).get('gnomad_vcf', REFERENCE_URLS['gnomad_vcf']).format(chrom=chrom)
        parallel_download(url, source, workers=cfg.get('download_workers', 4))
        download_reference(url + '.tbi', Path(str(source) + '.tbi'))

def validate_source(filename, key):
    """Validate downloaded payloads before making them available to the pipeline."""
    if not filename.stat().st_size:
        raise ValueError('server returned an empty file')
    if key == 'gwas_catalog':
        with zipfile.ZipFile(filename) as archive:
            catalog_member(archive)
            bad = archive.testzip()
            if bad:
                raise ValueError(f'Catalog ZIP checksum failed: {bad}')
    else:
        # The published RDS is gzip-compressed. Read through EOF to check its
        # trailer/CRC, rejecting HTML, truncated streams and compressed error pages.
        with gzip.open(filename, 'rb') as handle:
            if handle.read(2) not in (b'X\n', b'A\n', b'B\n'):
                raise ValueError('response is not a serialized R object')
            while handle.read(1024 * 1024):
                pass


def download_source(key, destination, urls, attempts=3):
    require_login_network()
    destination.parent.mkdir(parents=True, exist_ok=True)
    # A fresh temporary file per invocation avoids reusing an interrupted payload.
    temporary = destination.with_name(destination.name + '.partial')
    last_error = None
    for attempt in range(attempts):
        url = urls[attempt % len(urls)]
        print(f'DOWNLOAD {url} -> {destination} (attempt {attempt + 1}/{attempts})', flush=True)
        try:
            request = Request(url, headers={'User-Agent': 'harmonized-variant-set/0.1',
                                            'Accept-Encoding': 'identity'})
            with urlopen(request, timeout=120) as response, temporary.open('wb') as handle:
                content_type = response.headers.get('Content-Type', '').lower()
                if response.headers.get('x-amzn-waf-action') or 'text/html' in content_type:
                    raise ValueError('download endpoint returned a bot challenge/HTML instead of data')
                metadata = {'final_url': response.url,
                            'content_disposition': response.headers.get('Content-Disposition'),
                            'content_length': response.headers.get('Content-Length')}
                shutil.copyfileobj(response, handle, length=1024 * 1024)
            size = temporary.stat().st_size
            length = metadata['content_length']
            if length is not None and size != int(length):
                raise ValueError(f'incomplete response: expected {length} bytes, received {size}')
            validate_source(temporary, key)
            checksum = digest(temporary)
            temporary.replace(destination)
            receipt = {'url': url, 'downloaded_utc': datetime.now(timezone.utc).isoformat(),
                       **metadata, **identity(destination), 'sha256': checksum}
            destination.with_name(destination.name + '.download.json').write_text(json.dumps(receipt, indent=2) + '\n')
            return
        except (OSError, ValueError, EOFError, zipfile.BadZipFile, HTTPException) as exc:
            last_error = exc
            print(f'Download attempt failed for {key}: {exc}', flush=True)
            if attempt + 1 < attempts:
                sleep(attempt + 1)
    raise ValueError(f'Could not download {key} after {attempts} attempts: {last_error}. '
                     f'No completed source was installed. Check network access or stage the source at {destination}.') from last_error


def download_sources(cfg):
    """Stage discovery sources only; never download the large reference resources."""
    require_login_network()
    sources = {
        'hapmap3_plus': ('https://api.figshare.com/v2/file/download/37802721',
                        'https://ndownloader.figshare.com/files/37802721'),
        'gwas_catalog': ('https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip',),
    }
    for key, urls in sources.items():
        destination = local(cfg[key])
        if destination.is_file() and destination.stat().st_size:
            print(f'Reusing {destination}', flush=True)
            continue
        if destination.exists() and not destination.is_file():
            raise ValueError(f'Source path is not a file: {destination}')
        expected_suffix = '.rds' if key == 'hapmap3_plus' else '.zip'
        if destination.suffix.lower() != expected_suffix:
            raise ValueError(f'Automatic download requires {key} to use {expected_suffix}; '
                             f'stage your custom-format source manually at {destination}')
        download_source(key, destination, urls)
    resolved = resolve_catalog_release(dict(cfg))
    print('Sources staged. Catalog snapshot:', resolved['gwas_catalog_release'], flush=True)


REFERENCE_URLS = {
    'reference_grch38': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz',
    'reference_grch37': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz',
    'chain_38_to_37': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz',
    'chain_37_to_38': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz',
    'dbsnp_vcf': 'https://ftp.ncbi.nih.gov/snp/archive/b156/VCF/GCF_000001405.40.gz',
    'gnomad_vcf': 'https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz',
}


def nonempty(path):
    return path.is_file() and path.stat().st_size > 0


def reference_plan(cfg, chroms):
    urls = {**REFERENCE_URLS, **cfg.get('reference_urls', {})}
    plan = []
    for key in REFERENCE_URLS:
        if key == 'gnomad_vcf':
            continue  # Downloaded/extracted per chromosome on compute nodes.
        for chrom in chroms if key == 'gnomad_vcf' else (None,):
            target = local(cfg[key].format(chrom=chrom))
            fasta = key.startswith('reference_')
            if fasta and target.suffix not in ('.fa', '.fasta', '.fna'):
                raise ValueError(f'{key} must name an uncompressed FASTA (.fa/.fasta/.fna)')
            plan.append({'key': key, 'target': target,
                         'download': Path(str(target) + '.gz') if fasta else target,
                         'url': urls[key].format(chrom=chrom),
                         'vcf': key.endswith('_vcf')})
    return plan


def download_reference(url, destination):
    """Use curl's verified TLS, retry and HTTP Range support for large assets."""
    require_login_network()
    if nonempty(destination):
        print(f'Reusing {destination}', flush=True)
        return
    if not shutil.which('curl'):
        raise FileNotFoundError('curl is required for resumable reference downloads; install/load curl')
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + '.partial')
    print(f'DOWNLOAD/RESUME {url} -> {destination}', flush=True)
    base = ['curl', '--fail', '--location', '--silent', '--show-error', '--retry', '5',
            '--retry-delay', '2', '--connect-timeout', '30', '--speed-limit', '1024',
            '--speed-time', '120', '--output', str(partial)]
    result = subprocess.run([*base, '--continue-at', '-', url], check=False)
    if result.returncode == 33 or (result.returncode == 22 and nonempty(partial)):
        # Some servers do not support Range. Restart instead of appending a full response.
        result = subprocess.run([*base, url], check=False)
    if result.returncode:
        raise ValueError(f'Reference download failed for {destination} (curl exit {result.returncode}). '
                         'Rerun the same command to resume; the partial file is retained.')
    with partial.open('rb') as handle:
        if handle.read(2) != b'\x1f\x8b':
            raise ValueError(f'Reference server did not return gzip/BGZF data for {destination}')
    partial.replace(destination)
    receipt = {'url': url, 'downloaded_utc': datetime.now(timezone.utc).isoformat(),
               **identity(destination)}
    destination.with_name(destination.name + '.download.json').write_text(json.dumps(receipt, indent=2) + '\n')


def download_references(cfg, chroms):
    require_login_network()
    for item in reference_plan(cfg, chroms):
        if nonempty(item['target']):
            print(f'Reusing {item["target"]}', flush=True)
        else:
            download_reference(item['url'], item['download'])
        if item['vcf'] and not any(nonempty(Path(str(item['target']) + ext)) for ext in ('.tbi', '.csi')):
            receipt = item['download'].with_name(item['download'].name + '.download.json')
            # Never attach a public index to an independently staged/custom VCF.
            if receipt.exists() and json.loads(receipt.read_text()).get('url') == item['url']:
                download_reference(item['url'] + '.tbi', Path(str(item['target']) + '.tbi'))


def prepare_references(cfg, chroms):
    """Runs on a compute node: unpack FASTAs and build absent local indexes."""
    import pysam
    for key in ('reference_grch38', 'reference_grch37'):
        reference = local(cfg[key])
        if not nonempty(reference):
            partial = reference.with_name(reference.name + '.partial')
            with gzip.open(str(reference) + '.gz', 'rb') as incoming, partial.open('wb') as outgoing:
                if incoming.peek(1)[:1] != b'>':
                    raise ValueError(f'Not a FASTA download: {reference}.gz')
                shutil.copyfileobj(incoming, outgoing, length=1024 * 1024)
            partial.replace(reference)
        if not nonempty(Path(str(reference) + '.fai')):
            pysam.faidx(str(reference))
        dictionary = reference.with_suffix('.dict')
        if not nonempty(dictionary):
            partial = dictionary.with_name(dictionary.name + '.partial')
            if partial.exists():
                partial.unlink()
            run(*picard_command(cfg), 'CreateSequenceDictionary', f'R={reference}', f'O={partial}')
            partial.replace(dictionary)
    vcfs = {local(cfg['dbsnp_vcf'])} | {p for c in chroms if nonempty(p := local(cfg['gnomad_vcf'].format(chrom=c)))}
    for vcf in sorted(vcfs):
        if not any(nonempty(Path(str(vcf) + ext)) for ext in ('.tbi', '.csi')):
            index = Path(str(vcf) + '.tbi')
            partial = index.with_name(index.name + '.partial')
            run('bcftools', 'index', '--force', '--tbi', '--output', partial, vcf)
            partial.replace(index)
    snapshot = check_inputs(cfg, chroms)
    provenance = local(cfg['output_dir']) / 'provenance'
    partial = provenance / 'preflight.json.partial'
    partial.write_text(json.dumps(snapshot, indent=2) + '\n')
    partial.replace(provenance / 'preflight.json')
    return {'validated_inputs': len(snapshot['inputs'])}



def parallel_download(url, destination, workers=4, chunk_size=128 * 1024 * 1024):
    """Resume independently verified HTTP ranges in one preallocated scratch file."""
    require_login_network()
    if nonempty(destination):
        receipt = destination.with_name(destination.name + '.download.json')
        if not receipt.exists():
            raise ValueError(f'Missing download receipt for managed scratch file: {destination}')
        return json.loads(receipt.read_text())
    destination.parent.mkdir(parents=True, exist_ok=True)
    headers = {'Accept-Encoding': 'identity', 'User-Agent': 'harmonized-variant-set/0.1'}
    with urlopen(Request(url, headers=headers, method='HEAD'), timeout=120) as response:
        size = int(response.headers['Content-Length'])
        etag = response.headers.get('ETag')
        remote_hash = ', '.join(response.headers.get_all('x-goog-hash', []))
        modified = response.headers.get('Last-Modified')
    if size <= 0:
        raise ValueError(f'Empty remote object: {url}')
    source = {'url': url, 'size': size, 'etag': etag, 'last_modified': modified,
              'remote_hash': remote_hash, 'chunk_size': chunk_size}
    partial = destination.with_name(destination.name + '.partial')
    statefile = destination.with_name(destination.name + '.ranges.json')
    state = {'source': source, 'completed': []}
    if nonempty(partial) and partial.stat().st_size == size and statefile.exists() and (etag or modified):
        previous = json.loads(statefile.read_text())
        if previous.get('source') == source:
            state = previous
    chunks = (size + chunk_size - 1) // chunk_size
    done = {i for i in state['completed'] if isinstance(i, int) and 0 <= i < chunks}
    fd = os.open(partial, os.O_CREAT | os.O_RDWR, 0o600)
    try:
        os.ftruncate(fd, size)
        def fetch(index):
            start = index * chunk_size
            end = min(size, start + chunk_size) - 1
            request_headers = {**headers, 'Range': f'bytes={start}-{end}'}
            if etag:
                request_headers['If-Match'] = etag
            elif modified:
                request_headers['If-Unmodified-Since'] = modified
            for attempt in range(3):
                try:
                    with urlopen(Request(url, headers=request_headers), timeout=120) as response:
                        expected_range = f'bytes {start}-{end}/{size}'
                        if response.status != 206 or response.headers.get('Content-Range') != expected_range:
                            raise ValueError(f'Server did not honor byte range {expected_range}')
                        offset = start
                        while offset <= end:
                            data = response.read(min(1024 * 1024, end - offset + 1))
                            if not data:
                                raise ValueError(f'Truncated range at byte {offset}')
                            written = 0
                            while written < len(data):
                                n = os.pwrite(fd, data[written:], offset + written)
                                if not n:
                                    raise OSError('Scratch write made no progress')
                                written += n
                            offset += len(data)
                        if response.read(1):
                            raise ValueError('Server returned excess range data')
                    return index
                except (OSError, ValueError, HTTPException) as exc:
                    if attempt == 2:
                        raise ValueError(f'Range {start}-{end} failed for {url}: {exc}') from exc
                    sleep(attempt + 1)
        print(f'DOWNLOAD chr object: {size / 1e9:.2f} GB, {workers} connections, {len(done)}/{chunks} slices cached', flush=True)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(fetch, i) for i in range(chunks) if i not in done]
            for future in as_completed(futures):
                done.add(future.result())
                os.fsync(fd)  # Never mark a slice complete before its bytes are durable.
                state['completed'] = sorted(done)
                temporary = statefile.with_name(statefile.name + '.tmp')
                temporary.write_text(json.dumps(state))
                temporary.replace(statefile)
                if len(done) == chunks or len(done) % max(1, chunks // 20) == 0:
                    print(f'{destination.name}: {len(done)}/{chunks} slices complete', flush=True)
    finally:
        os.close(fd)
    # GCS supplies MD5 for the pinned gnomAD objects. Verify the entire object,
    # since a range response's checksum describes the whole object, not its slice.
    hashes = dict(item.strip().split('=', 1) for item in remote_hash.split(',') if '=' in item)
    if 'md5' in hashes:
        with partial.open('rb') as handle:
            actual = base64.b64encode(hashlib.file_digest(handle, 'md5').digest()).decode()
        if actual != hashes['md5']:
            statefile.unlink(missing_ok=True)  # Redownload all slices on the next attempt.
            raise ValueError(f'MD5 mismatch for {url}; partial object retained but its slices were invalidated')
    with partial.open('rb') as handle:
        if handle.read(2) != b'\x1f\x8b':
            raise ValueError(f'Not a BGZF/gzip object: {url}')
    receipt = {**source, 'downloaded_utc': datetime.now(timezone.utc).isoformat(),
               'md5_verified': 'md5' in hashes, 'download_workers': workers}
    # Publish receipt first, so a completed managed file always has its metadata.
    destination.with_name(destination.name + '.download.json').write_text(json.dumps(receipt, indent=2) + '\n')
    partial.replace(destination)
    statefile.unlink(missing_ok=True)
    return receipt


def candidate_counts(cfg, chrom, master, work, output):
    """Cache exact normalized candidate AC/AN records, then release managed VCFs."""
    import pysam
    groups = sorted({g for values in cfg['ancestries'].values() for g in values})
    fields = sorted(f'{metric}_{g}' for g in groups for metric in ('AC', 'AN'))
    source = local(cfg['gnomad_vcf'].format(chrom=chrom))
    owned = not nonempty(source)
    url = cfg.get('reference_urls', {}).get('gnomad_vcf', REFERENCE_URLS['gnomad_vcf']).format(chrom=chrom)
    signature = {'schema': 1, 'groups': groups, 'release': cfg.get('gnomad_release', '4.1'),
                 'candidates_sha256': hashlib.sha256(master[['VARIANT_KEY']].to_csv(index=False).encode()).hexdigest(),
                 'reference': identity(local(cfg['reference_grch38'])),
                 'source': {'url': url} if owned else identity(source)}
    managed_source = managed_gnomad(cfg, chrom, output)
    cache = output / 'gnomad_counts' / f'chr_{chrom}'
    cache.parent.mkdir(parents=True, exist_ok=True)
    with (cache.parent / f'.chr_{chrom}.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not cache.exists():
            receipt = None
            if owned:
                source = managed_source
                if not nonempty(source) or not nonempty(Path(str(source) + '.tbi')):
                    raise FileNotFoundError(f'Offline worker input was not staged: {source}. Resume the login-node controller.')
                receipt = json.loads(Path(str(source) + '.download.json').read_text())
            subset = work / 'gnomad.candidates.vcf.gz'
            subset_vcf(source, master[['CHROM', 'POS']].assign(CHROM=chrom), subset, work)
            normalized = work / 'gnomad.normalized.vcf.gz'
            normalize(subset, normalized, local(cfg['reference_grch38']), work)
            wanted = set(master.VARIANT_KEY)
            rows = []
            with pysam.VariantFile(str(normalized)) as handle:
                if not set(fields) <= set(handle.header.info):
                    raise ValueError(f'gnomAD header missing ancestry fields: {set(fields) - set(handle.header.info)}')
                for record in handle:
                    if not record.alts or len(record.alts) != 1:
                        continue
                    key = variant_key(*record_key(record))
                    if key in wanted:
                        rows.append([key, *record_key(record), ';'.join(record.filter),
                                     *(json.dumps(scalar(record.info.get(f))) for f in fields)])
            pending = cache.parent / f'.chr_{chrom}.partial'
            pending.mkdir(exist_ok=True)
            table = pending / 'counts.tsv.gz'
            pd.DataFrame(rows, columns=['VARIANT_KEY', *KEY, 'FILTER', *fields]).to_csv(table, sep='\t', index=False)
            manifest = {'signature': signature, 'ancestries': cfg['ancestries'], 'records': len(rows),
                        'counts_sha256': digest(table), 'download': receipt}
            (pending / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
            pending.replace(cache)
        manifest = json.loads((cache / 'manifest.json').read_text())
        if manifest['signature'] != signature or digest(cache / 'counts.tsv.gz') != manifest['counts_sha256']:
            raise ValueError(f'Candidate count cache does not match inputs or checksum: {cache}; use a fresh output version')
        table = pd.read_csv(cache / 'counts.tsv.gz', sep='\t', dtype=str, keep_default_na=False)
        frequencies, duplicate = {}, set()
        for row in table.to_dict('records'):
            key = row['VARIANT_KEY']
            if key in frequencies:
                duplicate.add(key)
            frequencies[key] = {'PASS': row['FILTER'] == 'PASS', **{f: json.loads(row[f]) for f in fields}}
        if owned:
            # Cleanup follows successful table checksum/read validation, including
            # recovery after a prior job committed its cache but died in cleanup.
            for suffix in ('', '.tbi', '.tbi.download.json', '.download.json'):
                Path(str(managed_source) + suffix).unlink(missing_ok=True)
        return frequencies, duplicate


def atomic_json(path, value):
    temporary = path.with_name(path.name + f'.{os.getpid()}.tmp')
    with temporary.open('w') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


def scientific_config(cfg):
    operational = {'slurm_account', 'scratch_dir', 'download_workers', 'wall_time', 'worker_network', 'gwas_catalog_snapshot'}
    return {key: value for key, value in cfg.items() if key not in operational}


def stage_complete(output, stage, cfg, chroms):
    directory = output / stage
    try:
        manifest = json.loads((directory / 'manifest.json').read_text())
        if scientific_config(manifest['config']) != scientific_config(cfg):
            return False
        if stage == 'candidates':
            files = [f'chr_{c}.tsv.gz' for c in range(1, 23)]
        elif stage == 'merge' or stage.startswith('chr_'):
            files = ['master.GRCh38.tsv.gz'] + [f'{a}.{build}.tsv.gz' for a in cfg['ancestries']
                      for build in ('GRCh38', 'GRCh37')]
            if stage.startswith('chr_'):
                files += [f'{a}.audit.tsv.gz' for a in cfg['ancestries']]
            elif manifest.get('chromosomes') != list(chroms):
                return False
        else:
            files = []
        return all(nonempty(directory / name) for name in files)
    except (OSError, ValueError, KeyError, TypeError):
        return False



def wall_time(value):
    match = re.fullmatch(r'(\d+):([0-5]\d):([0-5]\d)', value) if isinstance(value, str) else None
    if not match or not any(int(part) for part in match.groups()):
        raise ValueError('wall_time/--time must be a positive HH:MM:SS duration')
    return f'{int(match[1]):02d}:{match[2]}:{match[3]}'


def shorten_pending_jobs(active, limit):
    def seconds(value):
        if value == 'UNLIMITED':
            return float('inf')
        match = re.fullmatch(r'(?:(\d+)-)?(\d+):(\d+):(\d+)', value or '')
        if not match:
            return None
        days, hours, minutes, secs = match.groups()
        return int(days or 0) * 86400 + int(hours) * 3600 + int(minutes) * 60 + int(secs)
    updated = set()
    for job in active:
        current = seconds(job.get('wall_time'))
        if job['state'] != 'PENDING' or current is None or current <= seconds(limit):
            continue
        siblings = [j for j in active if j['root'] == job['root']]
        job_id = job['root'] if all(j['state'] == 'PENDING' for j in siblings) else job['id']
        if job_id in updated:
            continue
        subprocess.run(['scontrol', 'update', f'JobId={job_id}', f'TimeLimit={limit}'], check=True)
        updated.add(job_id)
        print(f'Reduced pending job {job_id} wall time to {limit}', flush=True)

def slurm_account(cfg):
    account = cfg.get('slurm_account') or os.environ.get('SBATCH_ACCOUNT')
    if not account:
        worker = ROOT / 'job.sbatch'
        if not worker.exists():
            worker = Path(__file__).with_name('job.sbatch')
        if worker.exists():
            match = re.search(r'^#SBATCH\s+(?:--account(?:=|\s+)|-A\s+)(\S+)', worker.read_text(), re.M)
            if match:
                account = match.group(1)
    if account and not re.fullmatch(r'[A-Za-z0-9_.-]+', account):
        raise ValueError('Invalid SLURM account name')
    return account


def scheduler_jobs(provenance, run_id, completed):
    """Discover live jobs, including submissions accepted before their ID was saved."""
    known = json.loads((provenance / 'jobs.json').read_text()) if (provenance / 'jobs.json').exists() else {}
    journal = json.loads((provenance / 'submissions.json').read_text()) if (provenance / 'submissions.json').exists() else []
    roots = {str(job): stage for stage, job in known.items()}
    roots.update({str(entry['job']): entry['stage'] for entry in journal if entry.get('job')})
    terminal_file = provenance / 'terminal_jobs.json'
    terminal_roots = set(json.loads(terminal_file.read_text())) if terminal_file.exists() else set()
    terminal = {'COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY',
                'NODE_FAIL', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE', 'REVOKED', 'SPECIAL_EXIT'}
    for attempt in range(3):
        try:
            raw = subprocess.check_output(['squeue', '--noheader', '--array', '--user', getpass.getuser(),
                                           '--format=%i|%T|%k|%R|%l'], text=True, stderr=subprocess.PIPE)
            active = []
            for line in raw.splitlines():
                fields = line.strip().split('|', 4)
                if len(fields) < 4:
                    continue
                job, state, comment, reason = (value.strip() for value in fields[:4])
                root = job.split('_')[0]
                label = comment.removeprefix(f'hvs:{run_id}:') if comment.startswith(f'hvs:{run_id}:') else roots.get(root, '')
                stage, _, protocol = label.partition('@')
                if stage in ('prepare', 'candidates', 'chromosome', 'merge'):
                    active.append({'id': job, 'root': root, 'stage': stage, 'state': state, 'reason': reason, 'protocol': protocol,
                                   'wall_time': fields[4].strip() if len(fields) == 5 else None})
            # Accounting confirms that a previously accepted job really ended;
            # an empty/broken queue response alone is never permission to duplicate it.
            inactive = {root: stage for root, stage in roots.items()
                        if root not in {job['root'] for job in active} and stage not in completed and root not in terminal_roots}
            if not inactive:
                return active
            raw = subprocess.check_output(['sacct', '--noheader', '--parsable2', '--jobs', ','.join(sorted(inactive)),
                                           '--format=JobID%64,State%40'], text=True, stderr=subprocess.PIPE)
            settled = set()
            unsettled = set()
            for line in raw.splitlines():
                fields = line.strip().split('|')
                if len(fields) < 2 or '.' in fields[0]:
                    continue
                root = fields[0].strip().split('_')[0]
                words = fields[1].strip().split()
                if not words:
                    continue
                state = words[0].rstrip('+')
                (settled if state in terminal else unsettled).add(root)
            terminal_roots.update(settled - unsettled)
            atomic_json(terminal_file, sorted(terminal_roots))
            if set(inactive) <= terminal_roots:
                return active
        except (OSError, subprocess.CalledProcessError) as exc:
            detail = getattr(exc, 'stderr', None) or str(exc)
            raise ValueError(f'Unable to reconcile SLURM state safely: {detail}. No duplicate jobs were submitted.') from exc
        if attempt < 2:
            sleep(attempt + 1)
    raise ValueError('SLURM accounting has not settled for a recorded submission. No jobs were duplicated; '
                     'rerun the same command after the scheduler catches up. Saved work is retained.')


def cancel_pending(job):
    if job['state'] != 'PENDING':
        raise ValueError(f'Cannot replace active {job["stage"]} job {job["id"]} while its prerequisite is incomplete')
    subprocess.run(['scancel', job['id']], check=True)
    print(f'Replacing pending {job["stage"]} job {job["id"]} with corrected dependencies', flush=True)


def submit(cfg, chroms, max_concurrent, retry_baseline=None):
    """Reconcile durable stage outputs and live SLURM jobs, then submit only missing work."""
    command = picard_command(cfg)
    for tool in ('sbatch', 'squeue', 'bcftools', command[0]):
        if not shutil.which(tool):
            raise FileNotFoundError(f'Missing executable: {tool}')
    output = local(cfg['output_dir'])
    output.mkdir(parents=True, exist_ok=True)
    # Serializes both first submission and all retries from concurrent launchers.
    with (output / '.submit.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        provenance = output / 'provenance'
        provenance.mkdir(exist_ok=True)
        run_file = provenance / 'run.json'
        run_id = hashlib.sha256(str(output).encode()).hexdigest()[:16]
        if run_file.exists():
            saved_run = json.loads(run_file.read_text())
            if saved_run['chromosomes'] != list(chroms):
                raise ValueError(f'This run uses chromosomes {saved_run["chromosomes"]}; use the same list or a new output_dir for a different panel')
        done_candidates = stage_complete(output, 'candidates', cfg, chroms)
        done_chroms = {c for c in chroms if stage_complete(output, f'chr_{c}', cfg, chroms)}
        done_merge = len(done_chroms) == len(chroms) and stage_complete(output, 'merge', cfg, chroms)
        completed = ({'candidates'} if done_candidates else set()) | ({'chromosome'} if len(done_chroms) == len(chroms) else set())
        if done_merge:
            print(f'Already complete. Final manifests: {output / "merge"}', flush=True)
            return {'complete': True, 'done_chroms': sorted(done_chroms), 'active_chroms': [], 'unstaged': []}
        try:
            snapshot = check_inputs(cfg, chroms)
            needs_prepare = False
            completed.add('prepare')
        except FileNotFoundError:
            snapshot = {'config': cfg, 'inputs': [], 'pending_reference_preparation': True}
            needs_prepare = True
        active = scheduler_jobs(provenance, run_id, completed)
        obsolete = [j for j in active if j['stage'] in ('chromosome', 'merge')
                    and j.get('protocol') != 'offline-v1' and j['state'] == 'PENDING']
        cancelled = set()
        for job in obsolete:
            siblings = [j for j in active if j['root'] == job['root']]
            target = job['root'] if all(j['state'] == 'PENDING' for j in siblings) else job['id']
            if target not in cancelled:
                cancel_pending({**job, 'id': target})
                cancelled.add(target)
        active = [j for j in active if j not in obsolete]
        limit = wall_time(cfg.get('wall_time', '01:00:00'))
        shorten_pending_jobs(active, limit)
        original = provenance / 'config.json'
        if original.exists():
            previous = json.loads(original.read_text())
            changed = scientific_config(previous) != scientific_config(cfg)
            has_results = any(output.glob('*/manifest.json'))
            if changed and (active or has_results):
                raise ValueError('Scientific configuration changed for an existing run; use a new output_dir for a different panel. '
                                 'Account, download workers and scratch settings can be changed in place.')
        account = slurm_account(cfg)
        preflight_file = provenance / 'preflight.json'
        if not preflight_file.exists() or (not needs_prepare and json.loads(preflight_file.read_text()).get('pending_reference_preparation')):
            atomic_json(preflight_file, snapshot)
        # Immutable revisions keep queued/running workers isolated from launcher updates.
        config_hash = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:16]
        config_file = provenance / f'config.{config_hash}.json'
        atomic_json(config_file, cfg)
        code_file = provenance / f'workflow.{digest(Path(__file__))[:16]}.py'
        if not code_file.exists():
            shutil.copyfile(Path(__file__), code_file)
        worker_source = ROOT / 'job.sbatch'
        worker_file = worker_source
        if worker_source.exists():
            worker_file = provenance / f'job.{digest(worker_source)[:16]}.sbatch'
            if not worker_file.exists():
                shutil.copyfile(worker_source, worker_file)
        if not original.exists() or (not active and not any(output.glob('*/manifest.json'))):
            atomic_json(original, cfg)
        if not (provenance / 'workflow.py').exists():
            shutil.copyfile(Path(__file__), provenance / 'workflow.py')
        atomic_json(run_file, {'chromosomes': list(chroms), 'run_id': run_id, 'account': account})
        jobs_file = provenance / 'jobs.json'
        jobs = json.loads(jobs_file.read_text()) if jobs_file.exists() else {}
        journal_file = provenance / 'submissions.json'
        journal = json.loads(journal_file.read_text()) if journal_file.exists() else []
        recorded = {entry.get('job') for entry in journal}
        for job in active:
            jobs[job['stage']] = job['root']
            if job['root'] not in recorded:
                journal.append({'stage': job['stage'], 'job': job['root'], 'state': 'adopted',
                                'created_utc': datetime.now(timezone.utc).isoformat()})
                recorded.add(job['root'])
        atomic_json(journal_file, journal)
        atomic_json(jobs_file, jobs)
        logs = ROOT / 'logs'
        logs.mkdir(exist_ok=True)
        spec = ','.join(map(str, chroms))
        def sbatch(stage, dependencies=(), selected=None, **resources):
            keys = [f'chr_{c}' for c in selected] if selected is not None else [stage]
            attempts = submission_counts(journal)
            if retry_baseline is not None and any(attempts.get(key, 0) > retry_baseline.get(key, 0) for key in keys):
                raise ValueError(f'{stage} failed during this controller session; stopping rather than repeatedly submitting the same error. '
                                 f'Worker logs: {logs}. Rerunning the launcher permits a new attempt using saved data.')
            resources.setdefault('time', limit)
            options = [f'--account={account}'] if account else []
            if dependencies:
                options += ['--dependency=afterok:' + ':'.join(dict.fromkeys(dependencies))]
            options += [f'--{key.replace("_", "-")}={value}' for key, value in resources.items()]
            if selected is not None:
                options += [f'--array={",".join(map(str, selected))}%{max_concurrent}']
            args = ['sbatch', '--parsable', '--kill-on-invalid-dep=yes', f'--comment=hvs:{run_id}:{stage}@offline-v1',
                    f'--job-name=hvs-{stage}', *options, '--output=' + str(logs / 'hvs-%x-%A_%a.out'),
                    str(worker_file), str(ROOT), str(config_file), str(code_file), stage, spec]
            intent = {'stage': stage, 'chromosomes': selected, 'state': 'submitting',
                      'created_utc': datetime.now(timezone.utc).isoformat(), 'command': args}
            journal.append(intent)
            atomic_json(journal_file, journal)
            print('RUN', ' '.join(args), flush=True)
            try:
                job = subprocess.check_output(args, text=True, stderr=subprocess.PIPE).strip().split(';')[0]
                if not job.isdigit():
                    raise ValueError(f'Unexpected sbatch response: {job!r}')
            except (OSError, ValueError, subprocess.CalledProcessError) as exc:
                detail = getattr(exc, 'stderr', None) or str(exc)
                intent.update(state='submission_failed', error=detail)
                atomic_json(journal_file, journal)
                raise ValueError(f'SLURM rejected/interrupted {stage} submission: {detail}\n'
                                 'The run is saved. Rerun the same command; completed stages and active jobs will be reused.') from exc
            intent.update(state='accepted', job=job)
            atomic_json(journal_file, journal)
            jobs[stage] = job
            atomic_json(jobs_file, jobs)
            return job
        def available(stage, upstream_replaced=False):
            matches = [job for job in active if job['stage'] == stage]
            keep = []
            for job in matches:
                if 'DependencyNeverSatisfied' in job['reason'] or (upstream_replaced and job['state'] == 'PENDING'):
                    cancel_pending(job)
                    for entry in journal:
                        if entry.get('job') == job['root'] and stage != 'chromosome':
                            entry['cancelled_pending'] = True
                    atomic_json(journal_file, journal)
                else:
                    keep.append(job)
            return keep
        preparation = available('prepare') if needs_prepare else []
        prep_new = needs_prepare and not preparation
        prep_ids = [j['id'] for j in preparation]
        if prep_new:
            prep_ids = [sbatch('prepare')]
        candidates = available('candidates', prep_new) if not done_candidates else []
        candidate_new = not done_candidates and not candidates
        candidate_ids = [j['id'] for j in candidates]
        if candidate_new:
            candidate_ids = [sbatch('candidates', prep_ids)]
        chromosomes = available('chromosome', candidate_new)
        active_chroms = set()
        dependencies = list(dict.fromkeys(prep_ids + candidate_ids))
        for job in chromosomes:
            task = job['id'].split('_')
            if len(task) == 2 and task[1].isdigit():
                if int(task[1]) in set(chroms) - done_chroms:
                    active_chroms.add(int(task[1]))
                    dependencies.append(job['id'])
            elif len(task) == 1:
                # Some schedulers render an array parent rather than expanded tasks.
                active_chroms.update(set(chroms) - done_chroms)
                dependencies.append(job['id'])
        unclaimed = sorted(set(chroms) - done_chroms - active_chroms)
        ready = staged_chromosomes(cfg, unclaimed)
        capacity = max(0, max_concurrent - len(active_chroms))
        missing = sorted(ready)[:capacity]
        if missing:
            dependencies.append(sbatch('chromosome', list(dict.fromkeys(prep_ids + candidate_ids)), selected=missing,
                                       cpus_per_task=1))
        covered = done_chroms | active_chroms | set(missing)
        merges = available('merge', bool(missing) or candidate_new or prep_new or covered != set(chroms))
        if not merges and covered == set(chroms):
            sbatch('merge', dependencies)
        print('Run reconciled:', json.dumps(jobs), '\nFinal manifests:', output / 'merge', flush=True)
        return {'complete': False, 'done_chroms': sorted(done_chroms),
                'active_chroms': sorted(active_chroms | set(missing)),
                'unstaged': sorted(set(chroms) - covered)}



def submission_counts(journal):
    counts = {}
    for entry in journal:
        if not entry.get('job') or entry.get('cancelled_pending'):
            continue
        keys = [f'chr_{c}' for c in entry.get('chromosomes') or []] if entry['stage'] == 'chromosome' else [entry['stage']]
        for key in keys:
            counts[key] = counts.get(key, 0) + 1
    return counts


def coordinate(cfg, chroms, max_concurrent):
    require_login_network()
    output = local(cfg['output_dir'])
    provenance = output / 'provenance'
    provenance.mkdir(parents=True, exist_ok=True)
    state_file = provenance / 'controller.json'
    journal_file = provenance / 'submissions.json'
    baseline = submission_counts(json.loads(journal_file.read_text()) if journal_file.exists() else [])
    state = {'pid': os.getpid(), 'host': __import__('socket').gethostname(), 'state': 'running',
             'started_utc': datetime.now(timezone.utc).isoformat(), 'network_location': 'login-node'}
    atomic_json(state_file, state)
    try:
        while True:
            controls = provenance / 'controller_options.json'
            if controls.exists():
                options = json.loads(controls.read_text())
                for key in ('slurm_account', 'scratch_dir', 'download_workers', 'wall_time'):
                    if key in options['config']:
                        cfg[key] = options['config'][key]
                max_concurrent = options['max_concurrent']
            progress = submit(cfg, chroms, max_concurrent, retry_baseline=baseline)
            state.update(progress, updated_utc=datetime.now(timezone.utc).isoformat())
            atomic_json(state_file, state)
            if progress['complete']:
                state['state'] = 'completed'
                atomic_json(state_file, state)
                return
            slots = max(0, max_concurrent - len(progress['active_chroms']))
            batch = progress['unstaged'][:slots]
            if batch:
                with ThreadPoolExecutor(max_workers=len(batch)) as pool:
                    futures = [pool.submit(stage_gnomad, cfg, chrom) for chrom in batch]
                    for future in as_completed(futures):
                        future.result()
                        submit(cfg, chroms, max_concurrent, retry_baseline=baseline)
            else:
                sleep(30)
    except BaseException as exc:
        detail = latest_worker_error(provenance) if 'failed during this controller session' in str(exc) else ''
        state.update(state='failed', error=str(exc), worker_log=detail, updated_utc=datetime.now(timezone.utc).isoformat())
        atomic_json(state_file, state)
        raise


def start_controller(cfg, chroms, max_concurrent):
    require_login_network()
    output = local(cfg['output_dir'])
    provenance = output / 'provenance'
    provenance.mkdir(parents=True, exist_ok=True)
    options = {'config': cfg, 'max_concurrent': max_concurrent, 'chromosomes': list(chroms)}
    lock = (output / '.controller.lock').open('a')
    try:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            previous_file = provenance / 'controller_options.json'
            previous = json.loads(previous_file.read_text()) if previous_file.exists() else options
            if scientific_config(previous['config']) != scientific_config(cfg) or previous.get('chromosomes', list(chroms)) != list(chroms):
                raise ValueError('A controller is active for different scientific settings; use a separate output_dir')
            atomic_json(previous_file, options)
            print(f'Login-node controller is already running; settings updated. Status: {provenance / "controller.json"}', flush=True)
            return
        atomic_json(provenance / 'controller_options.json', options)
        if all(stage_complete(output, f'chr_{c}', cfg, chroms) for c in chroms) and stage_complete(output, 'merge', cfg, chroms):
            print(f'Already complete. Final manifests: {output / "merge"}', flush=True)
            return
        config_file = provenance / f'controller-config.{hashlib.sha256(json.dumps(options, sort_keys=True).encode()).hexdigest()[:16]}.json'
        atomic_json(config_file, cfg)
        code_file = provenance / f'controller-workflow.{digest(Path(__file__))[:16]}.py'
        if not code_file.exists():
            shutil.copyfile(Path(__file__), code_file)
        logs = ROOT / 'logs'
        logs.mkdir(exist_ok=True)
        log_file = logs / f'hvs-controller-{hashlib.sha256(str(output).encode()).hexdigest()[:16]}.out'
        command = [__import__('sys').executable, '-u', str(code_file), '--stage', 'coordinate',
                   '--config', str(config_file), '--chromosomes', ','.join(map(str, chroms)),
                   '--max-concurrent', str(max_concurrent), '--controller-fd', str(lock.fileno())]
        environment = {**os.environ, 'HVS_ROOT': str(ROOT)}
        with log_file.open('a') as log:
            child = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                                     start_new_session=True, pass_fds=(lock.fileno(),), env=environment)
        print(f'Login-node controller started (PID {child.pid}). All SLURM workers are offline.\n'
              f'Progress/errors: {log_file}\nStatus: {provenance / "controller.json"}', flush=True)
    finally:
        lock.close()  # Child inherits the same open-file lock until it exits.


def latest_worker_error(provenance):
    journal_file = provenance / 'submissions.json'
    journal = json.loads(journal_file.read_text()) if journal_file.exists() else []
    for entry in reversed(journal):
        if not entry.get('job'):
            continue
        paths = sorted((ROOT / 'logs').glob(f'hvs-hvs-{entry["stage"]}-{entry["job"]}_*.out'),
                       key=lambda p: p.stat().st_mtime, reverse=True)
        for path in paths[:4]:
            with path.open('rb') as handle:
                handle.seek(max(0, path.stat().st_size - 16384))
                tail = handle.read().decode(errors='replace')
            if any(word in tail for word in ('Traceback', 'ERROR', 'Error', 'TIME LIMIT', 'OUT_OF_MEMORY')):
                return f'{path}\n{tail[-6000:]}'
    return ''


def show_status(cfg):
    output = local(cfg['output_dir'])
    provenance = output / 'provenance'
    state_file = provenance / 'controller.json'
    state = json.loads(state_file.read_text()) if state_file.exists() else {}
    run_file = provenance / 'run.json'
    run = json.loads(run_file.read_text()) if run_file.exists() else {}
    options_file = provenance / 'controller_options.json'
    options = json.loads(options_file.read_text()) if options_file.exists() else {}
    chroms = run.get('chromosomes', options.get('chromosomes', list(range(1, 23))))
    frozen = provenance / 'config.json'
    effective = json.loads(frozen.read_text()) if frozen.exists() else cfg
    done = [c for c in chroms if stage_complete(output, f'chr_{c}', effective, chroms)]
    merged = len(done) == len(chroms) and stage_complete(output, 'merge', effective, chroms)
    running = False
    lockfile = output / '.controller.lock'
    if lockfile.exists():
        with lockfile.open('r') as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                running = True
    print('Controller:', 'completed' if merged else 'running on login node' if running else state.get('state', 'not running'))
    print(f'Completed chromosomes: {len(done)}/{len(chroms)}')
    print('Final merge:', output / 'merge' if merged else 'not complete')
    if not running and not merged:
        print('Resume with: bash run.sh')
    if not running and state.get('error'):
        print('Last error:', state['error'])
        if state.get('worker_log'):
            print(state['worker_log'])

def merge(cfg, chroms, directory, output):
    names = ['master.GRCh38'] + [f'{a}.{build}' for a in cfg['ancestries'] for build in ('GRCh38', 'GRCh37')]
    summary = {}
    for name in names:
        first = True
        count = 0
        with gzip.open(directory / f'{name}.tsv.gz', 'wt') as destination:
            for chrom in chroms:
                with pd.read_csv(output / f'chr_{chrom}' / f'{name}.tsv.gz', sep='\t',
                                 dtype=str, keep_default_na=False, chunksize=50000) as chunks:
                    for chunk in chunks:
                        chunk.to_csv(destination, sep='\t', index=False, header=first)
                        first = False
                        count += len(chunk)
        summary[name] = count
    return summary


def execute_stage(cfg, stage, chrom, chroms):
    output = local(cfg['output_dir'])
    output.mkdir(parents=True, exist_ok=True)
    stage_name = f'chr_{chrom}' if stage == 'chromosome' else stage
    with (output / f'.{stage_name}.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        preflight = output / 'provenance' / 'preflight.json'
        if preflight.exists() and stage != 'prepare':
            snapshot = json.loads(preflight.read_text())
            if snapshot.get('pending_reference_preparation'):
                raise ValueError('Reference preparation has not completed; wait for the prepare job')
            for expected in snapshot['inputs']:
                if identity(expected['path']) != expected:
                    raise ValueError(f'Input changed after submission: {expected["path"]}')
        stage_name = f'chr_{chrom}' if stage == 'chromosome' else stage
        destination = output / stage_name
        if stage_complete(output, stage_name, cfg, chroms) and (stage != 'prepare' or not input_inventory(cfg, chroms)[1]):
            print(f'Reusing completed stage: {destination}', flush=True)
            return
        if destination.exists():
            recovery = output / 'recovery'
            recovery.mkdir(exist_ok=True)
            destination.replace(recovery / f'{stage_name}.{datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%f")}')
        temporary = output / f'.{stage_name}.{os.environ["SLURM_JOB_ID"]}.partial'
        # SLURM can requeue a job with the same ID; its unpublished stage directory
        # is safe to reuse, while the completed destination above remains protected.
        temporary.mkdir(exist_ok=True)
        if stage == 'prepare':
            summary = prepare_references(cfg, chroms)
        elif stage == 'candidates':
            summary = candidates(cfg, temporary)
        elif stage == 'chromosome':
            summary = process_chromosome(cfg, chrom, temporary, output)
        else:
            summary = merge(cfg, chroms, temporary, output)
        manifest = {'config': cfg, 'stage': stage_name, 'chromosomes': list(chroms), 'summary': summary,
                    'code_sha256': digest(Path(__file__))}
        (temporary / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
        temporary.replace(destination)
        print(json.dumps(manifest, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='config.json')
    parser.add_argument('--stage', default='submit', choices=['submit', 'coordinate', 'check', 'prepare', 'candidates', 'chromosome', 'merge'])
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--check-only', action='store_true')
    mode.add_argument('--status', action='store_true')
    mode.add_argument('--download-sources', action='store_true')
    mode.add_argument('--download-references', action='store_true')
    parser.add_argument('--max-concurrent', type=int, default=2)
    parser.add_argument('--download-workers', type=int)
    parser.add_argument('--scratch-dir')
    parser.add_argument('--account')
    parser.add_argument('--time')
    parser.add_argument('--controller-fd', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--chromosomes')
    parser.add_argument('--chrom', type=int, choices=range(1,23))
    args = parser.parse_args()
    if args.status:
        show_status(json.loads(local(args.config).read_text()))
        return
    if args.check_only:
        args.stage = 'check'
    if (args.download_sources or args.download_references or args.stage in ('submit', 'coordinate', 'check')) and os.environ.get('SLURM_JOB_ID'):
        parser.error('Run downloads/preflight/submission from the login node')
    if args.chromosomes is None:
        try:
            raw_cfg = json.loads(local(args.config).read_text())
            saved_run = local(raw_cfg['output_dir']) / 'provenance' / 'run.json'
            if not saved_run.exists():
                saved_run = saved_run.with_name('controller_options.json')
            saved = json.loads(saved_run.read_text()) if saved_run.exists() else {}
            args.chromosomes = ','.join(map(str, saved['chromosomes'])) if saved else '1-22'
        except (OSError, ValueError, KeyError) as exc:
            parser.error(str(exc))
    chroms = []
    for token in args.chromosomes.split(','):
        if not re.fullmatch(r'[0-9]+(?:-[0-9]+)?', token):
            parser.error('Use autosomes such as 1-22 or 1-3,22')
        parts = list(map(int, token.split('-')))
        if parts[-1] < parts[0]:
            parser.error('Chromosome ranges must ascend')
        chroms += list(range(parts[0], parts[-1]+1))
    chroms = sorted(set(chroms))
    if not chroms or min(chroms)<1 or max(chroms)>22:
        parser.error('Use autosomes 1-22')
    if args.max_concurrent < 1:
        parser.error('--max-concurrent must be positive')
    try:
        if args.download_sources or args.download_references or args.stage == 'submit':
            download_sources(json.loads(local(args.config).read_text()))
            if args.download_sources:
                return
        cfg = load_config(args.config)
        if args.download_workers is not None:
            if not 1 <= args.download_workers <= 32:
                parser.error('--download-workers must be between 1 and 32')
            cfg['download_workers'] = args.download_workers
        if args.scratch_dir is not None:
            cfg['scratch_dir'] = args.scratch_dir
        if args.account is not None:
            cfg['slurm_account'] = args.account
        if args.time is not None:
            cfg['wall_time'] = wall_time(args.time)
        if args.download_references or args.stage == 'submit':
            download_references(cfg, chroms)
            if args.download_references:
                return
    except (OSError, ValueError, EOFError, zipfile.BadZipFile, HTTPException) as exc:
        parser.error(str(exc))
    if args.stage == 'coordinate':
        if args.controller_fd is None:
            parser.error('The coordinator must be started by the login-node launcher')
        with os.fdopen(args.controller_fd, 'a'):
            coordinate(cfg, chroms, args.max_concurrent)
        return
    if args.stage in ('check', 'submit'):
        try:
            if args.stage == 'check':
                print(json.dumps(check_inputs(cfg, chroms), indent=2))
            else:
                start_controller(cfg, chroms, args.max_concurrent)
        except (OSError, ValueError, subprocess.CalledProcessError) as exc:
            parser.error(str(exc))
        return
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Run processing under SLURM; only --stage check belongs on login')
    if args.stage == 'chromosome' and args.chrom is None:
        parser.error('--chrom is required')
    execute_stage(cfg, args.stage, args.chrom, chroms)


if __name__ == '__main__':
    main()
