#!/usr/bin/env python3
"""Offline construction of a GRCh38 HapMap3+ / GWAS panel with ancestry and build exports."""
import argparse
import base64
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
from collections import defaultdict
from datetime import datetime, timezone
import gzip
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
    cfg.setdefault('scratch_dir', '/scratch/szabad')
    cfg.setdefault('download_workers', 4)
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
    namespace = hashlib.sha256(str(output.resolve()).encode()).hexdigest()[:16]
    managed_source = local(cfg.get('scratch_dir', '/scratch/szabad')) / 'harmonized-variant-set' / namespace / f'chr_{chrom}' / f'gnomad.chr{chrom}.vcf.bgz'
    cache = output / 'gnomad_counts' / f'chr_{chrom}'
    cache.parent.mkdir(parents=True, exist_ok=True)
    with (cache.parent / f'.chr_{chrom}.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not cache.exists():
            receipt = None
            if owned:
                source = managed_source
                source.parent.mkdir(parents=True, exist_ok=True)
                receipt = parallel_download(url, source, workers=cfg.get('download_workers', 4))
                download_reference(url + '.tbi', Path(str(source) + '.tbi'))
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


def submit(cfg, chroms, max_concurrent):
    """Small preflight and dependent submissions only; no reference data are parsed on login."""
    # Reference files may still be compressed and their indexes absent. Validate
    # those on the prepare job, then gate every analysis job on its success.
    command = picard_command(cfg)
    for tool in ('sbatch', 'bcftools', command[0]):
        if not shutil.which(tool):
            raise FileNotFoundError(f'Missing executable: {tool}')
    needs_prepare = False
    try:
        snapshot = check_inputs(cfg, chroms)
    except FileNotFoundError:
        needs_prepare = True
        snapshot = {'config': cfg, 'inputs': [], 'pending_reference_preparation': True}
    output = local(cfg['output_dir'])
    output.mkdir(parents=True, exist_ok=True)
    provenance = output / 'provenance'
    if provenance.exists():
        raise FileExistsError(f'{provenance} exists: inspect its jobs before using a fresh output_dir/version')
    provenance.mkdir()  # One submission per version; also prevents accidental duplicate arrays.
    (provenance / 'preflight.json').write_text(json.dumps(snapshot, indent=2)+'\n')
    config_file = provenance / 'config.json'
    config_file.write_text(json.dumps(cfg, indent=2)+'\n')
    code_file = provenance / 'workflow.py'
    shutil.copyfile(Path(__file__), code_file)
    jobs = {}
    spec = ','.join(map(str, chroms))
    logs = ROOT / 'logs'
    logs.mkdir(exist_ok=True)
    def sbatch(stage, *options):
        command = ['sbatch', '--parsable', '--kill-on-invalid-dep=yes', *options,
            '--output='+str(logs / 'hvs-%x-%A_%a.out'), str(ROOT / 'job.sbatch'),
            str(ROOT), str(config_file), str(code_file), stage, spec]
        print('RUN', ' '.join(command), flush=True)
        job = subprocess.check_output(command, text=True).strip().split(';')[0]
        if not job.isdigit():
            raise ValueError(f'Unexpected sbatch result: {job!r}')
        jobs[stage] = job
        (provenance / 'jobs.json').write_text(json.dumps(jobs, indent=2)+'\n')
        return job
    preparation = sbatch('prepare', '--job-name=hvs-prepare', '--time=48:00:00') if needs_prepare else None
    dependency = [f'--dependency=afterok:{preparation}'] if preparation else []
    first = sbatch('candidates', '--job-name=hvs-candidates', *dependency)
    array = sbatch('chromosome', '--job-name=hvs-chromosome', f'--dependency=afterok:{first}',
                   f'--array={spec}%{max_concurrent}', '--time=48:00:00',
                   f'--cpus-per-task={cfg.get("download_workers", 4)}')
    sbatch('merge', '--job-name=hvs-merge', f'--dependency=afterok:{array}')
    print('Submitted:', json.dumps(jobs), '\nFinal manifests:', output / 'merge', flush=True)


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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='config.json')
    parser.add_argument('--stage', default='submit', choices=['submit', 'check', 'prepare', 'candidates', 'chromosome', 'merge'])
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--check-only', action='store_true')
    mode.add_argument('--download-sources', action='store_true')
    mode.add_argument('--download-references', action='store_true')
    parser.add_argument('--max-concurrent', type=int, default=2)
    parser.add_argument('--download-workers', type=int)
    parser.add_argument('--scratch-dir')
    parser.add_argument('--chromosomes', default='1-22')
    parser.add_argument('--chrom', type=int, choices=range(1,23))
    args = parser.parse_args()
    if args.check_only:
        args.stage = 'check'
    if (args.download_sources or args.download_references or args.stage in ('submit', 'check')) and os.environ.get('SLURM_JOB_ID'):
        parser.error('Run downloads/preflight/submission from the login node')
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
        if args.download_references or args.stage == 'submit':
            download_references(cfg, chroms)
            if args.download_references:
                return
    except (OSError, ValueError, EOFError, zipfile.BadZipFile, HTTPException) as exc:
        parser.error(str(exc))
    if args.stage in ('check', 'submit'):
        try:
            if args.stage == 'check':
                print(json.dumps(check_inputs(cfg, chroms), indent=2))
            else:
                submit(cfg, chroms, args.max_concurrent)
        except (OSError, ValueError, subprocess.CalledProcessError) as exc:
            parser.error(str(exc))
        return
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Run processing under SLURM; only --stage check belongs on login')
    if args.stage == 'chromosome' and args.chrom is None:
        parser.error('--chrom is required')
    output = local(cfg['output_dir'])
    output.mkdir(parents=True, exist_ok=True)
    preflight = output / 'provenance' / 'preflight.json'
    if preflight.exists() and args.stage != 'prepare':
        snapshot = json.loads(preflight.read_text())
        if snapshot.get('pending_reference_preparation'):
            raise ValueError('Reference preparation has not completed; wait for the prepare job')
        for expected in snapshot['inputs']:
            if identity(expected['path']) != expected:
                raise ValueError(f'Input changed after submission: {expected["path"]}')
    stage = f'chr_{args.chrom}' if args.stage == 'chromosome' else args.stage
    destination = output / stage
    if destination.exists():
        raise FileExistsError(f'{destination} exists; use a new version/output_dir or remove only the failed stage')
    temporary = output / f'.{stage}.{os.environ["SLURM_JOB_ID"]}.partial'
    # SLURM can requeue a job with the same ID; its unpublished stage directory
    # is safe to reuse, while the completed destination above remains protected.
    temporary.mkdir(exist_ok=True)
    if args.stage == 'prepare':
        summary = prepare_references(cfg, chroms)
    elif args.stage == 'candidates':
        summary = candidates(cfg, temporary)
    elif args.stage == 'chromosome':
        summary = process_chromosome(cfg, args.chrom, temporary, output)
    else:
        summary = merge(cfg, chroms, temporary, output)
    manifest = {'config': cfg, 'stage': stage, 'summary': summary,
                'code_sha256': digest(Path(__file__))}
    (temporary / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    temporary.replace(destination)
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
