#!/usr/bin/env python3
"""Offline construction of a GRCh38 HapMap3+ / GWAS panel with ancestry and build exports."""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from urllib.request import urlopen
import zipfile

import pandas as pd

ROOT = Path(os.environ.get('HVS_ROOT', Path(__file__).resolve().parent)).resolve()
DNA = re.compile('[ACGT]+')
RSID = re.compile(r'rs[0-9]+')
KEY = ['CHROM', 'POS', 'REF', 'ALT']
CANDIDATES = ['CHROM', 'POS', 'RSID', 'A0', 'A1', 'SOURCE', 'POS37_NATIVE']


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

def load_config(filename):
    cfg = json.loads(local(filename).read_text())
    cfg['picard_command'] = picard_command(cfg)
    if cfg['gwas_catalog_build'] != 'GRCh38' or cfg['gnomad_data_type'] != 'genomes':
        raise ValueError('Initial workflow requires GRCh38 Catalog mappings and gnomAD genomes')
    if cfg['gnomad_release'] != '4.1':
        raise ValueError('This workflow is pinned to gnomAD genomes v4.1; audit headers before changing release')
    if not cfg['gwas_catalog_release'].strip() or 'REPLACE_' in cfg['gwas_catalog_release']:
        raise ValueError('Record the downloaded GWAS Catalog snapshot/release in config.json')
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
    return cfg


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


def read_hm3(filename):
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
    frame['CHROM'] = frame.CHROM.map(chrom_number)
    frame = frame.loc[frame.CHROM.notna(), CANDIDATES].copy()
    if frame.empty:
        raise ValueError('No autosomal HapMap3+ variants')
    frame['CHROM'] = frame.CHROM.astype(int)
    for column in ('POS', 'POS37_NATIVE'):
        number = pd.to_numeric(frame[column], errors='raise')
        if not (number.gt(0) & number.mod(1).eq(0)).all():
            raise ValueError(f'Invalid HapMap3+ coordinates: {column}')
        frame[column] = number.astype('int64')
    if not (frame.RSID.str.fullmatch(RSID.pattern) & frame.A0.str.fullmatch('[ACGT]')
            & frame.A1.str.fullmatch('[ACGT]') & frame.A0.ne(frame.A1)).fillna(False).all():
        raise ValueError('Invalid HapMap3+ SNP coordinates, IDs or alleles')
    return frame


def candidates(cfg, directory):
    frames = [read_hm3(local(cfg['hapmap3_plus']))]
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
    if filename.suffix == '.zip':
        with zipfile.ZipFile(filename) as archive:
            members = [m for m in archive.namelist() if not m.endswith('/') and not m.startswith('__MACOSX/')]
            if len(members) != 1:
                raise ValueError('Use the single full Catalog TSV ZIP, not the publication-year split ZIP')
            with archive.open(members[0]) as handle:
                consume(handle)
    else:
        consume(filename)
    table = pd.concat(frames, ignore_index=True).drop_duplicates()
    for chrom in range(1, 23):
        table.loc[table.CHROM.eq(chrom)].to_csv(directory / f'chr_{chrom}.tsv.gz', sep='\t', index=False)
    return {'candidate_rows': len(table), 'significant_catalog_rows': significant,
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
    run('bcftools', 'index', '--tbi', destination)
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
    run('bcftools', 'index', '--tbi', destination)


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
        # Genome-wide reference frequency query is restricted to canonical candidate positions.
        gnomad = tmp / 'gnomad.vcf.gz'
        subset_vcf(local(cfg['gnomad_vcf'].format(chrom=chrom)),
                   pd.DataFrame({'CHROM': [chrom]*len(master), 'POS': master.POS}), gnomad, tmp)
        normalized_gnomad = tmp / 'gnomad.normalized.vcf.gz'
        normalize(gnomad, normalized_gnomad, ref38, tmp)
        groups = sorted({g for values in cfg['ancestries'].values() for g in values})
        frequencies, duplicate = {}, set()
        with pysam.VariantFile(str(normalized_gnomad)) as handle:
            required = {f'{metric}_{g}' for g in groups for metric in ('AC', 'AN')}
            if not required <= set(handle.header.info):
                raise ValueError(f'gnomAD header is missing required ancestry fields: {required - set(handle.header.info)}')
            for record in handle:
                if not record.alts or len(record.alts) != 1:
                    continue
                key = variant_key(*record_key(record))
                values = {'PASS': list(record.filter) == ['PASS'],
                          **{field: record.info.get(field) for field in required}}
                if key in frequencies:
                    duplicate.add(key)
                frequencies[key] = values
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
        summary = {'candidates': len(table), 'resolved_source_rows': len(resolved), 'canonical_variants': len(master),
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


def check_inputs(cfg, chroms):
    required = [local(cfg[k]) for k in ('hapmap3_plus', 'gwas_catalog', 'dbsnp_vcf',
                'reference_grch38', 'reference_grch37', 'chain_38_to_37', 'chain_37_to_38')]
    for key in ('reference_grch38', 'reference_grch37'):
        reference = local(cfg[key])
        required += [Path(str(reference)+'.fai'), reference.with_suffix('.dict')]
    vcfs = [local(cfg['dbsnp_vcf'])] + [local(cfg['gnomad_vcf'].format(chrom=c)) for c in chroms]
    required += vcfs
    missing_indexes = []
    for vcf in vcfs:
        index = next((Path(str(vcf)+ext) for ext in ('.tbi', '.csi')
                      if Path(str(vcf)+ext).is_file() and Path(str(vcf)+ext).stat().st_size), None)
        if index is None:
            missing_indexes.append(str(vcf)+'.tbi/.csi')
        else:
            required.append(index)
    missing = [str(p) for p in required if not p.is_file() or not p.stat().st_size]
    missing += missing_indexes
    if missing:
        raise FileNotFoundError('Stage these local resources/indexes before submitting:\n'+'\n'.join(missing))
    command = picard_command(cfg)
    if '-jar' in command:
        required.append(Path(command[command.index('-jar') + 1]))
    for tool in ('bcftools', command[0]):
        if not shutil.which(tool):
            raise FileNotFoundError(f'Missing executable: {tool}')
    # Inspect small FASTA indexes, never the multi-GB FASTAs themselves on login.
    for key, chr1_length in (('reference_grch38', 248956422), ('reference_grch37', 249250621)):
        fai = pd.read_csv(str(local(cfg[key]))+'.fai', sep='\t', header=None, usecols=[0, 1])
        lengths = fai.set_index(0)[1].to_dict()
        if lengths.get('chr1') != chr1_length or not all(f'chr{c}' in lengths for c in chroms):
            raise ValueError(f'{key} must contain chr-prefixed primary contigs for its configured build')
    return {'config': cfg, 'inputs': [identity(p) for p in sorted(set(required))],
            'code_sha256': digest(Path(__file__)),
            'packages': {p: version(p) for p in ('pandas', 'numpy', 'pysam', 'pyreadr')},
            'bcftools_version': subprocess.check_output(['bcftools', '--version'], text=True).strip(),
            'picard_command': command}


def download_sources(cfg):
    """Download only the two small discovery sources; large indexed references are staged separately."""
    sources = {
        'hapmap3_plus': 'https://figshare.com/ndownloader/files/37802721',
        'gwas_catalog': 'https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip',
    }
    for key, url in sources.items():
        destination = local(cfg[key])
        if destination.exists():
            print(f'Reusing {destination}', flush=True)
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(destination.name+'.partial')
        print(f'DOWNLOAD {url} -> {destination}', flush=True)
        with urlopen(url, timeout=120) as response, temporary.open('wb') as handle:
            response_metadata = {'final_url': response.url,
                                 'content_disposition': response.headers.get('Content-Disposition')}
            shutil.copyfileobj(response, handle, length=1024*1024)
        temporary.replace(destination)
        receipt = {'url': url, 'downloaded_utc': datetime.now(timezone.utc).isoformat(),
                   **response_metadata, **identity(destination)}
        destination.with_name(destination.name+'.download.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print('Sources staged. Record the actual Catalog release/build in config.json before submission.')


def submit(cfg, chroms, max_concurrent):
    """Small preflight and dependent submissions only; no reference data are parsed on login."""
    snapshot = check_inputs(cfg, chroms)
    if not shutil.which('sbatch'):
        raise FileNotFoundError('Missing SLURM sbatch command')
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
    first = sbatch('candidates', '--job-name=hvs-candidates')
    array = sbatch('chromosome', '--job-name=hvs-chromosome', f'--dependency=afterok:{first}',
                   f'--array={spec}%{max_concurrent}')
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
    parser.add_argument('--stage', default='submit', choices=['submit', 'check', 'candidates', 'chromosome', 'merge'])
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--check-only', action='store_true')
    mode.add_argument('--download-sources', action='store_true')
    parser.add_argument('--max-concurrent', type=int, default=4)
    parser.add_argument('--chromosomes', default='1-22')
    parser.add_argument('--chrom', type=int, choices=range(1,23))
    args = parser.parse_args()
    if args.download_sources:
        if os.environ.get('SLURM_JOB_ID'):
            parser.error('Downloads belong on the login node')
        download_sources(json.loads(local(args.config).read_text()))
        return
    cfg = load_config(args.config)
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
    if args.check_only:
        args.stage = 'check'
    if args.stage in ('submit', 'check') and os.environ.get('SLURM_JOB_ID'):
        parser.error('Run preflight/submission from the login node')
    if args.stage == 'check':
        snapshot = check_inputs(cfg, chroms)
        print(json.dumps(snapshot, indent=2))
        return
    if args.stage == 'submit':
        submit(cfg, chroms, args.max_concurrent)
        return
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Run processing under SLURM; only --stage check belongs on login')
    if args.stage == 'chromosome' and args.chrom is None:
        parser.error('--chrom is required')
    output = local(cfg['output_dir'])
    output.mkdir(parents=True, exist_ok=True)
    preflight = output / 'provenance' / 'preflight.json'
    if preflight.exists():
        for expected in json.loads(preflight.read_text())['inputs']:
            if identity(expected['path']) != expected:
                raise ValueError(f'Input changed after submission: {expected["path"]}')
    stage = f'chr_{args.chrom}' if args.stage == 'chromosome' else args.stage
    destination = output / stage
    if destination.exists():
        raise FileExistsError(f'{destination} exists; use a new version/output_dir or remove only the failed stage')
    temporary = output / f'.{stage}.{os.environ["SLURM_JOB_ID"]}.partial'
    temporary.mkdir()
    if args.stage == 'candidates':
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
