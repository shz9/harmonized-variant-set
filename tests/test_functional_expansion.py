"""Functional expansion: real allele normalization, indexed annotations and bounded exports."""
from contextlib import ExitStack
import gzip
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
import pysam
import pysam.bcftools

from test_workflow import hvs


class FunctionalExpansionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.output = self.root / 'output'
        self.prepared = self.output / 'candidates'
        self.prepared.mkdir(parents=True)
        self.fasta = self.root / 'hg38.fa'
        self.fasta.write_text('>chr22\n' + 'A'*400 + '\n')
        pysam.faidx(str(self.fasta))
        self.gtf = self.root / 'basic.gtf'
        rows = []
        def feature(kind, start, end, transcript='T1', gene='G1', biotype='protein_coding', strand='+'):
            rows.append(f'chr22\ttest\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t'
                        f'gene_id "{gene}"; transcript_id "{transcript}"; gene_type "{biotype}"; tag "basic";\n')
        feature('transcript', 40, 120)
        feature('exon', 40, 60)
        feature('exon', 100, 120)
        feature('CDS', 40, 50)
        feature('UTR', 51, 60)
        feature('transcript', 200, 220, 'T2', 'G2', 'lncRNA', '-')
        feature('exon', 200, 220, 'T2', 'G2', 'lncRNA', '-')
        self.gtf.write_text(''.join(rows))
        self.bed = self.root / 'enhancer.bed'
        self.bed.write_text('chr22\t140\t150\tE1\n')
        self.scores = self.root / 'scores.tsv'
        self.scores.write_text('variant_id\tAVI\nchr22:250:A>C\t35\nchr22:260:A>G\t40\nchr22:251:A>C\t29\n')
        self.panel = self.root / 'panel.tsv'
        self.panel.write_text('VARIANT_KEY\nchr22:260:A:G\n')
        def resource(path, **extra):
            return {'path': str(path), 'release': 'test-v1', 'build': 'GRCh38', **extra}
        self.settings = {'enabled': True, 'min_variants': 3, 'target_variants': 5, 'max_variants': 5,
            'promoter_upstream_bp': 2, 'promoter_downstream_bp': 1, 'splice_flank_bp': 5,
            'gencode': resource(self.gtf),
            'beds': [resource(self.bed, label='ENCODE', category='ENHANCER')],
            'alphagenome1': [resource(self.scores, model='AlphaGenome1', kind='avi_phred', threshold=30)],
            'prs_panels': [resource(self.panel, name='SBayesRC_fixture', normalized=True)]}
        self.vcf = self.root / 'gnomad.vcf'
        self.vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr22,length=400>\n'
            '##FILTER=<ID=LowQual,Description="failed">\n'
            '##INFO=<ID=AC_nfe,Number=A,Type=Integer,Description="AC">\n'
            '##INFO=<ID=AN_nfe,Number=1,Type=Integer,Description="AN">\n'
            '##INFO=<ID=AC_afr,Number=A,Type=Integer,Description="AC">\n'
            '##INFO=<ID=AN_afr,Number=1,Type=Integer,Description="AN">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' + ''.join(
            f'chr22\t{pos}\trs{pos}\tA\t{alt}\t.\t{status}\tAC_nfe={ac};AN_nfe=1000;AC_afr={afr};AN_afr=1000\n'
            for pos, alt, status, ac, afr in [
                (45, 'G', 'PASS', 20, 0), (50, 'G', 'PASS', 30, 40), (55, 'G', 'PASS', 21, 0),
                (65, 'G', 'PASS', 30, 30), (80, 'G', 'PASS', 50, 50),
                (100, 'G', 'PASS', 30, 30), (145, 'AT', 'PASS', 30, 0),
                (170, 'G', 'LowQual', 30, 30), (205, 'G', 'PASS', 0, 30),
                (250, 'C', 'PASS', 30, 0), (251, 'C', 'PASS', 30, 30), (260, 'G', 'PASS', 30, 30)]))
        pysam.tabix_compress(str(self.vcf), str(self.vcf)+'.gz')
        pysam.tabix_index(str(self.vcf)+'.gz', preset='vcf')
        self.cfg = {'functional_expansion': self.settings,
            'reference_grch38': str(self.fasta), 'reference_grch37': str(self.fasta),
            'chain_38_to_37': 'unused', 'chain_37_to_38': 'unused',
            'gnomad_vcf': str(self.vcf)+'.gz', 'dbsnp_vcf': str(self.vcf)+'.gz',
            'ancestries': {'EUR': ['nfe'], 'AFR': ['afr']}, 'maf_gt': .005, 'mac_gt': 20}
        pd.DataFrame([[22, 100, 'rs100', 'A', 'G', 'HAPMAP3_PLUS', '']],
                     columns=hvs.CANDIDATES).to_csv(self.prepared / 'chr_22.tsv.gz', sep='\t', index=False)
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.stack.enter_context(patch.dict(os.environ, {'SLURM_TMPDIR': self.tmp.name}))
        self.stack.enter_context(patch.object(hvs, 'run', self.bcftools))
        self.stack.enter_context(patch.object(hvs, 'liftover', self.identity_liftover))

    def bcftools(self, *args):
        self.assertEqual(args[0], 'bcftools')
        getattr(pysam.bcftools, args[1])(*map(str, args[2:]), catch_stdout=False)

    def identity_liftover(self, cfg, source, destination, *unused):
        with pysam.VariantFile(str(source)) as incoming:
            with pysam.VariantFile(str(destination), 'w', header=incoming.header) as outgoing:
                for record in incoming:
                    outgoing.write(record)

    def build(self):
        hvs.validate_functional_config(self.cfg)
        hvs.functional_prepare(self.cfg, self.prepared)
        stage = self.output / 'chr_22'
        stage.mkdir()
        hvs.process_chromosome(self.cfg, 22, stage, self.output)
        return pd.read_csv(stage / 'master.GRCh38.tsv.gz', sep='\t', keep_default_na=False)

    def test_discovery_and_budget_preserve_ancestry_counts_and_backbone(self):
        master = self.build().set_index('POS')
        self.assertEqual(set(master.index), {50, 55, 65, 100, 145, 205, 250, 260})
        self.assertEqual(master.loc[145, 'VARIANT_TYPE'], 'INSERTION')
        self.assertEqual(master.loc[205, 'ELIGIBLE_ANCESTRIES'], 'AFR')
        self.assertEqual(master.loc[250, 'PRS_STATUS'], 'absent_from_supplied_panels')
        self.assertEqual(master.loc[260, 'PRS_STATUS'], 'present')
        self.assertLess(master.loc[250, 'PANEL_PRIORITY'], master.loc[260, 'PANEL_PRIORITY'])
        for repeat in ('merge1', 'merge2'):
            destination = self.output / repeat
            destination.mkdir()
            summary = hvs.merge(self.cfg, [22], destination, self.output)
            selected = pd.read_csv(destination / 'master.GRCh38.tsv.gz', sep='\t')
            self.assertEqual(set(selected.POS), {100, 250, 260, 50, 65})
            self.assertEqual(summary['master.GRCh38'], 5)
            self.assertEqual(summary['AFR.GRCh38'], 4)
            self.assertEqual(summary['EUR.GRCh38'], summary['EUR.GRCh37'])
            report = json.loads((destination / 'functional_expansion.json').read_text())
            self.assertEqual(report['backbone_preserved'], 1)
            self.assertEqual(report['size_status'], 'within_range')
            self.assertEqual(report['candidate_pool'], 8)

    def test_native_atlas_tabix_matches_exact_alt_outside_annotations(self):
        atlas = self.root / 'atlas.tsv'
        atlas.write_text('#CHROM\tPOS\tREF\tALT\traw_score\tPHRED\n'
                         'chr22\t250\tA\tC\t0.5\t35\n'
                         'chr22\t251\tA\tG\t0.9\t45\n'
                         'chr22\t260\tA\tG\t0.7\t40\n')
        pysam.tabix_compress(str(atlas), str(atlas)+'.gz')
        pysam.tabix_index(str(atlas)+'.gz', seq_col=0, start_col=1, end_col=1, zerobased=False)
        self.settings['alphagenome1'][0].update(path=str(atlas)+'.gz', format='atlas_tabix')
        master = self.build().set_index('POS')
        self.assertIn(250, master.index)
        self.assertIn(260, master.index)
        self.assertNotIn(251, master.index)  # High score belongs to another ALT.
        self.assertEqual(json.loads(master.loc[250, 'AG1_EVIDENCE'])[0]['score'], 35)

    def test_splice_and_promoter_boundaries_and_cache_validation(self):
        self.build()
        regions = pd.read_csv(self.prepared / 'functional_regions.chr_22.tsv.gz', sep='\t', keep_default_na=False)
        index = hvs.FunctionalIntervalIndex(regions)
        self.assertTrue(any(row[2] == 'PROMOTER' for row in index.overlaps(218, 219)))
        self.assertFalse(any(row[2] == 'PROMOTER' for row in index.overlaps(217, 218)))
        self.assertTrue(any(row[2] == 'SPLICE_FLANK' for row in index.overlaps(60, 61)))
        self.assertFalse(any(row[2] == 'SPLICE_FLANK' for row in index.overlaps(65, 66)))
        self.assertEqual(hvs.functional_span(140, 'AA', 'A'), (140, 141))
        self.assertEqual(hvs.functional_span(140, 'A', 'AT'), (139, 141))
        work = self.root / 'retry'
        work.mkdir()
        with patch.object(hvs, 'run', side_effect=AssertionError('A valid discovery cache must skip extraction')):
            self.assertEqual(len(hvs.functional_discover(self.cfg, 22, work, self.output)), 7)
        self.settings['splice_flank_bp'] = 6
        with self.assertRaisesRegex(ValueError, 'cache changed'):
            hvs.functional_discover(self.cfg, 22, work, self.output)

    def test_no_score_or_panel_does_not_claim_novelty_and_shortfall_is_reported(self):
        self.settings['alphagenome1'] = []
        self.settings['prs_panels'] = []
        self.settings.update(min_variants=20, target_variants=21, max_variants=22)
        master = self.build()
        added = master.loc[master.SOURCES.ne('HAPMAP3_PLUS')]
        self.assertEqual(set(added.PRS_STATUS), {'not_assessed'})
        destination = self.output / 'merge'
        destination.mkdir()
        hvs.merge(self.cfg, [22], destination, self.output)
        report = json.loads((destination / 'functional_expansion.json').read_text())
        self.assertEqual(report['size_status'], 'below_minimum')
        self.assertEqual(report['selected_master'], 6)

    def test_score_scales_and_resource_contracts_fail_explicitly(self):
        self.settings['alphagenome1'][0]['model'] = 'AlphaGenome2'
        with self.assertRaisesRegex(ValueError, 'AlphaGenome1'):
            hvs.validate_functional_config(self.cfg)
        self.settings['alphagenome1'][0].update(model='AlphaGenome1', kind='quantile', scorer='fixture', threshold=30)
        with self.assertRaisesRegex(ValueError, 'quantile threshold'):
            hvs.validate_functional_config(self.cfg)
        self.assertIsNone(hvs.functional_variant_type('AC', 'GT'))
        self.assertIsNone(hvs.functional_variant_type('A', 'A'+'T'*50))
        self.assertEqual(hvs.functional_variant_type('AT', 'A'), 'DELETION')
        self.scores.write_text('variant_id\tAVI\nchr22:250:A>C\tNaN\n')
        self.settings['alphagenome1'][0].update(kind='avi_phred', threshold=30)
        with self.assertRaisesRegex(ValueError, 'nonfinite'):
            hvs.functional_prepare(self.cfg, self.prepared)


if __name__ == '__main__':
    unittest.main()
