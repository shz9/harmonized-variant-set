import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

import pandas as pd
import pysam
import pysam.bcftools
import pyreadr

MODULE = Path(__file__).resolve().parents[1] / 'workflow.py'
spec = importlib.util.spec_from_file_location('hvs', MODULE)
hvs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hvs)


class StandaloneTests(unittest.TestCase):
    def test_default_resource_paths_stay_within_project(self):
        project = MODULE.parent
        cfg = json.loads((project / 'config.json').read_text())
        resources = ['output_dir', 'hapmap3_plus', 'gwas_catalog', 'gnomad_vcf',
            'dbsnp_vcf', 'reference_grch38', 'reference_grch37',
            'chain_38_to_37', 'chain_37_to_38']
        with patch.object(hvs, 'ROOT', project):
            for key in resources:
                with self.subTest(resource=key):
                    self.assertTrue(hvs.local(cfg[key]).is_relative_to(project))


class PicardCommandTests(unittest.TestCase):
    def test_easybuild_jar_expansion_and_execution(self):
        with tempfile.TemporaryDirectory(prefix='picard module ') as tmp:
            jar = (Path(tmp) / 'picard.jar').resolve()
            jar.write_bytes(b'test jar')
            cfg = {'picard_command': ['java', '-jar', '$EBROOTPICARD/picard.jar']}
            with patch.dict(os.environ, {'EBROOTPICARD': tmp}), patch.object(hvs, 'run') as run:
                self.assertEqual(hvs.picard_command(cfg), ['java', '-jar', str(jar)])
                hvs.liftover(cfg, 'in.vcf', 'out.vcf', 'chain', 'ref.fa', 'reject.vcf')
            self.assertEqual(run.call_args.args[:4], ('java', '-jar', str(jar), 'LiftoverVcf'))
            self.assertEqual(cfg['picard_command'][-1], '$EBROOTPICARD/picard.jar')

    def test_missing_module_variable_has_actionable_error(self):
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaisesRegex(ValueError, 'load the picard module'):
                hvs.picard_command({'picard_command': ['java', '-jar', '$EBROOTPICARD/picard.jar']})

    def test_missing_jar_fails_before_java_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(FileNotFoundError, 'Picard JAR'):
                hvs.picard_command({'picard_command': ['java', '-jar', str(Path(tmp) / 'missing.jar')]})

    def test_wrapper_and_malformed_commands(self):
        self.assertEqual(hvs.picard_command({'picard_command': ['picard']}), ['picard'])
        for command in ([], 'picard', ['java', '-jar'], ['']):
            with self.subTest(command=command), self.assertRaises(ValueError):
                hvs.picard_command({'picard_command': command})


class FrequencyTests(unittest.TestCase):
    def test_strict_minor_count_and_frequency(self):
        cases = [(20, 1000, 'below_maf_or_mac'), (21, 4200, 'below_maf_or_mac'),
                 (21, 4000, 'pass'), (979, 1000, 'pass'), (980, 1000, 'below_maf_or_mac')]
        for ac, an, status in cases:
            with self.subTest(ac=ac, an=an):
                result = hvs.frequency({'AC_nfe': ac, 'AN_nfe': an}, ['nfe'], .005, 20)
                self.assertEqual(result['STATUS'], status)
                self.assertEqual(result['MAC'], min(ac, an-ac))

    def test_pools_counts_instead_of_averaging_frequencies(self):
        result = hvs.frequency({'AC_nfe': 10, 'AN_nfe': 100, 'AC_fin': 90, 'AN_fin': 9000},
                               ['nfe', 'fin'], .005, 20)
        self.assertAlmostEqual(result['MAF'], 100/9100)

    def test_missing_invalid_and_zero_are_distinct(self):
        for counts, status in [({}, 'missing_counts'), ({'AC_nfe': -1, 'AN_nfe': 10}, 'invalid_counts'),
                ({'AC_nfe': 11, 'AN_nfe': 10}, 'invalid_counts'),
                ({'AC_nfe': .5, 'AN_nfe': 10}, 'invalid_counts'),
                ({'AC_nfe': 0, 'AN_nfe': 0}, 'no_called_alleles'),
                ({'AC_nfe': 0, 'AN_nfe': 1000}, 'below_maf_or_mac')]:
            with self.subTest(counts=counts):
                self.assertEqual(hvs.frequency(counts, ['nfe'], .005, 20)['STATUS'], status)


class CandidateTests(unittest.TestCase):
    def test_catalog_strict_p_and_unresolved_associations(self):
        frame = pd.DataFrame({
            'P-VALUE': ['1e-9', '5e-8', '1e-8', '1e-9', '1e-9', '-1'],
            'CHR_ID': ['22', '22', 'X', '22', '22', '22'],
            'CHR_POS': ['100', '200', '300', '400', '500', '600'],
            'SNP_ID_CURRENT': ['123', '124', '125', '126', '', '127'],
            'SNPS': ['rs123', 'rs124', 'rs125', 'rs126;rs999', 'unmapped', 'rs127']})
        accepted, rejected = hvs.parse_catalog(frame, 5e-8)
        self.assertEqual(accepted.RSID.tolist(), ['rs123'])
        self.assertEqual(len(rejected), 3)

    def test_hm3_uses_build38_and_preserves_native_build37(self):
        with tempfile.TemporaryDirectory() as tmp:
            filename = Path(tmp) / 'hm3.tsv'
            pd.DataFrame({'chr': [22], 'pos': [90], 'pos_hg38': [100],
                          'rsid': ['rs1'], 'a0': ['A'], 'a1': ['G']}).to_csv(filename, sep='\t', index=False)
            table = hvs.read_hm3(filename)
            self.assertEqual(table.POS.tolist(), [100])
            self.assertEqual(table.POS37_NATIVE.tolist(), [90])

    def test_reads_published_rds_format(self):
        with tempfile.TemporaryDirectory() as tmp:
            filename = Path(tmp) / 'map.rds'
            frame = pd.DataFrame({'chr':[22], 'pos':[90], 'pos_hg38':[100],
                'rsid':['rs1'], 'a0':['A'], 'a1':['G']})
            pyreadr.write_rds(str(filename), frame)
            self.assertEqual(hvs.read_hm3(filename).POS.tolist(), [100])

    def test_source_union_preserves_association_evidence_and_audits_unresolved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hm3 = root/'hm3.tsv'
            pd.DataFrame({'chr':[22], 'pos':[90], 'pos_hg38':[100], 'rsid':['rs1'],
                'a0':['A'], 'a1':['G']}).to_csv(hm3, sep='\t', index=False)
            catalog = root/'catalog.zip'
            frame = pd.DataFrame({'P-VALUE':['1e-9']*3, 'CHR_ID':['22']*3,
                'CHR_POS':['100', '100', ''], 'SNP_ID_CURRENT':['1', '1', ''],
                'SNPS':['rs1', 'rs1', 'unmapped'], 'STUDY ACCESSION':['GCST1', 'GCST2', 'GCST3']})
            with zipfile.ZipFile(catalog, 'w') as archive:
                archive.writestr('catalog.tsv', frame.to_csv(sep='\t', index=False))
            destination = root/'candidates'
            destination.mkdir()
            summary = hvs.candidates({'hapmap3_plus':str(hm3), 'gwas_catalog':str(catalog),
                'gwas_p_threshold':5e-8}, destination)
            self.assertEqual(summary['candidate_rows'], 2)  # Same allele, two source labels.
            self.assertEqual(summary['unresolved_catalog_rows'], 1)
            evidence = pd.read_csv(destination/'gwas_significant_associations.tsv.gz', sep='\t')
            self.assertEqual(evidence['STUDY ACCESSION'].tolist(), ['GCST1', 'GCST2', 'GCST3'])
            self.assertEqual(len(summary['source_files'][0]['sha256']), 64)


class MappingTests(unittest.TestCase):
    def setUp(self):
        self.original = pd.DataFrame([['chr22:100:A:G', 'chr22', 100, 'A', 'G', '90']],
            columns=['VARIANT_KEY', *hvs.KEY, 'POS37_NATIVE'])
        self.forward = pd.DataFrame([['chr22:100:A:G', 'chr22', 90, 'T', 'C']],
            columns=['VARIANT_KEY', *hvs.KEY])
        self.back = self.original.drop(columns='POS37_NATIVE').copy()

    def test_reverse_complement_keeps_stable_id(self):
        self.assertEqual(len(hvs.mapping_checks(self.original, self.forward, self.back)), 1)

    def test_wrong_allele_roundtrip_is_rejected(self):
        self.back.loc[0, 'ALT'] = 'T'
        self.assertTrue(hvs.mapping_checks(self.original, self.forward, self.back).empty)

    def test_collision_and_conflicting_native_positions_are_rejected(self):
        extra = self.forward.copy()
        extra['VARIANT_KEY'] = 'another_variant'
        self.assertTrue(hvs.mapping_checks(self.original, pd.concat([self.forward, extra]), self.back).empty)
        self.original.loc[0, 'POS37_NATIVE'] = '90;91'
        self.assertTrue(hvs.mapping_checks(self.original, self.forward, self.back).empty)

    def test_native_mismatch_is_rejected(self):
        self.forward.loc[0, 'POS'] = 91
        self.assertTrue(hvs.mapping_checks(self.original, self.forward, self.back).empty)

    def test_build37_export_uses_target_coordinates_and_labels_counts(self):
        selected = self.original.merge(self.forward.rename(columns={k:k+'37' for k in hvs.KEY}), on='VARIANT_KEY')
        selected['BUILD37_STATUS'] = 'roundtrip_pass'
        selected['AC'], selected['AN'], selected['MAF'], selected['MAC'] = 30, 1000, .03, 30
        result = hvs.build37_export(selected).iloc[0]
        self.assertEqual((result.CHROM, result.POS, result.REF, result.ALT), ('chr22', 90, 'T', 'C'))
        self.assertEqual(result.VARIANT_KEY, 'chr22:100:A:G')
        self.assertEqual(result.TARGET_VARIANT_KEY, 'chr22:90:T:C')
        self.assertEqual(result.AC_GRCH38_ALT, 30)
        selected['BUILD37_STATUS'] = 'failed_mapping'
        self.assertTrue(hvs.build37_export(selected).empty)
        self.assertEqual(selected.VARIANT_KEY.tolist(), ['chr22:100:A:G'])


class OfflineIntegrationTests(unittest.TestCase):
    """Real FASTA/VCF IO and bcftools via pysam; isolate Picard in a separate command-contract test."""
    def test_chromosome_selection_and_exports(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / 'hg38.fa'
            fasta.write_text('>chr22\n'+'A'*300+'\n')
            pysam.faidx(str(fasta))
            sources = root / 'output' / 'candidates'
            sources.mkdir(parents=True)
            pd.DataFrame([[22, p, 'rs'+str(p), 'A', 'G', 'HAPMAP3_PLUS', ''] for p in (100, 120, 140, 160)]
                + [[22, 180, 'rs180', '', '', 'GWAS_CATALOG', '']],
                columns=hvs.CANDIDATES).to_csv(sources/'chr_22.tsv.gz', sep='\t', index=False)
            vcf = root / 'gnomad.vcf'
            vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr22,length=300>\n'
                '##INFO=<ID=AC_nfe,Number=A,Type=Integer,Description="AC">\n'
                '##INFO=<ID=AN_nfe,Number=1,Type=Integer,Description="AN">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                'chr22\t100\trs100\tA\tG\t.\tPASS\tAC_nfe=21;AN_nfe=1000\n'
                'chr22\t120\trs120\tA\tG\t.\tPASS\tAC_nfe=20;AN_nfe=1000\n'
                'chr22\t140\trs140\tA\tG\t.\tPASS\tAC_nfe=21;AN_nfe=4200\n'
                'chr22\t180\trs180\tA\tAT\t.\tPASS\tAC_nfe=21;AN_nfe=1000\n')
            pysam.tabix_compress(str(vcf), str(vcf)+'.gz')
            pysam.tabix_index(str(vcf)+'.gz', preset='vcf')
            directory = root/'stage'
            directory.mkdir()
            cfg = {'reference_grch38': str(fasta), 'reference_grch37': str(fasta),
                'chain_38_to_37': 'unused', 'chain_37_to_38': 'unused',
                'gnomad_vcf': str(vcf)+'.gz', 'dbsnp_vcf': str(vcf)+'.gz',
                'ancestries': {'EUR':['nfe']}, 'maf_gt': .005, 'mac_gt': 20}
            def bcftools(*args):
                self.assertEqual(args[0], 'bcftools')
                getattr(pysam.bcftools, args[1])(*map(str, args[2:]), catch_stdout=False)
            def identity_liftover(config, source, destination, *unused):
                with pysam.VariantFile(str(source)) as incoming:
                    with pysam.VariantFile(str(destination), 'w', header=incoming.header) as outgoing:
                        for record in incoming:
                            outgoing.write(record)
            with patch.dict(os.environ, {'SLURM_TMPDIR': tmp}), patch.object(hvs, 'run', bcftools), \
                    patch.object(hvs, 'liftover', identity_liftover):
                summary = hvs.process_chromosome(cfg, 22, directory, root/'output')
            self.assertEqual(summary['ancestries']['EUR']['selected_grch38'], 2)
            selected = pd.read_csv(directory/'EUR.GRCh38.tsv.gz', sep='\t')
            self.assertEqual(selected.POS.tolist(), [100, 180])
            self.assertEqual(selected.ALT.tolist(), ['G', 'AT'])
            audit = pd.read_csv(directory/'EUR.audit.tsv.gz', sep='\t')
            self.assertEqual(audit.STATUS.value_counts().to_dict(),
                {'below_maf_or_mac': 2, 'pass': 2, 'no_frequency_record': 1})
            self.assertEqual(pd.read_csv(directory/'EUR.GRCh37.tsv.gz', sep='\t').POS.tolist(), [100, 180])

    def test_picard_command_checks_target_and_recovers_swapped_snps(self):
        with patch.object(hvs, 'run') as run:
            hvs.liftover({'picard_command':['picard']}, 'in.vcf', 'out.vcf', '38to37.chain', 'hg19.fa', 'reject.vcf')
        args = run.call_args.args
        self.assertIn('R=hg19.fa', args)
        self.assertIn('RECOVER_SWAPPED_REF_ALT=true', args)
        self.assertNotIn('ALLOW_MISSING_FIELDS_IN_HEADER=true', args)


class SubmissionTests(unittest.TestCase):
    def test_dependencies_snapshot_and_duplicate_submission_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = {'output_dir':str(Path(tmp)/'output')}
            with patch.object(hvs, 'ROOT', Path(tmp)), patch.object(hvs, 'check_inputs', return_value={}), \
                    patch.object(hvs.shutil, 'which', return_value='/bin/sbatch'), \
                    patch.object(hvs.subprocess, 'check_output', side_effect=['101\n', '102;cluster\n', '103\n']) as submit:
                hvs.submit(cfg, [21, 22], 2)
                calls = [c.args[0] for c in submit.call_args_list]
                self.assertIn('--dependency=afterok:101', calls[1])
                self.assertIn('--array=21,22%2', calls[1])
                self.assertIn('--dependency=afterok:102', calls[2])
                provenance = Path(tmp)/'output'/'provenance'
                self.assertEqual(json.loads((provenance/'config.json').read_text()), cfg)
                self.assertEqual((provenance/'workflow.py').read_bytes(), MODULE.read_bytes())
                with self.assertRaises(FileExistsError):
                    hvs.submit(cfg, [21, 22], 2)
                self.assertEqual(submit.call_count, 3)


if __name__ == '__main__':
    unittest.main()
