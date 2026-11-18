from contextlib import redirect_stdout
import gzip
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import importlib.util
import io
import json
import os
from pathlib import Path
import tempfile
from threading import Thread
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


class PreparationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name).resolve()
        self.root_patch = patch.object(hvs, 'ROOT', self.root)
        self.root_patch.start()
        self.addCleanup(self.root_patch.stop)
        pause = patch.object(hvs, 'sleep')
        pause.start()
        self.addCleanup(pause.stop)
        self.cfg = json.loads((MODULE.parent / 'config.json').read_text())
        self.cfg['picard_command'] = ['picard']
        self.config = self.root / 'config.json'
        self.catalog = self.root / self.cfg['gwas_catalog']
        self.catalog.parent.mkdir(parents=True)
        self.save_config()

    def save_config(self):
        self.config.write_text(json.dumps(self.cfg))

    def write_catalog(self, member='gwas_catalog_v1.0.2-associations_e114_r2026-10-30.tsv'):
        with zipfile.ZipFile(self.catalog, 'w') as archive:
            archive.writestr(member, 'P-VALUE\tCHR_ID\tCHR_POS\tSNP_ID_CURRENT\n')
        return self.catalog

    def test_auto_and_legacy_placeholders_resolve_existing_archive_offline(self):
        self.write_catalog()
        for release in ('auto', ' AUTO ', 'REPLACE_WITH_SNAPSHOT_DATE_AND_ENSEMBL_RELEASE', '', None):
            with self.subTest(release=release), patch.object(hvs, 'urlopen') as network:
                self.cfg['gwas_catalog_release'] = release
                self.save_config()
                loaded = hvs.load_config('config.json')
                self.assertEqual(loaded['gwas_catalog_release'], 'GWAS Catalog r2026-10-30 e114')
                self.assertEqual(loaded['gwas_catalog_snapshot']['identification'], 'versioned_filename')
                self.assertEqual(json.loads(self.config.read_text())['gwas_catalog_release'], release)
                network.assert_not_called()

    def test_missing_release_field_also_uses_auto(self):
        self.write_catalog()
        del self.cfg['gwas_catalog_release']
        self.save_config()
        self.assertEqual(hvs.load_config('config.json')['gwas_catalog_release'], 'GWAS Catalog r2026-10-30 e114')

    def test_generic_snapshot_has_stable_content_identity(self):
        self.write_catalog('catalog.tsv')
        first = hvs.load_config('config.json')
        checksum = hashlib.sha256(self.catalog.read_bytes()).hexdigest()
        self.assertEqual(first['gwas_catalog_release'], 'local-sha256:' + checksum)
        self.assertEqual(first['gwas_catalog_snapshot']['upstream_release'], 'unknown')
        self.assertEqual(hvs.load_config('config.json'), first)
        self.write_catalog('different.tsv')
        self.assertNotEqual(hvs.load_config('config.json')['gwas_catalog_release'], first['gwas_catalog_release'])

    def test_uncompressed_versioned_filename(self):
        catalog = self.root / 'gwas_catalog_e114_r2026-10-30.tsv'
        catalog.write_text('header\n')
        self.cfg['gwas_catalog'] = str(catalog)
        self.save_config()
        self.assertEqual(hvs.load_config('config.json')['gwas_catalog_release'], 'GWAS Catalog r2026-10-30 e114')

    def test_explicit_release_and_frozen_snapshot_are_preserved(self):
        self.write_catalog()
        resolved = hvs.load_config('config.json')
        self.config.write_text(json.dumps(resolved))
        with patch.object(hvs, 'catalog_snapshot', side_effect=AssertionError('Must use frozen metadata')):
            self.assertEqual(hvs.load_config('config.json'), resolved)
        self.cfg['gwas_catalog_release'] = 'manually verified Catalog snapshot'
        self.save_config()
        self.assertEqual(hvs.load_config('config.json')['gwas_catalog_release'], self.cfg['gwas_catalog_release'])

    def test_missing_empty_and_invalid_release_configuration(self):
        with self.assertRaisesRegex(FileNotFoundError, '--download-sources'):
            hvs.load_config('config.json')
        self.catalog.touch()
        with self.assertRaisesRegex(FileNotFoundError, 'empty'):
            hvs.load_config('config.json')
        self.cfg['gwas_catalog_release'] = 123
        self.save_config()
        with self.assertRaisesRegex(ValueError, 'must be a string'):
            hvs.load_config('config.json')

    def test_split_corrupt_and_empty_archives_fail_early(self):
        self.write_catalog()
        with zipfile.ZipFile(self.catalog, 'a') as archive:
            archive.writestr('another.tsv', 'header\n')
        with self.assertRaisesRegex(ValueError, 'single full Catalog'):
            hvs.load_config('config.json')
        self.catalog.write_bytes(b'not a ZIP archive')
        with self.assertRaises(zipfile.BadZipFile):
            hvs.load_config('config.json')
        with zipfile.ZipFile(self.catalog, 'w') as archive:
            archive.writestr('catalog.tsv', '')
        with self.assertRaisesRegex(ValueError, 'empty table'):
            hvs.load_config('config.json')

    def test_conflicting_release_metadata_is_rejected(self):
        self.catalog = self.catalog.with_name('catalog_r2026-10-01.zip')
        self.cfg['gwas_catalog'] = str(self.catalog)
        self.save_config()
        self.write_catalog()
        with self.assertRaisesRegex(ValueError, 'Conflicting'):
            hvs.load_config('config.json')

    def test_download_reuses_existing_sources_without_network_or_receipt(self):
        self.write_catalog()
        (self.root / self.cfg['hapmap3_plus']).write_bytes(b'existing RDS')
        with patch.object(hvs, 'urlopen') as network, redirect_stdout(io.StringIO()) as output:
            hvs.download_sources(self.cfg)
        network.assert_not_called()
        self.assertIn('GWAS Catalog r2026-10-30 e114', output.getvalue())

    def response(self, payload, length=None):
        response = io.BytesIO(payload)
        response.url = 'https://example.test/releases/latest/catalog.zip'
        response.headers = {} if length is None else {'Content-Length': str(length)}
        return response

    def test_new_download_records_checksum_and_is_reusable(self):
        self.write_catalog()
        payload = self.catalog.read_bytes()
        self.catalog.unlink()
        (self.root / self.cfg['hapmap3_plus']).write_bytes(b'existing RDS')
        with patch.object(hvs, 'urlopen', return_value=self.response(payload, len(payload))) as network:
            with redirect_stdout(io.StringIO()):
                hvs.download_sources(self.cfg)
        network.assert_called_once()
        self.assertEqual(self.catalog.read_bytes(), payload)
        receipt = json.loads(self.catalog.with_name(self.catalog.name + '.download.json').read_text())
        self.assertEqual(receipt['sha256'], hashlib.sha256(payload).hexdigest())
        self.assertEqual(hvs.load_config('config.json')['gwas_catalog_release'], 'GWAS Catalog r2026-10-30 e114')

    def test_failed_download_is_not_promoted_to_completed_source(self):
        (self.root / self.cfg['hapmap3_plus']).write_bytes(b'existing RDS')
        empty_archive = io.BytesIO()
        with zipfile.ZipFile(empty_archive, 'w') as archive:
            archive.writestr('catalog.tsv', '')
        for payload, length, error in ((b'', None, ValueError), (b'truncated', 100, ValueError),
                                       (b'not a zip', None, zipfile.BadZipFile),
                                       (empty_archive.getvalue(), None, ValueError)):
            with self.subTest(payload=payload):
                with patch.object(hvs, 'urlopen', side_effect=lambda *a, **k: self.response(payload, length)):
                    with redirect_stdout(io.StringIO()), self.assertRaisesRegex(ValueError, 'after 3 attempts'):
                        hvs.download_sources(self.cfg)
                self.assertFalse(self.catalog.exists())
                self.assertFalse(self.catalog.with_name(self.catalog.name + '.download.json').exists())

    def test_existing_empty_source_is_repaired_automatically(self):
        (self.root / self.cfg['hapmap3_plus']).touch()
        self.write_catalog()
        payload = gzip.compress(b'X\nserialized R fixture')
        with patch.object(hvs, 'urlopen', side_effect=lambda *a, **k: self.response(payload)) as network:
            with redirect_stdout(io.StringIO()):
                hvs.download_sources(self.cfg)
        network.assert_called_once()
        self.assertEqual((self.root / self.cfg['hapmap3_plus']).read_bytes(), payload)

    def test_submit_stages_sources_before_loading_configuration(self):
        def stage(cfg):
            self.write_catalog()
            (self.root / cfg['hapmap3_plus']).write_bytes(b'existing RDS')
        with patch.object(hvs, 'download_sources', side_effect=stage) as download:
            with patch.object(hvs, 'submit') as submit, patch.object(hvs, 'download_references'):
                with patch('sys.argv', ['workflow.py', '--chromosomes', '22']):
                    hvs.main()
        download.assert_called_once()
        self.assertEqual(submit.call_args.args[1], [22])
        self.assertIn('r2026-10-30', submit.call_args.args[0]['gwas_catalog_release'])

    def test_download_cli_reports_failure_without_traceback(self):
        with patch.object(hvs, 'download_sources', side_effect=ValueError('network unavailable')):
            with patch('sys.argv', ['workflow.py', '--download-sources']):
                with redirect_stdout(io.StringIO()), patch('sys.stderr', new_callable=io.StringIO) as stderr:
                    with self.assertRaises(SystemExit) as error:
                        hvs.main()
        self.assertEqual(error.exception.code, 2)
        self.assertIn('network unavailable', stderr.getvalue())
        self.assertNotIn('Traceback', stderr.getvalue())

    def test_preflight_and_submission_freeze_resolved_snapshot(self):
        self.write_catalog('catalog.tsv')
        for key in ('hapmap3_plus', 'dbsnp_vcf', 'reference_grch38', 'reference_grch37',
                    'chain_38_to_37', 'chain_37_to_38'):
            resource = self.root / self.cfg[key]
            resource.parent.mkdir(parents=True, exist_ok=True)
            resource.write_text('fixture\n')
        for key, length in (('reference_grch38', 248956422), ('reference_grch37', 249250621)):
            reference = self.root / self.cfg[key]
            Path(str(reference) + '.fai').write_text(f'chr1\t{length}\nchr22\t1000\n')
            reference.with_suffix('.dict').write_text('fixture dictionary\n')
        for filename in (self.cfg['dbsnp_vcf'], self.cfg['gnomad_vcf'].format(chrom=22)):
            resource = self.root / filename
            resource.parent.mkdir(parents=True, exist_ok=True)
            resource.write_text('fixture VCF\n')
            Path(str(resource) + '.tbi').write_text('fixture index\n')
        with patch.object(hvs.shutil, 'which', return_value='/bin/mock-tool'):
            with patch.object(hvs.subprocess, 'check_output', return_value='bcftools fixture'):
                with patch('sys.argv', ['workflow.py', '--check-only', '--chromosomes', '22']):
                    with redirect_stdout(io.StringIO()) as output:
                        hvs.main()
            preflight = json.loads(output.getvalue())
            self.assertTrue(preflight['config']['gwas_catalog_release'].startswith('local-sha256:'))
            self.assertFalse((self.root / self.cfg['output_dir']).exists())
            cfg = hvs.load_config('config.json')
            with patch.object(hvs.subprocess, 'check_output', side_effect=['bcftools fixture', '101', '102', '103']):
                with redirect_stdout(io.StringIO()):
                    hvs.submit(cfg, [22], 1)
        provenance = self.root / self.cfg['output_dir'] / 'provenance'
        frozen = json.loads((provenance / 'config.json').read_text())
        self.assertEqual(frozen['gwas_catalog_snapshot'], cfg['gwas_catalog_snapshot'])
        self.assertEqual(frozen['gwas_catalog_release'], cfg['gwas_catalog_release'])
        self.assertEqual(json.loads(self.config.read_text())['gwas_catalog_release'], 'auto')


class DownloadHTTPTests(unittest.TestCase):
    """Real urllib requests exercise challenges, truncation, fallback and retries."""
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.destination = Path(self.tmp.name) / 'map.rds'
        self.payload = gzip.compress(b'X\nserialized R fixture')
        calls = self.calls = []
        payload = self.payload

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def do_GET(self):
                calls.append(self.path)
                self.send_response(200)
                if self.path == '/challenge':
                    self.send_header('Content-Type', 'text/html')
                    self.send_header('x-amzn-waf-action', 'challenge')
                    body = b''
                elif self.path == '/truncated':
                    body = payload[:12]
                elif self.path == '/corrupt':
                    body = b'<html>error</html>'
                else:
                    body = payload
                self.send_header('Content-Length', str(len(payload) if self.path == '/truncated' else len(body)))
                self.end_headers()
                self.wfile.write(body)
                self.close_connection = True

        self.server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        thread = Thread(target=self.server.serve_forever, daemon=True)
        thread.start()
        self.addCleanup(self.server.server_close)
        self.addCleanup(thread.join, 5)
        self.addCleanup(self.server.shutdown)
        self.base = 'http://127.0.0.1:' + str(self.server.server_port)
        pause = patch.object(hvs, 'sleep')
        pause.start()
        self.addCleanup(pause.stop)

    def test_empty_bot_challenge_falls_back_to_working_endpoint(self):
        with redirect_stdout(io.StringIO()):
            hvs.download_source('hapmap3_plus', self.destination,
                                (self.base + '/challenge', self.base + '/good'))
        self.assertEqual(self.calls, ['/challenge', '/good'])
        self.assertEqual(self.destination.read_bytes(), self.payload)
        receipt = json.loads(self.destination.with_name('map.rds.download.json').read_text())
        self.assertEqual(receipt['sha256'], hashlib.sha256(self.payload).hexdigest())

    def test_truncated_download_fails_then_clean_rerun_recovers(self):
        with redirect_stdout(io.StringIO()), self.assertRaisesRegex(ValueError, 'after 3 attempts'):
            hvs.download_source('hapmap3_plus', self.destination, (self.base + '/truncated',))
        self.assertFalse(self.destination.exists())
        with redirect_stdout(io.StringIO()):
            hvs.download_source('hapmap3_plus', self.destination, (self.base + '/good',))
        self.assertEqual(self.destination.read_bytes(), self.payload)
        self.assertFalse(self.destination.with_name('map.rds.partial').exists())

    def test_html_with_success_status_cannot_be_installed_as_rds(self):
        with redirect_stdout(io.StringIO()), self.assertRaisesRegex(ValueError, 'after 3 attempts'):
            hvs.download_source('hapmap3_plus', self.destination, (self.base + '/corrupt',))
        self.assertFalse(self.destination.exists())


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

    def test_unmapped_coordinates_and_non_rsid_rows_are_audited(self):
        with tempfile.TemporaryDirectory() as tmp:
            filename = Path(tmp) / 'hm3.tsv'
            pd.DataFrame({'chr': [22] * 4, 'pos': [90, 91, 92, 93],
                          'pos_hg38': ['100', '', '102', '103'],
                          'rsid': ['rs1', 'rs2', '22:92:A:G', 'rs4'],
                          'a0': ['A', 'A', 'A', 'N'], 'a1': ['G'] * 4}).to_csv(filename, sep='\t', index=False)
            frame, rejected = hvs.read_hm3(filename, with_rejections=True)
            self.assertEqual(frame.RSID.tolist(), ['rs1'])
            self.assertEqual(rejected.REASON.tolist(), ['missing_or_invalid_grch38_coordinate',
                                                       'invalid_or_missing_rsid', 'invalid_snp_alleles'])

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
            cfg = {'output_dir':str(Path(tmp)/'output'), 'picard_command': ['picard']}
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
