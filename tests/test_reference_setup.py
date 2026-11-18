"""Reference staging and full orchestration from an empty data directory."""
from contextlib import redirect_stdout
import gzip
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
from threading import Thread
import unittest
from unittest.mock import patch
from urllib.request import Request, urlopen
import zipfile

import pandas as pd
import pyreadr
import pysam
import pysam.bcftools
from test_workflow import hvs, MODULE


class ReferenceSetupTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name).resolve()
        self.remote = self.root / 'remote'
        self.remote.mkdir()
        self.project = self.root / 'project'
        self.project.mkdir()
        root_patch = patch.object(hvs, 'ROOT', self.project)
        root_patch.start()
        self.addCleanup(root_patch.stop)
        self.cfg = json.loads((MODULE.parent / 'config.json').read_text())
        self.cfg.update(picard_command=['picard'], ancestries={'EUR': ['nfe']})
        pyreadr.write_rds(str(self.remote / 'hm3.rds'), pd.DataFrame({
            'chr': [22], 'pos': [100], 'pos_hg38': [100], 'rsid': ['rs1'], 'a0': ['A'], 'a1': ['G']}))
        rds = self.remote / 'hm3.rds'
        rds.write_bytes(gzip.compress(rds.read_bytes()))
        with zipfile.ZipFile(self.remote / 'catalog.zip', 'w') as archive:
            archive.writestr('catalog.tsv', 'P-VALUE\tCHR_ID\tCHR_POS\tSNP_ID_CURRENT\n1e-9\t22\t100\t1\n')
        for name in ('hg38.fa.gz', 'hg19.fa.gz'):
            with gzip.open(self.remote / name, 'wb') as handle:
                handle.write(b'>chr1\nAAAAAAAAAA\n>chr22\n' + b'A' * 300 + b'\n')
        for name in ('38to37.gz', '37to38.gz'):
            (self.remote / name).write_bytes(gzip.compress(b'fixture chain\n'))
        raw = self.remote / 'sites.vcf'
        raw.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr22,length=300>\n'
                       '##INFO=<ID=AC_nfe,Number=A,Type=Integer,Description="AC">\n'
                       '##INFO=<ID=AN_nfe,Number=1,Type=Integer,Description="AN">\n'
                       '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                       'chr22\t100\trs1\tA\tG\t.\tPASS\tAC_nfe=21;AN_nfe=1000\n')
        for name in ('dbsnp.gz', 'chr22.bgz'):
            pysam.tabix_compress(str(raw), str(self.remote / name))
            pysam.tabix_index(str(self.remote / name), preset='vcf')
        remote = self.remote
        requests = self.requests = []

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def do_GET(self):
                requests.append((self.path, self.headers.get('Range')))
                filename = remote / self.path.lstrip('/')
                if not filename.is_file():
                    self.send_error(404)
                    return
                data = filename.read_bytes()
                start = int(self.headers['Range'].split('=')[1].split('-')[0]) if 'Range' in self.headers else 0
                self.send_response(206 if start else 200)
                if start:
                    self.send_header('Content-Range', f'bytes {start}-{len(data)-1}/{len(data)}')
                self.send_header('Content-Length', str(len(data) - start))
                self.end_headers()
                self.wfile.write(data[start:])

        self.server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        thread = Thread(target=self.server.serve_forever, daemon=True)
        thread.start()
        self.addCleanup(self.server.server_close)
        self.addCleanup(thread.join, 5)
        self.addCleanup(self.server.shutdown)
        self.base = f'http://127.0.0.1:{self.server.server_port}'
        self.cfg['reference_urls'] = {key: self.base + '/' + name for key, name in {
            'reference_grch38': 'hg38.fa.gz', 'reference_grch37': 'hg19.fa.gz',
            'chain_38_to_37': '38to37.gz', 'chain_37_to_38': '37to38.gz',
            'dbsnp_vcf': 'dbsnp.gz', 'gnomad_vcf': 'chr{chrom}.bgz'}.items()}
        self.config = self.project / 'config.json'
        self.config.write_text(json.dumps(self.cfg))

    def stage_discovery_request(self, request, **kwargs):
        url = request.full_url
        name = 'hm3.rds' if 'figshare' in url else 'catalog.zip'
        return urlopen(Request(self.base + '/' + name), **kwargs)

    def native_tools(self, *args):
        if args[0] == 'bcftools':
            getattr(pysam.bcftools, args[1])(*map(str, args[2:]), catch_stdout=False)
            return
        self.assertEqual(args[0], 'picard')
        options = dict(str(arg).split('=', 1) for arg in args[2:])
        if args[1] == 'CreateSequenceDictionary':
            pysam.dict('-o', options['O'], options['R'])
        else:
            self.assertEqual(args[1], 'LiftoverVcf')
            with pysam.VariantFile(options['I']) as source:
                with pysam.VariantFile(options['O'], 'w', header=source.header) as target:
                    for record in source:
                        target.write(record)

    def test_empty_project_to_final_exports(self):
        jobs = []
        def external(command, **kwargs):
            if command[0] == 'bcftools':
                return 'bcftools fixture'
            self.assertEqual(command[0], 'sbatch')
            jobs.append(command)
            return str(100 + len(jobs))
        with patch.object(hvs, 'urlopen', side_effect=self.stage_discovery_request):
            with patch.object(hvs.shutil, 'which', return_value='/mock/tool'):
                with patch.object(hvs.subprocess, 'check_output', side_effect=external):
                    with patch.object(hvs, 'REFERENCE_CHR1_LENGTHS', {'reference_grch38': 10, 'reference_grch37': 10}):
                        with patch.object(hvs, 'run', side_effect=self.native_tools):
                            with patch('sys.argv', ['workflow.py', '--config', str(self.config), '--chromosomes', '22']):
                                with redirect_stdout(io.StringIO()):
                                    hvs.main()
                            self.assertEqual(len(jobs), 4)
                            self.assertIn('--time=48:00:00', jobs[0])
                            for index in range(1, 4):
                                self.assertIn(f'--dependency=afterok:{100 + index}', jobs[index])
                            output = self.project / self.cfg['output_dir']
                            provenance = output / 'provenance'
                            self.assertTrue(json.loads((provenance / 'preflight.json').read_text())['pending_reference_preparation'])
                            # Exercise building indexes for independently staged VCFs too.
                            for key in ('dbsnp_vcf', 'gnomad_vcf'):
                                Path(str(hvs.local(self.cfg[key].format(chrom=22))) + '.tbi').unlink()
                            for stage in ('prepare', 'candidates', 'chromosome', 'merge'):
                                args = ['workflow.py', '--config', str(provenance / 'config.json'), '--stage', stage,
                                        '--chromosomes', '22']
                                if stage == 'chromosome':
                                    args += ['--chrom', '22']
                                with patch.dict(os.environ, {'SLURM_JOB_ID': '123', 'SLURM_TMPDIR': str(self.root)}):
                                    with patch('sys.argv', args), redirect_stdout(io.StringIO()):
                                        hvs.main()
                            # Existing resources are reused without another HTTP request.
                            count = len(self.requests)
                            hvs.download_references(self.cfg, [22])
                            self.assertEqual(len(self.requests), count)
        for build in ('GRCh38', 'GRCh37'):
            exported = pd.read_csv(output / 'merge' / f'EUR.{build}.tsv.gz', sep='\t')
            self.assertEqual(exported.VARIANT_KEY.tolist(), ['chr22:100:A:G'])
        self.assertFalse(json.loads((provenance / 'preflight.json').read_text()).get('pending_reference_preparation'))
        self.assertTrue((output / 'merge' / 'manifest.json').exists())
        self.assertEqual(len(hvs.input_inventory(self.cfg, [22])[1]), 0)
        self.assertNotIn('/chr1.bgz', [name for name, _ in self.requests])

    def test_interrupted_reference_transfer_resumes_using_range(self):
        destination = self.project / 'hg38.fa.gz'
        partial = destination.with_name(destination.name + '.partial')
        payload = (self.remote / 'hg38.fa.gz').read_bytes()
        partial.write_bytes(payload[:12])
        hvs.download_reference(self.base + '/hg38.fa.gz', destination)
        self.assertEqual(destination.read_bytes(), payload)
        self.assertIn(('/hg38.fa.gz', 'bytes=12-'), self.requests)
        self.assertFalse(partial.exists())

    def test_custom_vcf_never_receives_unrelated_public_index(self):
        vcf = self.project / self.cfg['dbsnp_vcf']
        vcf.parent.mkdir(parents=True)
        vcf.write_bytes((self.remote / 'dbsnp.gz').read_bytes())
        with redirect_stdout(io.StringIO()):
            hvs.download_references(self.cfg, [22])
        self.assertNotIn('/dbsnp.gz.tbi', [name for name, _ in self.requests])
        self.assertIn('/chr22.bgz.tbi', [name for name, _ in self.requests])

    def test_failed_transfer_never_installs_completed_reference(self):
        destination = self.project / 'failed.gz'
        with patch.object(hvs.subprocess, 'run', return_value=subprocess.CompletedProcess([], 18)):
            with self.assertRaisesRegex(ValueError, 'resume'):
                hvs.download_reference(self.base + '/hg38.fa.gz', destination)
        self.assertFalse(destination.exists())

    def test_inventory_deduplicates_missing_resources(self):
        _, missing = hvs.input_inventory(self.cfg, [22])
        self.assertEqual(len(missing), len(set(missing)))
        self.assertEqual(missing.count(str(self.project / self.cfg['dbsnp_vcf'])), 1)
