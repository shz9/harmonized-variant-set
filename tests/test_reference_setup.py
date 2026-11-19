"""Reference staging and full orchestration from an empty data directory."""
from contextlib import redirect_stdout
import gzip
import base64
import hashlib
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
        self.cfg.update(picard_command=['picard'], ancestries={'EUR': ['nfe']},
                        dbsnp_vcf='data/reference/dbsnp.gz', scratch_dir=str(self.root / 'scratch'), download_workers=2)
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
        failed_starts = self.failed_starts = set()
        requests = self.requests = []

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def do_HEAD(self):
                self.do_GET(head=True)

            def do_GET(self, head=False):
                requests.append((self.path, self.headers.get('Range')))
                filename = remote / self.path.lstrip('/')
                if not filename.is_file():
                    self.send_error(404)
                    return
                data = filename.read_bytes()
                start = int(self.headers['Range'].split('=')[1].split('-')[0]) if 'Range' in self.headers else 0
                ranged = 'Range' in self.headers
                end_text = self.headers.get('Range', '=0-').split('-')[1]
                end = int(end_text) if end_text else len(data) - 1
                if ranged and start in failed_starts:
                    self.send_error(503)
                    return
                self.send_response(206 if ranged else 200)
                if ranged:
                    self.send_header('Content-Range', f'bytes {start}-{end}/{len(data)}')
                self.send_header('ETag', '"' + hashlib.md5(data).hexdigest() + '"')
                self.send_header('x-goog-hash', 'crc32c=fixture')
                self.send_header('x-goog-hash', 'md5=' + base64.b64encode(hashlib.md5(data).digest()).decode())
                self.send_header('Content-Length', str(end - start + 1))
                self.end_headers()
                if not head:
                    self.wfile.write(data[start:end + 1])

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
        if 'figshare' not in url and 'ebi.ac.uk' not in url:
            return urlopen(request, **kwargs)
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
            if command[0] in ('squeue', 'sacct'):
                return ''
            self.assertEqual(command[0], 'sbatch')
            jobs.append(command)
            return str(100 + len(jobs))
        def worker(stage, frozen, fail_once=False):
            args = ['workflow.py', '--config', str(frozen), '--stage', stage, '--chromosomes', '22']
            if stage == 'chromosome':
                args += ['--chrom', '22']
            with patch.dict(os.environ, {'SLURM_JOB_ID': '123', 'SLURM_TMPDIR': str(self.root)}):
                with patch.object(hvs, 'urlopen', side_effect=AssertionError('HTTP forbidden on compute node')):
                    with patch.object(hvs, 'parallel_download', side_effect=AssertionError('Download forbidden on compute node')):
                        with patch.object(hvs, 'download_reference', side_effect=AssertionError('curl forbidden on compute node')):
                            with patch('sys.argv', args), redirect_stdout(io.StringIO()):
                                if fail_once:
                                    with patch.object(hvs, 'candidate_counts', side_effect=ValueError('interrupted extraction')):
                                        with self.assertRaisesRegex(ValueError, 'interrupted extraction'):
                                            hvs.main()
                                hvs.main()
        with patch.object(hvs, 'urlopen', side_effect=self.stage_discovery_request):
            with patch.object(hvs.shutil, 'which', return_value='/mock/tool'):
                with patch.object(hvs.subprocess, 'check_output', side_effect=external):
                    with patch.object(hvs, 'REFERENCE_CHR1_LENGTHS', {'reference_grch38': 10, 'reference_grch37': 10}):
                        with patch.object(hvs, 'run', side_effect=self.native_tools):
                            # Execute controller cycles in-process so worker completions can be simulated.
                            with patch.object(hvs, 'start_controller', side_effect=hvs.submit):
                                with patch('sys.argv', ['workflow.py', '--config', str(self.config), '--chromosomes', '22']):
                                    with redirect_stdout(io.StringIO()):
                                        hvs.main()
                            self.assertEqual(len(jobs), 2)  # Only prepare and candidates: no VCF staged yet.
                            self.assertNotIn('/chr22.bgz', [name for name, _ in self.requests])
                            output = self.project / self.cfg['output_dir']
                            provenance = output / 'provenance'
                            frozen = provenance / 'config.json'
                            cfg = json.loads(frozen.read_text())
                            Path(str(hvs.local(cfg['dbsnp_vcf'])) + '.tbi').unlink()
                            worker('prepare', frozen)
                            worker('candidates', frozen)
                            # This transfer is explicitly outside the SLURM environment.
                            with redirect_stdout(io.StringIO()):
                                hvs.stage_gnomad(cfg, 22)
                                hvs.submit(cfg, [22], 2)
                            self.assertEqual(len(jobs), 4)
                            self.assertIn('--dependency=afterok:101', jobs[1])
                            self.assertIn('--dependency=afterok:103', jobs[3])
                            worker('chromosome', frozen, fail_once=True)
                            worker('merge', frozen)
                            count = len(self.requests)
                            submitted = len(jobs)
                            with patch.object(hvs, 'start_controller', side_effect=hvs.submit):
                                with patch('sys.argv', ['workflow.py', '--config', str(self.config), '--chromosomes', '22']):
                                    with redirect_stdout(io.StringIO()):
                                        hvs.main()
                            self.assertEqual(len(jobs), submitted)
                            self.assertEqual(len(self.requests), count)
        for build in ('GRCh38', 'GRCh37'):
            exported = pd.read_csv(output / 'merge' / f'EUR.{build}.tsv.gz', sep='\t')
            self.assertEqual(exported.VARIANT_KEY.tolist(), ['chr22:100:A:G'])
        counts = output / 'gnomad_counts' / 'chr_22' / 'counts.tsv.gz'
        self.assertEqual(pd.read_csv(counts, sep='\t').AC_nfe.tolist(), [21])
        self.assertFalse(list(Path(cfg['scratch_dir']).rglob('*.vcf.bgz')))
        self.assertIn('--cpus-per-task=1', jobs[2])
        self.assertIn('--array=22%2', jobs[2])
        self.assertEqual(len(hvs.input_inventory(cfg, [22])[1]), 0)
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
        self.assertNotIn('/chr22.bgz.tbi', [name for name, _ in self.requests])

    def test_failed_transfer_never_installs_completed_reference(self):
        destination = self.project / 'failed.gz'
        with patch.object(hvs.subprocess, 'run', return_value=subprocess.CompletedProcess([], 18)):
            with self.assertRaisesRegex(ValueError, 'resume'):
                hvs.download_reference(self.base + '/hg38.fa.gz', destination)
        self.assertFalse(destination.exists())

    def test_parallel_slices_resume_without_retransferring_completed_ranges(self):
        destination = self.root / 'parallel.bgz'
        self.failed_starts.add(32)
        with patch.object(hvs, 'sleep'), redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError, 'Range 32-63 failed'):
                hvs.parallel_download(self.base + '/chr22.bgz', destination, workers=3, chunk_size=32)
        statefile = destination.with_name(destination.name + '.ranges.json')
        state = json.loads(statefile.read_text())
        self.assertTrue(state['completed'])
        completed = set(state['completed'])
        self.failed_starts.clear()
        self.requests.clear()
        with redirect_stdout(io.StringIO()):
            receipt = hvs.parallel_download(self.base + '/chr22.bgz', destination, workers=3, chunk_size=32)
        self.assertTrue(receipt['md5_verified'])
        self.assertEqual(destination.read_bytes(), (self.remote / 'chr22.bgz').read_bytes())
        retried = {int(r.split('=')[1].split('-')[0]) // 32 for _, r in self.requests if r}
        self.assertTrue(completed.isdisjoint(retried))

    def test_count_cache_reuse_preserves_existing_local_vcf(self):
        source = self.project / self.cfg['gnomad_vcf'].format(chrom=22)
        source.parent.mkdir(parents=True)
        source.write_bytes((self.remote / 'chr22.bgz').read_bytes())
        Path(str(source) + '.tbi').write_bytes((self.remote / 'chr22.bgz.tbi').read_bytes())
        reference = self.project / self.cfg['reference_grch38']
        reference.parent.mkdir(parents=True)
        reference.write_bytes(gzip.decompress((self.remote / 'hg38.fa.gz').read_bytes()))
        pysam.faidx(str(reference))
        master = pd.DataFrame([['chr22:100:A:G', 'chr22', 100, 'A', 'G']], columns=['VARIANT_KEY', *hvs.KEY])
        work = self.root / 'work'
        work.mkdir()
        output = self.project / self.cfg['output_dir']
        with patch.object(hvs, 'run', side_effect=self.native_tools):
            first = hvs.candidate_counts(self.cfg, 22, master, work, output)
        with patch.object(hvs, 'run', side_effect=AssertionError('Cache must skip processing')):
            second = hvs.candidate_counts(self.cfg, 22, master, work, output)
        self.assertEqual(first, second)
        self.assertTrue(source.exists())
        self.assertEqual(second[0]['chr22:100:A:G']['AC_nfe'], 21)

    def test_extraction_failure_preserves_scratch_vcf_until_valid_cache_exists(self):
        reference = self.project / self.cfg['reference_grch38']
        reference.parent.mkdir(parents=True)
        reference.write_bytes(gzip.decompress((self.remote / 'hg38.fa.gz').read_bytes()))
        pysam.faidx(str(reference))
        master = pd.DataFrame([['chr22:100:A:G', 'chr22', 100, 'A', 'G']], columns=['VARIANT_KEY', *hvs.KEY])
        work = self.root / 'work'
        work.mkdir()
        output = self.project / self.cfg['output_dir']
        with redirect_stdout(io.StringIO()):
            hvs.stage_gnomad(self.cfg, 22)
        with patch.object(hvs, 'run', side_effect=self.native_tools):
            with patch.object(hvs, 'normalize', side_effect=ValueError('synthetic extraction failure')):
                with self.assertRaisesRegex(ValueError, 'synthetic extraction failure'):
                    hvs.candidate_counts(self.cfg, 22, master, work, output)
        self.assertTrue(list(Path(self.cfg['scratch_dir']).rglob('*.vcf.bgz')))
        self.assertFalse((output / 'gnomad_counts' / 'chr_22').exists())
        with patch.object(hvs, 'run', side_effect=self.native_tools):
            frequencies, _ = hvs.candidate_counts(self.cfg, 22, master, work, output)
        self.assertEqual(frequencies['chr22:100:A:G']['AN_nfe'], 1000)
        self.assertFalse(list(Path(self.cfg['scratch_dir']).rglob('*.vcf.bgz')))

    def test_parallel_download_rejects_incorrect_object_checksum(self):
        destination = self.root / 'bad.bgz'
        payload = (self.remote / 'chr22.bgz').read_bytes()
        # Preallocated bytes and a forged all-complete checkpoint must still pass checksum validation.
        url = self.base + '/chr22.bgz'
        source = {'url': url, 'size': len(payload), 'etag': '"' + hashlib.md5(payload).hexdigest() + '"',
                  'last_modified': None, 'remote_hash': 'crc32c=fixture, md5=' + base64.b64encode(hashlib.md5(payload).digest()).decode(),
                  'chunk_size': len(payload)}
        destination.with_name('bad.bgz.partial').write_bytes(b'X' * len(payload))
        state = destination.with_name('bad.bgz.ranges.json')
        state.write_text(json.dumps({'source': source, 'completed': [0]}))
        with self.assertRaisesRegex(ValueError, 'MD5 mismatch'):
            hvs.parallel_download(url, destination, workers=2, chunk_size=len(payload))
        self.assertFalse(destination.exists())
        self.assertFalse(state.exists())

    def test_inventory_deduplicates_missing_resources(self):
        _, missing = hvs.input_inventory(self.cfg, [22])
        self.assertEqual(len(missing), len(set(missing)))
        self.assertEqual(missing.count(str(self.project / self.cfg['dbsnp_vcf'])), 1)
