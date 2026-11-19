"""Fault injection for repeatable submission without a real SLURM installation."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from test_workflow import hvs, MODULE


class FakeSlurm:
    def __init__(self):
        self.jobs = {}
        self.commands = []
        self.cancellations = []
        self.reject_stage = None
        self.interrupt_stage = None
        self.query_failure = False

    def output(self, command, **kwargs):
        if command[0] == 'squeue':
            if self.query_failure:
                raise subprocess.CalledProcessError(1, command, stderr='scheduler unavailable')
            rows = []
            for root, job in self.jobs.items():
                for task, state in job['tasks'].items():
                    if state in ('PENDING', 'RUNNING', 'COMPLETING'):
                        jobid = root if task is None else f'{root}_{task}'
                        rows.append(f'{jobid}|{state}|{job["comment"]}|{job.get("reason", "Priority")}')
            return '\n'.join(rows)
        if command[0] == 'sacct':
            roots = command[command.index('--jobs') + 1].split(',')
            return '\n'.join(f'{root if task is None else root + "_" + str(task)}|{state}|'
                             for root in roots if root in self.jobs
                             for task, state in self.jobs[root]['tasks'].items())
        if command[0] != 'sbatch':
            raise AssertionError(command)
        self.commands.append(command)
        stage = command[-2]
        if stage == self.reject_stage:
            self.reject_stage = None
            raise subprocess.CalledProcessError(1, command, stderr='Please specify --account=desired_account')
        root = str(101 + len(self.jobs))
        array = next((arg.split('=', 1)[1].split('%')[0] for arg in command if arg.startswith('--array=')), None)
        tasks = {int(c): 'PENDING' for c in array.split(',')} if array else {None: 'PENDING'}
        self.jobs[root] = {'stage': stage, 'tasks': tasks, 'command': command,
                           'comment': next(arg.split('=', 1)[1] for arg in command if arg.startswith('--comment='))}
        if stage == self.interrupt_stage:
            self.interrupt_stage = None
            raise KeyboardInterrupt('Launcher terminated after scheduler accepted the job')
        return root + ';cluster\n'

    def run(self, command, **kwargs):
        if command[0] != 'scancel':
            raise AssertionError(command)
        jobid = command[1]
        self.cancellations.append(jobid)
        root, _, task = jobid.partition('_')
        tasks = self.jobs[root]['tasks']
        for c in list(tasks):
            if not task or c == int(task):
                tasks[c] = 'CANCELLED'
        return subprocess.CompletedProcess(command, 0)

    def by_stage(self, stage):
        return [(root, job) for root, job in self.jobs.items() if job['stage'] == stage]


class ResumeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name).resolve()
        shutil.copyfile(MODULE.parent / 'job.sbatch', self.root / 'job.sbatch')
        self.cfg = {'output_dir': 'output/v1', 'picard_command': ['picard'], 'ancestries': {'EUR': ['nfe']}}
        self.chroms = [21, 22]
        self.backend = FakeSlurm()
        self.ready = True
        for context in (patch.object(hvs, 'ROOT', self.root),
                        patch.object(hvs.shutil, 'which', return_value='/mock/tool'),
                        patch.object(hvs, 'check_inputs', side_effect=self.inputs),
                        patch.object(hvs.subprocess, 'check_output', side_effect=self.backend.output),
                        patch.object(hvs.subprocess, 'run', side_effect=self.backend.run),
                        patch.object(hvs, 'sleep')):
            context.start()
            self.addCleanup(context.stop)
        self.environment = patch.dict(hvs.os.environ, {'SBATCH_ACCOUNT': ''})
        self.environment.start()
        self.addCleanup(self.environment.stop)
        self.output = self.root / self.cfg['output_dir']
        self.provenance = self.output / 'provenance'

    def inputs(self, cfg, chroms):
        if not self.ready:
            raise FileNotFoundError('FASTA indexes are missing')
        return {'config': cfg, 'inputs': []}

    def launch(self):
        with redirect_stdout(io.StringIO()):
            hvs.submit(self.cfg, self.chroms, 2)

    def complete(self, stage):
        path = self.output / stage
        path.mkdir(parents=True, exist_ok=True)
        files = []
        if stage == 'candidates':
            files = [f'chr_{c}.tsv.gz' for c in range(1, 23)]
        elif stage == 'merge' or stage.startswith('chr_'):
            files = ['master.GRCh38.tsv.gz', 'EUR.GRCh38.tsv.gz', 'EUR.GRCh37.tsv.gz']
            if stage.startswith('chr_'):
                files += ['EUR.audit.tsv.gz']
        for name in files:
            (path / name).write_bytes(b'completed fixture')
        hvs.atomic_json(path / 'manifest.json', {'config': self.cfg, 'chromosomes': self.chroms, 'summary': {}})

    def test_initial_account_rejection_then_same_command_recovers(self):
        self.ready = False
        self.backend.reject_stage = 'prepare'
        with self.assertRaisesRegex(ValueError, 'run is saved'):
            self.launch()
        self.assertTrue(self.provenance.exists())
        self.assertFalse(self.backend.jobs)
        self.launch()
        self.assertEqual([job['stage'] for job in self.backend.jobs.values()],
                         ['prepare', 'candidates', 'chromosome', 'merge'])
        for command in self.backend.commands:
            self.assertIn('--account=def-sgravel', command)
        count = len(self.backend.commands)
        self.launch()
        self.assertEqual(len(self.backend.commands), count)

    def test_legacy_provenance_without_any_jobs_is_recoverable(self):
        self.provenance.mkdir(parents=True)
        hvs.atomic_json(self.provenance / 'config.json', self.cfg)
        hvs.atomic_json(self.provenance / 'preflight.json', {'inputs': [], 'pending_reference_preparation': True})
        self.launch()
        self.assertEqual(len(self.backend.jobs), 3)
        self.assertNotIn('pending_reference_preparation', json.loads((self.provenance / 'preflight.json').read_text()))

    def test_rejection_halfway_reuses_accepted_upstream_job(self):
        self.backend.reject_stage = 'chromosome'
        with self.assertRaises(ValueError):
            self.launch()
        self.assertEqual(len(self.backend.jobs), 1)
        self.launch()
        self.assertEqual(len(self.backend.by_stage('candidates')), 1)
        self.assertEqual(len(self.backend.by_stage('chromosome')), 1)
        command = self.backend.by_stage('chromosome')[0][1]['command']
        self.assertIn('--dependency=afterok:101', command)

    def test_active_jobs_are_reused_without_duplicate_submission(self):
        self.launch()
        count = len(self.backend.commands)
        for _ in range(3):
            self.launch()
        self.assertEqual(len(self.backend.commands), count)

    def test_accepted_job_without_saved_id_is_adopted_by_comment(self):
        self.backend.interrupt_stage = 'candidates'
        with self.assertRaises(KeyboardInterrupt):
            self.launch()
        self.assertEqual(json.loads((self.provenance / 'jobs.json').read_text()), {})
        self.launch()
        self.assertEqual(len(self.backend.by_stage('candidates')), 1)
        self.assertEqual(len(self.backend.jobs), 3)
        self.assertEqual(json.loads((self.provenance / 'jobs.json').read_text())['candidates'], '101')

    def test_completed_stages_and_failed_array_elements_resume_independently(self):
        self.chroms = [20, 21, 22]
        self.launch()
        self.complete('candidates')
        self.complete('chr_20')
        self.backend.jobs['101']['tasks'][None] = 'COMPLETED'
        self.backend.jobs['102']['tasks'] = {20: 'COMPLETED', 21: 'FAILED', 22: 'RUNNING'}
        self.backend.jobs['103']['reason'] = 'DependencyNeverSatisfied'
        self.launch()
        self.assertEqual(len(self.backend.by_stage('candidates')), 1)
        command = self.backend.by_stage('chromosome')[-1][1]['command']
        self.assertIn('--array=21%2', command)
        merge = self.backend.by_stage('merge')[-1][1]['command']
        dep = next(x for x in merge if x.startswith('--dependency='))
        self.assertIn('102_22', dep.split(':')[1:])
        self.assertIn('104', dep.split(':')[1:])
        self.assertNotIn('102', dep.split(':')[1:])
        self.assertEqual(self.backend.cancellations, ['103'])

    def test_failed_upstream_repairs_pending_dependency_chain(self):
        self.launch()
        self.backend.jobs['101']['tasks'][None] = 'FAILED'
        for root in ('102', '103'):
            self.backend.jobs[root]['reason'] = 'DependencyNeverSatisfied'
        self.launch()
        self.assertEqual(len(self.backend.by_stage('candidates')), 2)
        self.assertEqual(set(self.backend.cancellations), {'102_21', '102_22', '103'})
        new_chr = self.backend.by_stage('chromosome')[-1][1]['command']
        self.assertIn('--dependency=afterok:104', new_chr)

    def test_completed_full_run_is_a_noop(self):
        self.launch()
        for stage in ('candidates', 'chr_21', 'chr_22', 'merge'):
            self.complete(stage)
        count = len(self.backend.commands)
        self.backend.query_failure = True  # Completed outputs do not need the scheduler.
        self.launch()
        self.assertEqual(len(self.backend.commands), count)

    def test_scheduler_outage_never_causes_duplicate_jobs(self):
        self.launch()
        self.backend.query_failure = True
        count = len(self.backend.commands)
        with self.assertRaisesRegex(ValueError, 'Unable to reconcile'):
            self.launch()
        self.assertEqual(len(self.backend.commands), count)

    def test_immutable_worker_config_and_explicit_account_override(self):
        self.launch()
        first = self.backend.jobs['101']['command']
        frozen = Path(first[-4])
        before = frozen.read_bytes()
        self.cfg.update(slurm_account='def-liyue', download_workers=8)
        self.backend.jobs['101']['tasks'][None] = 'FAILED'
        for root in ('102', '103'):
            self.backend.jobs[root]['reason'] = 'DependencyNeverSatisfied'
        self.launch()
        self.assertEqual(frozen.read_bytes(), before)
        new = self.backend.by_stage('candidates')[-1][1]['command']
        self.assertIn('--account=def-liyue', new)
        self.assertNotEqual(new[-4], first[-4])

    def test_concurrent_launchers_do_not_duplicate_jobs(self):
        with ThreadPoolExecutor(max_workers=2) as pool, redirect_stdout(io.StringIO()):
            futures = [pool.submit(hvs.submit, self.cfg, self.chroms, 2) for _ in range(2)]
            for future in futures:
                future.result()
        self.assertEqual(len(self.backend.jobs), 3)

    def test_legacy_failed_submission_is_retried_from_accounting_state(self):
        self.provenance.mkdir(parents=True)
        hvs.atomic_json(self.provenance / 'config.json', self.cfg)
        hvs.atomic_json(self.provenance / 'jobs.json', {'candidates': '87'})
        self.backend.jobs['87'] = {'stage': 'candidates', 'tasks': {None: 'FAILED'}, 'comment': '', 'command': []}
        self.launch()
        self.assertEqual(len(self.backend.by_stage('candidates')), 2)
        self.assertTrue(any(command[0] == 'sbatch' for command in self.backend.commands))

    def test_cli_uses_saved_chromosomes_and_accepts_account_override(self):
        self.chroms = [22]
        self.launch()
        (self.root / 'config.json').write_text(json.dumps(self.cfg))
        with patch.object(hvs, 'load_config', return_value=dict(self.cfg)):
            with patch.object(hvs, 'download_sources'), patch.object(hvs, 'download_references'):
                with patch.object(hvs, 'submit') as plan:
                    with patch('sys.argv', ['workflow.py', '--account', 'def-liyue']):
                        hvs.main()
        self.assertEqual(plan.call_args.args[1], [22])
        self.assertEqual(plan.call_args.args[0]['slurm_account'], 'def-liyue')

    def test_completed_worker_is_idempotent(self):
        self.complete('chr_21')
        with patch.object(hvs, 'process_chromosome', side_effect=AssertionError('Completed work must be reused')):
            with redirect_stdout(io.StringIO()):
                hvs.execute_stage(self.cfg, 'chromosome', 21, self.chroms)
                hvs.execute_stage(self.cfg, 'chromosome', 21, self.chroms)

    def test_incomplete_stage_is_archived_and_rebuilt_automatically(self):
        broken = self.output / 'candidates'
        broken.mkdir(parents=True)
        (broken / 'unfinished.txt').write_text('retained for recovery')
        def rebuild(cfg, directory):
            for c in range(1, 23):
                (directory / f'chr_{c}.tsv.gz').write_bytes(b'completed fixture')
            return {'candidate_rows': 0}
        with patch.object(hvs, 'candidates', side_effect=rebuild):
            with patch.dict(hvs.os.environ, {'SLURM_JOB_ID': '123'}), redirect_stdout(io.StringIO()):
                hvs.execute_stage(self.cfg, 'candidates', None, self.chroms)
        self.assertTrue(hvs.stage_complete(self.output, 'candidates', self.cfg, self.chroms))
        self.assertEqual(len(list((self.output / 'recovery').glob('candidates.*/unfinished.txt'))), 1)

    def test_scientific_changes_are_rejected_without_modifying_results(self):
        self.launch()
        self.complete('candidates')
        self.cfg['maf_gt'] = 0.01
        count = len(self.backend.commands)
        with self.assertRaisesRegex(ValueError, 'Scientific configuration changed'):
            self.launch()
        self.assertEqual(len(self.backend.commands), count)
