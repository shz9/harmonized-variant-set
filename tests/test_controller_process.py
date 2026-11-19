"""Verify the actual detached-process lock handoff without network or SLURM."""
from contextlib import redirect_stdout
import fcntl
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest
from unittest.mock import patch

from test_workflow import hvs


class ControllerProcessTests(unittest.TestCase):
    def test_detached_child_inherits_lock_and_blocks_duplicate_launcher(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            marker = root / 'started'
            release = root / 'release'
            stub = root / 'stub.py'
            stub.write_text('''import argparse, json, os, pathlib, time
parser = argparse.ArgumentParser()
parser.add_argument('--controller-fd', type=int)
parser.add_argument('--config')
args, extra = parser.parse_known_args()
cfg = json.loads(pathlib.Path(args.config).read_text())
with os.fdopen(args.controller_fd, 'a'):
    pathlib.Path(cfg['marker']).write_text(str(os.getpid()))
    while not pathlib.Path(cfg['release']).exists():
        time.sleep(0.02)
''')
            cfg = {'output_dir': str(root / 'output'), 'marker': str(marker), 'release': str(release),
                   'ancestries': {'EUR': ['nfe']}}
            children = []
            real_popen = subprocess.Popen
            def spawn(*args, **kwargs):
                child = real_popen(*args, **kwargs)
                children.append(child)
                return child
            try:
                with patch.object(hvs, 'ROOT', root), patch.object(hvs, '__file__', str(stub)):
                    with patch.object(hvs.subprocess, 'Popen', side_effect=spawn), redirect_stdout(io.StringIO()):
                        hvs.start_controller(cfg, [22], 2)
                    deadline = time.monotonic() + 5
                    while not marker.exists() and time.monotonic() < deadline:
                        time.sleep(0.02)
                    self.assertTrue(marker.exists(), 'Detached controller did not start')
                    with patch.object(hvs.subprocess, 'Popen') as duplicate, redirect_stdout(io.StringIO()):
                        hvs.start_controller(cfg, [22], 2)
                    duplicate.assert_not_called()
                    with redirect_stdout(io.StringIO()) as status:
                        hvs.show_status(cfg)
                    self.assertIn('running on login node', status.getvalue())
            finally:
                release.touch()
                for child in children:
                    try:
                        child.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        child.terminate()
                        child.wait(timeout=5)
            with (root / 'output/.controller.lock').open('r') as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)

    def test_status_is_read_only_and_does_not_access_network(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = {'output_dir': str(Path(tmp) / 'absent'), 'ancestries': {'EUR': ['nfe']}}
            with patch.object(hvs, 'urlopen', side_effect=AssertionError('Status must be offline')):
                with redirect_stdout(io.StringIO()) as status:
                    hvs.show_status(cfg)
            self.assertIn('Completed chromosomes: 0/22', status.getvalue())
            self.assertFalse(Path(cfg['output_dir']).exists())
