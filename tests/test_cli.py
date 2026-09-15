import json
import subprocess
import sys
from mkin4py.cli import example


def test_cli_local_workflow():
    result = subprocess.run([sys.executable, '-m', 'mkin4py', 'solve', '-', '--json'],
        input=json.dumps(example()), capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)['success']


def test_cli_invalid_input():
    result = subprocess.run([sys.executable, '-m', 'mkin4py', 'solve', '-'],
        input='{}', capture_output=True, text=True)
    assert result.returncode == 2
    assert not json.loads(result.stderr)['success']
