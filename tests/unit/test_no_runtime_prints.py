"""Guardrail test to prevent executable print calls in runtime package modules."""

from pathlib import Path
import re

import pytest


PRINT_CALL_PATTERN = re.compile(r"^\s*print\(")


@pytest.mark.unit
def test_runtime_modules_do_not_use_executable_print_calls():
    repo_root = Path(__file__).resolve().parents[2]
    package_root = repo_root / "HydroEO"

    offenders = []
    for py_file in package_root.rglob("*.py"):
        rel = py_file.relative_to(repo_root)
        lines = py_file.read_text(encoding="utf-8").splitlines()
        for line_no, line in enumerate(lines, start=1):
            if PRINT_CALL_PATTERN.match(line):
                offenders.append(f"{rel}:{line_no}")

    assert not offenders, "Executable print() calls found:\n" + "\n".join(offenders)
