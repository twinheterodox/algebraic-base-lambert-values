"""Reproduce the two review programs and compare their output to the recorded runs."""

import difflib
import json
from pathlib import Path
import subprocess
import sys


def main():
    directory = Path(__file__).resolve().parent
    reference = json.loads((directory / "reference-results.json").read_text())
    for program, expected in reference["stdout_lines"].items():
        run = subprocess.run(
            [sys.executable, "-E", "-B", str(directory / program)],
            check=True, stdout=subprocess.PIPE, text=True,
        )
        print(run.stdout, end="", flush=True)
        actual = run.stdout.splitlines()
        if actual != expected:
            difference = difflib.unified_diff(
                expected, actual, fromfile="reference", tofile=program, lineterm="",
            )
            raise SystemExit("\n".join(difference))
    print("Both programs reproduced their reference results.")


if __name__ == "__main__":
    main()
