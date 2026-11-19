"""Rerun a paper's verification programs and compare output with the recorded runs."""

import argparse
import difflib
import json
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite", nargs="?", choices=["mersenne"],
                        help="run the reciprocal Mersenne cubes checks")
    args = parser.parse_args()
    directory = Path(__file__).resolve().parent
    if args.suite:
        directory = directory / args.suite
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
    if args.suite:
        print(f"All {len(reference['stdout_lines'])} programs reproduced their reference results.")
    else:
        print("Both programs reproduced their reference results.")


if __name__ == "__main__":
    main()
