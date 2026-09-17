#!/usr/bin/env python3
"""Validate metadata_pysradb_*.tsv files.

Exit code 0 if all valid, 1 if any invalid.
Prints one line per file: "OK <file>" or "BAD <file>: <reason>".

A file is VALID if:
  - non-empty,
  - header contains 'run_accession' and >1 tab-separated column,
  - the entire file parses as TSV without pandas raising,
  - no row has more fields than the header (catches truncation).
"""
import sys, glob, csv, os

def validate(path: str) -> tuple[bool, str]:
    # --- cheap check: empty file ---
    try:
        if os.path.getsize(path) == 0:
            return False, "empty file"
    except OSError as e:
        return False, f"cannot stat: {e}"

    # --- cheap check: header ---
    try:
        with open(path, "r", newline="") as fh:
            header_line = fh.readline()
    except OSError as e:
        return False, f"cannot read: {e}"

    if not header_line:
        return False, "empty file"
    if header_line.startswith("ValueError") or "Unable to parse xml" in header_line \
       or header_line.lstrip().startswith("<?xml"):
        return False, "XML error dump"

    header_cols = header_line.rstrip("\n").split("\t")
    if len(header_cols) < 2:
        return False, f"header has only {len(header_cols)} column(s)"
    if "run_accession" not in header_cols:
        return False, "no run_accession column in header"

    # --- full parse, streaming, no pandas ---
    # csv.reader handles quoted fields with embedded newlines correctly,
    # unlike awk/grep/wc. We only need the row lengths, not the values.
    expected = len(header_cols)
    n_rows = 0
    try:
        with open(path, "r", newline="") as fh:
            reader = csv.reader(fh, delimiter="\t", quotechar='"')
            try:
                next(reader)  # header
            except StopIteration:
                return False, "empty file"
            for row in reader:
                n_rows += 1
                if len(row) != expected:
                    return False, (
                        f"row {n_rows}: expected {expected} fields, "
                        f"got {len(row)} (likely truncated mid-quote)"
                    )
    except csv.Error as e:
        return False, f"csv error at row {n_rows}: {e}"
    except OSError as e:
        return False, f"read error: {e}"

    if n_rows == 0:
        return False, "header only, no data rows"
    return True, f"{n_rows} rows x {expected} cols"


def main():
    if len(sys.argv) > 1:
        files = sys.argv[1:]
    else:
        files = sorted(glob.glob("metadata_pysradb_*.tsv"))

    if not files:
        print("No files to validate.", file=sys.stderr)
        return 1

    bad = []
    for f in files:
        ok, reason = validate(f)
        if ok:
            print(f"OK  {f}  ({reason})")
        else:
            print(f"BAD {f}  ({reason})")
            bad.append(f)

    if bad:
        print(f"\n{len(bad)} invalid file(s):", file=sys.stderr)
        for f in bad:
            print(f"  {f}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())