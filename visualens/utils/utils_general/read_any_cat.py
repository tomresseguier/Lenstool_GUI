import re
import numpy as np
from astropy.table import Table, Column


def _is_skipped_line(line):
    """True for blank lines or comment lines starting with # or --."""
    stripped = line.strip()
    return not stripped or stripped.startswith('#') or stripped.startswith('--')


def _strip_leading_hashes(line):
    """Remove one or more leading '#' (may be glued to first column)."""
    stripped = line.strip()
    i = 0
    while i < len(stripped) and stripped[i] == '#':
        i += 1
    return stripped[i:].strip()


def _strip_inline_comment(line):
    """Remove trailing inline comments like '... 0  #2.2  #note'."""
    stripped = line.strip()
    if stripped.startswith('#') or stripped.startswith('--'):
        return stripped
    return re.split(r'\s+#', stripped)[0].strip()


def _detect_separator(line):
    """Pick whitespace vs comma from whichever yields more columns."""
    ws_cols = line.split()
    comma_cols = [c.strip() for c in line.split(',') if c.strip() != '']
    if len(ws_cols) >= len(comma_cols):
        return 'ws'
    return ','


def _split_fields(line, sep):
    line = _strip_inline_comment(line)
    if sep == 'ws':
        return line.split()
    return [c.strip() for c in line.split(',') if c.strip() != '']


def read_any_cat(path):
    """
    Read a whitespace- or comma-separated text file into an astropy Table.

    Header detection:
    - Find the first non-comment line (# or --).
    - Look backward among prior '#' lines for one with the same column count
      (after stripping leading '#' characters).
    - If found: that line is the header, first non-comment line starts data.
    - Otherwise: first non-comment line is the header.
    """
    with open(path, 'r') as f:
        lines = f.readlines()

    # 1) First non-comment line
    first_idx = None
    for i, line in enumerate(lines):
        if not _is_skipped_line(line):
            first_idx = i
            break

    if first_idx is None:
        return Table()

    first_content = _strip_inline_comment(lines[first_idx])
    sep = _detect_separator(first_content)
    ncols = len(_split_fields(first_content, sep))

    # 2) Look backward for a matching '#' header line
    header_idx = None
    for j in range(first_idx - 1, -1, -1):
        stripped = lines[j].strip()
        if not stripped.startswith('#'):
            continue
        content = _strip_leading_hashes(stripped)
        if not content:
            continue
        if len(_split_fields(content, sep)) == ncols:
            header_idx = j
            break  # closest match above first data line

    # 3) Decide header + data start
    if header_idx is not None:
        colnames = _split_fields(_strip_leading_hashes(lines[header_idx].strip()), sep)
        data_start = first_idx
    else:
        colnames = _split_fields(lines[first_idx], sep)
        data_start = first_idx + 1

    # 4) Parse data rows
    rows = []
    for i in range(data_start, len(lines)):
        stripped = lines[i].strip()
        if not stripped:
            continue
        if stripped.startswith('--'):
            continue

        if stripped.startswith('#'):
            content = _strip_leading_hashes(stripped)
            if not content:
                continue
            cols = _split_fields(content, sep)
        else:
            cols = _split_fields(stripped, sep)

        if len(cols) != ncols:
            continue
        rows.append(cols)

    # 5) Build astropy Table with numeric conversion where possible
    table = Table()
    for k, name in enumerate(colnames):
        col_vals = [row[k] for row in rows]
        try:
            col_array = np.array(col_vals, dtype=float)
        except ValueError:
            col_array = np.array(col_vals, dtype=str)
        table[name] = Column(col_array)
        
    return table
