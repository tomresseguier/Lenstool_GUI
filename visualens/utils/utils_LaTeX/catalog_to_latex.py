import os
import numpy as np



def catalog_to_latex(cat_instance, columns=[], file_path=None, use_best_fit_value=False, deluxetable=True, blank_str='...') :
    if len(columns)==0 :
        columns = list(cat_instance.cat.colnames)
    
    # Detect triplets/quadruplets involving _16/_50/_84_percentile columns and fold
    # them into a single column with asymmetric uncertainties as super/subscript.
    #
    # use_best_fit_value=True  → value source is the base column (name);
    #                            _50_percentile columns are suppressed entirely.
    # use_best_fit_value=False → value source is name_50_percentile when available,
    #                            falling back to the base column otherwise.
    columns_set = set(columns)
    percentile_map = {}   # base_col -> (col_16, col_value, col_84)
    percentile_cols = set()
    for col in columns :
        if col.endswith('_16_percentile') :
            base   = col[:-len('_16_percentile')]
            col_84 = base + '_84_percentile'
            col_50 = base + '_50_percentile'
            if col_84 not in columns_set :
                continue
            # Require the base column when using best-fit values; otherwise it is
            # optional because _50_percentile can stand in as the display column.
            if use_best_fit_value and base not in columns_set :
                continue
            if not use_best_fit_value and base not in columns_set and col_50 not in columns_set :
                continue

            col_value = base if (use_best_fit_value or col_50 not in columns_set) else col_50
            percentile_map[base] = (col, col_value, col_84)
            percentile_cols.add(col)
            percentile_cols.add(col_84)
            if col_50 in columns_set :
                percentile_cols.add(col_50)

    # When using the median as central value, the base column is not necessarily
    # in the catalog; add a sentinel so display_columns still contains it.
    if not use_best_fit_value :
        for base, (col_16, col_value, col_84) in percentile_map.items() :
            if base not in columns_set :
                columns_set.add(base)
                # Insert the synthetic base name right before its _16_percentile column
                idx = columns.index(col_16)
                columns.insert(idx, base)

    # Merge z_in / z_opt into a single 'z' column.
    # z_in takes priority; z_opt is used (marked with *) when z_in is absent or invalid.
    merge_z = 'z_in' in columns_set and 'z_opt' in columns_set
    suppressed_cols = percentile_cols.copy()
    if merge_z :
        suppressed_cols.add('z_opt')

    # Only display the base columns; percentile partners are merged into them.
    display_columns = [col for col in columns if col not in suppressed_cols]

    def _colhead(col) :
        return 'z' if (merge_z and col == 'z_in') else col

    _DISPLAY_NAMES = {
        'magnification'          : r"$\mu$",
        'shear'                  : r"$\gamma$",
        'tangential_magnification': r"$\mu_\parallel$",
        'radial_magnification'   : r"$\mu_\perp$",
        'f200w_mag'              : r"$\rm{mag}_{F200W}$",
        'f115w_mag'              : r"$\rm{mag}_{F115W}$",
        'mag_F200W'              : r"$\rm{mag}_{F200W}$",
        'mag_F115W'              : r"$\rm{mag}_{F115W}$",
        'ra'                     : r"R.A.",
        'dec'                    : r"Dec"
    }

    def _display_name(col) :
        logical = _colhead(col)
        return _DISPLAY_NAMES.get(logical, logical)

    _UNITS_MAP = {
        'ra'    : 'deg',
        'dec'   : 'deg',
        'theta' : 'deg',
        'time'  : 'days',
        'x'     : 'pix',
        'y'     : 'pix',
        'a'     : 'pix',
        'b'     : 'pix',
    }
    def _units(col) :
        name = _colhead(col)
        if name in _UNITS_MAP :
            return _UNITS_MAP[name]
        return ''

    units_entries = [f"({_units(col)})" if _units(col) else "" for col in display_columns]
    has_units     = any(u for u in units_entries)

    col_fmt = "l" + " c" * len(display_columns)

    def _format_uncertainty(val, upper, lower) :
        """Return a LaTeX math string for val with uncertainties.

        If both uncertainties are zero, the value is displayed without any
        uncertainty notation.  Uses ``\\pm`` when the upper and lower
        uncertainties are equal at the displayed precision; otherwise uses
        super/subscript asymmetric form.
        """
        upper_str = f"{upper:.3g}"
        lower_str = f"{lower:.3g}"
        if upper_str == "0" and lower_str == "0" :
            return f"${val:.3g}$"
        if upper_str == lower_str :
            return f"${val:.3g} \\pm {upper_str}$"
        return f"${val:.3g}^{{+{upper_str}}}_{{-{lower_str}}}$"

    def _build_data_rows() :
        rows = ""
        for src in cat_instance.cat :
            for col in display_columns :
                if merge_z and col == 'z_in' :
                    z_in_val = src['z_in']
                    if not (np.isnan(float(z_in_val)) or float(z_in_val) <= 0) :
                        to_parse = f"{z_in_val:.3g}"
                    else :
                        if 'z_opt' in percentile_map :
                            col_16, col_value, col_84 = percentile_map['z_opt']
                            val    = src[col_value]
                            upper  = src[col_84] - val
                            lower  = val - src[col_16]
                            to_parse = _format_uncertainty(val, upper, lower)
                        else :
                            z_opt_val = src['z_opt']
                            to_parse = f"{z_opt_val:.3g}*"
                elif col in percentile_map :
                    col_16, col_value, col_84 = percentile_map[col]
                    val    = src[col_value]
                    upper  = src[col_84] - val
                    lower  = val - src[col_16]
                    to_parse = _format_uncertainty(val, upper, lower)
                elif type(src[col]) in [np.str_, str] :
                    to_parse = src[col]
                elif col in ['ra', 'dec'] :
                    to_parse = f"{src[col]:#.8g}"
                else :
                    to_parse = f"{src[col]:.3g}"
                if 'np.nan' in to_parse :
                    to_parse = blank_str
                rows += to_parse + " &\n"
            rows = rows[:-3] + "\\\\\n"
        return rows

    if deluxetable :
        # Build both header lines; only emit the units line when at least one
        # column has a known unit (otherwise the second line is all empty).
        col_names_row = " &\n".join(f"\\colhead{{{_display_name(col)}}}" for col in display_columns)
        units_row     = " &\n".join(f"\\colhead{{{u}}}" for u in units_entries)

        table_str = ("\\begin{deluxetable*}{" + col_fmt + "}\n"
                     "\\tablecaption{\\label{tab:}}\n"
                     "\\tablewidth{\\columnwidth}\n"
                     "\\tablehead{\n")
        if has_units :
            table_str += col_names_row + " \\\\\n"
            table_str += units_row + "\n"
        else :
            table_str += col_names_row + "\n"
        table_str += ("}\n"
                      "\\startdata\n")
        table_str += _build_data_rows()
        table_str += ("\\enddata\n"
                      "\\tablecomments{}\n"
                      "\\end{deluxetable*}")
    else :
        col_names_row = " & ".join(_display_name(col) for col in display_columns)
        units_row     = " & ".join(units_entries)

        table_str = ("\\begin{table}\n"
                     "\\centering\n"
                     "\\caption{\\label{tab:}}\n"
                     "\\begin{tabular}{" + col_fmt + "}\n"
                     "\\hline\\hline\n")
        table_str += col_names_row + " \\\\\n"
        if has_units :
            table_str += units_row + " \\\\\n"
        table_str += "\\hline\n"
        table_str += _build_data_rows()
        table_str += ("\\hline\n"
                      "\\end{tabular}\n"
                      "\\tablecomments{}\n"
                      "\\end{table}")
    cat_instance._vprint(table_str)
    
    if file_path is None and cat_instance.path is not None :
        file_path = os.path.join(os.path.dirname(cat_instance.path), 'catalog_latex_table.txt')
    if file_path is not None :
        with open(file_path, 'w') as f :
            f.write(table_str)
        cat_instance._vprint('LaTeX table saved at ' + file_path)
    return table_str
    





import re
import math
from decimal import Decimal, ROUND_HALF_UP


def _round_unc_1sig(unc_str: str):
    """Round a positive uncertainty string to 1 significant figure.
    Returns (Decimal rounded_value, int n_decimal_places)."""
    x = float(unc_str)
    if x == 0:
        return Decimal('0'), 0

    magnitude = math.floor(math.log10(x))
    n_dec = max(0, -magnitude)
    quant = Decimal('1') if n_dec == 0 else Decimal('0.' + '0' * n_dec)
    rounded = Decimal(unc_str).quantize(quant, rounding=ROUND_HALF_UP)

    # Re-derive n_dec from the rounded value: handles cases like 0.097 → 0.10 → 0.1
    if rounded > 0:
        new_magnitude = math.floor(math.log10(float(rounded)))
        n_dec = max(0, -new_magnitude)

    return rounded, n_dec


def _fmt_decimal(value_str: str, n_dec: int) -> str:
    """Round value_str to n_dec decimal places and return as string."""
    quant = Decimal('1') if n_dec == 0 else Decimal('0.' + '0' * n_dec)
    rounded = Decimal(value_str).quantize(quant, rounding=ROUND_HALF_UP)
    if rounded == 0:
        rounded = abs(rounded)  # avoid "-0"
    return str(rounded)


def reformat_latex_uncertainties(table_str: str) -> str:
    """
    Reformat a LaTeX string so that values with uncertainties use the minimal
    number of digits:
      - Uncertainties are rounded to 1 significant figure.
      - Central values are rounded to match that decimal place.
      - Asymmetric $v^{+u}_{-l}$ is collapsed to $v \\pm u$ when u == l after rounding.
      - Already-symmetric $v \\pm u$ entries are also reformatted.
    """

    def _replace_asym(match):
        val_str = match.group(1)
        up_str  = match.group(2).lstrip('+')
        lo_str  = match.group(3).lstrip('-')
        try:
            up_rounded, up_dec = _round_unc_1sig(up_str)
            lo_rounded, lo_dec = _round_unc_1sig(lo_str)
        except Exception:
            return match.group(0)  # leave unchanged on error

        n_dec = max(up_dec, lo_dec)
        quant = Decimal('1') if n_dec == 0 else Decimal('0.' + '0' * n_dec)

        val_fmt = _fmt_decimal(val_str, n_dec)
        up_fmt  = str(up_rounded.quantize(quant))
        lo_fmt  = str(lo_rounded.quantize(quant))

        if up_fmt == lo_fmt:
            return f'${val_fmt} \\pm {up_fmt}$'
        return f'${val_fmt}^{{+{up_fmt}}}_{{-{lo_fmt}}}$'

    def _replace_sym(match):
        val_str = match.group(1)
        unc_str = match.group(2)
        try:
            unc_rounded, n_dec = _round_unc_1sig(unc_str)
        except Exception:
            return match.group(0)

        quant = Decimal('1') if n_dec == 0 else Decimal('0.' + '0' * n_dec)
        val_fmt = _fmt_decimal(val_str, n_dec)
        unc_fmt = str(unc_rounded.quantize(quant))
        return f'${val_fmt} \\pm {unc_fmt}$'

    # Pass 1: asymmetric  $value^{+upper}_{-lower}$
    result = re.sub(
        r'\$(-?[0-9.]+)\^\{(\+[0-9.]+)\}_\{(-[0-9.]+)\}\$',
        _replace_asym,
        table_str,
    )
    # Pass 2: symmetric  $value \pm uncertainty$
    result = re.sub(
        r'\$(-?[0-9.]+)\s*\\pm\s*([0-9.]+)\$',
        _replace_sym,
        result,
    )
    return result

