import os
import numpy as np



def catalog_to_latex(cat_instance, columns=[], file_path=None, use_best_fit_value=False, deluxetable=True) :
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
    