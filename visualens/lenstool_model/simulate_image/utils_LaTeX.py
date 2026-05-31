import numpy as np
import math
from astropy.table import Table



kwargs_name_to_model_name_dict = {
    'kwargs_lens': 'lens_model_list',
    'kwargs_source': 'source_light_model_list',
    'kwargs_source_ps': 'point_source_model_list',
}


def _fmt(v):
    if isinstance(v, (float, np.floating)):
        return f'{float(v):.4g}'
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return str(v).replace('_', r'\_')

def _fmt_unc(self, val, lo, hi):
    ref_unc = min(lo, hi)
    if ref_unc <= 0:
        return self._fmt(val)
    mag = math.floor(math.log10(ref_unc))
    decimals = max(0, 1 - mag)
    fs = f'.{decimals}f'
    return rf'${val:{fs}}^{{+{hi:{fs}}}}_{{-{lo:{fs}}}}$'

def _hpd(samples, credible_level=0.6827):
    s = np.sort(samples)
    n = len(s)
    n_in = max(1, int(np.floor(credible_level * n)))
    if n_in >= n:
        return s[0], s[-1]
    widths = s[n_in:] - s[:n - n_in]
    i_min = int(np.argmin(widths))
    return float(s[i_min]), float(s[i_min + n_in])


def make_lm_table(lm_instance, kwargs_name_list=['kwargs_source', 'kwargs_source_ps'], columns=None) :
    params_all = []
    for kwargs_name1 in ['kwargs_MCMC', 'kwargs_fixed'] :
        for kwargs_name2 in kwargs_name_list :
            for i, kwargs in enumerate(lm_instance.local[kwargs_name1][kwargs_name2]) :
                if kwargs_name2=='kwargs_lens' and i==0 : # The INTERPOL base model is not included in the MCMC or fixed parameters
                    continue
                for param in kwargs.keys() :
                    params_all.append(param)
    params_all = list(set(params_all))

    def _param_names(colname) :
        if colname=='x' :
            if np.array([p in params_all for p in ['center_x', 'ra_source']]).all() :
                return ['center_x', 'ra_source']
            elif 'center_x' in params_all :
                return ['center_x']
            elif 'ra_source' in params_all :
                return ['ra_source']
            else :
                return [colname]
        elif colname=='y' :
            if np.array([p in params_all for p in ['center_y', 'dec_source']]).all() :
                return ['center_y', 'dec_source']
            elif 'center_y' in params_all :
                return ['center_y']
            elif 'dec_source' in params_all :
                return ['dec_source']
            else :
                return [colname]
        elif colname=='size' :
            if np.array([p in params_all for p in ['R_sersic', 'sigma']]).all() :
                return ['R_sersic', 'sigma']
            elif 'R_sersic' in params_all :
                return ['R_sersic']
            elif 'sigma' in params_all :
                return ['sigma']
            else :
                return [colname]
        elif colname=='amp' :
            if np.array([p in params_all for p in ['amp', 'source_amp']]).all() :
                return ['amp', 'source_amp']
            elif 'amp' in params_all :
                return ['amp']
            elif 'source_amp' in params_all :
                return ['source_amp']
            else :
                return [colname]
        else :
            return [colname]

    if columns is None :
        params_to_include = params_all
    else :
        # This is to handle the case where the user gives standard column names like 'x', 'y', 'size', that are not in the param_list
        params_to_include = []
        for colname in columns :
            params = _param_names(colname)
            for param in params :
                if param in params_all :
                    params_to_include.append(param)
                else :
                    print(f"Invalid column: {colname}")
    params_to_include = list(set(params_to_include))

    def _standard_colname(param) :
        if param in ['center_x', 'ra_source'] :
            return 'x'
        elif param in ['center_y', 'dec_source'] :
            return 'y'
        elif param in ['R_sersic', 'sigma'] and np.array([p in params_to_include for p in ['R_sersic', 'sigma']]).all() :
            return 'size'
        elif param in ['amp', 'source_amp'] :
            return 'amp'
        else :
            return param

    standard_colnames = list(set([_standard_colname(param) for param in params_to_include]))

    nlines = 0
    for kwargs_name2 in kwargs_name_list :
        nlines += len(lm_instance.world['kwargs'][kwargs_name2])
    
    lm_table = Table()
    for standard_colname in standard_colnames :
        lm_table.add_column(np.full(nlines, np.nan), name=standard_colname)
        lm_table.add_column(np.full(nlines, np.nan), name=f'{standard_colname}_16_percentile')
        lm_table.add_column(np.full(nlines, np.nan), name=f'{standard_colname}_50_percentile')
        lm_table.add_column(np.full(nlines, np.nan), name=f'{standard_colname}_84_percentile')
        lm_table.add_column(np.full(nlines, np.nan), name=f'{standard_colname}_hpd_low')
        lm_table.add_column(np.full(nlines, np.nan), name=f'{standard_colname}_hpd_high')

    id_col = []
    i = 0
    for kwargs_name2 in kwargs_name_list :
        model_name = kwargs_name_to_model_name_dict[kwargs_name2]
        for j, model in enumerate(lm_instance.local['models'][model_name]) :
            if model=='INTERPOL' :
                continue
            id_col.append(model)

            param_kwargs_fixed = lm_instance.local['kwargs_fixed'][kwargs_name2][j]
            for param in param_kwargs_fixed.keys() :
                standard_colname = _standard_colname(param)
                if param in params_to_include :
                    lm_table[standard_colname][i] = param_kwargs_fixed[param]

            param_kwargs_MCMC = lm_instance.local['kwargs_MCMC'][kwargs_name2][j]
            for param in param_kwargs_MCMC.keys() :
                standard_colname = _standard_colname(param)
                if param in params_to_include :
                    lm_table[standard_colname][i] = lm_instance.local['kwargs_opt'][kwargs_name2][j][param]
                    p16, p50, p84 = np.percentile(param_kwargs_MCMC[param], [16, 50, 84])
                    lm_table[f'{standard_colname}_16_percentile'][i] = p16
                    lm_table[f'{standard_colname}_50_percentile'][i] = p50
                    lm_table[f'{standard_colname}_84_percentile'][i] = p84
                    hpd_low, hpd_high = _hpd(param_kwargs_MCMC[param])
                    lm_table[f'{standard_colname}_hpd_low'][i] = hpd_low
                    lm_table[f'{standard_colname}_hpd_high'][i] = hpd_high
            
            i += 1
    lm_table.add_column(id_col, name='id', index=0)
    columns = ['id'] + standard_colnames
    return lm_table, columns


def table_to_LaTeX_str(table, columns=None, uncertainty='best', deluxetable=True) :
    if columns is None :
        columns = ['id']
        for col in table.columns :
            has_suffix = any(col.endswith(s) for s in
                             ['_16_percentile', '_50_percentile', '_84_percentile',
                              '_hpd_low', '_hpd_high'])
            if not has_suffix :
                columns.append(col)

    def _format_uncertainty(val, upper, lower) :
        """Return a LaTeX math string for val with asymmetric uncertainties.

        Uses \\pm when upper == lower at displayed precision; otherwise
        super/subscript form.  Precision is driven by the smaller uncertainty.
        """
        ref_unc = min(abs(upper), abs(lower))
        if ref_unc <= 0 :
            return f'${_fmt(val)}$'
        mag = math.floor(math.log10(ref_unc))
        decimals = max(0, 1 - mag)
        fs = f'.{decimals}f'
        val_str   = f'{val:{fs}}'
        upper_str = f'{upper:{fs}}'
        lower_str = f'{lower:{fs}}'
        if upper_str == lower_str :
            return rf'${val_str} \pm {upper_str}$'
        return rf'${val_str}^{{+{upper_str}}}_{{-{lower_str}}}$'

    def _format_cell(row, col) :
        """Format a single table cell as a LaTeX string.

        Fixed parameters are wrapped in brackets.  Sampled parameters show
        the central value with uncertainties as super/subscript (or \\pm).
        """
        # String columns (e.g. 'id') — no numeric formatting needed.
        raw = row[col]
        try :
            val = float(raw)
        except (ValueError, TypeError) :
            return str(raw).replace('_', r'\_')

        if np.isnan(val) :
            return '--'

        p16_col     = f'{col}_16_percentile'
        p50_col     = f'{col}_50_percentile'
        p84_col     = f'{col}_84_percentile'
        hpd_low_col = f'{col}_hpd_low'
        hpd_high_col = f'{col}_hpd_high'

        # Fixed parameter: no percentile columns, or percentile value is NaN.
        if p16_col not in table.colnames or np.isnan(float(row[p16_col])) :
            return f'[{_fmt(val)}]'

        # Sampled parameter — derive centre and uncertainties per mode.
        if uncertainty == 'hpd' :
            center = val
            lo = center - float(row[hpd_low_col])
            hi = float(row[hpd_high_col]) - center
        elif uncertainty == 'median' :
            center = float(row[p50_col])
            lo = center - float(row[p16_col])
            hi = float(row[p84_col]) - center
        else :  # 'best' (default)
            center = val
            lo = center - float(row[p16_col])
            hi = float(row[p84_col]) - center

        return _format_uncertainty(center, hi, lo)

    def _build_data_rows() :
        rows = ''
        for row in table :
            cells = [_format_cell(row, col) for col in columns]
            rows += ' & '.join(cells) + ' \\\\\n'
        return rows

    display_names = [col.replace('_', r'\_') for col in columns]
    col_fmt_str   = 'l' + ' c' * (len(columns) - 1)

    if deluxetable :
        col_names_row = ' &\n'.join(f'\\colhead{{{name}}}' for name in display_names)
        latex_str = (
            f'\\begin{{deluxetable*}}{{{col_fmt_str}}}\n'
            '\\tablecaption{\\label{tab:}}\n'
            '\\tablewidth{\\columnwidth}\n'
            '\\tablehead{\n'
            + col_names_row + '\n'
            + '}\n'
            '\\startdata\n'
        )
        latex_str += _build_data_rows()
        latex_str += (
            '\\enddata\n'
            '\\tablecomments{}\n'
            '\\end{deluxetable*}'
        )
    else :
        col_names_row = ' & '.join(display_names)
        latex_str = (
            '\\begin{table}\n'
            '\\centering\n'
            '\\caption{\\label{tab:}}\n'
            f'\\begin{{tabular}}{{{col_fmt_str}}}\n'
            '\\hline\\hline\n'
            + col_names_row + ' \\\\\n'
            + '\\hline\n'
        )
        latex_str += _build_data_rows()
        latex_str += (
            '\\hline\n'
            '\\end{tabular}\n'
            '\\tablecomments{}\n'
            '\\end{table}'
        )

    return latex_str


def make_lm_LaTeX(lm_instance, kwargs_name_list=['kwargs_source', 'kwargs_source_ps'], columns=None, uncertainty='best', deluxetable=True) :
    """
    Args:
        lm_instance : lenstronomy_model instance
        kwargs_name_list : list of strings, the names of the kwargs to include in the table
        columns : list of strings, the columns to include in the table
        uncertainty : string, the uncertainty to use in the table
                    'hpd' : 68.3% HPD intervals
                    'best' : best-fit value and 16th/84th percentiles
                    'median' : median and 16th/84th percentiles
    Returns:
        latex_str : string, the LaTeX string of the table
    """
    lm_table, columns = make_lm_table(lm_instance, kwargs_name_list=kwargs_name_list, columns=columns)
    latex_str = table_to_LaTeX_str(lm_table, columns=columns, uncertainty=uncertainty, deluxetable=deluxetable)
    return latex_str



