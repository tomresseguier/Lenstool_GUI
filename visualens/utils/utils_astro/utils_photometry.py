from astropy.io import ascii
import os


default_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'NRC_ZPs_1126pmap.txt')
def load_jwst_phot_table(file_path=default_path):
    # We use guess=True but specifically suggest the delimiter
    # 'header_start' and 'data_start' help skip potential comment lines
    table = ascii.read(file_path, 
                       format='fixed_width',
                       delimiter='|',
                       header_start=0)
    
    # Strip any extra whitespace from string columns (like 'pupil+filter')
    for col in table.colnames:
        if table[col].dtype.kind in ['U', 'S']:
            table[col] = [val.strip() for val in table[col]]
    
    return table

