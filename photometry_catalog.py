from astropy.table import Table
import pandas as pd



class photometry_catalog() :
    def __init__(self, photometry_catalog_path) :
        self.photometry_catalog_path = photometry_catalog_path
        def open_cat(photometry_catalog_path) :
            with open(photometry_catalog_path, 'r') as raw_cat :
                first_line = raw_cat.readlines()[0]
            if len(first_line.split()) > len(first_line.split(',')) :
                cat_df = pd.read_csv(photometry_catalog_path, delim_whitespace=True)[1:].apply(pd.to_numeric, errors='coerce')
            else :
                cat_df = pd.read_csv(photometry_catalog_path)[1:].apply(pd.to_numeric, errors='coerce')
            cat = Table.from_pandas(cat_df)
            return cat
        self.photometry_catalog = open_cat(photometry_catalog_path)
        




