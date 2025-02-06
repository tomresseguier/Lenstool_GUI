import os
import glob
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse, Polygon, Circle, Rectangle
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import astropy.constants as c
from reproject import reproject_interp
#from astropy.visualization.wcsaxes import *
from astropy.coordinates import SkyCoord
from astropy.visualization import ZScaleInterval, ImageNormalize
from matplotlib.patches import Ellipse, Polygon, Circle
from tqdm import tqdm
import requests
import types

import PyQt5
from PyQt5.QtWidgets import QMainWindow, QWidget, QHBoxLayout
import pyqtgraph as pg
from PyQt5.QtWidgets import QGraphicsEllipseItem, QGraphicsScene, QGraphicsView, QGraphicsSceneMouseEvent, QApplication
#from PyQt5.QtCore import Qt
from pyqtgraph.Qt import QtCore
from astropy.table import Table
import pandas as pd

###############################################################################
import sys
#module_dir = os.path.dirname(os.path.abspath(__file__))
#sys.path.append(module_dir)

from .source_extraction.source_extract import source_extract, source_extract_DIM
from .source_extraction.match_cat import run_match
#sys.path.append(module_dir + "/utils")
from .utils.utils_plots.plot_utils_general import *
from .utils.utils_Qt.selectable_classes import *
from .utils.utils_Qt.utils_general import *
from .utils.utils_general.utils_general import find_close_coord, make_colnames_dict
from .utils.utils_plots.plt_framework import plt_framework
from .utils.utils_Lenstool.redshift_extractors import make_source_z_dict, find_param_file
from .utils.utils_Lenstool.param_extractors import read_potfile, read_bayes_file, make_param_latex_table
#from utils_general.utils import flux_muJy_to_magAB
from .utils.utils_astro.set_cosmology import set_cosmo
cosmo = set_cosmo()









class lenstool_model :
    def __init__(self, model_path, fits_image) :
        self.model_dir = model_path if os.path.isdir(model_path) else os.path.dirname(model_path)
        all_par_file_paths = glob.glob(os.path.join(self.model_dir, "*.par"))
        if os.path.isfile(model_path) :
            self.param_file_path = model_path
        else :
            for file_path in all_par_file_paths :
                if not os.path.basename(file_path).startswith('best') :
                    with open(file_path, 'r') as file :
                        for line in file :
                            stripped_line = line.strip()
                            if stripped_line and not stripped_line.startswith('#') :  # Skip empty and comment lines
                                if stripped_line.startswith('runmode') :
                                    self.param_file_path = file_path
        all_par_file_names = [ os.path.basename(file_path) for file_path in all_par_file_paths ]
        self.has_run = 'best.par' in all_par_file_names
        self.best_file_path = os.path.join(self.model_dir, 'best.par') if 'best.par' in all_par_file_names else None
        self.bestopt_file_path = os.path.join(self.model_dir, 'bestopt.par') if 'bestopt.par' in all_par_file_names else None
        potfile_paths_list = glob.glob(os.path.join(self.model_dir, "*potfile*.lenstool"))
        self.potfile_path = potfile_paths_list[0] if len(potfile_paths_list)>=1 else None
        
        potfile_Table = read_potfile(self.potfile_path)
        self.potfile = fits_image.make_catalog(potfile_Table, color=[1.,0.,0.], unit_is_pixel=True)
        
        
        
        
        