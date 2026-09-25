import matplotlib
from PyQt5.QtWidgets import QApplication
import sys
import os
import numpy as np


#from visualens import Visualens
sys.path.append( os.path.join(os.path.expanduser("~"), 'Library/Mobile Documents/com~apple~CloudDocs/RESEARCH/PROCESS/') )
from Lenstool_GUI.visualens import Visualens

#DATA_dir = #os.path.join( os.getcwd(), 'DATA' )
DATA_dir = os.path.join( os.path.expanduser("~"), 'RESEARCH_DATA/MACS0308/DATA/' )
Lenstool_dir = os.path.join(os.path.expanduser("~"), 'Library/Mobile Documents/com~apple~CloudDocs/RESEARCH/PROCESS/Lenstool_runs/')

vl = Visualens()
vl.import_image(DATA_dir + "/RGB_cropped.fits")

vl.image.load_filters()

vl.import_catalog(DATA_dir + '/phot-eazy_magRS.fits')

Lenstool_dir = os.path.join(os.path.expanduser("~"), 'Library/Mobile Documents/com~apple~CloudDocs/RESEARCH/PROCESS/Lenstool_runs/')
vl.import_lens_model(Lenstool_dir + "/MACS0308/v2/" + "r06_1H4G_shear")





vl.import_image(DATA_dir + 'macs0308_rgb.fits')
vl.image.load_filters()

vl.import_catalog(os.path.join(os.path.expanduser("~"), 'Library/Mobile Documents/com~apple~CloudDocs/RESEARCH/PROCESS/Lenstool_GUI/examples/DATA/' + 'phot-eazy_magRS.fits'))
vl.import_lens_model(Lenstool_dir + "/MACS0308/v1/" + "RUN_031_MCMC_without_perturber")



vl.lens_model.start_simulate_image(which_filter="F200W")

vl.lens_model.imsim.load()
vl.lens_model.imsim.lm_imported.send_to_imsim()



vl.lens_model.plot_bayes()
vl.lens_model.plot_burnin()





for z in np.arange(1, 10, 0.5) :
    vl.lens_model.set_lt_z(z)


