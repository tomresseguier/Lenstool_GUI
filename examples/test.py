from PyQt5.QtWidgets import QApplication
import sys

app = QApplication.instance()
if app is None:
    app = QApplication(sys.argv)







import os
import sys

def initialize_workspace() :
    workdir = os.path.join(os.path.expanduser("~"), 'Library/Mobile Documents/com~apple~CloudDocs/RESEARCH/PROCESS/MACS0308/MACS0308_process')
    sys.path.append(workdir)
    os.chdir(workdir)
initialize_workspace()






import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy import constants as astro_constants
from scipy.integrate import quad
from scipy.stats import mode


DATA_dir = os.path.join( os.path.expanduser("~"), 'RESEARCH_DATA/MACS0308/DATA/' )
Lenstool_dir = os.path.join(os.path.expanduser("~"), 'Library/Mobile Documents/com~apple~CloudDocs/RESEARCH/PROCESS/Lenstool_runs')


sys.path.append( os.path.join(os.path.expanduser("~"), 'Library/Mobile Documents/com~apple~CloudDocs/RESEARCH/PROCESS/') )
from Lenstool_GUI.visualens.fits_image import fits_image










im = fits_image(DATA_dir + "macs0308_rgb.fits")
#im = fits_image(DATA_dir + "macs0308_lw_rgb.fits")
#im = fits_image(os.path.join(os.path.expanduser("~"), 'Library/Mobile Documents/com~apple~CloudDocs/RESEARCH/TRANSPORT') + "/trilogy_full_image_rgb.fits")
#im.plot_image()
im.boost()
im.load_filters()

model_path = Lenstool_dir + "/MACS0308/RUN_031/arc_optimization_v5_TEST/"
model_path = Lenstool_dir + "/MACS0308/RUN_031/arc_optimization_v6_TEST/"
model_path = Lenstool_dir + "/MACS0308/RUN_031/arc_optimization_v4/"


#im.import_lenstool(model_path, use_wrapper=True)
im.import_lenstool(model_path, use_wrapper=False)
im.lt.set_lt_z(6.2)

from PyQt5.QtCore import QPointF
im.image_widget.current_ROI.setState({'pos': QPointF(12011.579155, 20853.592623),
                                      'size': QPointF(-249.243113, 238.214657),
                                      'angle': 0.0})

im.lt.start_simulate_image()

im.lt.imsim.load()
im.lt.imsim.lm_imported.send_to_imsim()



if __name__ == '__main__':
    #app = QApplication(sys.argv)
    window = im.lt.imsim.window
    window.show()
    sys.exit(app.exec())




