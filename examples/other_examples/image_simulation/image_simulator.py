from PyQt5.QtWidgets import QApplication
import sys
import os

app = QApplication.instance()
if app is None:
    app = QApplication(sys.argv)
    
from visualens.fits_image import fits_image

sys.path.append( os.path.join(os.path.expanduser("~"), 'Library/Mobile Documents/com~apple~CloudDocs/RESEARCH/PROCESS/') )
from Lenstool_GUI.visualens.fits_image import fits_image

DATA_dir = os.path.join( os.path.dirname(os.getcwd()), 'DATA' )

im = fits_image(DATA_dir + "/RGB_cropped.fits")
im.boost()

im.load_filters()

model_path = DATA_dir + "/lens_model/arc_optimization/"
im.import_lenstool(model_path, use_wrapper=True)
im.lt.set_lt_z(6.2)

im.lt.plot()

im.lt.start_simulate_image()