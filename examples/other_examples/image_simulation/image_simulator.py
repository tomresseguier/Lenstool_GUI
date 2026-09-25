from PyQt5.QtWidgets import QApplication
import sys
import os

app = QApplication.instance()
if app is None:
    app = QApplication(sys.argv)
    
from visualens import Visualens

sys.path.append( os.path.join(os.path.expanduser("~"), 'Library/Mobile Documents/com~apple~CloudDocs/RESEARCH/PROCESS/') )
from Lenstool_GUI.visualens import Visualens

DATA_dir = os.path.join( os.path.dirname(os.getcwd()), 'DATA' )

vl = Visualens()
vl.import_image(DATA_dir + "/RGB_cropped.fits")
vl.image.boost()

vl.image.load_filters()

model_path = DATA_dir + "/lens_model/arc_optimization/"
vl.import_lenstool(model_path)
vl.lens_model.set_lt_z(6.2)

vl.lens_model.plot()

vl.lens_model.start_simulate_image()
