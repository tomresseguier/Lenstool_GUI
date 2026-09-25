from PyQt5.QtWidgets import QApplication
import sys
import os
from visualens import Visualens

app = QApplication.instance()
if app is None:
    app = QApplication(sys.argv)


DATA_dir = os.path.join( os.path.realpath(__file__), 'DATA')

vl = Visualens()
vl.import_image(DATA_dir + "/RGB_cropped.fits")
vl.image.boost()
vl.image.load_filters()

model_path = DATA_dir + "/lens_model/"
vl.import_lenstool(model_path)
vl.lens_model.set_lt_z(6.2)


#if __name__ == '__main__':
#    #app = QApplication(sys.argv)
#    window = vl.lens_model.imsim.window
#    window.show()
#    sys.exit(app.exec())
