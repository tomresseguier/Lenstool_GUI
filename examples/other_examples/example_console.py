from PyQt5.QtWidgets import QApplication
import sys
import os
from visualens.fits_image import fits_image

app = QApplication.instance()
if app is None:
    app = QApplication(sys.argv)


DATA_dir = os.path.join( os.path.realpath(__file__), 'DATA')

im = fits_image(DATA_dir + "/RGB_cropped.fits")
im.boost()
im.load_filters()

model_path = DATA_dir + "/lens_model/"
im.import_lenstool(model_path, use_wrapper=True)
im.lt.set_lt_z(6.2)


#if __name__ == '__main__':
#    #app = QApplication(sys.argv)
#    window = im.lt.imsim.window
#    window.show()
#    sys.exit(app.exec())




