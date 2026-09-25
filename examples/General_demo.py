import matplotlib
from PyQt5.QtWidgets import QApplication
import sys
import os


#from visualens import Visualens
sys.path.append( os.path.join(os.path.expanduser("~"), 'Library/Mobile Documents/com~apple~CloudDocs/RESEARCH/PROCESS/') )
from Lenstool_GUI.visualens import Visualens

#DATA_dir = #os.path.join( os.getcwd(), 'DATA' )
DATA_dir = '/Users/tomresseguier/Library/Mobile Documents/com~apple~CloudDocs/RESEARCH/PROCESS/Lenstool_GUI/examples/DATA'

vl = Visualens()
vl.import_image(DATA_dir + "/RGB_cropped.fits")

vl.import_catalog(DATA_dir + '/phot-eazy_magRS.fits')

vl.catalog.cat

vl.catalog.plot(scale=1., color=[0,1,1], text_column=None, linewidth=3, marker=None)

mask = vl.catalog.selection_mask
print(vl.catalog.cat[mask])

vl.catalog.plot_column('z_phot')

vl.catalog.clear()

vl.catalog.export_to_mult_file()

# Create a new column for a specific quantity we want to look at. Here, color.
vl.catalog.cat['f115w_mag-f200w_mag'] = vl.catalog.cat['f115w_mag'] - vl.catalog.cat['f200w_mag']

# Create the selection panel
vl.catalog.make_selection_panel(xy_axes=['f200w_mag', 'f115w_mag-f200w_mag'])

vl.catalog.export_to_potfile()

vl.catalog.remove_selection_panel()
vl.catalog.clear()


vl.toggle_split4()
