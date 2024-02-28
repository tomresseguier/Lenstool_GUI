import numpy as np
import matplotlib.pyplot as plt
from fits_image import fits_image


fig, ax = plt.subplots()
ax.scatter(np.linspace(0,1,100), np.linspace(0,1,100))
plt.show()


image_path = '../../spt0615/DATA/v7-wisps/color/' + 'spt0615_color_rgb.fits'
mult_file_path = '../../spt0615/spt0615_process/SL_models/' + 'spt0615_mult_images.lenstool'



image = fits_image(image_path)
#image.extract_sources()
image.import_multiple_images(mult_file_path)
image.plot_multiple_images()

plt.show()
