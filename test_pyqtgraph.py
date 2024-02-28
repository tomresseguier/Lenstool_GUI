#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 23:10:04 2023

@author: Tom
"""






import pyqtgraph as pg
import numpy as np
from tqdm import tqdm



pg.setConfigOptions(imageAxisOrder='row-major')

image_data = np.transpose(image.image, axes=[1,0,2])
image_data = np.flip(image.image, axis=0)
image_data = image.image
PlotItem = pg.image(image_data)



semi_major = image.sources['A_IMAGE']
semi_minor = image.sources['B_IMAGE']
angle = image.sources['THETA_IMAGE']

def generate_ellipse_points(x_center, y_center, semi_major, semi_minor, angle, n_points=20):
    t = np.linspace(0, 2 * np.pi, n_points)
    x = x_center + semi_major * np.cos(t) * np.cos(np.radians(angle)) - semi_minor * np.sin(t) * np.sin(np.radians(angle))
    y = y_center + semi_major * np.cos(t) * np.sin(np.radians(angle)) + semi_minor * np.sin(t) * np.cos(np.radians(angle))
    return x, y


x = np.zeros((100,20))
y = np.zeros((100,20))
for i in tqdm(range(100)) :
    x[i], y[i] = generate_ellipse_points(image.sources['RA'][i], image.sources['DEC'][i], semi_major[i], semi_minor[i], angle[i], n_points=20)




scatter = pg.ScatterPlotItem(size=10, brush=pg.mkBrush(255, 255, 255, 120))
pos = np.array([np.ndarray.flatten(x), np.ndarray.flatten(y)])
spots = [{'pos': pos[:, i], 'data': 1}
         for i in range(len(np.ndarray.flatten(x)))] + [{'pos': [0, 0], 'data': 1}]
scatter.addPoints(spots)
PlotItem.addItem(scatter)


p_ellipse = PyQt5.QtWidgets.QGraphicsEllipseItem(20, 20, 10, 10)  # x, y, width, height
p_ellipse.setPen(pg.mkPen((0, 0, 0, 50)))
p_ellipse.setBrush(pg.mkBrush((50, 50, 50, 50)))
PlotItem.addItem(p_ellipse)






polygon_item = pg.ScatterPlotItem(size=10)

(polygon_points, brush=(255, 0, 0, 100))



# importing Qt widgets
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *


class Window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PyQtGraph")
        self.setGeometry(100, 100, 600, 500)
        icon = QIcon("skin.png")
        self.setWindowIcon(icon)
        self.UiComponents()
        self.show()
    
    def UiComponents(self):
        widget = QWidget()
        label = QLabel("Geeksforgeeks Scatter Plot")
        label.setWordWrap(True)
        plot = pg.plot()
        
        n = 20
        scatter = pg.ScatterPlotItem(size=10, brush=pg.mkBrush(255, 255, 255, 120))
        pos = np.array([x[0], y[0]])
        spots = [{'pos': pos[:, i], 'data': 1}
                 for i in range(n)] + [{'pos': [0, 0], 'data': 1}]
        scatter.addPoints(spots)
        plot.addItem(scatter)
        
        layout = QGridLayout()
        label.setMinimumWidth(130)
        widget.setLayout(layout)
        layout.addWidget(label, 1, 0)
        layout.addWidget(plot, 0, 1, 3, 1)
        self.setCentralWidget(widget)
 
 
# create pyqt5 app
App = QApplication(sys.argv)
 
# create the instance of our Window
window = Window()
 
# start the app
sys.exit(App.exec())
