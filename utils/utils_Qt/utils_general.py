import pyqtgraph as pg
import numpy as np


def make_handles(selection_ROI) :
    selection_ROI.addScaleHandle([0.5,0], [0.5,1])
    selection_ROI.addScaleHandle([1,0.5], [0,0.5])
    selection_ROI.addScaleHandle([0.5,1], [0.5,0])
    selection_ROI.addScaleHandle([0,0.5], [1,0.5])
    
    selection_ROI.addScaleHandle([0,0], [1,1])
    selection_ROI.addScaleHandle([0,1], [1,0])
    selection_ROI.addScaleHandle([1,1], [0,0])
    selection_ROI.addScaleHandle([1,0], [0,1])
    
    selection_ROI.addScaleRotateHandle([0, 0.125], [1,1])
    selection_ROI.addScaleRotateHandle([0.125, 0], [1,1])
    selection_ROI.addScaleRotateHandle([0.875, 0], [0,1])
    selection_ROI.addScaleRotateHandle([1, 0.125], [0,1])
    selection_ROI.addScaleRotateHandle([0, 0.875], [1,0])
    selection_ROI.addScaleRotateHandle([0.125, 1], [1,0])
    selection_ROI.addScaleRotateHandle([0.875, 1], [0,0])
    selection_ROI.addScaleRotateHandle([1, 0.875], [0,0])
    
    selection_ROI.addScaleRotateHandle([0.5, 1.125], [0.5,0.5])
    selection_ROI.addScaleRotateHandle([1.125, 0.5], [0.5,0.5])
    selection_ROI.addScaleRotateHandle([0.5, -0.125], [0.5,0.5])
    selection_ROI.addScaleRotateHandle([ -0.125, 0.5], [0.5,0.5])



def plot_panel(x, y, image_widget_layout, qt_plot) :
    image_widget_layout.setStretchFactor(qt_plot, 6)
    RS_widget = pg.PlotWidget()
    RS_widget.setTitle('Red sequence')
    
    RS_widget.plot(x, y, pen=None, symbol='o', symbolBrush='g', symbolSize=2)
    RS_widget.setAspectLocked(lock=True, ratio=1)
    RS_widget.autoRange()
    #RS_widget.setSizePolicy(pg.QtWidgets.QSizePolicy.Fixed, pg.QtWidgets.QSizePolicy.Expanding)
    image_widget_layout.addWidget(RS_widget)
    image_widget_layout.setStretchFactor(RS_widget, 4)
    
    center_x = np.mean(x)
    center_y = np.mean(y)
    selection_ROI = pg.ROI([center_x-2, center_y-1], [4, 2], removable=True)
    make_handles(selection_ROI)
    RS_widget.addItem(selection_ROI)
    return RS_widget, selection_ROI





import sys
from PyQt5.QtWidgets import QApplication, QMainWindow
import pyqtgraph as pg
import numpy as np

class PlotWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Plot with ROI')
        self.setGeometry(100, 100, 800, 600)

        self.central_widget = pg.GraphicsLayoutWidget()
        self.setCentralWidget(self.central_widget)

        self.plot_widget = self.central_widget.addPlot()
        self.plot_widget.setLabel('left', 'Y Axis')
        self.plot_widget.setLabel('bottom', 'X Axis')

        # Add some data to the plot
        x = np.linspace(0, 10, 100)
        y = np.sin(x)
        self.plot_widget.plot(x, y, pen='b')

        # Create an ROI object
        self.roi = pg.RectROI([2, -0.5], [3, 1], pen=(0,9))
        self.plot_widget.addItem(self.roi)

        # Connect the signal
        self.roi.sigRegionChangeFinished.connect(self.roi_changed)

    def roi_changed(self):
        print("ROI changed. New state:", self.roi.saveState())

window = PlotWindow()
window.show()
