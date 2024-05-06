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



