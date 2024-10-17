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



def transform_rectangle(x0, y0, a, b, angle) :
    '''
    Transforms the params from the selection rectangle (base corner is where user clicked first)
    to standardized params to be used to make the selection and be saved.
    '''
    angle = angle%(2*np.pi)
    
    if a < 0. :
        a = -a
        x0 = x0 - a*np.cos(-angle)
        y0 = y0 + a*np.sin(-angle)
    if b < 0. :
        b = -b
        x0 = x0 - b*np.sin(-angle)
        y0 = y0 - b*np.cos(-angle)
        
    if angle > np.pi/2 and angle < 3*np.pi/2 :
        angle = angle - np.pi
        x0 = x0 - (a*np.cos(angle) - b*np.sin(angle))
        y0 = y0 - (a*np.sin(angle) + b*np.cos(angle))
    return x0, y0, a, b, angle



def InRectangle(x_array, y_array, rect_params) :
    x0, y0, a, b, angle = rect_params
    angle_bis = (np.pi/2-angle)#%(2*np.pi)
    
    mask_x = (x_array > x0 - (y_array-y0)/np.tan(angle_bis)) & (x_array < x0 - (y_array-y0)/np.tan(angle_bis) + a/np.cos(angle))
    mask_y = (y_array > y0 + (x_array-x0)*np.tan(angle)) & (y_array < y0 + (x_array-x0)*np.tan(angle) + b/np.cos(angle))
    full_mask = mask_x & mask_y
    return full_mask






