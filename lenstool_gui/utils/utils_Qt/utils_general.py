import pyqtgraph as pg
import numpy as np



def make_handles(selection_ROI, make_rotation_handles=True, fixed_ratio=False) :
    while len(selection_ROI.handles)!=0 :
        for handle in selection_ROI.handles :
            selection_ROI.removeHandle(handle['item'])
    
    selection_ROI.addScaleHandle([0.5,0], [0.5,1])
    selection_ROI.addScaleHandle([1,0.5], [0,0.5])
    selection_ROI.addScaleHandle([0.5,1], [0.5,0])
    selection_ROI.addScaleHandle([0,0.5], [1,0.5])
    
    selection_ROI.addScaleHandle([0,0], [1,1])
    selection_ROI.addScaleHandle([0,1], [1,0])
    selection_ROI.addScaleHandle([1,1], [0,0])
    selection_ROI.addScaleHandle([1,0], [0,1])
    
    if make_rotation_handles :
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


def transform_rectangle(x0, y0, a, b, angle) :
    '''
    Transforms the params from the selection rectangle (base corner is where user clicked first)
    to standardized params to be used to make the selection and be saved (base corner is at the top left).
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


def transform_ROI_params(roi) :
    a = roi.size()[0]
    b = roi.size()[1]
    angle = roi.angle() * np.pi/180
    if a>0 :
        beta = (np.arctan(b/a) + angle)
    else :
        beta = (np.arctan(b/a) + angle - np.pi)
    d = (a**2 + b**2)**0.5 / 2
    x_center = roi.x() + d*np.cos(beta)
    y_center = roi.y() + d*np.sin(beta)
    
    if a>b :
        semi_major = abs(a)/2
        semi_minor = abs(b)/2
    else :
        semi_major = abs(b)/2
        semi_minor = abs(a)/2
        angle+=np.pi/2 %np.pi
    return x_center, y_center, semi_major, semi_minor, angle


def InRectangle(x_array, y_array, rect_params) :
    x0, y0, a, b, angle = rect_params
    angle_bis = (np.pi/2-angle)#%(2*np.pi)
    
    mask_x = (x_array > x0 - (y_array-y0)/np.tan(angle_bis)) & (x_array < x0 - (y_array-y0)/np.tan(angle_bis) + a/np.cos(angle))
    mask_y = (y_array > y0 + (x_array-x0)*np.tan(angle)) & (y_array < y0 + (x_array-x0)*np.tan(angle) + b/np.cos(angle))
    full_mask = mask_x & mask_y
    return full_mask


def make_full_color(color) :
    full_color = color.copy()
    if len(full_color)==3 :
        full_color = full_color + [0]
    if len(full_color)==4 :
        edge_color = 0.5 + round(full_color[-1]/2)
        full_color = full_color + [edge_color]
    return list(np.array(full_color)*255)



