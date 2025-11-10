from PyQt5.QtWidgets import QGraphicsEllipseItem, QGraphicsSceneMouseEvent, QApplication
from PyQt5.QtCore import Qt, QObject, QEvent
from pyqtgraph.Qt import QtCore
import PyQt5
import pyqtgraph as pg
import numpy as np
from tqdm import tqdm

from .utils_general import transform_rectangle, InRectangle
from .drag_widgets import DragPlotWidget


class SelectableEllipse(QGraphicsEllipseItem) :
    def __init__(self, x, y, width, height, idx, selection_mask, qtItems, initial_color, selection_color=[255, 255, 255], 
                 scatter_pos=None, Scatter_widget=None, alpha=127, linewidth=3):
        #super(SelectableEllipse, self).__init__(x, y, width, height)
        super().__init__(x, y, width, height)
        self.idx = idx
        self.selection_mask = selection_mask
        self.qtItems = qtItems
        self.linewidth = linewidth
        self.setFlag(QGraphicsEllipseItem.ItemIsSelectable, True)
        
        self.alpha = alpha
        if len(initial_color)==4 :
            self.alpha = initial_color[-1]
            initial_color = initial_color[0:3]
        
        self.selection_color = selection_color
        self.initial_color = initial_color
        self.setPen( pg.mkPen(initial_color + [255], width=linewidth) )
        self.setBrush( pg.mkBrush(initial_color + [self.alpha]) )
        
        self.scatter_pos = scatter_pos
        self.Scatter_widget = Scatter_widget
        self.selection_scatter = None

    def mousePressEvent(self, event: QGraphicsSceneMouseEvent):
        if event.button() == Qt.LeftButton:
            print(f"Ellipse selected: {self.rect().x()}, {self.rect().y()}")
            print('\n########\n Index: ' + str(self.idx) + '\n ########\n')
            
            self.selection_mask[self.idx] = not self.selection_mask[self.idx]
            if self.selection_mask[self.idx] :
                self.qtItems[self.idx].setPen( pg.mkPen(self.selection_color + [255], width=self.linewidth) )
                self.qtItems[self.idx].setBrush( pg.mkBrush(self.selection_color + [self.alpha]) )
            else :
                self.qtItems[self.idx].setPen( pg.mkPen(self.initial_color + [255], width=self.linewidth) )
                self.qtItems[self.idx].setBrush( pg.mkBrush(self.initial_color + [self.alpha]) )
            
            
            if self.Scatter_widget is not None :
                if self.selection_scatter is None :
                    self.selection_scatter = pg.ScatterPlotItem()
                    self.Scatter_widget.addItem(self.selection_scatter)
                if self.selection_mask[self.idx] :
                    self.selection_scatter.setData([self.scatter_pos[0]], [self.scatter_pos[1]], pen=(255, 0, 0), size=10)
                else :
                    self.selection_scatter.setData([], [])
            
            
        
class SelectableScatter(DragPlotWidget) :
    def __init__(self, xy, selection_mask, qtItems=None, color=None, selection_color=[255, 0, 0]) :
        super().__init__()
        self.xy = xy
        self.plot(xy[0], xy[1], pen=None, symbol='o', symbolBrush='g', symbolSize=2)
        
        self.full_mask = selection_mask
        self.qtItems = qtItems
        self.initial_color = color
        self.selection_color = selection_color
        
        self.selection_scatter = pg.ScatterPlotItem()
        self.selection_scatter.setData([],[], pen=(255, 0, 0), size=2)
        self.addItem(self.selection_scatter)
    
    def checkLongPress(self) :
        super().checkLongPress()
        if not (QApplication.mouseButtons() & Qt.LeftButton) and not (QApplication.mouseButtons() & Qt.ControlModifier):
            current_roi = self.last_roi
            def ROI_changed() :
                x0 = current_roi.getState()['pos'][0]
                y0 = current_roi.getState()['pos'][1]
                a = current_roi.getState()['size'][0]
                b = current_roi.getState()['size'][1]
                angle = current_roi.getState()['angle']
                
                rect_params = transform_rectangle(x0, y0, a, b, angle*np.pi/180)
                current_roi.selection_mask = InRectangle(self.xy[0], self.xy[1], rect_params)
                
                mask_matrix = np.full((len(self.roi_list), len(self.xy[0])), False)
                for i, roi in enumerate(self.roi_list) :
                    mask_matrix[i] = roi.selection_mask
                    
                self.full_mask_previous = self.full_mask.copy()
                self.full_mask = np.any(mask_matrix, axis=0)
                to_plot_x = self.xy[0][self.full_mask]
                to_plot_y = self.xy[1][self.full_mask]
                self.selection_scatter.setData(to_plot_x, to_plot_y)                
                
                to_change = np.where( np.logical_xor(self.full_mask_previous, self.full_mask) )[0]
                for i in tqdm(to_change) :
                    if self.full_mask[i] :
                        self.qtItems[i].setPen( pg.mkPen(self.selection_color + [255]) )
                        self.qtItems[i].setBrush( pg.mkBrush(self.selection_color + [127]) )
                    else :
                        self.qtItems[i].setPen( pg.mkPen(self.initial_color + [255]) )
                        self.qtItems[i].setBrush( pg.mkBrush(self.initial_color + [127]) )
                
            ROI_changed()
            current_roi.sigRegionChangeFinished.connect(ROI_changed)
    
    def keyPressEvent(self, event) : #To insure galaxies get deselected when ROIs get deleted
        if event.key() in (QtCore.Qt.Key_Backspace, QtCore.Qt.Key_Delete) :
            for roi in self.roi_list :
                if roi.selected :
                    state = roi.getState()
                    state['size'] = (0,0)
                    roi.setState(state)
        if event.key() == Qt.Key_Escape :
            for roi in self.roi_list :
                state = roi.getState()
                state['size'] = (0,0)
                roi.setState(state)
        super().keyPressEvent(event)
        
        
class SelectSources() : #pg.PlotWidget()
    def __init__(self, cat, qt_image, selection_ROI, selection_mask, selection_regions, window=None, qtItems=None, color=None, selection_color=[255, 0, 0]) :
        self.cat = cat
        self.selection_ROI = selection_ROI
        self.selection_mask = selection_mask
        self.selection_mask_temp = np.full(len(cat), False)
        self.selection_ROI.sigRegionChangeFinished.connect(self.selection_ROI_changed)
        self.selection_scatter = None
        self.qt_image = qt_image
        self.qtItems = qtItems
        self.initial_color = color
        self.selection_color = selection_color
        self.selection_regions = selection_regions
        self.confirm_selection_filter = SelectSources_KeyPressFilter(self.selection_ROI, self.selection_mask, self.selection_mask_temp, \
                                                                     self.qt_image, self.qtItems, self.initial_color, self.selection_regions)
        #self.window = window
        window.installEventFilter(self.confirm_selection_filter)
        self.make_selection()
        
        
    def selection_ROI_changed(self) :
        self.make_selection()
    
        
    def make_selection(self) :
        size_y = self.qt_image.image.shape[0]
        
        x = self.cat['x']
        y = size_y - self.cat['y']
        
        x0 = self.selection_ROI.getState()['pos'][0]
        y0 = self.selection_ROI.getState()['pos'][1]
        a = self.selection_ROI.getState()['size'][0]
        b = self.selection_ROI.getState()['size'][1]
        angle = self.selection_ROI.getState()['angle']
        
        rect_params = transform_rectangle(x0, y0, a, b, angle*np.pi/180)
        full_mask = InRectangle(x, y, rect_params)
        
        self.selection_mask_temp[np.where(full_mask)] = True
        self.selection_mask_temp[np.where( np.logical_not(full_mask) )] = False
        
        
        for qtItem in self.qtItems[~self.selection_mask] :
            qtItem.setPen( pg.mkPen(self.initial_color + [255]) )
            qtItem.setBrush( pg.mkBrush(self.initial_color + [127]) )
        
        for qtItem in self.qtItems[self.selection_mask_temp | self.selection_mask] :
            qtItem.setPen( pg.mkPen(self.selection_color + [255]) )
            qtItem.setBrush( pg.mkBrush(self.selection_color + [127]) )
        
        #for i in np.where(self.selection_mask_temp)[0] :
        #    self.qtItems[i].setPen( pg.mkPen(self.selection_color + [255]) )
        #    self.qtItems[i].setBrush( pg.mkBrush(self.selection_color + [127]) )
    

class SelectSources_KeyPressFilter(QObject) :
    def __init__(self, ROI, selection_mask, selection_mask_temp, qt_image, qtItems, initial_color, selection_regions) :
        super().__init__()
        self.ROI = ROI
        self.selection_mask = selection_mask
        self.selection_mask_temp = selection_mask_temp
        self.qt_image = qt_image
        self.qtItems = qtItems
        self.initial_color = initial_color
        self.selection_regions = selection_regions

    def eventFilter(self, obj, event) :
        if event.type() == QEvent.KeyPress :
            key = event.key()
            if key in [Qt.Key_Enter, Qt.Key_Return, Qt.Key_Space] :
                self.selection_mask[self.selection_mask_temp] = True
                
                x0 = self.ROI.getState()['pos'][0]
                y0 = self.ROI.getState()['pos'][1]
                a = self.ROI.getState()['size'][0]
                b = self.ROI.getState()['size'][1]
                angle = self.ROI.getState()['angle']
                rect_params = transform_rectangle(x0, y0, a, b, angle*np.pi/180)
                self.selection_regions.append(rect_params)
                
            if key in [Qt.Key_Backspace, Qt.Key_Escape, Qt.Key_D] and self.ROI in self.qt_image.getView().allChildren() :
                self.qt_image.removeItem(self.ROI)
                for qtItem in self.qtItems[~self.selection_mask] :
                    qtItem.setPen( pg.mkPen(self.initial_color + [255]) )
                    qtItem.setBrush( pg.mkBrush(self.initial_color + [127]) )
            return True  # Event has been handled
        return False  # Pass the event to the parent









class ellipse_maker_ROI(pg.EllipseROI) :
    def __init__(self, pos, size, qt_image, window, cat, color=[1, 1, 0]) :
        super(ellipse_maker_ROI, self).__init__(pos, size, removable=True)
        self.qt_image = qt_image
        self.color = color
        self.cat = cat
        window.keyPressEvent = self.keyPressEvent.__get__(window, window)
        
    def keyPressEvent(self, event) :
        print(self.pos())
        if event.key() == Qt.Key_Return or event.key() == Qt.Key_Enter :
            x, y = self.pos()
            a, b = self.size()
            theta = self.angle()
            
            color=list(np.array(self.color)*255)
            
            ellipse = QGraphicsEllipseItem(x, y, a, b)
            ellipse.setTransformOriginPoint( PyQt5.QtCore.QPointF(x, y) ) 
            
            ellipse.setRotation(theta)
            ellipse.setPen( pg.mkPen(color + [255]) )
            ellipse.setBrush( pg.mkBrush(color + [127]) )
            
            self.qt_image.addItem(ellipse)
            
            size_y = self.qt_image.image.shape[0]
            y = size_y-y
            theta = -theta
            
            alpha = np.arctan(b/a)
            beta = theta*np.pi/180 - alpha
            r = (a**2 + b**2)**0.5/2
            
            x = x + r*np.cos(beta)
            y = y + r*np.sin(beta)
            
            self.cat.add_row([x, y, a, b, theta])
    
    
    
    
    
            
            