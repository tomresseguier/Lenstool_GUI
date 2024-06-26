from PyQt5.QtWidgets import QGraphicsEllipseItem, QGraphicsScene, QGraphicsView, QGraphicsSceneMouseEvent, QApplication
from PyQt5.QtCore import Qt, QTimer, QPointF, QObject, QEvent
from PyQt5.QtGui import QKeyEvent
from PyQt5.QtWidgets import QWidget
import PyQt5
import pyqtgraph as pg
import numpy as np
from tqdm import tqdm

import os
import sys
module_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(module_dir)

sys.path.append(module_dir + "/utils")
from utils_Qt.utils_general import *



class SelectableEllipse(QGraphicsEllipseItem) :
    def __init__(self, x, y, width, height, idx, selection_mask, qtItems, initial_color, selection_color=[255, 255, 255], scatter_pos=None, RS_widget=None):
        super(SelectableEllipse, self).__init__(x, y, width, height)
        self.idx = idx
        self.selection_mask = selection_mask
        self.qtItems = qtItems
        self.selection_color = selection_color
        self.setFlag(QGraphicsEllipseItem.ItemIsSelectable, True)
        self.initial_color = initial_color
        self.setPen( pg.mkPen(initial_color + [255]) )
        self.setBrush( pg.mkBrush(initial_color + [127]) )
        
        self.scatter_pos = scatter_pos
        self.RS_widget = RS_widget
        self.selection_scatter = None

    def mousePressEvent(self, event: QGraphicsSceneMouseEvent):
        if event.button() == Qt.LeftButton:
            print(f"Ellipse selected: {self.rect().x()}, {self.rect().y()}")
            print('\n########\n Index: ' + str(self.idx) + '\n ########\n')
            
            self.selection_mask[self.idx] = not self.selection_mask[self.idx]
            if self.selection_mask[self.idx] :
                self.qtItems[self.idx].setPen( pg.mkPen(self.selection_color + [255]) )
                self.qtItems[self.idx].setBrush( pg.mkBrush(self.selection_color + [127]) )
            else :
                self.qtItems[self.idx].setPen( pg.mkPen(self.initial_color + [255]) )
                self.qtItems[self.idx].setBrush( pg.mkBrush(self.initial_color + [127]) )
            
            
            if self.RS_widget is not None :
                if self.selection_scatter is None :
                    self.selection_scatter = pg.ScatterPlotItem()
                    self.RS_widget.addItem(self.selection_scatter)
                if self.selection_mask[self.idx] :
                    self.selection_scatter.setData([self.scatter_pos[0]], [self.scatter_pos[1]], pen=(255, 0, 0), size=10)
                else :
                    self.selection_scatter.setData([], [])
            
            
            
        
class SelectableScatter() : #pg.PlotWidget()
    def __init__(self, RS_widget, selection_ROI, data, selection_mask, qtItems=None, color=None, selection_color=[255, 0, 0]) :
        self.selection_ROI = selection_ROI
        self.selection_mask = selection_mask
        self.RS_widget = RS_widget
        self.data = data
        self.selection_ROI.sigRegionChangeFinished.connect(self.selection_ROI_changed)
        self.selection_scatter = None
        
        self.qtItems = qtItems
        self.initial_color = color
        self.selection_color = selection_color
    
    def selection_ROI_changed(self) :
        if self.selection_scatter is not None :
            self.selection_scatter.setData([],[])
        
        x = self.data[0]
        y = self.data[1]
        
        x0 = self.selection_ROI.getState()['pos'][0]
        y0 = self.selection_ROI.getState()['pos'][1]
        a = self.selection_ROI.getState()['size'][0]
        b = self.selection_ROI.getState()['size'][1]
        angle = self.selection_ROI.getState()['angle']
        
        rect_params = transform_rectangle(x0, y0, a, b, angle)
        full_mask = InRectangle(x, y, rect_params)
        
        
        #angle = ((self.selection_ROI.getState()['angle'])*np.pi/180)%(2*np.pi)
        #angle_bis = (np.pi/2-angle)#%(2*np.pi)
        #
        #if angle > np.pi/2 and angle < 3*np.pi/2 :
        #    angle = angle - np.pi
        #    angle_bis = (np.pi/2-angle)#%(2*np.pi)
        #    x0 = x0 - (a*np.cos(angle) - b*np.sin(angle))
        #    y0 = y0 - (a*np.sin(angle) + b*np.cos(angle))
        #mask_x = (x > x0 - (y-y0)/np.tan(angle_bis)) & (x < x0 - (y-y0)/np.tan(angle_bis) + a/np.cos(angle))
        #mask_y = (y > y0 + (x-x0)*np.tan(angle)) & (y < y0 + (x-x0)*np.tan(angle) + b/np.cos(angle))
        #full_mask = mask_x & mask_y
        
        
        self.selection_mask[np.where(full_mask)] = True
        self.selection_mask[np.where( np.logical_not(full_mask) )] = False
        
        self.selection_scatter = pg.ScatterPlotItem()
        to_plot_x = self.data[0][self.selection_mask]
        to_plot_y = self.data[1][self.selection_mask]
        self.selection_scatter.setData(to_plot_x, to_plot_y, pen=(255, 0, 0), size=2)
        
        self.RS_widget.addItem(self.selection_scatter)
        
        
        for i in tqdm(range(len(self.qtItems))) :
            self.qtItems[i].setPen( pg.mkPen(self.initial_color + [255]) )
            self.qtItems[i].setBrush( pg.mkBrush(self.initial_color + [127]) )
        
        for i in np.where(self.selection_mask)[0] :
            self.qtItems[i].setPen( pg.mkPen(self.selection_color + [255]) )
            self.qtItems[i].setBrush( pg.mkBrush(self.selection_color + [127]) )
        
        
        
        
class SelectSources() : #pg.PlotWidget()
    def __init__(self, cat, qt_plot, selection_ROI, selection_mask, selection_regions, window=None, qtItems=None, color=None, selection_color=[255, 0, 0]) :
        self.cat = cat
        self.selection_ROI = selection_ROI
        self.selection_mask = selection_mask
        self.selection_mask_temp = np.full(len(cat), False)
        self.selection_ROI.sigRegionChangeFinished.connect(self.selection_ROI_changed)
        self.selection_scatter = None
        self.qt_plot = qt_plot
        self.qtItems = qtItems
        self.initial_color = color
        self.selection_color = selection_color
        self.selection_regions = selection_regions
        self.confirm_selection_filter = KeyPressFilter(self.selection_ROI, self.selection_mask, self.selection_mask_temp, \
                                                       self.qt_plot, self.qtItems, self.initial_color, self.selection_regions)
        #self.window = window
        window.installEventFilter(self.confirm_selection_filter)
        self.make_selection()
        
        
    def selection_ROI_changed(self) :
        self.make_selection()
    
        
    def make_selection(self) :
        size_y = self.qt_plot.image.shape[0]
        
        x = self.cat['x']
        y = size_y - self.cat['y']
        
        x0 = self.selection_ROI.getState()['pos'][0]
        y0 = self.selection_ROI.getState()['pos'][1]
        a = self.selection_ROI.getState()['size'][0]
        b = self.selection_ROI.getState()['size'][1]
        angle = self.selection_ROI.getState()['angle']
        
        rect_params = transform_rectangle(x0, y0, a, b, angle)
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
    
            
    

class KeyPressFilter(QObject) :
    def __init__(self, ROI, selection_mask, selection_mask_temp, qt_plot, qtItems, initial_color, selection_regions) :
        super().__init__()
        self.ROI = ROI
        self.selection_mask = selection_mask
        self.selection_mask_temp = selection_mask_temp
        self.qt_plot = qt_plot
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
                angle = ((self.ROI.getState()['angle'])*np.pi/180)%(2*np.pi)
                rect_params = transform_rectangle(x0, y0, a, b, angle)
                self.selection_regions.append(rect_params)
                
            if key in [Qt.Key_Backspace, Qt.Key_Escape, Qt.Key_D] and self.ROI in self.qt_plot.getView().allChildren() :
                self.qt_plot.removeItem(self.ROI)
                for qtItem in self.qtItems[~self.selection_mask] :
                    qtItem.setPen( pg.mkPen(self.initial_color + [255]) )
                    qtItem.setBrush( pg.mkBrush(self.initial_color + [127]) )
            return True  # Event has been handled
        return False  # Pass the event to the parent



            
    

class DragWidget(QWidget):
    def __init__(self, qt_plot):#, cat, selection_mask):
        super().__init__()
        self.initUI()
        self.qt_plot = qt_plot
        self.qt_plot.scene.sigMouseMoved.connect(self.mouse_moved)
        self.cat = None
        self.drawing = False
        self.current_ROI = pg.RectROI([-100, -100], [0, 0], pen='r', invertible=True)
        #self.selection_mask = selection_mask
        #self.cat = cat

    def initUI(self):
        self.timer = QTimer()
        self.timer.setInterval(1)  # Check every 10 ms
        self.timer.timeout.connect(self.checkLongPress)

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton and event.modifiers() == Qt.ShiftModifier :
            self.qt_plot.removeItem(self.current_ROI)
            self.qt_plot.view.setMouseEnabled(x=False, y=False)
            self.timer.start()
            ###################################################################
            self.start_pos = self.qt_plot.view.mapToView(event.pos() + QPointF(-13, -14))
            self.current_ROI = pg.RectROI([self.start_pos.x(), self.start_pos.y()], [0, 0], pen='r', invertible=True)
            for handle in self.current_ROI.handles :
                self.current_ROI.removeHandle(handle['item'])
            self.qt_plot.addItem(self.current_ROI)
            self.drawing = True

    def checkLongPress(self):
        if not (QApplication.mouseButtons() & Qt.LeftButton) :
            self.drawing = False
            if self.timer.isActive() :
                self.timer.stop()
            make_handles(self.current_ROI)
            self.qt_plot.view.setMouseEnabled(x=True, y=True)
            if self.cat is not None :
                self.cat.make_cleaner_ROI()
            
    def mouse_moved(self, pos):
        if self.drawing and self.current_ROI is not None:
            current_pos = self.qt_plot.view.mapToView(pos)
            width = current_pos.x() - self.start_pos.x()
            height = current_pos.y() - self.start_pos.y()
            self.current_ROI.setSize([width, height])
            
    
            
    
    
    
class ellipse_maker_ROI(pg.EllipseROI) :
    def __init__(self, pos, size, qt_plot, window, cat, color=[1, 1, 0]) :
        super(ellipse_maker_ROI, self).__init__(pos, size, removable=True)
        self.qt_plot = qt_plot
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
            
            self.qt_plot.addItem(ellipse)
            
            size_y = self.qt_plot.image.shape[0]
            y = size_y-y
            theta = -theta
            
            alpha = np.arctan(b/a)
            beta = theta*np.pi/180 - alpha
            r = (a**2 + b**2)**0.5/2
            
            x = x + r*np.cos(beta)
            y = y + r*np.sin(beta)
            
            self.cat.add_row([x, y, a, b, theta])
    
    
    
    
    
    
    
            
"""
class DragWidget(QWidget):
    def __init__(self, qt_plot):
        super().__init__()
        self.show()
        self.qt_plot = qt_plot        

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_R :
            print('R pressed')

    def keyReleaseEvent(self, event):
        if event.key() == Qt.Key_R :
            print('R released')
"""


    
    
    
    
"""
def make_image_ROI(self) :
    center_y = self.image_data.shape[0]/2
    center_x = self.image_data.shape[1]/2
    self.selection_ROI_image = pg.EllipseROI([center_x-200, center_y-100], [400, 200], removable=True)            
    self.qt_plot.addItem(self.selection_ROI_image)
        
        
        def selection_ROI_image_changed() :
            x, y = self.selection_ROI_image.pos()
            a, b = self.selection_ROI_image.size()
            theta = self.selection_ROI_image.angle()
            
            color=list(np.array(self.color)*255)
            color=list(np.array([1., 0., 1.])*255)
            
            test_ellipse = QGraphicsEllipseItem(x, y, a, b)
            test_ellipse.setTransformOriginPoint( PyQt5.QtCore.QPointF(x, y) ) 
            
            test_ellipse.setRotation(theta)
            test_ellipse.setPen( pg.mkPen(color + [255]) )
            test_ellipse.setBrush( pg.mkBrush(color + [127]) )
            
            self.qt_plot.addItem(test_ellipse)
            
            
            def mousePressEvent(self, event: QGraphicsSceneMouseEvent):
                if event.button() == Qt.LeftButton:
                    
            #Draw ellipse
            #Add ellipse to catalog when enter is pressed or double click in the ellipse (then can use native sellectable ellipse function)
            
            
        self.selection_ROI_image.sigRegionChangeFinished.connect(selection_ROI_image_changed)
    """
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            