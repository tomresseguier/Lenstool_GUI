from PyQt5.QtWidgets import QWidget, QApplication, QHBoxLayout, QVBoxLayout, QLabel, QGraphicsProxyWidget
from PyQt5.QtCore import Qt, QTimer, QPointF
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore
import numpy as np

from .utils_general import make_handles
from .sliders import TripleSlider
from ...lenstool_model.simulate_image.utils import get_light_model_ranges




class DragWidget(QWidget):
    """
    This class implements the option to create a ROI anywhere by just pressing SHIFT and clicking.
    """
    def __init__(self, qt_image):
        super().__init__()
        self.initUI()
        self.qt_image = qt_image
        self.qt_image.scene.sigMouseMoved.connect(self.mouse_moved)
        
        layout = QHBoxLayout(self) #These lines are necessary for the image to actually display in the QMainWindow later
        layout.addWidget(self.qt_image)
        self.setLayout(layout)
        
        self.cat = None
        self.drawing = False
        self.current_ROI = pg.RectROI([-100, -100], [0, 0], pen='r', invertible=True) #Try initializing with a None
        
    def initUI(self):
        self.timer = QTimer()
        self.timer.setInterval(1)  # Check every 10 ms
        self.timer.timeout.connect(self.checkLongPress)

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton and event.modifiers() == Qt.ShiftModifier :
            try:
                self.qt_image.removeItem(self.current_ROI)
            except Exception:
                pass
            self.qt_image.view.setMouseEnabled(x=False, y=False)
            self.timer.start()
            ###################################################################
            self.start_pos = self.qt_image.view.mapToView(event.pos() + QPointF(-13, -14)) #the QPointF(-13, -14) offset corrects a bug where the anchor point would move.
            self.current_ROI = pg.RectROI([self.start_pos.x(), self.start_pos.y()], [0, 0], pen='r', invertible=True)
            for handle in self.current_ROI.handles :
                self.current_ROI.removeHandle(handle['item'])
            self.qt_image.addItem(self.current_ROI)
            self.drawing = True
    
    def keyPressEvent(self, event) :
        if event.key() == Qt.Key_Escape :
            try:
                self.qt_image.removeItem(self.current_ROI)
            except Exception:
                pass
    
    def checkLongPress(self):
        if not (QApplication.mouseButtons() & Qt.LeftButton) :
            self.drawing = False
            if self.timer.isActive() :
                self.timer.stop()
            make_handles(self.current_ROI)
            self.qt_image.view.setMouseEnabled(x=True, y=True)
            if self.cat is not None :
                self.cat.make_selection_ROI()
            
    def mouse_moved(self, pos):
        if self.drawing and self.current_ROI is not None:
            current_pos = self.qt_image.view.mapToView(pos)
            width = current_pos.x() - self.start_pos.x()
            height = current_pos.y() - self.start_pos.y()
            self.current_ROI.setSize([width, height])
            

class DragPlotWidget(pg.PlotWidget):
    """
    Similar to DragWidget but for pyqtgraph's PlotWidget.
    Let's the user create two types of selectable ROI.
    """
    def __init__(self):
        super().__init__()
        self.scene().sigMouseMoved.connect(self.mouse_moved)
        self.view = self.getViewBox()
        self.view.setAspectLocked(True)
        self.initUI()
        self.last_roi = None
        self.roi_list = []
        self.drawing1 = False
        self.drawing2 = False
        
    def initUI(self):
        self.timer = QTimer()
        self.timer.setInterval(1) #Check every 10 ms
        self.timer.timeout.connect(self.checkLongPress)

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton and event.modifiers() == Qt.ShiftModifier:
            self.view.setMouseEnabled(x=False, y=False)
            self.timer.start()
            mouse_point = self.view.mapSceneToView(event.pos()) # Convert click position to data coordinates
            self.start_pos = QPointF(mouse_point.x(), mouse_point.y())
            self.last_roi = SelectableRectangleROI([self.start_pos.x(), self.start_pos.y()], [1e-9, 1e-9], pen='r', invertible=True)
            self.roi_list.append(self.last_roi)
            for handle in self.last_roi.handles:
                self.last_roi.removeHandle(handle['item'])
            self.addItem(self.last_roi)
            self.drawing1 = True
        elif event.button() == Qt.LeftButton and event.modifiers() == Qt.ControlModifier:
            self.view.setMouseEnabled(x=False, y=False)
            self.timer.start()
            mouse_point = self.view.mapSceneToView(event.pos()) # Convert click position to data coordinates
            self.start_pos = QPointF(mouse_point.x(), mouse_point.y())
            self.last_roi = SelectableEllipseROI([self.start_pos.x(), self.start_pos.y()], radius=1e-9, pen='g', invertible=True)
            self.roi_list.append(self.last_roi)
            for handle in self.last_roi.handles:
                self.last_roi.removeHandle(handle['item'])
            self.addItem(self.last_roi)
            self.drawing2 = True
        else:
            super().mousePressEvent(event) #Default PlotWidget behavior, so that the user is still able to move around by click & drag

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Escape:
            for roi in self.roi_list :
                self.removeItem(roi)
                del roi
            self.last_roi = None
            self.roi_list = []
            event.accept()
        elif event.key() in (QtCore.Qt.Key_Backspace, QtCore.Qt.Key_Delete) :
            roi_list_new = []
            for roi in self.roi_list :
                if roi.selected :
                    self.removeItem(roi)
                    del roi
                else :
                    roi_list_new.append(roi)
            self.roi_list = roi_list_new
            if len(self.roi_list)>0:
                self.last_roi = self.roi_list[-1]
            event.accept()
            
    def checkLongPress(self):
        if not (QApplication.mouseButtons() & Qt.LeftButton) and not (QApplication.mouseButtons() & Qt.ControlModifier):
            self.drawing1 = False
            self.drawing2 = False
            if self.timer.isActive():
                self.timer.stop()
            if len(self.roi_list)>0:
                make_handles(self.last_roi)
            self.view.setMouseEnabled(x=True, y=True)

    def mouse_moved(self, pos):
        if self.drawing1 or self.drawing2:
            mouse_point = self.view.mapSceneToView(pos)
            width = max( mouse_point.x() - self.start_pos.x(), 1e-9, key=abs )
            height = max( mouse_point.y() - self.start_pos.y(), 1e-9, key=abs )
            self.last_roi.setSize([width, height])


class DragPlotWidget_special(pg.PlotWidget):
    """
    Same as DragPlotWidget but with selecatble elliptical and two different types of circular ROI.
    """
    def __init__(self):
        super().__init__()
        self.scene().sigMouseMoved.connect(self.mouse_moved)
        self.view = self.getViewBox()
        self.view.setAspectLocked(True)
        self.initUI()
        self.last_roi = None
        self.roi_list = []
        self.drawing1 = False
        self.drawing2 = False
        self.drawing3 = False
        self.ranges_dict = get_light_model_ranges()
        
    def initUI(self):
        self.timer = QTimer()
        self.timer.setInterval(1) #Check every 10 ms
        self.timer.timeout.connect(self.checkLongPress)

    def mousePressEvent(self, event):
        # Circular Sérsic
        if event.button() == Qt.LeftButton and (event.modifiers() & Qt.ShiftModifier) and (event.modifiers() & Qt.ControlModifier):
            self.view.setMouseEnabled(x=False, y=False)
            self.timer.start()
            mouse_point = self.view.mapSceneToView(event.pos()) # Convert click position to data coordinates
            self.start_pos = QPointF(mouse_point.x(), mouse_point.y())
            self.last_roi = SelectableCircleROI([self.start_pos.x(), self.start_pos.y()], PlotWidget=self, radius=1e-9, pen='b', invertible=True)
            self.last_roi.type = 3
            self.roi_list.append(self.last_roi)
            for handle in self.last_roi.handles:
                self.last_roi.removeHandle(handle['item'])
            self.addItem(self.last_roi)
            self.drawing3 = True
            
            # Add sliders to control light source parameters
            params = ['amp', 'n_sersic']
            for i, param in enumerate(params) : #, 'R_sersic'
                min_range, max_range = self.ranges_dict[param][0], self.ranges_dict[param][1]
                if not hasattr(self.last_roi, 'sliders') :
                    self.last_roi.sliders = {param : TripleSlider(min_range, max_range, PlotWidget=self, label=param, roi=self.last_roi)}
                else :
                    offset = self.last_roi.sliders[params[i-1]].offset + self.last_roi.sliders[params[i-1]].bounding_height + self.last_roi.sliders[params[i-1]].offset_text
                    self.last_roi.sliders[param] = TripleSlider(min_range, max_range, PlotWidget=self, label=param, roi=self.last_roi, offset=offset)
                self.last_roi.sliders[param].setParentItem(self.last_roi)
        # 
        elif event.button() == Qt.LeftButton and event.modifiers() == Qt.AltModifier:
            print('Alt modifier function not yet implemented')
        # Elliptical Sérsic
        elif event.button() == Qt.LeftButton and event.modifiers() == Qt.ShiftModifier:
            self.view.setMouseEnabled(x=False, y=False)
            self.timer.start()
            mouse_point = self.view.mapSceneToView(event.pos()) # Convert click position to data coordinates
            self.start_pos = QPointF(mouse_point.x(), mouse_point.y())
            self.last_roi = SelectableEllipseROI([self.start_pos.x(), self.start_pos.y()], [1e-9, 1e-9], pen='r', invertible=True)
            self.last_roi.type = 1
            self.roi_list.append(self.last_roi)
            for handle in self.last_roi.handles:
                self.last_roi.removeHandle(handle['item'])
            self.addItem(self.last_roi)
            self.drawing1 = True
            
            # Add sliders to control light source parameters
            params = ['amp', 'n_sersic']
            for i, param in enumerate(params) : #, 'R_sersic'
                min_range, max_range = self.ranges_dict[param][0], self.ranges_dict[param][1]
                if not hasattr(self.last_roi, 'sliders') :
                    self.last_roi.sliders = {param : TripleSlider(min_range, max_range, PlotWidget=self, label=param, roi=self.last_roi)}
                else :
                    offset = self.last_roi.sliders[params[i-1]].offset + self.last_roi.sliders[params[i-1]].bounding_height + self.last_roi.sliders[params[i-1]].offset_text
                    self.last_roi.sliders[param] = TripleSlider(min_range, max_range, PlotWidget=self, label=param, roi=self.last_roi, offset=offset)
                self.last_roi.sliders[param].setParentItem(self.last_roi)
        # Gaussian light source
        elif event.button() == Qt.LeftButton and event.modifiers() == Qt.ControlModifier:
            self.view.setMouseEnabled(x=False, y=False)
            self.timer.start()
            mouse_point = self.view.mapSceneToView(event.pos()) # Convert click position to data coordinates
            self.start_pos = QPointF(mouse_point.x(), mouse_point.y())
            self.last_roi = SelectableCircleROI([self.start_pos.x(), self.start_pos.y()], PlotWidget=self, radius=1e-9, pen='g', invertible=True)
            self.last_roi.type = 2
            self.roi_list.append(self.last_roi)
            for handle in self.last_roi.handles:
                self.last_roi.removeHandle(handle['item'])
            self.addItem(self.last_roi)
            self.drawing2 = True
            
            # Add sliders to control light source parameters
            param = 'amp'
            min_range, max_range = self.ranges_dict[param][0], self.ranges_dict[param][1]
            self.last_roi.sliders = {param : TripleSlider(min_range, max_range, PlotWidget=self, label=param, roi=self.last_roi)}
            self.last_roi.sliders[param].setParentItem(self.last_roi)
        else:
            super().mousePressEvent(event) #Default PlotWidget behavior, so that the user is still able to move around by click & drag

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Escape:
            for roi in self.roi_list :
                self.removeItem(roi)
                del roi
            self.last_roi = None
            self.roi_list = []
            event.accept()
        elif event.key() in (QtCore.Qt.Key_Backspace, QtCore.Qt.Key_Delete) :
            roi_list_new = []
            for roi in self.roi_list :
                if roi.selected :
                    self.removeItem(roi)
                    del roi
                else :
                    roi_list_new.append(roi)
            self.roi_list = roi_list_new
            if len(self.roi_list)>0:
                self.last_roi = self.roi_list[-1]
            event.accept()
            
    def checkLongPress(self):
        if not (QApplication.mouseButtons() & Qt.LeftButton) and not (QApplication.mouseButtons() & Qt.ControlModifier):
            self.drawing1 = False
            self.drawing2 = False
            self.drawing3 = False
            if self.timer.isActive():
                self.timer.stop()
            if len(self.roi_list)>0:
                if type(self.last_roi).__name__ == 'SelectableEllipseROI' :
                    make_handles(self.last_roi)
                elif type(self.last_roi).__name__ == 'SelectableCircleROI' :
                    #make_handles(self.roi2, make_rotation_handles=False, fixed_ratio=True)
                    self.last_roi.addScaleHandle([0.5+2**-1.5, 0.5+2**-1.5], [0.5, 0.5])
            self.view.setMouseEnabled(x=True, y=True)

    def mouse_moved(self, pos):
        if self.drawing1 :
            mouse_point = self.view.mapSceneToView(pos)
            width = max( mouse_point.x() - self.start_pos.x(), 1e-9, key=abs )
            height = max( mouse_point.y() - self.start_pos.y(), 1e-9, key=abs )
            self.last_roi.setSize([width, height])
        if self.drawing2 :
            mouse_point = self.view.mapSceneToView(pos)
            w2 = abs(mouse_point.x() - self.start_pos.x())
            h2 = abs(mouse_point.y() - self.start_pos.y())
            radius = max( (w2**2 + h2**2)**0.5, 1e-9 )
            new_corner_x = self.start_pos.x() - radius
            new_corner_y = self.start_pos.y() - radius
            self.last_roi.setPos([new_corner_x, new_corner_y])
            self.last_roi.setSize([2*radius, 2*radius])
        if self.drawing3 :
            mouse_point = self.view.mapSceneToView(pos)
            w2 = abs(mouse_point.x() - self.start_pos.x())
            h2 = abs(mouse_point.y() - self.start_pos.y())
            radius = max( (w2**2 + h2**2)**0.5, 1e-9 )
            new_corner_x = self.start_pos.x() - radius
            new_corner_y = self.start_pos.y() - radius
            self.last_roi.setPos([new_corner_x, new_corner_y])
            self.last_roi.setSize([2*radius, 2*radius])


class SelectableRectangleROI(pg.ROI):
    def __init__(self, pos, size, **kwargs):
        super().__init__(pos, size, **kwargs)
        self.selected = False
        self.normal_pen = self.pen #pg.mkPen('r', width=1)
        self.selected_pen = pg.mkPen('y', width=3)
        self.setPen(self.normal_pen)

    def mouseClickEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton:
            self.selected = not self.selected
            self.setPen(self.selected_pen if self.selected else self.normal_pen)
            event.accept()  # Prevent propagation
        else:
            super().mouseClickEvent(event)

class SelectableEllipseROI(pg.EllipseROI):
    def __init__(self, pos, size, **kwargs):
        super().__init__(pos, size, **kwargs)
        self.selected = False
        self.normal_pen = self.pen #pg.mkPen('r', width=1)
        self.selected_pen = pg.mkPen('y', width=3)
        self.setPen(self.normal_pen)

    def mouseClickEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton:
            self.selected = not self.selected
            self.setPen(self.selected_pen if self.selected else self.normal_pen)
            event.accept()  # Prevent propagation
        else:
            super().mouseClickEvent(event)

class SelectableCircleROI(pg.CircleROI):
    def __init__(self, pos, PlotWidget=None, **kwargs):
        super().__init__(pos, **kwargs)
        self.selected = False
        self.normal_pen = self.pen #pg.mkPen('g', width=1)
        self.selected_pen = pg.mkPen('y', width=3)
        self.setPen(self.normal_pen)
        
        if False :
            #self.viewbox = PlotWidget.getViewBox()
            self.slider = TripleSlider(PlotWidget=PlotWidget, label='param1')
            self.slider.setParentItem(self)
            self.slider2 = TripleSlider(PlotWidget=PlotWidget, label='param2', offset=self.slider.height + self.slider.offset_text)
            self.slider2.setParentItem(self)
        
        
    def _updateLabel(self, vmin, vmid, vmax):
        self.slider_label.setText(f"n_Sersic: [{vmin:.2f}, {vmid:.2f}, {vmax:.2f}]")

    def _updatePanelPos(self):
        """
        Position floating widget BELOW the ROI, but keep it horizontal.
        """
        br = self.boundingRect()
        local_anchor = QPointF(br.left(), br.top())
        scene_pos = self.pos() + local_anchor
        self.proxy.setPos(scene_pos)
        #######################################################################

    def mouseClickEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton:
            self.selected = not self.selected
            self.setPen(self.selected_pen if self.selected else self.normal_pen)
            event.accept()  # Prevent propagation
        else:
            super().mouseClickEvent(event)






