from PyQt5.QtCore import Qt, pyqtSignal, QRectF, QPointF
from PyQt5.QtGui import QPainter, QColor, QPen, QBrush, QTransform, QFont, QFontMetricsF
from PyQt5.QtWidgets import QWidget, QGraphicsItem
import pyqtgraph as pg
import numpy as np




class TripleSlider(pg.GraphicsObject):
    valuesChanged = pyqtSignal(float, float, float)  # min, mid, max
    def __init__(self, min_range=0, max_range=10, label=None, above_text=None, PlotWidget=None, roi=None, offset=0.) :
        super().__init__()
        #self.roi = roi
        self.offset = offset
        self.offset_text = 0. if above_text is None else 20.
        self.slider_width, self.slider_height = 100, 20
        self._bottom_margin = 9
        self.bounding_width, self.bounding_height = self.slider_width, self.slider_height + self._bottom_margin
        self.total_offset = -self.offset - self.offset_text - self.bounding_height
        
        self.min_range = min_range
        self.max_range = max_range

        self.vmin = min_range
        self.vmid = (min_range + max_range) / 2
        self.vmax = max_range

        self.dragging = None

        # fraction of width/height for bar/handles
        self._left_pad_frac = 0.1   # 5% padding
        self._right_pad_frac = 0.1
        self._bar_height_frac = 0.3
        self._handle_width_frac = 0.06
        self._handle_height_frac = 0.6
        self._handle_color = QColor(60, 60, 60)
        self._mid_color = QColor(120, 120, 120)
        
        #self.setFlag(QGraphicsItem.ItemIgnoresTransformations)
        self.setAcceptedMouseButtons(Qt.LeftButton)
        self.setAcceptHoverEvents(True)
        self.setFlag(QGraphicsItem.ItemIsFocusable)
        
        if PlotWidget is not None :
            self.PlotWidget = PlotWidget
            self.viewbox = PlotWidget.getViewBox()
            def data_to_pixel_func() :
                x0 = self.viewbox.mapViewToScene(pg.Point(0, 0)).x()
                x1 = self.viewbox.mapViewToScene(pg.Point(1, 0)).x()
                return x1-x0
            self.data_to_pixel = data_to_pixel_func
            
            # Ensure the anchor point scales with the ROI when handles are dragged
            if roi is not None :
                self.roi = roi
                def anchor_func() :
                    a, b = self.roi.size()/2
                    return (1 + 2**-0.5)*a, (1 - 2**-0.5)*b
                self.get_anchor = anchor_func
                self.roi.sigRegionChanged.connect(self.adjust_all)
            
            # Ensure size is constant
            PlotWidget.getViewBox().sigRangeChanged.connect(self.adjust_all)
            self.adjust_all()
        
        if above_text is not None :
            self.add_above_text(above_text)
            
        if label is not None :
            self.label_str = label + ": "
            self._label_font = QFont()
            self._label_font.setPointSizeF(12)
            self._label_color = QColor(0, 0, 0)
            metrics = QFontMetricsF(self._label_font)
            self._label_width = metrics.horizontalAdvance(self.label_str)
            self.bounding_width += self._label_width
            
            
        self._value_font = QFont()
        self._value_font.setPointSizeF(9)
        self._value_color = QColor(20, 20, 20)
        self._value_fmt = lambda v: f"{v:.2g}".replace("e+0", "e")
        self._value_offset = -1       # pixels above handle
        
    
    def boundingRect(self):
        return QRectF(0, 0, self.bounding_width, self.bounding_height)
    
    def add_above_text(self, text_str) :
        self.above_text = pg.TextItem(text=text_str, anchor=(0, 0), color=(255, 255, 255))
        self.above_text.setParentItem(self)
        self.above_text.setPos( 0, self.bounding_height + self.offset_text )
    
    def update_tranform(self):
        """Update slider rotation to maintain constant horizontal orientation"""
        self.Transform = QTransform.fromScale(1.0/self.data_to_pixel(), 1.0/self.data_to_pixel())
        if self.roi is not None :
            self.Transform.rotate( -self.roi.angle() )
        self.setTransform(self.Transform)
    
    def update_position(self) :
        anchor = self.get_anchor() if self.roi is not None else (0, 0)
        offset = self.total_offset / self.data_to_pixel()
        offset_x = np.sin( np.deg2rad(self.roi.angle()) )*offset
        offset_y = np.cos( np.deg2rad(self.roi.angle()) )*offset
        self.setPos( anchor[0] + offset_x, anchor[1] + offset_y )
    
    def adjust_all(self) :
        """Update slider scale to maintain constant size regardless of zoom level"""
        # Update transform to account for scale and rotation
        self.update_tranform()
        # Move the slider so that it's top matches the anchor point
        self.update_position()
    
    def val_to_x(self, v):
        w = self.slider_width * (1 - self._left_pad_frac - self._right_pad_frac)
        return self._label_width + self.slider_width * self._left_pad_frac + (v - self.min_range) / (self.max_range - self.min_range) * w

    def x_to_val(self, x):
        w = self.slider_width * (1 - self._left_pad_frac - self._right_pad_frac)
        x = min(max(x - self._label_width - self.slider_width * self._left_pad_frac, 0), w)
        return self.min_range + x / w * (self.max_range - self.min_range)

    def paint(self, p, *_):
        # Draw background
        p.setBrush(QColor(200, 200, 255))
        #p.drawRect(0, 0, self.slider_width, self.bounding_height)
        p.drawRoundedRect(self.boundingRect(), self.bounding_height/4, self.bounding_height/4)
        
        # Draw background bar
        bar_left = self._label_width + self.slider_width * self._left_pad_frac
        bar_width = self.slider_width * (1 - self._left_pad_frac - self._right_pad_frac)
        bar_bottom = self.bounding_height - self.slider_height * (1 + self._bar_height_frac)/2
        bar_height = self.slider_height * self._bar_height_frac
        p.setPen(QPen(QColor(0, 0, 120), 1))
        p.setBrush(QBrush(QColor(0, 0, 120)))
        p.drawRect(int(bar_left), int(bar_bottom), int(bar_width), int(bar_height))
        
        # Draw handles
        handle_width = self.slider_width * self._handle_width_frac
        handle_height = self.slider_height * self._handle_height_frac
        handle_bottom = self.bounding_height - self.slider_height/2 - handle_height/2

        p.setPen(Qt.NoPen)
        for v in (self.vmin, self.vmid, self.vmax):
            x = self.val_to_x(v)
            p.setBrush(QBrush(self._mid_color if v == self.vmid else self._handle_color))
            p.drawRect(int(x - handle_width/2), int(handle_bottom), int(handle_width), int(handle_height))
        
        # Draw values under handles (correct text orientation)
        p.save()
        p.scale(1, -1) # Flip text back upright
        
        p.setFont(self._value_font)
        p.setPen(self._value_color)
        metrics = QFontMetricsF(self._value_font)
        
        for v in (self.vmin, self.vmid, self.vmax):
            x = self.val_to_x(v)
            text = self._value_fmt(v)
            text_width = metrics.horizontalAdvance(text)
            text_height = metrics.height()
        
            text_x = x - text_width / 2
            text_y = -(handle_bottom - self._value_offset - text_height)
        
            p.drawText(QPointF(text_x, text_y), text)
        
        # Draw label
        p.setFont(self._label_font)
        p.setPen(self._label_color)
        text_x = self._left_pad_frac*self.slider_width/2
        text_y = -bar_bottom
        p.drawText(QPointF(text_x, text_y), self.label_str)
        
        p.restore()
        
        
    def mousePressEvent(self, e) :
        pos = e.pos()
        dist = {"min": abs(self.val_to_x(self.vmin) - pos.x()),
                "mid": abs(self.val_to_x(self.vmid) - pos.x()),
                "max": abs(self.val_to_x(self.vmax) - pos.x())}
        self.dragging = min(dist, key=dist.get)

    def mouseMoveEvent(self, e) :
        if self.dragging is None :
            return
        val = self.x_to_val(e.pos().x())
        
        eps = 1e-6
        if self.dragging == "min":
            self.vmin = min(val, self.vmid - eps)
        elif self.dragging == "mid":
            self.vmid = max(min(val, self.vmax - eps), self.vmin + eps)
        elif self.dragging == "max":
            self.vmax = max(val, self.vmid + eps)
            
        self.valuesChanged.emit(self.vmin, self.vmid, self.vmax)
        self.update()
        
    def mouseReleaseEvent(self, e):
        self.dragging = None











class TripleSlider_ProportionalSize(QWidget):
    valuesChanged = pyqtSignal(float, float, float)  # min, mid, max
    def __init__(self, min_range, max_range, roi=None, panel=None, viewbox=None, proxy=None) :
        super().__init__()
        self.roi = roi
        self.panel = panel
        self.viewbox = viewbox
        self.proxy = proxy
        
        def data_to_pixel_func() :
            x0 = self.viewbox.mapViewToScene(pg.Point(0, 0)).x()
            x1 = self.viewbox.mapViewToScene(pg.Point(1, 0)).x()
            return x1-x0
        self.data_to_pixel = data_to_pixel_func

        self.min_range = min_range
        self.max_range = max_range

        self.vmin = min_range
        self.vmid = (min_range + max_range) / 2
        self.vmax = max_range

        self.dragging = None

        # fraction of width/height for bar/handles
        self._left_pad_frac = 0.05   # 5% padding
        self._right_pad_frac = 0.05
        self._bar_height_frac = 0.1  # fraction of ROI height
        self._handle_width_frac = 0.08
        self._handle_height_frac = 0.3
        self._handle_color = QColor(50, 50, 50)
        self._mid_color = QColor(120, 120, 120)

        # Connect ROI signal
        if self.roi is not None:
            self.update_geometry()
            self.roi.sigRegionChanged.connect(self.update_geometry)

    def update_geometry(self):
        """Update slider width to match ROI width"""
        if self.roi is not None :
            npix = self.data_to_pixel()
            self.proxy.setScale(1/npix)
            
            br = self.roi.boundingRect()
            self._roi_width = br.width() * npix
            self._roi_height = br.height() * npix
            
            self.panel.setFixedWidth(round(self._roi_width /2))
            self.panel.setFixedHeight(round(self._roi_width /4))
            self.panel.update()

    def val_to_x(self, v):
        w = self.width() * (1 - self._left_pad_frac - self._right_pad_frac)
        return self.width() * self._left_pad_frac + (v - self.min_range) / (self.max_range - self.min_range) * w

    def x_to_val(self, x):
        w = self.width() * (1 - self._left_pad_frac - self._right_pad_frac)
        x = min(max(x - self.width() * self._left_pad_frac, 0), w)
        return self.min_range + x / w * (self.max_range - self.min_range)

    def paintEvent(self, event):
        p = QPainter(self)

        w = self.width()
        h = self.height()
        
        # compute bar dimensions proportional
        bar_left = w * self._left_pad_frac
        bar_width = w * (1 - self._left_pad_frac - self._right_pad_frac)
        bar_top = h * 0.4  # center vertically
        bar_height = h * self._bar_height_frac

        # draw bar
        p.setPen(QPen(Qt.black, 1))
        p.drawRect(int(bar_left), int(bar_top), int(bar_width), int(bar_height))

        # draw handles
        handle_width = w * self._handle_width_frac
        handle_height = h * self._handle_height_frac
        handle_top = (h - handle_height)/2

        p.setPen(Qt.NoPen)
        for v in (self.vmin, self.vmid, self.vmax):
            x = self.val_to_x(v) - handle_width/2
            p.setBrush(QBrush(self._mid_color if v == self.vmid else self._handle_color))
            p.drawRect(int(x), int(handle_top), int(handle_width), int(handle_height))

    def mousePressEvent(self, e):
        pos = e.pos()
        dist = {
            "min": abs(self.val_to_x(self.vmin) - pos.x()),
            "mid": abs(self.val_to_x(self.vmid) - pos.x()),
            "max": abs(self.val_to_x(self.vmax) - pos.x()),
        }
        self.dragging = min(dist, key=dist.get)

    def mouseMoveEvent(self, e):
        if not self.dragging:
            return

        val = self.x_to_val(e.pos().x())

        if self.dragging == "min":
            self.vmin = min(val, self.vmid)
        elif self.dragging == "mid":
            self.vmid = max(min(val, self.vmax), self.vmin)
        elif self.dragging == "max":
            self.vmax = max(val, self.vmid)

        self.valuesChanged.emit(self.vmin, self.vmid, self.vmax)
        self.update()

    def mouseReleaseEvent(self, e):
        self.dragging = None




