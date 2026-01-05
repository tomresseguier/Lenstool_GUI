import pyqtgraph as pg
from PyQt5.QtCore import Qt, pyqtSignal, QRectF, QPointF
from PyQt5.QtGui import QPainter, QColor, QPen, QBrush


class TripleSlider(pg.GraphicsObject):
    valuesChanged = pyqtSignal(float, float, float)  # min, mid, max
    
    def __init__(self, min_range=0, max_range=10, roi=None):
        super().__init__()
        self.width, self.height = 40, 10
        self.min_range = min_range
        self.max_range = max_range
        
        self.vmin = min_range
        self.vmid = (min_range + max_range) / 2
        self.vmax = max_range
        
        self._value = 0.5
        self.setAcceptHoverEvents(True)
        self.setAcceptedMouseButtons(Qt.LeftButton)
        
        # fraction of width/height for bar/handles
        self._left_pad_frac = 0.05   # 5% padding
        self._right_pad_frac = 0.05
        self._bar_height_frac = 0.1  # fraction of ROI height
        self._handle_width_frac = 0.08
        self._handle_height_frac = 0.3
        self._handle_color = QColor(50, 50, 50)
        self._mid_color = QColor(120, 120, 120)
        
        self.dragging = None
        
        
    def add_label(self) :
        label = pg.TextItem(text="Parameter", anchor=(0.5, 1.0), color=(220, 220, 220))
        label.setParentItem(self)
        label.setPos(slider.width / 2, -4)

    def boundingRect(self):
        return QRectF(0, 0, self.width, self.height)
    
    def val_to_x(self, v):
        w = self.width * (1 - self._left_pad_frac - self._right_pad_frac)
        return self.width * self._left_pad_frac + (v - self.min_range) / (self.max_range - self.min_range) * w

    def x_to_val(self, x):
        w = self.width * (1 - self._left_pad_frac - self._right_pad_frac)
        x = min(max(x - self.width * self._left_pad_frac, 0), w)
        return self.min_range + x / w * (self.max_range - self.min_range)
    
    def paint(self, p, *_):
        # background
        #p.setBrush(QColor(60, 60, 60))
        p.setBrush(QColor(200, 200, 255))
        #p.drawRect(0, 0, self.width, self.height)
        p.drawRoundedRect(self.boundingRect(), 0.2, 0.2)

        # handle
        w = self.width
        h = self.height
        
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

    def _setValue(self, v):
        self._value = max(0.0, min(1.0, v))
        self.update()
        self.valueChanged.emit(self._value)



plot = pg.PlotWidget()
plot.show()

roi = pg.RectROI([1, 1], [3, 2])
plot.addItem(roi)

slider = TripleSlider_ConstantSize(PlotWidget=plot, label='blabla')
slider.setParentItem(roi)
#slider.setZValue(1000)
#slider.setPos(0, 0)
def on_view_changed():
    slider.update_scale()
plot.getViewBox().sigRangeChanged.connect(on_view_changed)
slider.update_scale()


def updateScale():
    vb = plot.getViewBox()
    t = vb.transform()
    inv = QtGui.QTransform.fromScale(1 / t.m11(), 1 / t.m22())
    slider.setTransform(inv)

plot.getViewBox().sigTransformChanged.connect(updateScale)
updateScale()
