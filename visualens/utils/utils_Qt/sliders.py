from PyQt5.QtCore import Qt, pyqtSignal, QRectF, QPointF
from PyQt5.QtGui import QPainter, QColor, QPen, QBrush, QTransform, QFont, QFontMetricsF
from PyQt5.QtWidgets import QWidget, QGraphicsItem
import pyqtgraph as pg
import numpy as np

from .utils_general import transform_ROI_params, transform_ROI_params_inverse




class TripleSlider(pg.GraphicsObject):
    valuesChanged = pyqtSignal(float, float, float)  # min, mid, max
    def __init__(self, min_range=0, max_range=10, label=None, above_text=None,
                 PlotWidget=None, roi=None, offset=0., log_scale=False,
                 single_handle=False) :
        super().__init__()
        #self.debug = 0
        self.param_name = label
        #self.roi = roi
        self.offset = offset
        self.offset_text = 0. if above_text is None else 20.
        self.slider_width, self.slider_height = 100, 20
        self._bottom_margin = 9
        self.eps = 1e-8
        self.bounding_width, self.bounding_height = self.slider_width, self.slider_height + self._bottom_margin
        self.total_offset = -self.offset - self.offset_text - self.bounding_height
        
        self.min_range = min_range
        self.max_range = max_range
        self.log_scale = log_scale
        # If True, only the middle handle (vmid) is shown and draggable
        self.single_handle = single_handle

        self.vmin = min_range
        self.vmid = (min_range + max_range) / 2 if not log_scale else 10**((np.log10(min_range) + np.log10(max_range)) / 2)
        self.vmax = max_range

        self.dragging = None
        # When True, keep the slider widget position fixed while dragging a handle.
        # This prevents jitter when ROI geometry changes during slider drags.
        self._freeze_position_while_dragging = True
        self._detached_from_roi = False
        self._syncing_from_roi = False
        self._syncing_from_slider = False
        self._siblings_hidden = False
        self._hidden_sibling_sliders = []

        # fraction of width/height for bar/handles
        self._left_pad_frac = 0.1   # 5% padding
        self._right_pad_frac = 0.1
        self.slider_width_effective = self.slider_width * (1 - self._left_pad_frac - self._right_pad_frac)
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
                    a, b = self.roi.size() / 2.0
                    return (1 + 1/np.sqrt(2))*a, (1 - 1/np.sqrt(2))*b
                self.get_anchor = anchor_func
                self.roi.sigRegionChanged.connect(self.adjust_all)
                # Optional 2-way binding between ROI size and vmid for sigma/R_sersic
                self.roi.sigRegionChanged.connect(self._update_from_roi)
            
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
        self._value_font_bold = QFont()
        self._value_font_bold.setPointSizeF(9)
        self._value_font_bold.setBold(True)
        self._value_color = QColor(20, 20, 20)
        self._value_fmt = lambda v: f"{v:.2g}".replace("e+0", "e")
        self._value_offset = -1       # pixels above handle

        # Initialize vmid from ROI geometry if applicable
        self._update_from_roi()


    def _update_from_roi(self):
        """If bound sigma, R_sersic, q or phi, update vmid from ROI geometry."""
        #print('Debug: running _update_from_roi')
        if self._syncing_from_slider:
            return
        if getattr(self, "roi", None) is None or getattr(self, "param_name", None) is None:
            return
        
        roi_type = getattr(getattr(self, "roi", None), "type_str", None)
        param = self.param_name

        if param == "sigma" and roi_type == "GAUSSIAN":
            new_vmid = abs(self.roi.size()[0])
        elif param == "R_sersic" and roi_type == "SERSIC" :
            new_vmid = abs(self.roi.size()[0]) / 2.0
        elif param == "R_sersic" and roi_type == "SERSIC_ELLIPSE_Q_PHI" :
            _, _, semi_major, semi_minor, _ = transform_ROI_params(self.roi)
            # keep consistent with event_filters: R_sersic=(semi_major+semi_minor)/2
            new_vmid = np.sqrt(semi_major * semi_minor)
        elif param == "q" and roi_type == "SERSIC_ELLIPSE_Q_PHI":
            _, _, semi_major, semi_minor, _ = transform_ROI_params(self.roi)
            if semi_major <= 0:
                return
            new_vmid = semi_minor / semi_major
        elif param == "phi" and roi_type == "SERSIC_ELLIPSE_Q_PHI":
            _, _, _, _, angle = transform_ROI_params(self.roi)
            new_vmid = float(np.arctan2(np.sin(angle), np.cos(angle)))
        else:
            return

        self._syncing_from_roi = True
        try:
            new_vmid = max(min(new_vmid, self.vmax - self.eps), self.vmin + self.eps)
            self.vmid = new_vmid
            self.update()
        finally:
            self._syncing_from_roi = False

    def _apply_to_roi(self):
        """If bound sigma/R_sersic/q/phi, apply vmid to ROI geometry."""
        if self._syncing_from_roi:
            return
        if getattr(self, "roi", None) is None or getattr(self, "param_name", None) is None:
            return

        roi_type = getattr(getattr(self, "roi", None), "type_str", None)
        param = self.param_name

        if param == "sigma" and roi_type == "GAUSSIAN":
            new_radius = float(self.vmid) / 2.0
            #mode = "circle_diameter"
        elif param == "R_sersic" and roi_type in ("SERSIC", "SERSIC_ELLIPSE_Q_PHI"):
            new_radius = float(self.vmid)
            #mode = "sersic_radius"
        elif roi_type == "SERSIC_ELLIPSE_Q_PHI" and param in ("q", "phi"):
            pass
        else:
            return

        self._syncing_from_slider = True
        try:
            if type(self.roi).__name__ == "SelectableCircleROI":
                # CircleROI pos is top-left; size is [diameter, diameter]
                current_radius = float(abs(self.roi.size()[0])) / 2.0
                center_x = self.roi.pos()[0] + current_radius
                center_y = self.roi.pos()[1] + current_radius
                #if mode == "sersic_radius":
                #    #current code uses size()[0] as R_sersic for circle Sersic
                #    new_d = new_radius
                #else:
                #    new_d = new_radius
                #new_r = new_d / 2.0
                new_diameter = new_radius * 2.0
                self.roi.setPos([center_x - new_radius, center_y - new_radius])
                self.roi.setSize([new_diameter, new_diameter])
            else:
                # EllipseROI: scale semi-major/minor proportionally to change the mean radius
                x_center, y_center, semi_major, semi_minor, angle = transform_ROI_params(self.roi)
                if param == "R_sersic":
                    current_R = np.sqrt(semi_major * semi_minor)
                    if current_R <= 0:
                        return
                    scale = new_radius / current_R
                    semi_major_new = semi_major * scale
                    semi_minor_new = semi_minor * scale
                    angle_new = angle
                elif param == "q":
                    q_new = max(min(float(self.vmid), 1.0), self.eps)
                    current_R = np.sqrt(semi_major * semi_minor)
                    if current_R <= 0:
                        return
                    semi_major_new = current_R / np.sqrt(q_new)
                    semi_minor_new = current_R * np.sqrt(q_new)
                    angle_new = angle
                elif param == "phi":
                    semi_major_new = semi_major
                    semi_minor_new = semi_minor
                    angle_new = float(np.arctan2(np.sin(self.vmid), np.cos(self.vmid)))
                else:
                    return
                x, y, a, b, angle_deg = transform_ROI_params_inverse(
                    x_center, y_center, semi_major_new, semi_minor_new, angle_new
                )
                self.roi.setPos([x, y])
                self.roi.setSize([a, b])
                self.roi.setAngle(angle_deg)
        finally:
            self._syncing_from_slider = False
        
    
    def boundingRect(self):
        return QRectF(0, 0, self.bounding_width, self.bounding_height)
    
    def add_above_text(self, text_str) :
        self.above_text = pg.TextItem(text=text_str, anchor=(0, 0), color=(255, 255, 255))
        self.above_text.setParentItem(self)
        self.above_text.setPos( 0, self.bounding_height + self.offset_text )

    def hide_slider(self):
        """Hide this slider (and its optional above_text) and disable interaction."""
        self.setVisible(False)
        self.setAcceptedMouseButtons(Qt.NoButton)
        if hasattr(self, "above_text"):
            self.above_text.setVisible(False)

    def show_slider(self):
        """Show this slider (and its optional above_text) and re-enable interaction."""
        self.setVisible(True)
        self.setAcceptedMouseButtons(Qt.LeftButton)
        if hasattr(self, "above_text"):
            self.above_text.setVisible(True)

    def _is_roi_resizing_slider(self) :
        """Return True if dragging this slider can resize the ROI."""
        roi_type = getattr(getattr(self, "roi", None), "type_str", None)
        if self.param_name == "sigma" and roi_type == "GAUSSIAN":
            return True
        if self.param_name == "R_sersic" and roi_type in ("SERSIC", "SERSIC_ELLIPSE_Q_PHI"):
            return True
        if self.param_name in ("q", "phi") and roi_type == "SERSIC_ELLIPSE_Q_PHI":
            return True
        return False

    def _hide_sibling_sliders(self):
        """Hide all sliders attached to the ROI except this one."""
        if self._siblings_hidden:
            return
        roi = getattr(self, "roi", None)
        if roi is None:
            return
        sliders = getattr(roi, "sliders", None)
        if not isinstance(sliders, dict):
            return

        self._hidden_sibling_sliders = []
        for _, s in sliders.items():
            if s is self:
                continue
            # Only hide if currently visible (so we don't override external state)
            try:
                if s.isVisible():
                    s.hide_slider()
                    self._hidden_sibling_sliders.append(s)
            except Exception:
                continue
        self._siblings_hidden = True

    def _show_sibling_sliders(self):
        """Restore the sliders that were hidden by _hide_sibling_sliders()."""
        if not self._siblings_hidden:
            return
        for s in self._hidden_sibling_sliders:
            try:
                s.show_slider()
            except Exception:
                pass
        self._hidden_sibling_sliders = []
        self._siblings_hidden = False
    
    def update_tranform(self):
        """Update slider rotation to maintain constant horizontal orientation"""
        self.Transform = QTransform.fromScale(1.0/self.data_to_pixel(), 1.0/self.data_to_pixel())
        # When detached from ROI we want the slider unrotated (angle = 0)
        if self.roi is not None and not getattr(self, "_detached_from_roi", False) :
            self.Transform.rotate( -self.roi.angle() )
        
        #I'm trying this line to see if it prevents the random rare crashes:
        self.prepareGeometryChange()
        self.setTransform(self.Transform) # This might be the line that causes the annoying crash sometimes
    
    def update_position(self) :
        anchor = self.get_anchor() if self.roi is not None else (0, 0)
        offset = self.total_offset / self.data_to_pixel()
        offset_x = np.sin( np.deg2rad(self.roi.angle()) )*offset
        offset_y = np.cos( np.deg2rad(self.roi.angle()) )*offset
        self.setPos( anchor[0] + offset_x, anchor[1] + offset_y )
    
    def adjust_all(self) :
        """Update slider scale to maintain constant size regardless of zoom level"""
        #print('Debug: running adjust_all', self.roi, self.param_name, self.debug)
        #self.debug += 1
        # Update transform to account for scale and rotation
        self.update_tranform()
        # Move the slider so that it's top matches the anchor point
        if not (self._freeze_position_while_dragging and self.dragging is not None):
            self.update_position()

    def _detach_from_roi_keep_view_pos(self):
        """
        Temporarily detach from ROI so ROI transforms (move/scale/rotate) don't move the slider.
        Keeps the slider at its current ViewBox (data) position.
        """
        if self._detached_from_roi:
            return
        if getattr(self, "PlotWidget", None) is None or getattr(self, "viewbox", None) is None:
            return
        if getattr(self, "roi", None) is None:
            return

        # Current origin of this item in scene coords -> convert to view (data) coords
        scene_pt = self.mapToScene(QPointF(0, 0))
        view_pt = self.viewbox.mapSceneToView(scene_pt)

        # Detach and add to the plot so it's independent of ROI parenting
        self.setParentItem(None)
        try:
            self.PlotWidget.addItem(self)
        except Exception:
            pass
        self.setPos(view_pt)
        self._detached_from_roi = True
        # Reset slider angle to zero while detached
        self.update_tranform()

    def _reattach_to_roi(self):
        """Re-attach to ROI after drag ends and snap back to anchored position."""
        if not self._detached_from_roi:
            return
        if getattr(self, "PlotWidget", None) is None or getattr(self, "roi", None) is None:
            return

        try:
            self.PlotWidget.removeItem(self)
        except Exception:
            pass
        self.setParentItem(self.roi)
        self._detached_from_roi = False
        self.adjust_all()

    def val_to_x(self, v):
        # Map a value to a normalized [0, 1] coordinate, optionally in log space.
        if not self.log_scale:
            v_clipped = min(max(v, self.min_range), self.max_range)
            x_normalized = (v_clipped - self.min_range) / (self.max_range - self.min_range)
        else :
            v_pos = max(v, self.eps, self.min_range if self.min_range > 0 else self.eps)
            v_log = np.log10(v_pos)
            min_log = np.log10(max(self.min_range, self.eps))
            max_log = np.log10(max(self.max_range, self.min_range + self.eps))
            x_normalized = (v_log - min_log) / (max_log - min_log)
        
        # Convert fraction to full x coordinate with respect to slider box
        x = self._label_width + self.slider_width * self._left_pad_frac + x_normalized * self.slider_width_effective
        return x

    def x_to_val(self, x):
        # Convert full x coordinate with respect to slider box to normalized [0, 1] coordinate
        x_relative = min(max(x - self._label_width - self.slider_width * self._left_pad_frac, 0), self.slider_width_effective)
        x_normalized = x_relative / self.slider_width_effective

        # Map normalized [0, 1] coordinate to value.
        v_clipped = min(max(x_normalized, 0.0), 1.0)
        if not self.log_scale:
            v = self.min_range + v_clipped * (self.max_range - self.min_range)
        else :
            min_log = np.log10(max(self.min_range, self.eps))
            max_log = np.log10(max(self.max_range, self.min_range + self.eps))
            v_log = min_log + v_clipped * (max_log - min_log)
            v = 10 ** v_log
        
        return v

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
        handles = (self.vmid,) if self.single_handle else (self.vmin, self.vmid, self.vmax)
        for v in handles:
            x = self.val_to_x(v)
            p.setBrush(QBrush(self._mid_color if v == self.vmid else self._handle_color))
            p.drawRect(int(x - handle_width/2), int(handle_bottom), int(handle_width), int(handle_height))
        
        # Draw values under handles (correct text orientation)
        p.save()
        p.scale(1, -1) # Flip text back upright
        
        p.setFont(self._value_font)
        p.setPen(self._value_color)
        metrics = QFontMetricsF(self._value_font)
        
        for i, v in enumerate(handles) :
            show_value_mid = (i==1 or self.single_handle) and not (self.dragging == "min" or self.dragging == "max")
            show_value_minmax = (i==0 and self.dragging == "min") or (i==2 and self.dragging == "max")
            show_value = show_value_mid or show_value_minmax
            if show_value :
                x = self.val_to_x(v)
                text = self._value_fmt(v)
                text_width = metrics.horizontalAdvance(text)
                text_height = metrics.height()
            
                text_x = x - text_width / 2
                text_y = -(handle_bottom - self._value_offset - text_height)
                
                p.setFont(self._value_font_bold) if show_value_mid else None
                p.drawText(QPointF(text_x, text_y), text)
                p.setFont(self._value_font) if show_value_mid else None
        
        # Draw label
        p.setFont(self._label_font)
        p.setPen(self._label_color)
        text_x = self._left_pad_frac*self.slider_width/2
        text_y = -bar_bottom
        p.drawText(QPointF(text_x, text_y), self.label_str)
        
        p.restore()
        
        
    def mousePressEvent(self, e) :
        pos = e.pos()
        if self.single_handle:
            self.dragging = "mid"
        else:
            dist = {"min": abs(self.val_to_x(self.vmin) - pos.x()),
                    "mid": abs(self.val_to_x(self.vmid) - pos.x()),
                    "max": abs(self.val_to_x(self.vmax) - pos.x())}
            self.dragging = min(dist, key=dist.get)
        if self._freeze_position_while_dragging:
            self._detach_from_roi_keep_view_pos()
        # When dragging a slider that resizes the ROI, hide siblings for clarity.
        if self.dragging == "mid" and self._is_roi_resizing_slider():
            self._hide_sibling_sliders()

    def mouseMoveEvent(self, e) :
        if self.dragging is None :
            return
        val = self.x_to_val(e.pos().x())
        
        if self.dragging == "min":
            self.vmin = min(val, self.vmid - self.eps)
        elif self.dragging == "mid":
            self.vmid = max(min(val, self.vmax - self.eps), self.vmin + self.eps)
        elif self.dragging == "max":
            self.vmax = max(val, self.vmid + self.eps)
            
        self.valuesChanged.emit(self.vmin, self.vmid, self.vmax)
        self.update()

        # If vmid is bound to ROI geometry (sigma/R_sersic), apply it while dragging.
        if self.dragging == "mid":
            self._apply_to_roi()

        # Sync this parameter to other selected ROIs that have the same slider
        if self.PlotWidget is not None and self.param_name is not None and self.roi.selected :
            for roi in self.PlotWidget.roi_list:
                if roi is self.roi or not roi.selected :
                    continue
                other_slider = getattr(roi, "sliders", {}).get(self.param_name)
                if other_slider is None:
                    continue
                if self.dragging == "min" :
                    other_slider.vmin = min(self.vmin, other_slider.vmid - other_slider.eps)
                elif self.dragging == "mid" :
                    other_slider.vmid = max(min(self.vmid, other_slider.vmax - other_slider.eps), other_slider.vmin + other_slider.eps)
                elif self.dragging == "max" :
                    other_slider.vmax = max(self.vmax, other_slider.vmid + other_slider.eps)
                if self.dragging == "mid":
                    other_slider._apply_to_roi()
                other_slider.update()
        
    def mouseReleaseEvent(self, e):
        self.dragging = None
        if self._freeze_position_while_dragging:
            self._reattach_to_roi()
        else:
            self.adjust_all()
        self._show_sibling_sliders()











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
        w = self.slider_width_effective
        return self.width() * self._left_pad_frac + (v - self.min_range) / (self.max_range - self.min_range) * w

    def x_to_val(self, x):
        w = self.slider_width_effective
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




