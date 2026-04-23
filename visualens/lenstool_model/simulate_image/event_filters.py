import numpy as np
from PyQt5.QtCore import Qt, QObject, QEvent
import pyqtgraph as pg

from ...utils.utils_astro.utils_general import world_to_relative


def _get_view(iv) :
    """Return the ViewBox for both pg.ImageView and PlotWidget-based widgets."""
    if hasattr(iv, 'getViewBox'):
        return iv.getViewBox()
    return iv.getView()



def _simulate_optimize_filter_function(filter_object, obj, event):
    if event.type() == QEvent.KeyPress and event.key() in (Qt.Key_Return, Qt.Key_Enter) :
        if event.modifiers() == Qt.ShiftModifier :
            filter_object.imsim.optimize()
        else :
            filter_object.imsim.simulate()
        return True  # Stop propagation
    return False     # Let other events pass through  


class SourceFilter(QObject) :
    def __init__(self, imsim) :
        super().__init__()
        self.imsim = imsim
        
    def eventFilter(self, obj, event):
        return _simulate_optimize_filter_function(self, obj, event)


class ImageFilter(QObject) :
    def __init__(self, imsim) :
        super().__init__()
        self.imsim = imsim
        
        self._hover_text_items = {}
        for name, iv in getattr(self.imsim, 'image_plane_views', {}).items():
            txt = pg.TextItem(anchor=(0, 1), color=(0, 255, 255))
            txt.hide()
            iv.addItem(txt)
            self._hover_text_items[name] = txt
            _get_view(iv).scene().sigMouseMoved.connect(
                lambda pos, _name=name, _iv=iv: self._update_hover_text(_name, _iv, pos)
            )

            # Use pyqtgraph's mouse click signal (scene coordinates) for accurate mapping.
            _get_view(iv).scene().sigMouseClicked.connect(
                lambda ev, _name=name, _iv=iv: self._handle_image_double_click(_name, _iv, ev)
            )

    def _handle_image_double_click(self, panel_name, image_view, ev):
        """Map a double-click in an image panel to a source-plane point."""
        try:
            if not ev.double() or ev.button() != Qt.LeftButton:
                return
        except Exception:
            return

        vb = _get_view(image_view)
        scene_pos = ev.scenePos()
        data_pos = vb.mapSceneToView(scene_pos)
        x, y = data_pos.x(), data_pos.y()
        
        # Convert from local pixel to full-image pixel coordinates (bottom-origin convention).
        x_im_full = self.imsim._crop_x0 + x
        y_im_full = self.imsim._crop_y0 + y - self.imsim._crop_npix

        ra, dec = self.imsim.fits_image.image_to_world(x_im_full, y_im_full)
        xr, yr = world_to_relative(ra, dec, self.imsim.center_world)
        src_xr, src_yr = self.imsim.LensModel.ray_shooting(xr, yr, self.imsim.LensModel_kwargs)
        src_xr = src_xr - self.imsim.source_center_coordinates[0]
        src_yr = src_yr - self.imsim.source_center_coordinates[1]

        if not hasattr(self.imsim, 'interactive_source_plot') :
            self.imsim.interactive_source_plot = pg.ScatterPlotItem()
            self.imsim.source_plane_widget.addItem(self.imsim.interactive_source_plot)
        self.imsim.interactive_source_plot.setData([src_xr], [src_yr], symbol='x')

    def _get_panel_data(self, panel_name):
        if panel_name == 'observed':
            return self.imsim.image_data
        if panel_name == 'simulated':
            return self.imsim.simulated_image
        if panel_name == 'residual':
            return self.imsim.residual_image
        if panel_name == 'rgb':
            return self.imsim.rgb_image_data
        return None

    def _update_hover_text(self, panel_name, image_view, scene_pos):
        text_item = self._hover_text_items.get(panel_name)
        if text_item is None:
            return

        data = self._get_panel_data(panel_name)
        if data is None:
            text_item.hide()
            return

        vb = _get_view(image_view)
        data_pos = vb.mapSceneToView(scene_pos)
        x = int(np.floor(data_pos.x()))
        y_plot = int(np.floor(data_pos.y()))
        h, w = data.shape[:2]
        if x < 0 or x >= w or y_plot < 0 or y_plot >= h:
            text_item.hide()
            return

        # Images are displayed flipped in y (data[::-1, :]); convert plotted y to array row.
        y = h - 1 - y_plot
        val = data[y, x]
        if np.ndim(val) == 0:
            value_str = f"{float(val):.4g}"
        else:
            value_str = "(" + ", ".join(f"{float(v):.4g}" for v in np.asarray(val).ravel()) + ")"

        text_item.setText(f"{value_str}")
        text_item.setPos(data_pos.x() + 1, data_pos.y() + 1)
        text_item.show()

    def _hide_hover_text(self, panel_name=None):
        if panel_name is None:
            for txt in self._hover_text_items.values():
                txt.hide()
            return
        txt = self._hover_text_items.get(panel_name)
        if txt is not None:
            txt.hide()

    def eventFilter(self, obj, event) :
        _simulate_optimize_filter_function(self, obj, event)
        if event.type() == QEvent.KeyPress and event.key() == Qt.Key_Space :
            print("Spacebar function not yet implemented.")
            return True  # Stop propagation
        if event.type() == QEvent.Leave:
            for _name, _iv in getattr(self.imsim, 'image_plane_views', {}).items():
                if obj is _iv:
                    self._hide_hover_text(_name)
                    break
            else:
                self._hide_hover_text()
            return False
        # Double-click mapping is handled via sigMouseClicked (scene coords) in __init__.
        #return super().eventFilter(obj, event)
        return False     # Let other events pass through







