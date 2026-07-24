import pyqtgraph as pg
from PyQt5.QtCore import Qt
from astropy.table import Table


class ImageView_custom_selector(pg.ImageView) :
    def __init__(self, wcs, *args, **kwargs) :
        super().__init__(*args, **kwargs)
        self.wcs = wcs
        self.prefix = '1'
        self.catalog = Table(names=['id','ra','dec'], dtype=['str', 'float64', 'float64'])
        self._suffix = 0
        
        self._pending_x = []
        self._pending_y = []
        self._confirmed_x = []
        self._confirmed_y = []
        
        self._plus_scatter = pg.ScatterPlotItem(
            size=18, symbol='+', pen=pg.mkPen('w', width=2), brush=pg.mkBrush('w')
        )
        self._circle_scatter = pg.ScatterPlotItem(
            size=14, symbol='o', pen=pg.mkPen('w', width=2), brush=pg.mkBrush(0, 0, 0, 0)
        )
        self.addItem(self._plus_scatter)
        self.addItem(self._circle_scatter)
        
        self.setFocusPolicy(Qt.StrongFocus)
        self.scene.sigMouseClicked.connect(self._on_mouse_clicked)
    
    
    def _on_mouse_clicked(self, evt) :
        if not evt.double() :
            return
        pos = evt.scenePos()
        if not self.getView().sceneBoundingRect().contains(pos) :
            return
        
        mouse_point = self.getView().mapSceneToView(pos)
        x, y_display = mouse_point.x(), mouse_point.y()
        
        # Display y -> FITS/WCS pixel y (image is typically flipped for pyqtgraph)
        if self.image is not None :
            y_wcs = self.image.shape[0] - y_display
        else :
            y_wcs = y_display
        
        world = self.wcs.pixel_to_world(x, y_wcs)
        ra, dec = float(world.ra.deg), float(world.dec.deg)
        
        self._suffix += 1
        obj_id = self.prefix + '.' + str(self._suffix)
        self.catalog.add_row({'id': obj_id, 'ra': ra, 'dec': dec})
        
        self._pending_x.append(x)
        self._pending_y.append(y_display)
        self._plus_scatter.setData(self._pending_x, self._pending_y)
    
    
    def keyPressEvent(self, event) :
        if event.key() in (Qt.Key_Return, Qt.Key_Enter) :
            self._confirm_markers()
        else :
            super().keyPressEvent(event)
    
    
    def _confirm_markers(self) :
        self._confirmed_x.extend(self._pending_x)
        self._confirmed_y.extend(self._pending_y)
        self._pending_x.clear()
        self._pending_y.clear()
        self._plus_scatter.setData([], [])
        self._circle_scatter.setData(self._confirmed_x, self._confirmed_y)
        self.prefix = str(int(self.prefix) + 1)
        self._suffix = 0
    
    
    def clear_selection(self) :
        self.prefix = '1'
        self._suffix = 0
        self._pending_x.clear()
        self._pending_y.clear()
        self._confirmed_x.clear()
        self._confirmed_y.clear()
        self._plus_scatter.setData([], [])
        self._circle_scatter.setData([], [])
        if len(self.catalog) > 0 :
            self.catalog.remove_rows(slice(0, len(self.catalog)))
