import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
import numpy as np


class ParameterSliderWidget(QtWidgets.QWidget):
    """Widget containing min, value, max sliders for a single parameter"""
    valueChanged = QtCore.Signal(str, dict)  # param_name, {min, value, max}
    
    def __init__(self, param_name, initial_min=0, initial_val=1, initial_max=2):
        super().__init__()
        self.param_name = param_name
        
        # Make widget frameless and translucent
        self.setWindowFlags(QtCore.Qt.FramelessWindowHint | QtCore.Qt.WindowStaysOnTopHint)
        self.setAttribute(QtCore.Qt.WA_TranslucentBackground)
        
        layout = QtWidgets.QVBoxLayout()
        layout.setContentsMargins(5, 5, 5, 5)
        
        # Label
        label = QtWidgets.QLabel(param_name)
        label.setStyleSheet("background-color: rgba(50, 50, 50, 200); color: white; padding: 2px;")
        layout.addWidget(label)
        
        # Create slider container with background
        slider_container = QtWidgets.QWidget()
        slider_container.setStyleSheet("background-color: rgba(50, 50, 50, 200); border-radius: 5px;")
        slider_layout = QtWidgets.QVBoxLayout(slider_container)
        
        # Three sliders
        self.min_slider = self._create_slider("Min", initial_min)
        self.val_slider = self._create_slider("Val", initial_val)
        self.max_slider = self._create_slider("Max", initial_max)
        
        slider_layout.addLayout(self.min_slider['layout'])
        slider_layout.addLayout(self.val_slider['layout'])
        slider_layout.addLayout(self.max_slider['layout'])
        
        layout.addWidget(slider_container)
        self.setLayout(layout)
        
        # Connect signals
        for slider_dict in [self.min_slider, self.val_slider, self.max_slider]:
            slider_dict['slider'].valueChanged.connect(self._on_value_changed)
    
    def _create_slider(self, label, initial_value):
        """Create a labeled slider"""
        layout = QtWidgets.QHBoxLayout()
        
        lbl = QtWidgets.QLabel(f"{label}:")
        lbl.setStyleSheet("color: white;")
        lbl.setFixedWidth(30)
        
        slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        slider.setMinimum(0)
        slider.setMaximum(100)
        slider.setValue(int(initial_value * 10))
        slider.setStyleSheet("QSlider::groove:horizontal { background: #555; }")
        
        value_label = QtWidgets.QLabel(f"{initial_value:.1f}")
        value_label.setStyleSheet("color: white;")
        value_label.setFixedWidth(35)
        
        layout.addWidget(lbl)
        layout.addWidget(slider)
        layout.addWidget(value_label)
        
        return {'layout': layout, 'slider': slider, 'label': value_label}
    
    def _on_value_changed(self):
        """Emit current values when any slider changes"""
        values = {
            'min': self.min_slider['slider'].value() / 10.0,
            'value': self.val_slider['slider'].value() / 10.0,
            'max': self.max_slider['slider'].value() / 10.0
        }
        
        # Update labels
        self.min_slider['label'].setText(f"{values['min']:.1f}")
        self.val_slider['label'].setText(f"{values['value']:.1f}")
        self.max_slider['label'].setText(f"{values['max']:.1f}")
        
        self.valueChanged.emit(self.param_name, values)
    
    def get_values(self):
        """Get current parameter values"""
        return {
            'min': self.min_slider['slider'].value() / 10.0,
            'value': self.val_slider['slider'].value() / 10.0,
            'max': self.max_slider['slider'].value() / 10.0
        }


class SourceROI(pg.EllipseROI):
    """Extended EllipseROI with parameter controls"""
    
    def __init__(self, pos, size, parameters=None, **kwargs):
        super().__init__(pos, size, **kwargs)
        
        # Store parameters and their widgets
        self.parameters = parameters or {'n_Sersic': {'min': 0.5, 'value': 4.0, 'max': 8.0}}
        self.param_widgets = {}
        
        # Connect to ROI movement
        self.sigRegionChanged.connect(self._update_widget_positions)
        
        # Create parameter widgets
        self._create_parameter_widgets()
    
    def _create_parameter_widgets(self):
        """Create slider widgets for each parameter"""
        for i, (param_name, param_vals) in enumerate(self.parameters.items()):
            widget = ParameterSliderWidget(
                param_name,
                param_vals.get('min', 0),
                param_vals.get('value', 1),
                param_vals.get('max', 2)
            )
            
            widget.valueChanged.connect(lambda name, vals, pn=param_name: self._on_param_changed(pn, vals))
            self.param_widgets[param_name] = widget
            
            # Show widget
            widget.show()
        
        self._update_widget_positions()
    
    def _update_widget_positions(self):
        """Update widget positions to follow ROI"""
        if not self.param_widgets:
            return
        
        # Get ROI position in scene coordinates
        roi_pos = self.pos()
        roi_size = self.size()
        
        # Get the view widget to convert scene to screen coordinates
        if self.scene() and self.scene().views():
            view = self.scene().views()[0]
            
            # Position to the right of the ROI
            offset_x = roi_pos[0] + roi_size[0] + 10
            offset_y = roi_pos[1]
            
            # Convert to screen coordinates
            screen_pos = view.mapToGlobal(view.mapFromScene(QtCore.QPointF(offset_x, offset_y)))
            
            # Stack widgets vertically
            current_y = screen_pos.y()
            for widget in self.param_widgets.values():
                widget.move(screen_pos.x(), current_y)
                current_y += widget.height() + 5
    
    def _on_param_changed(self, param_name, values):
        """Handle parameter value changes"""
        self.parameters[param_name] = values
        print(f"Parameter {param_name} changed: {values}")
    
    def get_source_parameters(self):
        """Get all source parameters including position, size, and custom params"""
        return {
            'position': self.pos(),
            'size': self.size(),
            'angle': self.angle(),
            'parameters': {name: widget.get_values() 
                          for name, widget in self.param_widgets.items()}
        }
    
    def remove(self):
        """Clean up widgets when ROI is removed"""
        for widget in self.param_widgets.values():
            widget.close()
        super().removeROI()


# Example usage

# Create plot widget
win = pg.GraphicsLayoutWidget()
win.show()
win.setWindowTitle('Source ROI with Parameters')

view = win.addViewBox()
view.setAspectLocked()

# Add some background image
img = np.random.normal(size=(100, 100))
img_item = pg.ImageItem(img)
view.addItem(img_item)

# Create source ROI with parameters
roi = SourceROI(
    pos=[30, 30],
    size=[20, 15],
    parameters={
        'n_Sersic': {'min': 0.5, 'value': 4.0, 'max': 8.0},
        'ellipticity': {'min': 0.0, 'value': 0.3, 'max': 1.0}
    },
    pen='r'
)
view.addItem(roi)

# Button to print current parameters
def print_params():
    params = roi.get_source_parameters()
    print("\nCurrent Source Parameters:")
    print(f"Position: {params['position']}")
    print(f"Size: {params['size']}")
    print(f"Parameters: {params['parameters']}")

QtWidgets.QShortcut(QtGui.QKeySequence("Space"), win, print_params)

