import math

from PyQt5.QtCore import QObject, QPointF, QSize, Qt, pyqtSignal
from PyQt5.QtGui import QColor, QImage, QPainter, QPen
from PyQt5.QtWidgets import QHBoxLayout, QLabel, QSlider, QWidget


SLIDER_MAX = 255


class ColorDisk(QWidget) :
    """
    Hue/saturation disk: hue is the angle, saturation the distance from the
    centre (white at the centre, fully saturated colours on the rim).
    """
    changed = pyqtSignal(float, float)   # hue, saturation (both in [0, 1])
    released = pyqtSignal()

    DIAMETER = 140

    def __init__(self, hue=0., saturation=0., parent=None) :
        super().__init__(parent)
        self.hue = hue
        self.saturation = saturation
        self.setFixedSize(self.DIAMETER, self.DIAMETER)
        self.setCursor(Qt.CrossCursor)
        self._disk = self._render_disk(self.DIAMETER)

    def sizeHint(self) :
        return QSize(self.DIAMETER, self.DIAMETER)

    @staticmethod
    def _render_disk(diameter) :
        image = QImage(diameter, diameter, QImage.Format_ARGB32)
        image.fill(Qt.transparent)
        radius = diameter / 2.
        for y in range(diameter) :
            for x in range(diameter) :
                dx = (x + 0.5 - radius) / radius
                dy = (y + 0.5 - radius) / radius
                dist = math.hypot(dx, dy)
                alpha = min(1., max(0., (1. - dist) * radius))   # anti-aliased rim
                if alpha <= 0. :
                    continue
                hue = (math.atan2(-dy, dx) / (2 * math.pi)) % 1.
                color = QColor.fromHsvF(hue, min(dist, 1.), 1., alpha)
                image.setPixelColor(x, y, color)
        return image

    def set_hs(self, hue, saturation) :
        self.hue = hue
        self.saturation = saturation
        self.update()

    def paintEvent(self, event) :
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        if not self.isEnabled() :
            painter.setOpacity(0.35)
        painter.drawImage(0, 0, self._disk)

        radius = self.DIAMETER / 2.
        angle = self.hue * 2 * math.pi
        marker = QPointF(
            radius + self.saturation * radius * math.cos(angle),
            radius - self.saturation * radius * math.sin(angle),
        )
        painter.setBrush(Qt.NoBrush)
        painter.setPen(QPen(Qt.black, 3))
        painter.drawEllipse(marker, 5, 5)
        painter.setPen(QPen(Qt.white, 1.5))
        painter.drawEllipse(marker, 5, 5)

    def mousePressEvent(self, event) :
        if event.button() == Qt.LeftButton :
            self._pick(event.pos())

    def mouseMoveEvent(self, event) :
        if event.buttons() & Qt.LeftButton :
            self._pick(event.pos())

    def mouseReleaseEvent(self, event) :
        if event.button() == Qt.LeftButton :
            self.released.emit()

    def _pick(self, pos) :
        radius = self.DIAMETER / 2.
        dx = (pos.x() - radius) / radius
        dy = (pos.y() - radius) / radius
        saturation = min(math.hypot(dx, dy), 1.)
        # At the exact centre the angle is undefined: keep the current hue.
        hue = (math.atan2(-dy, dx) / (2 * math.pi)) % 1. if saturation > 1e-6 else self.hue
        self.set_hs(hue, saturation)
        self.changed.emit(hue, saturation)


class ColorControls(QObject) :
    """
    Colour disk + brightness and opacity sliders, exposed as separate widgets /
    layouts (``disk``, ``brightness_row``, ``opacity_row``) so that the caller
    decides how to arrange them. The caller must keep a reference to this object.
    """
    changed = pyqtSignal(float, float, float, float)   # hue, saturation, value, opacity (all in [0, 1])
    released = pyqtSignal()

    def __init__(self, hue=0., saturation=0., value=1., opacity=0., parent=None) :
        super().__init__(parent)
        self.value = value
        self.opacity = opacity

        self.disk = ColorDisk(hue, saturation)
        self.disk.changed.connect(self._on_disk_changed)
        self.disk.released.connect(self.released)

        self.brightness_slider = self._make_slider(value, self._on_brightness_changed)
        self.opacity_slider = self._make_slider(opacity, self._on_opacity_changed)

        self.swatch = QLabel()
        self.swatch.setFixedSize(28, 18)

        brightness_label = QLabel('Brightness')
        opacity_label = QLabel('Opacity')
        opacity_label.setFixedWidth(brightness_label.sizeHint().width())

        # Same width as the swatch so that both sliders have the same length.
        swatch_spacer = QWidget()
        swatch_spacer.setFixedSize(self.swatch.size())

        self.brightness_row = QHBoxLayout()
        self.brightness_row.addWidget(brightness_label)
        self.brightness_row.addWidget(self.brightness_slider, 1)
        self.brightness_row.addWidget(self.swatch)

        self.opacity_row = QHBoxLayout()
        self.opacity_row.addWidget(opacity_label)
        self.opacity_row.addWidget(self.opacity_slider, 1)
        self.opacity_row.addWidget(swatch_spacer)

        self._update_swatch()

    def _make_slider(self, initial, on_value_changed) :
        slider = QSlider(Qt.Horizontal)
        slider.setRange(0, SLIDER_MAX)
        slider.setSingleStep(1)
        slider.setValue(round(initial * SLIDER_MAX))
        slider.valueChanged.connect(lambda raw : on_value_changed(raw / SLIDER_MAX, slider))
        slider.sliderReleased.connect(self.released)
        return slider

    def hsv(self) :
        return self.disk.hue, self.disk.saturation, self.value

    def set_enabled(self, enabled) :
        for widget in (self.disk, self.brightness_slider, self.opacity_slider, self.swatch) :
            widget.setEnabled(enabled)

    def _emit_changed(self) :
        self._update_swatch()
        self.changed.emit(self.disk.hue, self.disk.saturation, self.value, self.opacity)

    def _on_disk_changed(self, hue, saturation) :
        self._emit_changed()

    def _on_brightness_changed(self, value, slider) :
        self.value = value
        self._emit_changed()
        self._release_if_not_dragging(slider)

    def _on_opacity_changed(self, opacity, slider) :
        self.opacity = opacity
        self._emit_changed()
        self._release_if_not_dragging(slider)

    def _release_if_not_dragging(self, slider) :
        # Keyboard / wheel / track clicks never fire sliderReleased.
        if not slider.isSliderDown() :
            self.released.emit()

    def _update_swatch(self) :
        hue, saturation, value = self.hsv()
        color = QColor.fromHsvF(hue, saturation, value)
        self.swatch.setStyleSheet(
            f"background-color: {color.name()}; border: 1px solid palette(mid); border-radius: 3px;"
        )
