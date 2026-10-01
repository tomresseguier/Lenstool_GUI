from contextlib import contextmanager

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QFrame, QGridLayout, QLabel, QListWidget, QMessageBox,
    QPushButton, QScrollArea, QSizePolicy, QVBoxLayout, QWidget,
)


class ControlBuilder :
    """Factory helpers for the small widgets that the side panel builds over and over."""

    # Both states must share the same box model: if only ``:checked`` is styled,
    # Qt draws checked buttons via the stylesheet and unchecked ones natively,
    # and the two end up with slightly different sizes/offsets in a row.
    TOGGLE_STYLE = (
        "QPushButton {"
        "  padding: 4px 10px;"
        "  border: 1px solid palette(mid);"
        "  border-radius: 4px;"
        "  background-color: palette(button);"
        "  color: palette(button-text);"
        "}"
        "QPushButton:pressed {"
        "  background-color: palette(midlight);"
        "}"
        "QPushButton:checked {"
        "  background-color: #3d7a4a;"
        "  color: white;"
        "}"
    )

    @staticmethod
    @contextmanager
    def signals_blocked(*widgets) :
        previous = [w.blockSignals(True) for w in widgets]
        try :
            yield
        finally :
            for w, was_blocked in zip(widgets, previous) :
                w.blockSignals(was_blocked)

    @staticmethod
    def set_checked_silently(btn, checked) :
        with ControlBuilder.signals_blocked(btn) :
            btn.setChecked(checked)

    @staticmethod
    def scroll_area(content_widget) :
        """Wrap ``content_widget`` in a vertical scroll area for constrained layouts."""
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll.setWidget(content_widget)
        return scroll

    @staticmethod
    def separator() :
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        return line

    @staticmethod
    def toggle_button(label, callback=None, checked=False) :
        btn = QPushButton(label)
        btn.setCheckable(True)
        btn.setStyleSheet(ControlBuilder.TOGGLE_STYLE)
        ControlBuilder.set_checked_silently(btn, checked)
        if callback is not None :
            btn.toggled.connect(callback)
        return btn

    @staticmethod
    def list_widget(items, height, current=None) :
        list_widget = QListWidget()
        list_widget.setFixedHeight(height)
        list_widget.addItems([str(item) for item in items])
        if current is not None :
            matches = list_widget.findItems(current, Qt.MatchExactly)
            if matches :
                with ControlBuilder.signals_blocked(list_widget) :
                    list_widget.setCurrentItem(matches[0])
        return list_widget

    @staticmethod
    def labeled_list(label_text, items, height, current=None) :
        """Return ``(layout, list_widget)``: a centred caption above a fixed-height list."""
        v_layout = QVBoxLayout()
        v_layout.setContentsMargins(0, 0, 0, 0)
        v_layout.setSpacing(2)

        label = QLabel(label_text)
        label.setAlignment(Qt.AlignCenter)
        # The list below has a fixed height, so the label is the only "flexible"
        # item in this mini layout; without pinning it, it soaks up any extra
        # vertical space handed to the row, pushing it away from the list.
        label.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        v_layout.addWidget(label)

        list_widget = ControlBuilder.list_widget(items, height, current)
        v_layout.addWidget(list_widget)
        return v_layout, list_widget

    @staticmethod
    def reset_grid(container) :
        """Delete everything in ``container``'s layout and install a fresh, empty grid layout."""
        old_layout = container.layout()
        if old_layout is not None :
            while old_layout.count() > 0 :
                widget = old_layout.takeAt(0).widget()
                if widget is not None :
                    widget.deleteLater()
            # Reparent away from the container before deleting, otherwise Qt
            # keeps the (now childless) old layout installed on the widget.
            QWidget().setLayout(old_layout)

        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        container.setLayout(grid)
        return grid

    @staticmethod
    def clear_tabs(tab_widget) :
        while tab_widget.count() > 0 :
            widget = tab_widget.widget(0)
            tab_widget.removeTab(0)
            widget.deleteLater()

    @staticmethod
    def warn(parent, title, message) :
        QMessageBox.warning(parent, title, message)

    @staticmethod
    def error(parent, title, exc) :
        QMessageBox.critical(parent, title, str(exc))
