import os

from PyQt5.QtWidgets import (
    QTabWidget, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QFrame, QLabel, QFileDialog,
)




class ImageSidePanel(QTabWidget) :
    """
    Collapsible left-side panel shown next to the ImageView, with one tab
    each for Image, Catalog and Lens model controls.
    """
    def __init__(self, workspace, *args, **kwargs) :
        super().__init__(*args, **kwargs)
        self.workspace = workspace
        
        self.image_tab = self._make_image_tab()
        self.catalog_tab = self._make_tab_with_open_button(workspace.open_catalog_dialog)
        self.lens_model_tab = self._make_tab_with_open_button(workspace.open_lens_model_dialog)
        
        self.addTab(self.image_tab, 'Image')
        self.addTab(self.catalog_tab, 'Catalog')
        self.addTab(self.lens_model_tab, 'Lens model')

    def _make_image_tab(self) :
        tab = QWidget()
        layout = QVBoxLayout(tab)

        top_row = QHBoxLayout()
        open_btn = QPushButton('Open…')
        open_btn.clicked.connect(self.workspace.open_image_dialog)
        top_row.addWidget(open_btn, 1)

        second_window_btn = QPushButton('2nd window')
        second_window_btn.clicked.connect(self.workspace.plot_secondary_image)
        top_row.addWidget(second_window_btn, 2)

        split_btn = QPushButton('Split')
        split_btn.clicked.connect(self.workspace.toggle_split4)
        top_row.addWidget(split_btn, 2)

        layout.addLayout(top_row)

        link_row = QHBoxLayout()
        self.zoom_link_btn = self._make_link_toggle('Zoom link', self._toggle_zoom_link)
        self.pan_link_btn = self._make_link_toggle('Pan link', self._toggle_pan_link)
        link_row.addWidget(self.zoom_link_btn)
        link_row.addWidget(self.pan_link_btn)
        layout.addLayout(link_row)

        layout.addWidget(self._make_separator())

        hand_select_label = QLabel('Hand selected multiple images')
        layout.addWidget(hand_select_label)

        hand_select_row = QHBoxLayout()
        export_btn = QPushButton('export selection…')
        export_btn.clicked.connect(self._export_hand_selection)
        hand_select_row.addWidget(export_btn)

        clear_btn = QPushButton('clear_selection')
        clear_btn.clicked.connect(self._clear_hand_selection)
        hand_select_row.addWidget(clear_btn)
        layout.addLayout(hand_select_row)

        layout.addStretch()
        self.sync_link_button_states()
        return tab

    def _make_separator(self) :
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        return line

    def _export_hand_selection(self) :
        img = self.workspace.image
        if img is None :
            return
        base_dir = os.path.dirname(img.image_path) if img.image_path is not None else os.getcwd()
        default_path = os.path.join(base_dir, 'hand_selected_catalog.lenstool')
        path_str, _ = QFileDialog.getSaveFileName(
            self.workspace._dialog_parent(),
            'Export hand selection',
            default_path,
            'Lenstool mult files (*.lenstool);;All files (*)',
        )
        if not path_str :
            return
        img.export_hand_selection_as_mult_file(path_str)

    def _clear_hand_selection(self) :
        img = self.workspace.image
        if img is None :
            return
        img.clear_hand_selection()

    def _make_link_toggle(self, label, callback) :
        btn = QPushButton(label)
        btn.setCheckable(True)
        btn.toggled.connect(callback)
        btn.setStyleSheet(
            "QPushButton:checked {"
            "  background-color: #3d7a4a;"
            "  color: white;"
            "}"
        )
        return btn

    def _toggle_zoom_link(self, linked) :
        img = self.workspace.image
        if img is None :
            return
        if linked :
            img.link_zoom()
        else :
            img.unlink_zoom()

    def _toggle_pan_link(self, linked) :
        img = self.workspace.image
        if img is None :
            return
        if linked :
            img.link_pan()
        else :
            img.unlink_pan()

    def sync_link_button_states(self) :
        img = self.workspace.image
        if img is None :
            self.zoom_link_btn.setChecked(False)
            self.pan_link_btn.setChecked(False)
            return
        sync = img._view_range_sync_state()
        self.zoom_link_btn.blockSignals(True)
        self.pan_link_btn.blockSignals(True)
        self.zoom_link_btn.setChecked(sync['zoom'])
        self.pan_link_btn.setChecked(sync['pan'])
        self.zoom_link_btn.blockSignals(False)
        self.pan_link_btn.blockSignals(False)

    def _make_tab_with_open_button(self, callback) :
        tab = QWidget()
        layout = QVBoxLayout(tab)
        open_btn = QPushButton('Open…')
        open_btn.clicked.connect(callback)
        layout.addWidget(open_btn)
        layout.addStretch()
        return tab
