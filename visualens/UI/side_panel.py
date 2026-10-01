import os
from functools import partial

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QTabWidget, QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QLabel, QFileDialog, QCheckBox, QSlider,
)

from .catalog_controls import CatalogControlsTab
from .controls import ControlBuilder
from .ui_state import CatalogUIState, LensModelUIState


LT_Z_STEP = 0.5   # redshift slider step; 30 steps -> redshift range [0, 15]
LT_Z_N_STEPS = 30
BROAD_FAMILY_COLUMNS = 5


class ImageSidePanel(QTabWidget) :
    """
    Collapsible left-side panel shown next to the ImageView, with one tab
    each for Image, Catalog and Lens model controls.

    The panel subscribes to ``workspace.signals`` to keep itself in sync with
    the workspace; nothing needs to call its ``refresh_*`` methods by hand.
    """
    def __init__(self, workspace, *args, **kwargs) :
        super().__init__(*args, **kwargs)
        self.workspace = workspace

        self._critical_curve_toggle = None
        self._caustic_curve_toggle = None

        self.image_tab = self._make_image_tab()
        self.catalog_tab, self.catalog_subtabs = self._make_catalog_tab()
        self.lens_model_tab = self._make_lens_model_tab()

        self.addTab(ControlBuilder.scroll_area(self.image_tab), 'Image')
        self.addTab(ControlBuilder.scroll_area(self.catalog_tab), 'Catalog')
        self.addTab(ControlBuilder.scroll_area(self.lens_model_tab), 'Lens model')

        self._connect_workspace_signals()
        self.refresh_all()

    def _connect_workspace_signals(self) :
        signals = self.workspace.signals
        signals.image_changed.connect(self.sync_link_button_states)
        signals.catalogs_changed.connect(self.refresh_catalog_tabs)
        signals.lens_model_changed.connect(self.refresh_lens_model_tab)
        signals.window_attached.connect(self.refresh_all)

    def refresh_all(self) :
        self.sync_link_button_states()
        self.refresh_catalog_tabs()
        self.refresh_lens_model_tab()

    def _dialog_parent(self) :
        return self.workspace._dialog_parent()

    # ------------------------------------------------------------------
    # Image tab
    # ------------------------------------------------------------------
    def _make_image_tab(self) :
        tab = QWidget()
        layout = QVBoxLayout(tab)

        top_row = QHBoxLayout()
        for label, callback, stretch in (
            ('Open…', self.workspace.open_image_dialog, 1),
            ('Show viewer', self.workspace.reopen_main_window, 2),
            ('2nd window', self.workspace.plot_secondary_image, 2),
            ('Split', self.workspace.toggle_split4, 2),
        ) :
            btn = QPushButton(label)
            btn.clicked.connect(callback)
            top_row.addWidget(btn, stretch)
        layout.addLayout(top_row)

        link_row = QHBoxLayout()
        self.zoom_link_btn = ControlBuilder.toggle_button('Zoom link', self._toggle_zoom_link)
        self.pan_link_btn = ControlBuilder.toggle_button('Pan link', self._toggle_pan_link)
        link_row.addWidget(self.zoom_link_btn)
        link_row.addWidget(self.pan_link_btn)
        layout.addLayout(link_row)

        layout.addWidget(ControlBuilder.separator())

        layout.addWidget(QLabel('Hand selected multiple images'))

        hand_select_row = QHBoxLayout()
        export_btn = QPushButton('export selection…')
        export_btn.clicked.connect(self._export_hand_selection)
        hand_select_row.addWidget(export_btn)

        clear_btn = QPushButton('clear_selection')
        clear_btn.clicked.connect(self._clear_hand_selection)
        hand_select_row.addWidget(clear_btn)
        layout.addLayout(hand_select_row)

        layout.addStretch()
        return tab

    def _export_hand_selection(self) :
        img = self.workspace.image
        if img is None :
            return
        base_dir = os.path.dirname(img.image_path) if img.image_path is not None else os.getcwd()
        path_str, _ = QFileDialog.getSaveFileName(
            self._dialog_parent(),
            'Export hand selection',
            os.path.join(base_dir, 'hand_selected_catalog.lenstool'),
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
        sync = img._view_range_sync_state() if img is not None else {'zoom' : False, 'pan' : False}
        with ControlBuilder.signals_blocked(self.zoom_link_btn, self.pan_link_btn) :
            self.zoom_link_btn.setChecked(sync['zoom'])
            self.pan_link_btn.setChecked(sync['pan'])

    # ------------------------------------------------------------------
    # Lens model tab
    # ------------------------------------------------------------------
    def _make_lens_model_tab(self) :
        tab = QWidget()
        layout = QVBoxLayout(tab)

        open_btn = QPushButton('Open…')
        open_btn.clicked.connect(self.workspace.open_lens_model_dialog)
        layout.addWidget(open_btn)

        header_row = QHBoxLayout()
        header_row.addWidget(QLabel('Set multiple image system(s):'))
        all_btn = QPushButton('All')
        all_btn.clicked.connect(self._select_all_broad_families)
        header_row.addWidget(all_btn)
        none_btn = QPushButton('None')
        none_btn.clicked.connect(self._select_no_broad_families)
        header_row.addWidget(none_btn)
        header_row.addStretch()
        layout.addLayout(header_row)

        self._lens_model_families_container = QWidget()
        layout.addWidget(self._lens_model_families_container)

        layout.addWidget(ControlBuilder.separator())

        z_row = QHBoxLayout()
        z_row.addWidget(QLabel('Source plane redshift'))
        self._lt_z_value_label = QLabel('0.0')
        z_row.addWidget(self._lt_z_value_label)
        z_row.addStretch()
        layout.addLayout(z_row)

        self.lt_z_slider = QSlider(Qt.Horizontal)
        self.lt_z_slider.setMinimum(0)
        self.lt_z_slider.setMaximum(LT_Z_N_STEPS)
        self.lt_z_slider.setSingleStep(1)
        self.lt_z_slider.setValue(0)
        self.lt_z_slider.valueChanged.connect(self._on_lt_z_slider_value_changed)
        self.lt_z_slider.sliderReleased.connect(self._on_lt_z_slider_released)
        layout.addWidget(self.lt_z_slider)

        layout.addWidget(ControlBuilder.separator())

        self.lens_model_subtabs = QTabWidget()
        layout.addWidget(self.lens_model_subtabs, 1)

        layout.addWidget(ControlBuilder.separator())

        on_the_fly_btn = QPushButton('Start on the fly predicted images')
        on_the_fly_btn.clicked.connect(self._start_on_the_fly_predicted_images)
        layout.addWidget(on_the_fly_btn)

        layout.addWidget(ControlBuilder.separator())

        open_imsim_btn = QPushButton('Open image simulator')
        open_imsim_btn.clicked.connect(self._open_image_simulator)
        layout.addWidget(open_imsim_btn)

        return tab

    def refresh_lens_model_tab(self) :
        """Rebuild the dynamic parts of the Lens model tab to match ``workspace.lens_model``."""
        self._rebuild_broad_family_checkboxes()
        self._rebuild_lens_model_subtabs()

        lens_model = self.workspace.lens_model
        if lens_model is not None and getattr(lens_model, 'lt_z', None) is not None :
            with ControlBuilder.signals_blocked(self.lt_z_slider) :
                self.lt_z_slider.setValue(round(lens_model.lt_z / LT_Z_STEP))
            self._lt_z_value_label.setText(f"{lens_model.lt_z:.1f}")

    def _rebuild_broad_family_checkboxes(self) :
        grid = ControlBuilder.reset_grid(self._lens_model_families_container)

        lens_model = self.workspace.lens_model
        if lens_model is None :
            return

        broad_families = lens_model.broad_families
        if not isinstance(broad_families, list) :
            broad_families = broad_families.tolist()

        state = LensModelUIState.of(lens_model)
        state.sync_with(broad_families)

        for i, name in enumerate(broad_families) :
            checkbox = QCheckBox(str(name))
            checkbox.setChecked(name in state.checked_broad_families)
            checkbox.toggled.connect(partial(self._on_broad_family_toggled, name))
            row, col = divmod(i, BROAD_FAMILY_COLUMNS)
            grid.addWidget(checkbox, row, col)

    def _on_broad_family_toggled(self, name, checked) :
        lens_model = self.workspace.lens_model
        if lens_model is None :
            return
        state = LensModelUIState.of(lens_model)
        state.set_checked(name, checked)

        lens_model.set_which(state.checked_in_order(lens_model.broad_families))
        self._replot_lens_model_which(lens_model)

    def _select_all_broad_families(self) :
        lens_model = self.workspace.lens_model
        if lens_model is None :
            return
        broad_families = lens_model.broad_families
        if not isinstance(broad_families, list) :
            broad_families = broad_families.tolist()
        state = LensModelUIState.of(lens_model)
        state.checked_broad_families = set(broad_families)
        lens_model.set_which(broad_families)
        self._replot_lens_model_which(lens_model)
        self._rebuild_broad_family_checkboxes()

    def _select_no_broad_families(self) :
        lens_model = self.workspace.lens_model
        if lens_model is None :
            return
        state = LensModelUIState.of(lens_model)
        state.checked_broad_families = set()
        lens_model.set_which([])
        self._replot_lens_model_which(lens_model)
        self._rebuild_broad_family_checkboxes()

    def _replot_lens_model_which(self, lens_model) :
        for cat in (lens_model.mult, lens_model.source, lens_model.images_filtered) :
            if cat is None :
                continue
            cat_state = CatalogUIState.existing(cat)
            if cat_state is not None :
                cat_state.replot(cat)

    def _on_lt_z_slider_value_changed(self, value) :
        self._lt_z_value_label.setText(f"{value * LT_Z_STEP:.1f}")

    def _on_lt_z_slider_released(self) :
        lens_model = self.workspace.lens_model
        if lens_model is None :
            return
        z = self.lt_z_slider.value() * LT_Z_STEP
        try :
            lens_model.set_lt_z(
                z,
                plot_critical=self._is_checked(self._critical_curve_toggle),
                plot_caustic=self._is_checked(self._caustic_curve_toggle),
            )
        except Exception as exc :
            ControlBuilder.error(self._dialog_parent(), 'Failed to set redshift', exc)

    @staticmethod
    def _is_checked(btn) :
        return btn is not None and btn.isChecked()

    def _rebuild_lens_model_subtabs(self) :
        subtabs = self.lens_model_subtabs
        current_index = subtabs.currentIndex()
        ControlBuilder.clear_tabs(subtabs)

        lens_model = self.workspace.lens_model

        for title, attr in (
            ('Multiple images', 'mult'),
            ('Sources', 'source'),
            ('Predicted images', 'images_filtered'),
        ) :
            cat = getattr(lens_model, attr, None) if lens_model is not None else None
            if cat is not None and attr == 'mult' :
                CatalogUIState.of(cat, plot_column=True, colname='id')
            subtabs.addTab(self._make_catalog_controls_tab(cat) if cat is not None else QWidget(), title)

        subtabs.addTab(self._make_curves_tab(lens_model), 'Curves')

        if 0 <= current_index < subtabs.count() :
            subtabs.setCurrentIndex(current_index)

    def _make_curves_tab(self, lens_model) :
        tab = QWidget()
        layout = QVBoxLayout(tab)

        critical_on = lens_model is not None and getattr(lens_model, 'critical_curve_plot', None) is not None
        caustic_on = lens_model is not None and getattr(lens_model, 'caustic_curve_plot', None) is not None

        self._critical_curve_toggle = ControlBuilder.toggle_button(
            'Plot critical curve', partial(self._toggle_lt_curve, which='critical'), critical_on,
        )
        self._caustic_curve_toggle = ControlBuilder.toggle_button(
            'Plot caustic curve', partial(self._toggle_lt_curve, which='caustic'), caustic_on,
        )
        layout.addWidget(self._critical_curve_toggle)
        layout.addWidget(self._caustic_curve_toggle)
        layout.addStretch()
        return tab

    def _toggle_lt_curve(self, checked, which) :
        lens_model = self.workspace.lens_model
        if lens_model is None :
            return
        if checked :
            if getattr(lens_model, 'lt_z', None) is None :
                ControlBuilder.warn(
                    self._dialog_parent(), 'No redshift set',
                    'Please set the source plane redshift with the slider first.',
                )
                toggle = self._critical_curve_toggle if which == 'critical' else self._caustic_curve_toggle
                ControlBuilder.set_checked_silently(toggle, False)
                return
            try :
                lens_model.plot_lt_curve(which=which)
            except Exception as exc :
                ControlBuilder.error(self._dialog_parent(), 'Failed to plot curve', exc)
        else :
            lens_model.clear_lt_curve(which=which)

    def _start_on_the_fly_predicted_images(self) :
        lens_model = self._require_lens_model(need_redshift=True)
        if lens_model is None :
            return
        try :
            lens_model.start_im2source()
        except Exception as exc :
            ControlBuilder.error(self._dialog_parent(), 'Failed to start on the fly predicted images', exc)

    def _open_image_simulator(self) :
        lens_model = self._require_lens_model()
        if lens_model is None :
            return
        try :
            lens_model.start_simulate_image()
        except Exception as exc :
            ControlBuilder.error(self._dialog_parent(), 'Failed to open image simulator', exc)

    def _require_lens_model(self, need_redshift=False) :
        """Return the workspace lens model, or warn the user and return None if it is not usable."""
        lens_model = self.workspace.lens_model
        if lens_model is None :
            ControlBuilder.warn(self._dialog_parent(), 'No model', 'Please import a Lenstool model first.')
            return None
        if need_redshift and getattr(lens_model, 'lt_z', None) is None :
            ControlBuilder.warn(
                self._dialog_parent(), 'No redshift set',
                'Please set the source plane redshift with the slider first.',
            )
            return None
        return lens_model

    # ------------------------------------------------------------------
    # Catalog tab: one sub-tab per catalog in workspace.catalogs
    # ------------------------------------------------------------------
    def _make_catalog_tab(self) :
        tab = QWidget()
        layout = QVBoxLayout(tab)

        open_btn = QPushButton('Open…')
        open_btn.clicked.connect(self.workspace.open_catalog_dialog)
        layout.addWidget(open_btn)

        subtabs = QTabWidget()
        layout.addWidget(subtabs, 1)
        return tab, subtabs

    def refresh_catalog_tabs(self) :
        """Rebuild the per-catalog sub-tabs to match ``workspace.catalogs``."""
        current_index = self.catalog_subtabs.currentIndex()
        ControlBuilder.clear_tabs(self.catalog_subtabs)

        for i, cat in enumerate(self.workspace.catalogs) :
            self.catalog_subtabs.addTab(self._make_catalog_controls_tab(cat), f'Cat {i + 1}')

        if 0 <= current_index < self.catalog_subtabs.count() :
            self.catalog_subtabs.setCurrentIndex(current_index)

    def _make_catalog_controls_tab(self, cat) :
        return CatalogControlsTab(self.workspace, cat)
