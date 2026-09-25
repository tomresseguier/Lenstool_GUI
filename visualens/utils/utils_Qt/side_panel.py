import os

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QTabWidget, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QPushButton,
    QFrame, QLabel, QFileDialog, QListWidget, QMessageBox, QInputDialog,
    QSizePolicy, QCheckBox, QSlider,
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
        self.catalog_tab, self.catalog_subtabs = self._make_catalog_tab()
        self.lens_model_tab = self._make_lens_model_tab()
        
        self.addTab(self.image_tab, 'Image')
        self.addTab(self.catalog_tab, 'Catalog')
        self.addTab(self.lens_model_tab, 'Lens model')
        
        self.refresh_catalog_tabs()
        self.refresh_lens_model_tab()

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

    # ------------------------------------------------------------------
    # Lens model tab
    # ------------------------------------------------------------------
    def _make_lens_model_tab(self) :
        tab = QWidget()
        layout = QVBoxLayout(tab)

        open_btn = QPushButton('Open…')
        open_btn.clicked.connect(self.workspace.open_lens_model_dialog)
        layout.addWidget(open_btn)

        layout.addWidget(QLabel('Set multiple image system:'))
        self._lens_model_families_container = QWidget()
        layout.addWidget(self._lens_model_families_container)

        layout.addWidget(self._make_separator())

        z_row = QHBoxLayout()
        z_row.addWidget(QLabel('Source plane redshift'))
        self._lt_z_value_label = QLabel('0.0')
        z_row.addWidget(self._lt_z_value_label)
        z_row.addStretch()
        layout.addLayout(z_row)

        self.lt_z_slider = QSlider(Qt.Horizontal)
        self.lt_z_slider.setMinimum(0)
        self.lt_z_slider.setMaximum(30)  # 30 steps of 0.5 -> redshift range [0, 15]
        self.lt_z_slider.setSingleStep(1)
        self.lt_z_slider.setValue(0)
        self.lt_z_slider.valueChanged.connect(self._on_lt_z_slider_value_changed)
        self.lt_z_slider.sliderReleased.connect(self._on_lt_z_slider_released)
        layout.addWidget(self.lt_z_slider)

        self.lens_model_subtabs = QTabWidget()
        layout.addWidget(self.lens_model_subtabs, 1)

        layout.addWidget(self._make_separator())

        on_the_fly_btn = QPushButton('Start on the fly predicted images')
        on_the_fly_btn.clicked.connect(self._start_on_the_fly_predicted_images)
        layout.addWidget(on_the_fly_btn)

        layout.addWidget(self._make_separator())

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
            self.lt_z_slider.blockSignals(True)
            self.lt_z_slider.setValue(round(lens_model.lt_z / 0.5))
            self.lt_z_slider.blockSignals(False)
            self._lt_z_value_label.setText(f"{lens_model.lt_z:.1f}")

    def _rebuild_broad_family_checkboxes(self) :
        container = self._lens_model_families_container
        old_layout = container.layout()
        if old_layout is not None :
            while old_layout.count() > 0 :
                item = old_layout.takeAt(0)
                widget = item.widget()
                if widget is not None :
                    widget.deleteLater()
            # Reparent away from the container before deleting, otherwise Qt
            # keeps the (now childless) old layout installed on the widget.
            QWidget().setLayout(old_layout)

        grid = QGridLayout()
        grid.setContentsMargins(0, 0, 0, 0)
        container.setLayout(grid)

        lens_model = self.workspace.lens_model
        if lens_model is None :
            return

        broad_families = lens_model.broad_families
        if not isinstance(broad_families, list) :
            broad_families = broad_families.tolist()

        checked = getattr(lens_model, '_checked_broad_families', None)
        if checked is None :
            checked = set(broad_families)
            lens_model._checked_broad_families = checked
        seen = getattr(lens_model, '_seen_broad_families', set())
        for name in broad_families :
            if name not in seen :
                # Newly discovered broad family: checked by default.
                checked.add(name)
        lens_model._seen_broad_families = set(broad_families)

        n_columns = 5
        for i, name in enumerate(broad_families) :
            checkbox = QCheckBox(str(name))
            self._set_checked_silently(checkbox, name in checked)
            checkbox.stateChanged.connect(
                lambda state, name=name : self._on_broad_family_toggled(lens_model, name, state)
            )
            row, col = divmod(i, n_columns)
            grid.addWidget(checkbox, row, col)

    def _on_broad_family_toggled(self, lens_model, name, state) :
        checked = lens_model._checked_broad_families
        if state == Qt.Checked :
            checked.add(name)
        else :
            checked.discard(name)

        lens_model.set_which(list(checked))
        self._replot_lens_model_which(lens_model)

    def _replot_lens_model_which(self, lens_model) :
        for cat in (lens_model.mult, lens_model.images_filtered) :
            if cat is None :
                continue
            state = getattr(cat, '_plot_ui_state', None)
            if state is not None and state.get('plot') :
                cat.plot(marker='o', filled_markers=False)
                if state.get('plot_column') and state.get('colname') is not None :
                    cat.plot_column(state['colname'])

    def _on_lt_z_slider_value_changed(self, value) :
        self._lt_z_value_label.setText(f"{value * 0.5:.1f}")

    def _on_lt_z_slider_released(self) :
        lens_model = self.workspace.lens_model
        if lens_model is None :
            return
        z = self.lt_z_slider.value() * 0.5
        plot_critical = getattr(self, '_critical_curve_toggle', None) is not None and self._critical_curve_toggle.isChecked()
        plot_caustic = getattr(self, '_caustic_curve_toggle', None) is not None and self._caustic_curve_toggle.isChecked()
        try :
            lens_model.set_lt_z(z, plot_critical=plot_critical, plot_caustic=plot_caustic)
        except Exception as exc :
            QMessageBox.critical(self.workspace._dialog_parent(), 'Failed to set redshift', str(exc))

    def _rebuild_lens_model_subtabs(self) :
        subtabs = self.lens_model_subtabs
        current_index = subtabs.currentIndex()
        while subtabs.count() > 0 :
            widget = subtabs.widget(0)
            subtabs.removeTab(0)
            widget.deleteLater()

        lens_model = self.workspace.lens_model

        if lens_model is not None and lens_model.mult is not None :
            mult_tab = self._make_catalog_controls_tab(lens_model.mult)
        else :
            mult_tab = QWidget()
        subtabs.addTab(mult_tab, 'Multiple images')

        if lens_model is not None and lens_model.source is not None :
            source_tab = self._make_catalog_controls_tab(lens_model.source)
        else :
            source_tab = QWidget()
        subtabs.addTab(source_tab, 'Sources')

        if lens_model is not None and lens_model.images_filtered is not None :
            predicted_tab = self._make_catalog_controls_tab(lens_model.images_filtered)
        else :
            predicted_tab = QWidget()
        subtabs.addTab(predicted_tab, 'Predicted images')

        subtabs.addTab(self._make_curves_tab(lens_model), 'Curves')

        if 0 <= current_index < subtabs.count() :
            subtabs.setCurrentIndex(current_index)

    def _make_curves_tab(self, lens_model) :
        tab = QWidget()
        layout = QVBoxLayout(tab)

        self._critical_curve_toggle = self._make_link_toggle(
            'Plot critical curve',
            lambda checked : self._toggle_lt_curve(checked, 'critical'),
        )
        self._caustic_curve_toggle = self._make_link_toggle(
            'Plot caustic curve',
            lambda checked : self._toggle_lt_curve(checked, 'caustic'),
        )
        layout.addWidget(self._critical_curve_toggle)
        layout.addWidget(self._caustic_curve_toggle)
        layout.addStretch()

        critical_on = lens_model is not None and getattr(lens_model, 'critical_curve_plot', None) is not None
        caustic_on = lens_model is not None and getattr(lens_model, 'caustic_curve_plot', None) is not None
        self._set_checked_silently(self._critical_curve_toggle, critical_on)
        self._set_checked_silently(self._caustic_curve_toggle, caustic_on)

        return tab

    def _toggle_lt_curve(self, checked, which) :
        lens_model = self.workspace.lens_model
        if lens_model is None :
            return
        if checked :
            if getattr(lens_model, 'lt_z', None) is None :
                QMessageBox.warning(
                    self.workspace._dialog_parent(),
                    'No redshift set',
                    'Please set the source plane redshift with the slider first.',
                )
                toggle = self._critical_curve_toggle if which=='critical' else self._caustic_curve_toggle
                self._set_checked_silently(toggle, False)
                return
            try :
                lens_model.plot_lt_curve(which=which)
            except Exception as exc :
                QMessageBox.critical(self.workspace._dialog_parent(), 'Failed to plot curve', str(exc))
        else :
            lens_model.clear_lt_curve(which=which)

    def _start_on_the_fly_predicted_images(self) :
        lens_model = self.workspace.lens_model
        if lens_model is None :
            QMessageBox.warning(self.workspace._dialog_parent(), 'No model', 'Please import a Lenstool model first.')
            return
        if getattr(lens_model, 'lt_z', None) is None :
            QMessageBox.warning(
                self.workspace._dialog_parent(),
                'No redshift set',
                'Please set the source plane redshift with the slider first.',
            )
            return
        try :
            lens_model.start_im2source()
        except Exception as exc :
            QMessageBox.critical(self.workspace._dialog_parent(), 'Failed to start on the fly predicted images', str(exc))

    def _open_image_simulator(self) :
        lens_model = self.workspace.lens_model
        if lens_model is None :
            QMessageBox.warning(self.workspace._dialog_parent(), 'No model', 'Please import a Lenstool model first.')
            return
        try :
            lens_model.start_simulate_image()
        except Exception as exc :
            QMessageBox.critical(self.workspace._dialog_parent(), 'Failed to open image simulator', str(exc))

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
        while self.catalog_subtabs.count() > 0 :
            widget = self.catalog_subtabs.widget(0)
            self.catalog_subtabs.removeTab(0)
            widget.deleteLater()

        for i, cat in enumerate(self.workspace.catalogs) :
            tab = self._make_catalog_controls_tab(cat)
            self.catalog_subtabs.addTab(tab, f'Cat {i + 1}')

        if 0 <= current_index < self.catalog_subtabs.count() :
            self.catalog_subtabs.setCurrentIndex(current_index)

    def _make_catalog_controls_tab(self, cat) :
        tab = QWidget()
        layout = QVBoxLayout(tab)

        # Persist the toggle/selection state on the Catalog instance itself, so
        # it survives the sub-tabs being rebuilt whenever a catalog is added.
        state = getattr(cat, '_plot_ui_state', None)
        if state is None :
            state = {
                'plot' : True, 'plot_column' : False, 'colname' : None,
                'selection_panel' : False, 'x_colname' : None, 'y_colname' : None,
                'r' : round(cat.color[0] * 255),
                'g' : round(cat.color[1] * 255),
                'b' : round(cat.color[2] * 255),
            }
            cat._plot_ui_state = state
        else :
            for channel in ('r', 'g', 'b') :
                if channel not in state :
                    idx = {'r' : 0, 'g' : 1, 'b' : 2}[channel]
                    state[channel] = round(cat.color[idx] * 255)

        layout.addWidget(self._make_separator())
        layout.addWidget(QLabel('Plotting options'))

        toggle_row = QHBoxLayout()
        plot_btn = self._make_link_toggle('Plot', lambda checked : self._toggle_catalog_plot(cat, checked))
        plot_column_btn = self._make_link_toggle('Plot column', lambda checked : self._toggle_catalog_plot_column(cat, checked))
        self._set_checked_silently(plot_btn, state['plot'])
        self._set_checked_silently(plot_column_btn, state['plot_column'])
        toggle_row.addWidget(plot_btn)
        toggle_row.addWidget(plot_column_btn)
        layout.addLayout(toggle_row)

        for channel, label_text in (('r', 'R'), ('g', 'G'), ('b', 'B')) :
            slider_row, slider, value_label = self._make_color_slider_row(label_text, state[channel])
            slider.valueChanged.connect(
                lambda value, channel=channel, value_label=value_label :
                self._on_catalog_color_slider_value_changed(cat, channel, value, value_label)
            )
            slider.sliderReleased.connect(lambda channel=channel : self._on_catalog_color_slider_released(cat))
            layout.addLayout(slider_row)

        colname_list = QListWidget()
        colname_list.setFixedHeight(150)
        colname_list.addItems([str(c) for c in cat.cat.colnames])
        if state['colname'] is not None :
            matches = colname_list.findItems(state['colname'], Qt.MatchExactly)
            if matches :
                colname_list.blockSignals(True)
                colname_list.setCurrentItem(matches[0])
                colname_list.blockSignals(False)
        colname_list.currentTextChanged.connect(lambda text : self._select_catalog_colname(cat, text))
        layout.addWidget(colname_list)

        export_mult_btn = QPushButton('Export to multiple image file…')
        export_mult_btn.clicked.connect(lambda : self._export_catalog_mult_file(cat))
        layout.addWidget(export_mult_btn)

        layout.addWidget(self._make_separator())
        layout.addWidget(QLabel('Interactive selection panel'))

        xy_row = QHBoxLayout()
        x_layout, x_list = self._make_labeled_colname_list('x-axis', cat, state['x_colname'])
        y_layout, y_list = self._make_labeled_colname_list('y-axis', cat, state['y_colname'])
        xy_row.addLayout(x_layout)
        xy_row.addLayout(y_layout)
        layout.addLayout(xy_row)

        selection_row = QHBoxLayout()
        selection_btn = self._make_link_toggle(
            'Make selection panel',
            lambda checked : self._toggle_catalog_selection_panel(cat, checked, selection_btn),
        )
        self._set_checked_silently(selection_btn, state['selection_panel'])
        selection_row.addWidget(selection_btn)

        export_potfile_btn = QPushButton('Export to potfile…')
        export_potfile_btn.clicked.connect(lambda : self._export_catalog_potfile(cat))
        selection_row.addWidget(export_potfile_btn)
        layout.addLayout(selection_row)

        x_list.currentTextChanged.connect(lambda text : self._select_catalog_x_colname(cat, text))
        y_list.currentTextChanged.connect(lambda text : self._select_catalog_y_colname(cat, text))

        layout.addStretch()
        return tab

    def _make_color_slider_row(self, label_text, initial_value) :
        row = QHBoxLayout()
        row.addWidget(QLabel(label_text))
        value_label = QLabel(str(initial_value))
        value_label.setMinimumWidth(28)
        row.addWidget(value_label)

        slider = QSlider(Qt.Horizontal)
        slider.setMinimum(0)
        slider.setMaximum(255)
        slider.setSingleStep(1)
        slider.setValue(initial_value)
        row.addWidget(slider, 1)
        return row, slider, value_label

    def _on_catalog_color_slider_value_changed(self, cat, channel, value, value_label) :
        value_label.setText(str(value))
        state = cat._plot_ui_state
        state[channel] = value

    def _on_catalog_color_slider_released(self, cat) :
        state = cat._plot_ui_state
        cat.color[0] = state['r'] / 255.
        cat.color[1] = state['g'] / 255.
        cat.color[2] = state['b'] / 255.
        if state['plot'] :
            cat.plot()
            if state['plot_column'] and state['colname'] is not None :
                cat.plot_column(state['colname'])

    def _make_labeled_colname_list(self, label_text, cat, initial_colname) :
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

        list_widget = QListWidget()
        list_widget.setFixedHeight(120)
        list_widget.addItems([str(c) for c in cat.cat.colnames])
        if initial_colname is not None :
            matches = list_widget.findItems(initial_colname, Qt.MatchExactly)
            if matches :
                list_widget.blockSignals(True)
                list_widget.setCurrentItem(matches[0])
                list_widget.blockSignals(False)
        v_layout.addWidget(list_widget)
        return v_layout, list_widget

    def _set_checked_silently(self, btn, checked) :
        btn.blockSignals(True)
        btn.setChecked(checked)
        btn.blockSignals(False)

    def _toggle_catalog_plot(self, cat, checked) :
        state = cat._plot_ui_state
        state['plot'] = checked
        if checked :
            cat.plot()
            if state['plot_column'] and state['colname'] is not None :
                cat.plot_column(state['colname'])
        else :
            cat.clear_ellipses()

    def _toggle_catalog_plot_column(self, cat, checked) :
        state = cat._plot_ui_state
        state['plot_column'] = checked
        if checked and state['colname'] is not None :
            cat.plot_column(state['colname'])
        else :
            cat.clear_column()

    def _select_catalog_colname(self, cat, colname) :
        if not colname :
            return
        state = cat._plot_ui_state
        state['colname'] = colname
        if state['plot_column'] :
            cat.plot_column(colname)

    def _toggle_catalog_selection_panel(self, cat, checked, btn) :
        state = cat._plot_ui_state
        if checked :
            x, y = state.get('x_colname'), state.get('y_colname')
            if x is None or y is None :
                QMessageBox.warning(
                    self.workspace._dialog_parent(),
                    'Missing axis',
                    'Please select both an x and y column first.',
                )
                self._set_checked_silently(btn, False)
                return
            if cat.Scatter_widget is not None :
                cat.remove_selection_panel()
            cat.make_selection_panel(xy_axes=[x, y])
            state['selection_panel'] = True
        else :
            if cat.Scatter_widget is not None :
                cat.remove_selection_panel()
            state['selection_panel'] = False

    def _select_catalog_x_colname(self, cat, colname) :
        if not colname :
            return
        state = cat._plot_ui_state
        state['x_colname'] = colname
        if state.get('selection_panel') and state.get('y_colname') is not None :
            cat.remove_selection_panel()
            cat.make_selection_panel(xy_axes=[colname, state['y_colname']])

    def _select_catalog_y_colname(self, cat, colname) :
        if not colname :
            return
        state = cat._plot_ui_state
        state['y_colname'] = colname
        if state.get('selection_panel') and state.get('x_colname') is not None :
            cat.remove_selection_panel()
            cat.make_selection_panel(xy_axes=[state['x_colname'], colname])

    def _export_catalog_mult_file(self, cat) :
        if cat.image is None or cat.image.image_path is None :
            base_dir = os.path.dirname(cat.path) if cat.path else os.getcwd()
        else :
            base_dir = os.path.dirname(cat.image.image_path)
        default_path = os.path.join(base_dir, 'mult.lenstool')
        path_str, _ = QFileDialog.getSaveFileName(
            self.workspace._dialog_parent(),
            'Export to multiple image file',
            default_path,
            'Lenstool mult files (*.lenstool);;All files (*)',
        )
        if not path_str :
            return

        try :
            cat.export_to_mult_file(file_path=path_str)
        except Exception as exc :
            QMessageBox.critical(
                self.workspace._dialog_parent(),
                'Failed to export multiple image file',
                str(exc),
            )

    def _export_catalog_potfile(self, cat) :
        colnames = [str(c) for c in cat.cat.colnames]
        state = cat._plot_ui_state
        default_idx = colnames.index(state['x_colname']) if state.get('x_colname') in colnames else 0

        mag_col, ok = QInputDialog.getItem(
            self.workspace._dialog_parent(),
            'Export to potfile',
            'Select magnitude column',
            colnames,
            default_idx,
            False,
        )
        if not ok :
            return

        base_dir = os.path.dirname(cat.path) if cat.path else os.getcwd()
        default_path = os.path.join(base_dir, 'exported_potfile.lenstool')
        path_str, _ = QFileDialog.getSaveFileName(
            self.workspace._dialog_parent(),
            'Export to potfile',
            default_path,
            'Lenstool potfiles (*.lenstool);;All files (*)',
        )
        if not path_str :
            return

        try :
            cat.export_to_potfile(file_path=path_str, mag_col=mag_col)
        except Exception as exc :
            QMessageBox.critical(self.workspace._dialog_parent(), 'Failed to export potfile', str(exc))
