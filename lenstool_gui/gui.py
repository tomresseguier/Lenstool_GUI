from __future__ import annotations

import sys
from pathlib import Path

from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QFileDialog,
    QAction,
    QMessageBox,
)

from PyQt5.QtCore import Qt

from .fits_image import fits_image


class LensToolMainWindow(QMainWindow):
    """Main application window for the Lenstool GUI.

    The window currently provides a minimal File ▸ Open action that lets the
    user select a FITS file from disk. The selected file is used to create a
    :class:`lenstool_gui.fits_image.fits_image` instance, which is stored in
    :pyattr:`~LensToolMainWindow.image`.  Future menu actions and widgets can
    operate on this attribute to expose more of the *fits_image* API.
    """

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Lenstool GUI")
        self.resize(800, 600)

        # The currently loaded fits_image object (None until a file is opened)
        self.image: fits_image | None = None

        # Build menus / actions
        self._create_menus()

        # Placeholder central widget — replace with actual visualisation later
        # For now we just keep the default empty central widget created by Qt.

    # ---------------------------------------------------------------------
    # UI creation helpers
    # ---------------------------------------------------------------------
    def _create_menus(self) -> None:
        """Create the main menu bar and its actions."""
        menubar = self.menuBar()

        # ----- Image menu ---------------------------------------------------
        image_menu = menubar.addMenu("&Image")

        open_img_action = QAction("&Open…", self)
        open_img_action.setShortcut("Ctrl+O")
        open_img_action.triggered.connect(self._open_fits)  # type: ignore[arg-type]
        image_menu.addAction(open_img_action)

        boost_action = QAction("&Boost", self)
        boost_action.setShortcut("Ctrl+B")
        boost_action.triggered.connect(self._boost_image)  # type: ignore[arg-type]
        image_menu.addAction(boost_action)

        image_menu.addSeparator()

        exit_action = QAction("E&xit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)  # type: ignore[arg-type]
        image_menu.addAction(exit_action)

        # ----- Catalog menu -------------------------------------------------
        catalog_menu = menubar.addMenu("&Catalog")

        open_cat_action = QAction("&Open…", self)
        open_cat_action.triggered.connect(self._open_catalog)  # type: ignore[arg-type]
        catalog_menu.addAction(open_cat_action)

        plot_cat_action = QAction("&Plot", self)
        plot_cat_action.triggered.connect(self._plot_catalog)  # type: ignore[arg-type]
        catalog_menu.addAction(plot_cat_action)

        plot_column_action = QAction("Plot &Column…", self)
        plot_column_action.triggered.connect(self._plot_catalog_column)  # type: ignore[arg-type]
        catalog_menu.addAction(plot_column_action)

        selection_panel_action = QAction("Plot &Selection Panel…", self)
        selection_panel_action.triggered.connect(self._plot_selection_panel)  # type: ignore[arg-type]
        catalog_menu.addAction(selection_panel_action)

        clear_cat_action = QAction("&Clear", self)
        clear_cat_action.triggered.connect(self._clear_catalog)  # type: ignore[arg-type]
        catalog_menu.addAction(clear_cat_action)

        # ----- Lenstool menu -----------------------------------------------
        lt_menu = menubar.addMenu("&Lenstool")

        open_lt_action = QAction("&Open…", self)
        open_lt_action.triggered.connect(self._open_lenstool)  # type: ignore[arg-type]
        lt_menu.addAction(open_lt_action)

        plot_lt_action = QAction("&Plot", self)
        plot_lt_action.triggered.connect(self._plot_lenstool)  # type: ignore[arg-type]
        lt_menu.addAction(plot_lt_action)

        im2src_action = QAction("Start &im2source", self)
        im2src_action.triggered.connect(self._start_im2source)  # type: ignore[arg-type]
        lt_menu.addAction(im2src_action)

        bayes_action = QAction("Plot &Bayes", self)
        bayes_action.triggered.connect(self._plot_bayes)  # type: ignore[arg-type]
        lt_menu.addAction(bayes_action)

        clear_lt_action = QAction("&Clear", self)
        clear_lt_action.triggered.connect(self._clear_lenstool)  # type: ignore[arg-type]
        lt_menu.addAction(clear_lt_action)

    # ------------------------------------------------------------------
    # Slots / callbacks
    # ------------------------------------------------------------------
    def _open_fits(self) -> None:
        """Open a FITS image and create a :class:`fits_image` instance."""

        path_str, _ = QFileDialog.getOpenFileName(
            self,
            "Open FITS file",
            "",  # Start directory: current working directory
            "FITS images (*.fits *.fit *.fits.gz *.fit.gz);;All files (*)",
        )

        if not path_str:
            # User cancelled the dialog.
            return

        path = Path(path_str)

        try:
            self.image = fits_image(str(path), main_window=self)
        except Exception as exc:  # noqa: BLE001  (broad but user-facing)
            QMessageBox.critical(
                self,
                "Failed to load FITS",
                f"Could not load '{path.name}':\n{exc}",
                QMessageBox.StandardButton.Ok,
            )
            return

        # Success: inform the user (for now).
        # TODO: Update central widget once visualisation components are added.

    def _boost_image(self) -> None:
        """Toggle boost on currently loaded image."""
        if self.image is None:
            QMessageBox.warning(self, "No image", "Please open an image first.")
            return
        # If already boosted, unboost; else boost.
        if getattr(self.image, "boosted", False):
            self.image.unboost()
        else:
            self.image.boost()

    def _open_catalog(self) -> None:
        """Open a catalog file and import it into the current image."""
        if self.image is None:
            QMessageBox.warning(self, "No image", "Please open an image before importing a catalog.")
            return

        path_str, _ = QFileDialog.getOpenFileName(
            self,
            "Open catalog (FITS or ASCII)",
            "",
            "Catalog files (*.fits *.cat *.txt *.csv);;All files (*)",
        )
        if not path_str:
            return

        try:
            from .catalog import open_cat  # local import to avoid cost on startup

            cat = open_cat(path_str)
            self.image.import_catalog(cat)
            self._plot_catalog()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Failed to import catalog", str(exc))

    def _plot_catalog(self) -> None:
        """Plot the currently imported catalog."""
        if getattr(self.image, "imported_cat", None) is not None:
            self.image.imported_cat.plot()

    def _plot_catalog_column(self) -> None:
        """Plot a specific column from the currently imported catalog."""
        if getattr(self.image, "imported_cat", None) is None:
            QMessageBox.warning(self, "No catalog", "Please import a catalog first.")
            return

        cat = self.image.imported_cat
        colnames = [str(c) for c in cat.cat.colnames]

        from PyQt5.QtWidgets import QInputDialog

        col, ok = QInputDialog.getItem(self, "Plot column", "Select column", colnames, 0, False)
        if ok and col:
            cat.plot_column(col)

    def _plot_selection_panel(self) -> None:
        """Open selection panel after choosing x and y columns."""
        if getattr(self.image, "imported_cat", None) is None:
            QMessageBox.warning(self, "No catalog", "Please import a catalog first.")
            return

        cat_obj = self.image.imported_cat
        colnames = [str(c) for c in cat_obj.cat.colnames]

        from PyQt5.QtWidgets import QInputDialog

        x_col, ok1 = QInputDialog.getItem(self, "Select X axis", "X column", colnames, 0, False)
        if not ok1 or not x_col:
            return
        y_col, ok2 = QInputDialog.getItem(self, "Select Y axis", "Y column", colnames, 0, False)
        if not ok2 or not y_col:
            return

        # Clear existing catalog plots
        cat_obj.clear()

        # Re-import catalog with chosen columns as mag_colnames
        # Need original table; assume stored in cat_obj.cat
        self.image.import_catalog(cat_obj.cat, mag_colnames=[x_col, y_col])
        new_cat = self.image.imported_cat
        if new_cat is not None:
            new_cat.plot_selection_panel(xy_axes=[x_col, y_col])

    def _clear_catalog(self) -> None:
        """Clear the currently imported catalog."""
        if getattr(self.image, "imported_cat", None) is not None:
            self.image.imported_cat.clear()

    def _open_lenstool(self) -> None:
        """Open a Lenstool model directory and import it."""
        if self.image is None:
            QMessageBox.warning(self, "No image", "Please open an image before importing a Lenstool model.")
            return

        dir_path = QFileDialog.getExistingDirectory(self, "Select Lenstool model directory")
        if not dir_path:
            return
        try:
            self.image.import_lenstool(dir_path)
            self._plot_lenstool()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Failed to import model", str(exc))

    def _plot_lenstool(self) -> None:
        """Plot the currently imported Lenstool model."""
        if getattr(self.image, "lt", None) is not None:
            self.image.lt.plot()

    def _clear_lenstool(self) -> None:
        """Clear the currently imported Lenstool model."""
        if getattr(self.image, "lt", None) is not None:
            self.image.lt.clear()

    def _plot_bayes(self) -> None:
        lt_obj = getattr(self.image, "lt", None)
        if lt_obj is None:
            QMessageBox.warning(self, "No model", "Please import a Lenstool model first.")
            return
        try:
            lt_obj.plot_bayes()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Plot Bayes failed", str(exc))

    def _start_im2source(self) -> None:
        """Run SafeMode → set_lt_z → start_im2source sequence."""
        lt_obj = getattr(self.image, "lt", None)
        if lt_obj is None:
            QMessageBox.warning(self, "No model", "Please import a Lenstool model first.")
            return

        # Ask user for redshift value
        from PyQt5.QtWidgets import QInputDialog

        z_val, ok = QInputDialog.getDouble(self, "Set redshift", "Enter redshift (z):", decimals=4, min=0.0)
        if not ok:
            return

        try:
            lt_obj.SafeMode()
            lt_obj.set_lt_z(z_val)
            lt_obj.start_im2source()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "im2source failed", str(exc))


# -------------------------------------------------------------------------
# Entry point
# -------------------------------------------------------------------------

def main() -> None:  # pragma: no cover – CLI entry point
    """Run the Lenstool GUI application."""

    # A single QApplication is allowed per process.
    app = QApplication.instance() or QApplication(sys.argv)
    app.setAttribute(Qt.ApplicationAttribute.AA_ShareOpenGLContexts)

    window = LensToolMainWindow()
    window.show()

    # Start the Qt event loop.
    sys.exit(app.exec())


if __name__ == "__main__":  # pragma: no cover
    main()
