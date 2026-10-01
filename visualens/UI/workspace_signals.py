from PyQt5.QtCore import QObject, pyqtSignal


class WorkspaceSignals(QObject) :
    """
    Change notifications emitted by the ``Visualens`` workspace.

    Views (e.g. the side panel) subscribe to these instead of being called
    directly by the workspace, so the workspace does not need to know which
    views exist, or whether any has been created yet.
    """
    image_changed = pyqtSignal()
    catalogs_changed = pyqtSignal()
    lens_model_changed = pyqtSignal()
    # Emitted when an existing side panel is moved into another image window,
    # meaning every part of it has to be re-synchronised with the workspace.
    window_attached = pyqtSignal()
