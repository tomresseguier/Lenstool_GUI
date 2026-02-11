import sys
from PyQt5.QtWidgets import QApplication, QWidget
from PyQt5.QtCore import Qt

class ScrollDetector(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Scroll Detector")
        self.resize(400, 300)

    def wheelEvent(self, event):
        # Angle delta: traditional mouse wheel
        angle_delta = event.angleDelta().y()

        # Pixel delta: trackpads (usually)
        pixel_delta = event.pixelDelta().y()

        if pixel_delta != 0:
            print(f"Trackpad scroll detected: pixel delta = {pixel_delta}")
        else:
            print(f"Mouse wheel scroll detected: angle delta = {angle_delta}")

        event.accept()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ScrollDetector()
    window.show()
    sys.exit(app.exec_())
