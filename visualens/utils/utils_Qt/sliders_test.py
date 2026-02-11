from PyQt5.QtWidgets import QMainWindow
import os
import sys

sys.path.append( os.path.abspath("/Users/tomresseguier/Library/Mobile Documents/com~apple~CloudDocs/RESEARCH/PROCESS/") )
from Lenstool_GUI.visualens.utils.utils_Qt.drag_widgets import DragPlotWidget_special

w = QMainWindow()
p = DragPlotWidget_special()
w.setCentralWidget(p)
w.show()

