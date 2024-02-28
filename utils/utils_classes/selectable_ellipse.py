#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 00:33:34 2024

@author: Tom
"""


from PyQt5.QtWidgets import QGraphicsEllipseItem, QGraphicsScene, QGraphicsView, QGraphicsSceneMouseEvent, QApplication
from PyQt5.QtCore import Qt
import pyqtgraph as pg


class SelectableEllipse(QGraphicsEllipseItem):
    def __init__(self, x, y, width, height, idx, selection_mask, qtItems, selection_color=[255, 255, 255]):
        super(SelectableEllipse, self).__init__(x, y, width, height)
        self.idx = idx
        self.selection_mask = selection_mask
        self.qtItems = qtItems
        self.selection_color = selection_color
        self.setFlag(QGraphicsEllipseItem.ItemIsSelectable, True)

    def mousePressEvent(self, event: QGraphicsSceneMouseEvent):
        if event.button() == Qt.LeftButton:
            print(f"Ellipse selected: {self.rect().x()}, {self.rect().y()}")
            print('\n########\n Index: ' + str(self.idx) + '\n ########\n')
            self.selection_mask[self.idx] = not self.selection_mask[self.idx]
            self.qtItems[self.idx].setPen( pg.mkPen(self.selection_color + [255]) )
            self.qtItems[self.idx].setBrush( pg.mkBrush(self.selection_color + [127]) )
            
            