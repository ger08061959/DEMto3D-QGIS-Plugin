# -*- coding: utf-8 -*-
"""
/***************************************************************************
 AppONCE
                                 A QGIS plugin
 Creación de mapas en 3D
                              -------------------
        begin                : 2015-03-17
        git sha              : $Format:%H$
        copyright            : (C) 2015 by Francisco Javier Venceslá Simón
        email                : jawensi@gmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
from __future__ import absolute_import

from qgis.PyQt.QtCore import Qt
from qgis.PyQt.QtWidgets import QDialog, QListWidgetItem
from .SelectLayer_dialog_base import Ui_SelectLayer_dialog_base

""" 
a. Never used; 
b. QString is not contained in QtCore. 
c. QStrings are UTF8, lambda expression doesn't change a thing
try:
    _fromUtf8 = QtCore.QString().fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s
"""


class Dialog(QDialog, Ui_SelectLayer_dialog_base):
    def __init__(self, layers):
        """Constructor for the dialog."""
        QDialog.__init__(self, None, Qt.WindowStaysOnTopHint)
        self.ui = Ui_SelectLayer_dialog_base()
        self.ui.setupUi(self)

        self.ui.LayerList.clear()
        for layer in layers:
            item = QListWidgetItem()
            item.setText(layer.name())
            self.ui.LayerList.addItem(item)

        self.ui.buttonBox.accepted.connect(self.accept)
        self.ui.buttonBox.rejected.connect(self.reject)

    def get_layer(self):
        try:
            return self.ui.LayerList.currentItem().text()
        except AttributeError:
            pass
