# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'DEMto3D_dialog_base.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_DEMto3DDialogBase(object):
    def setupUi(self, DEMto3DDialogBase):
        DEMto3DDialogBase.setObjectName("DEMto3DDialogBase")
        DEMto3DDialogBase.setWindowModality(QtCore.Qt.WindowModal)
        DEMto3DDialogBase.resize(639, 740)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/plugins/DEMto3D/icons/demto3d.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        DEMto3DDialogBase.setWindowIcon(icon)
        self.verticalLayout = QtWidgets.QVBoxLayout(DEMto3DDialogBase)
        self.verticalLayout.setObjectName("verticalLayout")
        self.groupBox = QtWidgets.QGroupBox(DEMto3DDialogBase)
        self.groupBox.setMinimumSize(QtCore.QSize(0, 53))
        self.groupBox.setMaximumSize(QtCore.QSize(16777215, 53))
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.LayerComboBox = QtWidgets.QComboBox(self.groupBox)
        self.LayerComboBox.setEnabled(True)
        self.LayerComboBox.setMinimumSize(QtCore.QSize(0, 23))
        self.LayerComboBox.setObjectName("LayerComboBox")
        self.verticalLayout_2.addWidget(self.LayerComboBox)
        self.verticalLayout.addWidget(self.groupBox)
        self.groupBox_2 = QtWidgets.QGroupBox(DEMto3DDialogBase)
        self.groupBox_2.setMinimumSize(QtCore.QSize(0, 118))
        self.groupBox_2.setMaximumSize(QtCore.QSize(16777215, 118))
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(self.groupBox_2)
        self.label.setMinimumSize(QtCore.QSize(0, 25))
        self.label.setMaximumSize(QtCore.QSize(16777215, 25))
        self.label.setPixmap(QtGui.QPixmap(":/plugins/DEMto3D/icons/upleft.png"))
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.XMinLineEdit = QtWidgets.QLineEdit(self.groupBox_2)
        self.XMinLineEdit.setObjectName("XMinLineEdit")
        self.gridLayout.addWidget(self.XMinLineEdit, 0, 2, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.groupBox_2)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 0, 1, 1, 1)
        self.YMaxLineEdit = QtWidgets.QLineEdit(self.groupBox_2)
        self.YMaxLineEdit.setObjectName("YMaxLineEdit")
        self.gridLayout.addWidget(self.YMaxLineEdit, 0, 4, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.groupBox_2)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 0, 3, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.groupBox_2)
        self.label_4.setText("")
        self.label_4.setPixmap(QtGui.QPixmap(":/plugins/DEMto3D/icons/downright.png"))
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 1, 0, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.groupBox_2)
        self.label_5.setObjectName("label_5")
        self.gridLayout.addWidget(self.label_5, 1, 1, 1, 1)
        self.XMaxLineEdit = QtWidgets.QLineEdit(self.groupBox_2)
        self.XMaxLineEdit.setObjectName("XMaxLineEdit")
        self.gridLayout.addWidget(self.XMaxLineEdit, 1, 2, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.groupBox_2)
        self.label_6.setObjectName("label_6")
        self.gridLayout.addWidget(self.label_6, 1, 3, 1, 1)
        self.YMinLineEdit = QtWidgets.QLineEdit(self.groupBox_2)
        self.YMinLineEdit.setObjectName("YMinLineEdit")
        self.gridLayout.addWidget(self.YMinLineEdit, 1, 4, 1, 1)
        self.verticalLayout_3.addLayout(self.gridLayout)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setContentsMargins(330, -1, -1, -1)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.FullExtToolButton = QtWidgets.QToolButton(self.groupBox_2)
        self.FullExtToolButton.setMinimumSize(QtCore.QSize(20, 20))
        self.FullExtToolButton.setMaximumSize(QtCore.QSize(20, 20))
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(":/plugins/DEMto3D/icons/full_extension.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.FullExtToolButton.setIcon(icon1)
        self.FullExtToolButton.setObjectName("FullExtToolButton")
        self.horizontalLayout.addWidget(self.FullExtToolButton)
        self.LayerExtToolButton = QtWidgets.QToolButton(self.groupBox_2)
        self.LayerExtToolButton.setMinimumSize(QtCore.QSize(20, 20))
        self.LayerExtToolButton.setMaximumSize(QtCore.QSize(20, 20))
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(":/plugins/DEMto3D/icons/layer_extension.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.LayerExtToolButton.setIcon(icon2)
        self.LayerExtToolButton.setObjectName("LayerExtToolButton")
        self.horizontalLayout.addWidget(self.LayerExtToolButton)
        self.CustomExtToolButton = QtWidgets.QToolButton(self.groupBox_2)
        self.CustomExtToolButton.setMinimumSize(QtCore.QSize(20, 20))
        self.CustomExtToolButton.setMaximumSize(QtCore.QSize(20, 20))
        self.CustomExtToolButton.setText("")
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(":/plugins/DEMto3D/icons/cursor_extension.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.CustomExtToolButton.setIcon(icon3)
        self.CustomExtToolButton.setObjectName("CustomExtToolButton")
        self.horizontalLayout.addWidget(self.CustomExtToolButton)
        self.verticalLayout_3.addLayout(self.horizontalLayout)
        self.verticalLayout.addWidget(self.groupBox_2)
        self.groupBox_3 = QtWidgets.QGroupBox(DEMto3DDialogBase)
        self.groupBox_3.setObjectName("groupBox_3")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.groupBox_3)
        self.verticalLayout_4.setContentsMargins(15, -1, -1, -1)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.label_9 = QtWidgets.QLabel(self.groupBox_3)
        self.label_9.setObjectName("label_9")
        self.horizontalLayout_3.addWidget(self.label_9)
        self.SpacingLineEdit = QtWidgets.QLineEdit(self.groupBox_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.SpacingLineEdit.sizePolicy().hasHeightForWidth())
        self.SpacingLineEdit.setSizePolicy(sizePolicy)
        self.SpacingLineEdit.setMinimumSize(QtCore.QSize(85, 20))
        self.SpacingLineEdit.setMaximumSize(QtCore.QSize(85, 20))
        self.SpacingLineEdit.setObjectName("SpacingLineEdit")
        self.horizontalLayout_3.addWidget(self.SpacingLineEdit)
        self.label18 = QtWidgets.QLabel(self.groupBox_3)
        self.label18.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label18.setObjectName("label18")
        self.horizontalLayout_3.addWidget(self.label18)
        self.RecomSpacinglabel = QtWidgets.QLabel(self.groupBox_3)
        self.RecomSpacinglabel.setObjectName("RecomSpacinglabel")
        self.horizontalLayout_3.addWidget(self.RecomSpacinglabel)
        self.verticalLayout_4.addLayout(self.horizontalLayout_3)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setContentsMargins(-1, -1, 9, -1)
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.label_11 = QtWidgets.QLabel(self.groupBox_3)
        self.label_11.setObjectName("label_11")
        self.horizontalLayout_4.addWidget(self.label_11)
        self.label_12 = QtWidgets.QLabel(self.groupBox_3)
        self.label_12.setObjectName("label_12")
        self.horizontalLayout_4.addWidget(self.label_12)
        self.HeightLineEdit = QtWidgets.QLineEdit(self.groupBox_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.HeightLineEdit.sizePolicy().hasHeightForWidth())
        self.HeightLineEdit.setSizePolicy(sizePolicy)
        self.HeightLineEdit.setMinimumSize(QtCore.QSize(85, 20))
        self.HeightLineEdit.setMaximumSize(QtCore.QSize(85, 20))
        self.HeightLineEdit.setObjectName("HeightLineEdit")
        self.horizontalLayout_4.addWidget(self.HeightLineEdit)
        self.label_13 = QtWidgets.QLabel(self.groupBox_3)
        self.label_13.setObjectName("label_13")
        self.horizontalLayout_4.addWidget(self.label_13)
        self.WidthLineEdit = QtWidgets.QLineEdit(self.groupBox_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.WidthLineEdit.sizePolicy().hasHeightForWidth())
        self.WidthLineEdit.setSizePolicy(sizePolicy)
        self.WidthLineEdit.setMinimumSize(QtCore.QSize(85, 20))
        self.WidthLineEdit.setMaximumSize(QtCore.QSize(85, 20))
        self.WidthLineEdit.setObjectName("WidthLineEdit")
        self.horizontalLayout_4.addWidget(self.WidthLineEdit)
        self.verticalLayout_4.addLayout(self.horizontalLayout_4)
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_5.setContentsMargins(-1, -1, 230, -1)
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.label_15 = QtWidgets.QLabel(self.groupBox_3)
        self.label_15.setObjectName("label_15")
        self.horizontalLayout_5.addWidget(self.label_15)
        self.label_14 = QtWidgets.QLabel(self.groupBox_3)
        self.label_14.setMinimumSize(QtCore.QSize(12, 20))
        self.label_14.setMaximumSize(QtCore.QSize(20, 20))
        self.label_14.setObjectName("label_14")
        self.horizontalLayout_5.addWidget(self.label_14)
        self.ScaleLineEdit = QtWidgets.QLineEdit(self.groupBox_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ScaleLineEdit.sizePolicy().hasHeightForWidth())
        self.ScaleLineEdit.setSizePolicy(sizePolicy)
        self.ScaleLineEdit.setMinimumSize(QtCore.QSize(85, 20))
        self.ScaleLineEdit.setMaximumSize(QtCore.QSize(85, 20))
        self.ScaleLineEdit.setObjectName("ScaleLineEdit")
        self.horizontalLayout_5.addWidget(self.ScaleLineEdit)
        self.verticalLayout_4.addLayout(self.horizontalLayout_5)
        self.verticalLayout.addWidget(self.groupBox_3)
        self.groupBox_4 = QtWidgets.QGroupBox(DEMto3DDialogBase)
        self.groupBox_4.setMinimumSize(QtCore.QSize(0, 53))
        self.groupBox_4.setMaximumSize(QtCore.QSize(16777215, 53))
        self.groupBox_4.setObjectName("groupBox_4")
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout(self.groupBox_4)
        self.horizontalLayout_6.setContentsMargins(15, -1, -1, -1)
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.label_16 = QtWidgets.QLabel(self.groupBox_4)
        self.label_16.setMinimumSize(QtCore.QSize(0, 23))
        self.label_16.setObjectName("label_16")
        self.horizontalLayout_6.addWidget(self.label_16)
        self.ZScaleDoubleSpinBox = QtWidgets.QDoubleSpinBox(self.groupBox_4)
        self.ZScaleDoubleSpinBox.setMinimumSize(QtCore.QSize(0, 23))
        self.ZScaleDoubleSpinBox.setDecimals(1)
        self.ZScaleDoubleSpinBox.setMaximum(10.0)
        self.ZScaleDoubleSpinBox.setSingleStep(0.5)
        self.ZScaleDoubleSpinBox.setProperty("value", 1.0)
        self.ZScaleDoubleSpinBox.setObjectName("ZScaleDoubleSpinBox")
        self.horizontalLayout_6.addWidget(self.ZScaleDoubleSpinBox)
        self.verticalLayout.addWidget(self.groupBox_4)
        self.groupBox_5 = QtWidgets.QGroupBox(DEMto3DDialogBase)
        self.groupBox_5.setMinimumSize(QtCore.QSize(0, 81))
        self.groupBox_5.setMaximumSize(QtCore.QSize(16777215, 81))
        self.groupBox_5.setObjectName("groupBox_5")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.groupBox_5)
        self.verticalLayout_5.setContentsMargins(15, -1, -1, -1)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.BaseHeightLineEdit = QtWidgets.QLineEdit(self.groupBox_5)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.BaseHeightLineEdit.sizePolicy().hasHeightForWidth())
        self.BaseHeightLineEdit.setSizePolicy(sizePolicy)
        self.BaseHeightLineEdit.setObjectName("BaseHeightLineEdit")
        self.gridLayout_2.addWidget(self.BaseHeightLineEdit, 0, 1, 1, 1)
        self.ZMinLabel = QtWidgets.QLabel(self.groupBox_5)
        self.ZMinLabel.setMinimumSize(QtCore.QSize(75, 20))
        self.ZMinLabel.setMaximumSize(QtCore.QSize(75, 20))
        self.ZMinLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.ZMinLabel.setObjectName("ZMinLabel")
        self.gridLayout_2.addWidget(self.ZMinLabel, 0, 3, 1, 1)
        self.label_17 = QtWidgets.QLabel(self.groupBox_5)
        self.label_17.setObjectName("label_17")
        self.gridLayout_2.addWidget(self.label_17, 0, 0, 1, 1)
        self.ZMaxLabel = QtWidgets.QLabel(self.groupBox_5)
        self.ZMaxLabel.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.ZMaxLabel.setObjectName("ZMaxLabel")
        self.gridLayout_2.addWidget(self.ZMaxLabel, 1, 3, 1, 1)
        self.label_8 = QtWidgets.QLabel(self.groupBox_5)
        self.label_8.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_8.setObjectName("label_8")
        self.gridLayout_2.addWidget(self.label_8, 0, 2, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.groupBox_5)
        self.label_10.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_10.setObjectName("label_10")
        self.gridLayout_2.addWidget(self.label_10, 1, 2, 1, 1)
        self.HeightModelLabel = QtWidgets.QLabel(self.groupBox_5)
        self.HeightModelLabel.setObjectName("HeightModelLabel")
        self.gridLayout_2.addWidget(self.HeightModelLabel, 1, 1, 1, 1)
        self.label_21 = QtWidgets.QLabel(self.groupBox_5)
        self.label_21.setObjectName("label_21")
        self.gridLayout_2.addWidget(self.label_21, 1, 0, 1, 1)
        self.verticalLayout_5.addLayout(self.gridLayout_2)
        self.verticalLayout.addWidget(self.groupBox_5)
        self.groupBox_6 = QtWidgets.QGroupBox(DEMto3DDialogBase)
        self.groupBox_6.setMinimumSize(QtCore.QSize(0, 53))
        self.groupBox_6.setMaximumSize(QtCore.QSize(16777215, 53))
        self.groupBox_6.setObjectName("groupBox_6")
        self.formLayout = QtWidgets.QFormLayout(self.groupBox_6)
        self.formLayout.setObjectName("formLayout")
        self.RevereseZCheckBox = QtWidgets.QCheckBox(self.groupBox_6)
        self.RevereseZCheckBox.setMinimumSize(QtCore.QSize(0, 20))
        self.RevereseZCheckBox.setMaximumSize(QtCore.QSize(16777215, 20))
        self.RevereseZCheckBox.setObjectName("RevereseZCheckBox")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.RevereseZCheckBox)
        self.verticalLayout.addWidget(self.groupBox_6)
        self.horizontalLayout_9 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_9.setContentsMargins(220, -1, -1, -1)
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.STLToolButton = QtWidgets.QToolButton(DEMto3DDialogBase)
        self.STLToolButton.setMinimumSize(QtCore.QSize(150, 25))
        self.STLToolButton.setMaximumSize(QtCore.QSize(100, 25))
        self.STLToolButton.setObjectName("STLToolButton")
        self.horizontalLayout_9.addWidget(self.STLToolButton)
        self.CancelToolButton = QtWidgets.QToolButton(DEMto3DDialogBase)
        self.CancelToolButton.setMinimumSize(QtCore.QSize(100, 25))
        self.CancelToolButton.setMaximumSize(QtCore.QSize(100, 25))
        self.CancelToolButton.setObjectName("CancelToolButton")
        self.horizontalLayout_9.addWidget(self.CancelToolButton)
        self.verticalLayout.addLayout(self.horizontalLayout_9)
        self.label_2.setBuddy(self.XMinLineEdit)
        self.label_3.setBuddy(self.YMaxLineEdit)
        self.label_5.setBuddy(self.XMaxLineEdit)
        self.label_6.setBuddy(self.YMinLineEdit)
        self.label_9.setBuddy(self.SpacingLineEdit)
        self.label_12.setBuddy(self.HeightLineEdit)
        self.label_13.setBuddy(self.WidthLineEdit)
        self.label_15.setBuddy(self.ScaleLineEdit)
        self.label_16.setBuddy(self.ZScaleDoubleSpinBox)
        self.label_17.setBuddy(self.BaseHeightLineEdit)

        self.retranslateUi(DEMto3DDialogBase)
        QtCore.QMetaObject.connectSlotsByName(DEMto3DDialogBase)
        DEMto3DDialogBase.setTabOrder(self.CancelToolButton, self.STLToolButton)
        DEMto3DDialogBase.setTabOrder(self.STLToolButton, self.LayerComboBox)
        DEMto3DDialogBase.setTabOrder(self.LayerComboBox, self.XMinLineEdit)
        DEMto3DDialogBase.setTabOrder(self.XMinLineEdit, self.YMaxLineEdit)
        DEMto3DDialogBase.setTabOrder(self.YMaxLineEdit, self.XMaxLineEdit)
        DEMto3DDialogBase.setTabOrder(self.XMaxLineEdit, self.YMinLineEdit)
        DEMto3DDialogBase.setTabOrder(self.YMinLineEdit, self.FullExtToolButton)
        DEMto3DDialogBase.setTabOrder(self.FullExtToolButton, self.LayerExtToolButton)
        DEMto3DDialogBase.setTabOrder(self.LayerExtToolButton, self.CustomExtToolButton)
        DEMto3DDialogBase.setTabOrder(self.CustomExtToolButton, self.SpacingLineEdit)
        DEMto3DDialogBase.setTabOrder(self.SpacingLineEdit, self.HeightLineEdit)
        DEMto3DDialogBase.setTabOrder(self.HeightLineEdit, self.WidthLineEdit)
        DEMto3DDialogBase.setTabOrder(self.WidthLineEdit, self.ScaleLineEdit)
        DEMto3DDialogBase.setTabOrder(self.ScaleLineEdit, self.ZScaleDoubleSpinBox)
        DEMto3DDialogBase.setTabOrder(self.ZScaleDoubleSpinBox, self.BaseHeightLineEdit)
        DEMto3DDialogBase.setTabOrder(self.BaseHeightLineEdit, self.RevereseZCheckBox)

    def retranslateUi(self, DEMto3DDialogBase):
        _translate = QtCore.QCoreApplication.translate
        DEMto3DDialogBase.setWindowTitle(_translate("DEMto3DDialogBase", "DEM 3D printing"))
        self.groupBox.setTitle(_translate("DEMto3DDialogBase", "Layer to print"))
        self.groupBox_2.setTitle(_translate("DEMto3DDialogBase", "Print extent"))
        self.label_2.setText(_translate("DEMto3DDialogBase", "&X:"))
        self.label_3.setText(_translate("DEMto3DDialogBase", "&Y:"))
        self.label_5.setText(_translate("DEMto3DDialogBase", "X:"))
        self.label_6.setText(_translate("DEMto3DDialogBase", "Y:"))
        self.FullExtToolButton.setToolTip(_translate("DEMto3DDialogBase", "Select full extent"))
        self.FullExtToolButton.setStatusTip(_translate("DEMto3DDialogBase", "Select full extent"))
        self.FullExtToolButton.setWhatsThis(_translate("DEMto3DDialogBase", "Select full extent"))
        self.FullExtToolButton.setAccessibleName(_translate("DEMto3DDialogBase", "Select full extent"))
        self.LayerExtToolButton.setToolTip(_translate("DEMto3DDialogBase", "Select layer extent"))
        self.LayerExtToolButton.setStatusTip(_translate("DEMto3DDialogBase", "Select layer extent"))
        self.LayerExtToolButton.setWhatsThis(_translate("DEMto3DDialogBase", "Select layer extent"))
        self.LayerExtToolButton.setAccessibleName(_translate("DEMto3DDialogBase", "Select layer extent"))
        self.CustomExtToolButton.setToolTip(_translate("DEMto3DDialogBase", "Draw extent"))
        self.CustomExtToolButton.setStatusTip(_translate("DEMto3DDialogBase", "Draw extent"))
        self.CustomExtToolButton.setWhatsThis(_translate("DEMto3DDialogBase", "Draw extent"))
        self.CustomExtToolButton.setAccessibleName(_translate("DEMto3DDialogBase", "Draw extent"))
        self.groupBox_3.setTitle(_translate("DEMto3DDialogBase", "Model size"))
        self.label_9.setText(_translate("DEMto3DDialogBase", "Spaci&ng (mm):"))
        self.label18.setText(_translate("DEMto3DDialogBase", "Minimum recommended"))
        self.RecomSpacinglabel.setText(_translate("DEMto3DDialogBase", "0.2 mm"))
        self.label_11.setText(_translate("DEMto3DDialogBase", "Size:"))
        self.label_12.setText(_translate("DEMto3DDialogBase", "Wid&th (mm):"))
        self.label_13.setText(_translate("DEMto3DDialogBase", "&Lenght (mm):"))
        self.label_15.setText(_translate("DEMto3DDialogBase", "Scale:"))
        self.label_14.setText(_translate("DEMto3DDialogBase", "1:"))
        self.groupBox_4.setTitle(_translate("DEMto3DDialogBase", "Exaggeration terrain"))
        self.label_16.setText(_translate("DEMto3DDialogBase", "Exa&ggeration factor:"))
        self.groupBox_5.setTitle(_translate("DEMto3DDialogBase", "Height base"))
        self.ZMinLabel.setText(_translate("DEMto3DDialogBase", "0 m"))
        self.label_17.setText(_translate("DEMto3DDialogBase", "Height (m):"))
        self.ZMaxLabel.setText(_translate("DEMto3DDialogBase", "0 m"))
        self.label_8.setText(_translate("DEMto3DDialogBase", "Lowest point:"))
        self.label_10.setText(_translate("DEMto3DDialogBase", "Highest point:"))
        self.HeightModelLabel.setText(_translate("DEMto3DDialogBase", "0 mm"))
        self.label_21.setText(_translate("DEMto3DDialogBase", "Model height:"))
        self.groupBox_6.setTitle(_translate("DEMto3DDialogBase", "Other parameters"))
        self.RevereseZCheckBox.setText(_translate("DEMto3DDialogBase", "Terrain inversion"))
        self.STLToolButton.setText(_translate("DEMto3DDialogBase", "Export to STL"))
        self.CancelToolButton.setText(_translate("DEMto3DDialogBase", "Cancel"))