# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file './canvas/TransformGuiTemplate.ui'
#
# Created: Wed Apr  4 12:25:29 2012
#      by: pyside-uic 0.2.11 running on PySide 1.0.6
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(169, 82)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        self.verticalLayout = QtGui.QVBoxLayout(Form)
        self.verticalLayout.setSpacing(1)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.translateLabel = QtGui.QLabel(Form)
        self.translateLabel.setObjectName("translateLabel")
        self.verticalLayout.addWidget(self.translateLabel)
        self.rotateLabel = QtGui.QLabel(Form)
        self.rotateLabel.setObjectName("rotateLabel")
        self.verticalLayout.addWidget(self.rotateLabel)
        self.scaleLabel = QtGui.QLabel(Form)
        self.scaleLabel.setObjectName("scaleLabel")
        self.verticalLayout.addWidget(self.scaleLabel)
        self.mirrorImageBtn = QtGui.QPushButton(Form)
        self.mirrorImageBtn.setToolTip("")
        self.mirrorImageBtn.setObjectName("mirrorImageBtn")
        self.verticalLayout.addWidget(self.mirrorImageBtn)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(QtGui.QApplication.translate("Form", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.translateLabel.setText(QtGui.QApplication.translate("Form", "Translate:", None, QtGui.QApplication.UnicodeUTF8))
        self.rotateLabel.setText(QtGui.QApplication.translate("Form", "Rotate:", None, QtGui.QApplication.UnicodeUTF8))
        self.scaleLabel.setText(QtGui.QApplication.translate("Form", "Scale:", None, QtGui.QApplication.UnicodeUTF8))
        self.mirrorImageBtn.setText(QtGui.QApplication.translate("Form", "Mirror", None, QtGui.QApplication.UnicodeUTF8))

