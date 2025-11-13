"""
MorphoScope - Main Window Implementation

This module contains the main window class for the MorphoScope application,
a tool for quantifying structural plasticity in 3D microscopy images of
neuronal projections.

Features:
- Multi-format image loading (CZI, TIF, LSM)
- Interactive ROI selection
- Image filtering and preprocessing
- Structural plasticity quantification
- Batch processing and CSV export

Author: Francisco Tassara
Date: 2025-11-12
Based on: Petsakou, Sapsis & Blau, Cell 2015
"""


from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor,
    QFont, QFontDatabase, QGradient, QIcon,
    QImage, QKeySequence, QLinearGradient, QPainter,
    QPalette, QPixmap, QRadialGradient, QTransform)
from PySide6.QtWidgets import (QApplication, QComboBox, QGroupBox, QHBoxLayout,
    QLabel, QLayout, QLineEdit, QListWidget, QCheckBox,
    QListWidgetItem, QMainWindow, QPushButton, QRadioButton,
    QSizePolicy, QSpacerItem, QSpinBox, QTextEdit,
    QVBoxLayout, QWidget, QProgressDialog)

from pyqtgraph import ImageView
import pyqtgraph

from PySide6.QtWidgets import QMessageBox, QFileDialog, QInputDialog
import sys
import os
import numpy as np
from matplotlib.path import Path as MPLPath
from scipy import ndimage as ndi
from shapely.geometry import Polygon

import csv
import matplotlib.pyplot as plt

import tifffile
from pylibCZIrw import czi

from scipy.ndimage import median_filter
from skimage.filters import gaussian

from config import Config, setup_logging
from image_processor import ImageProcessor, validate_parameters

from typing import Optional, Tuple, List

# Configurar logging
import logging
setup_logging()
logger = logging.getLogger(__name__)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(1316, 815)
        icon = QIcon()
        icon.addFile(u"logo.ico", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        MainWindow.setWindowIcon(icon)
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        self.horizontalLayout_10 = QHBoxLayout(self.centralwidget)
        self.horizontalLayout_10.setObjectName(u"horizontalLayout_10")
        self.verticalLayout_11 = QVBoxLayout()
        self.verticalLayout_11.setObjectName(u"verticalLayout_11")
        self.horizontalLayout_7 = QHBoxLayout()
        self.horizontalLayout_7.setObjectName(u"horizontalLayout_7")
        self.pushButton_loadImage = QPushButton(self.centralwidget)
        self.pushButton_loadImage.setObjectName(u"pushButton_loadImage")
        font = QFont()
        font.setPointSize(10)
        font.setBold(True)
        self.pushButton_loadImage.setFont(font)

        self.horizontalLayout_7.addWidget(self.pushButton_loadImage)

        self.pushButton_clean = QPushButton(self.centralwidget)
        self.pushButton_clean.setObjectName(u"pushButton_clean")
        font1 = QFont()
        font1.setPointSize(10)
        self.pushButton_clean.setFont(font1)

        self.horizontalLayout_7.addWidget(self.pushButton_clean)


        self.verticalLayout_11.addLayout(self.horizontalLayout_7)

        self.listWidget_images = QListWidget(self.centralwidget)
        self.listWidget_images.setObjectName(u"listWidget_images")

        self.verticalLayout_11.addWidget(self.listWidget_images)

        self.groupBox_properties = QGroupBox(self.centralwidget)
        self.groupBox_properties.setObjectName(u"groupBox_properties")
        self.groupBox_properties.setFont(font1)
        self.horizontalLayout = QHBoxLayout(self.groupBox_properties)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.verticalLayout_10 = QVBoxLayout()
        self.verticalLayout_10.setObjectName(u"verticalLayout_10")
        self.label_5 = QLabel(self.groupBox_properties)
        self.label_5.setObjectName(u"label_5")

        self.verticalLayout_10.addWidget(self.label_5)

        self.label_2 = QLabel(self.groupBox_properties)
        self.label_2.setObjectName(u"label_2")

        self.verticalLayout_10.addWidget(self.label_2)


        self.horizontalLayout.addLayout(self.verticalLayout_10)

        self.verticalLayout_3 = QVBoxLayout()
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.lineEdit_image_size_X = QLineEdit(self.groupBox_properties)
        self.lineEdit_image_size_X.setObjectName(u"lineEdit_image_size_X")

        self.verticalLayout_3.addWidget(self.lineEdit_image_size_X)

        self.lineEdit_pixel_size_X = QLineEdit(self.groupBox_properties)
        self.lineEdit_pixel_size_X.setObjectName(u"lineEdit_pixel_size_X")

        self.verticalLayout_3.addWidget(self.lineEdit_pixel_size_X)


        self.horizontalLayout.addLayout(self.verticalLayout_3)

        self.verticalLayout_5 = QVBoxLayout()
        self.verticalLayout_5.setObjectName(u"verticalLayout_5")
        self.label_6 = QLabel(self.groupBox_properties)
        self.label_6.setObjectName(u"label_6")

        self.verticalLayout_5.addWidget(self.label_6)

        self.label_9 = QLabel(self.groupBox_properties)
        self.label_9.setObjectName(u"label_9")

        self.verticalLayout_5.addWidget(self.label_9)


        self.horizontalLayout.addLayout(self.verticalLayout_5)

        self.verticalLayout_6 = QVBoxLayout()
        self.verticalLayout_6.setObjectName(u"verticalLayout_6")
        self.lineEdit_image_size_Y = QLineEdit(self.groupBox_properties)
        self.lineEdit_image_size_Y.setObjectName(u"lineEdit_image_size_Y")

        self.verticalLayout_6.addWidget(self.lineEdit_image_size_Y)

        self.lineEdit_pixel_size_Y = QLineEdit(self.groupBox_properties)
        self.lineEdit_pixel_size_Y.setObjectName(u"lineEdit_pixel_size_Y")

        self.verticalLayout_6.addWidget(self.lineEdit_pixel_size_Y)


        self.horizontalLayout.addLayout(self.verticalLayout_6)

        self.verticalLayout_7 = QVBoxLayout()
        self.verticalLayout_7.setObjectName(u"verticalLayout_7")
        self.label_7 = QLabel(self.groupBox_properties)
        self.label_7.setObjectName(u"label_7")

        self.verticalLayout_7.addWidget(self.label_7)

        self.label_10 = QLabel(self.groupBox_properties)
        self.label_10.setObjectName(u"label_10")

        self.verticalLayout_7.addWidget(self.label_10)


        self.horizontalLayout.addLayout(self.verticalLayout_7)

        self.verticalLayout_8 = QVBoxLayout()
        self.verticalLayout_8.setObjectName(u"verticalLayout_8")
        self.lineEdit_image_size_Z = QLineEdit(self.groupBox_properties)
        self.lineEdit_image_size_Z.setObjectName(u"lineEdit_image_size_Z")

        self.verticalLayout_8.addWidget(self.lineEdit_image_size_Z)

        self.lineEdit_pixel_size_Z = QLineEdit(self.groupBox_properties)
        self.lineEdit_pixel_size_Z.setObjectName(u"lineEdit_pixel_size_Z")

        self.verticalLayout_8.addWidget(self.lineEdit_pixel_size_Z)


        self.horizontalLayout.addLayout(self.verticalLayout_8)

        self.verticalLayout_9 = QVBoxLayout()
        self.verticalLayout_9.setObjectName(u"verticalLayout_9")
        self.label_8 = QLabel(self.groupBox_properties)
        self.label_8.setObjectName(u"label_8")

        self.verticalLayout_9.addWidget(self.label_8)

        self.label_11 = QLabel(self.groupBox_properties)
        self.label_11.setObjectName(u"label_11")

        self.verticalLayout_9.addWidget(self.label_11)


        self.horizontalLayout.addLayout(self.verticalLayout_9)


        self.verticalLayout_11.addWidget(self.groupBox_properties)

        self.groupBox_visualizer = QGroupBox(self.centralwidget)
        self.groupBox_visualizer.setObjectName(u"groupBox_visualizer")
        self.groupBox_visualizer.setFont(font1)
        self.verticalLayout = QVBoxLayout(self.groupBox_visualizer)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.horizontalLayout_4 = QHBoxLayout()
        self.horizontalLayout_4.setObjectName(u"horizontalLayout_4")
        self.label_4 = QLabel(self.groupBox_visualizer)
        self.label_4.setObjectName(u"label_4")

        self.horizontalLayout_4.addWidget(self.label_4)

        self.comboBox_channel_selector = QComboBox(self.groupBox_visualizer)
        self.comboBox_channel_selector.setObjectName(u"comboBox_channel_selector")

        self.horizontalLayout_4.addWidget(self.comboBox_channel_selector)


        self.verticalLayout.addLayout(self.horizontalLayout_4)

        self.horizontalLayout_6 = QHBoxLayout()
        self.horizontalLayout_6.setObjectName(u"horizontalLayout_6")
        self.label_12 = QLabel(self.groupBox_visualizer)
        self.label_12.setObjectName(u"label_12")

        self.horizontalLayout_6.addWidget(self.label_12)

        self.radioButton_plotZproject = QRadioButton(self.groupBox_visualizer)
        self.radioButton_plotZproject.setObjectName(u"radioButton_plotZproject")
        self.radioButton_plotZproject.setChecked(True)

        self.horizontalLayout_6.addWidget(self.radioButton_plotZproject)

        self.radioButton_plotStack = QRadioButton(self.groupBox_visualizer)
        self.radioButton_plotStack.setObjectName(u"radioButton_plotStack")

        self.horizontalLayout_6.addWidget(self.radioButton_plotStack)


        self.verticalLayout.addLayout(self.horizontalLayout_6)


        self.verticalLayout_11.addWidget(self.groupBox_visualizer)

        self.groupBox_channels = QGroupBox(self.centralwidget)
        self.groupBox_channels.setObjectName(u"groupBox_channels")
        font2 = QFont()
        font2.setPointSize(10)
        font2.setBold(False)
        self.groupBox_channels.setFont(font2)
        self.verticalLayout_2 = QVBoxLayout(self.groupBox_channels)
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.horizontalLayout_2 = QHBoxLayout()
        self.horizontalLayout_2.setSpacing(1)
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.label = QLabel(self.groupBox_channels)
        self.label.setObjectName(u"label")

        self.horizontalLayout_2.addWidget(self.label)

        self.comboBox_plasticityChannel = QComboBox(self.groupBox_channels)
        self.comboBox_plasticityChannel.setObjectName(u"comboBox_plasticityChannel")

        self.horizontalLayout_2.addWidget(self.comboBox_plasticityChannel)

        self.comboBox_filter_type_chP = QComboBox(self.groupBox_channels)
        self.comboBox_filter_type_chP.addItem("")
        self.comboBox_filter_type_chP.addItem("")
        self.comboBox_filter_type_chP.addItem("")
        self.comboBox_filter_type_chP.addItem("")
        self.comboBox_filter_type_chP.addItem("")
        self.comboBox_filter_type_chP.setObjectName(u"comboBox_filter_type_chP")

        self.horizontalLayout_2.addWidget(self.comboBox_filter_type_chP)

        self.pushButton_applyfilter_chP = QPushButton(self.groupBox_channels)
        self.pushButton_applyfilter_chP.setObjectName(u"pushButton_applyfilter_chP")

        self.horizontalLayout_2.addWidget(self.pushButton_applyfilter_chP)

        self.pushButton_undofilter_chP = QPushButton(self.groupBox_channels)
        self.pushButton_undofilter_chP.setObjectName(u"pushButton_undofilter_chP")

        self.horizontalLayout_2.addWidget(self.pushButton_undofilter_chP)

        self.horizontalLayout_2.setStretch(0, 1)
        self.horizontalLayout_2.setStretch(1, 5)
        self.horizontalLayout_2.setStretch(2, 1)
        self.horizontalLayout_2.setStretch(3, 1)
        self.horizontalLayout_2.setStretch(4, 1)

        self.verticalLayout_2.addLayout(self.horizontalLayout_2)

        self.horizontalLayout_3 = QHBoxLayout()
        self.horizontalLayout_3.setSpacing(1)
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.horizontalLayout_3.setSizeConstraint(QLayout.SizeConstraint.SetDefaultConstraint)
        self.label_3 = QLabel(self.groupBox_channels)
        self.label_3.setObjectName(u"label_3")

        self.horizontalLayout_3.addWidget(self.label_3)

        self.comboBox_fluoChannel = QComboBox(self.groupBox_channels)
        self.comboBox_fluoChannel.setObjectName(u"comboBox_fluoChannel")

        self.horizontalLayout_3.addWidget(self.comboBox_fluoChannel)

        self.comboBox_filter_type_chF = QComboBox(self.groupBox_channels)
        self.comboBox_filter_type_chF.addItem("")
        self.comboBox_filter_type_chF.addItem("")
        self.comboBox_filter_type_chF.addItem("")
        self.comboBox_filter_type_chF.addItem("")
        self.comboBox_filter_type_chF.addItem("")
        self.comboBox_filter_type_chF.setObjectName(u"comboBox_filter_type_chF")

        self.horizontalLayout_3.addWidget(self.comboBox_filter_type_chF)

        self.pushButton_applyfilter_chF = QPushButton(self.groupBox_channels)
        self.pushButton_applyfilter_chF.setObjectName(u"pushButton_applyfilter_chF")

        self.horizontalLayout_3.addWidget(self.pushButton_applyfilter_chF)

        self.pushButton_undofilter_chF = QPushButton(self.groupBox_channels)
        self.pushButton_undofilter_chF.setObjectName(u"pushButton_undofilter_chF")

        self.horizontalLayout_3.addWidget(self.pushButton_undofilter_chF)

        self.horizontalLayout_3.setStretch(0, 1)
        self.horizontalLayout_3.setStretch(1, 5)
        self.horizontalLayout_3.setStretch(2, 1)
        self.horizontalLayout_3.setStretch(3, 1)
        self.horizontalLayout_3.setStretch(4, 1)

        self.verticalLayout_2.addLayout(self.horizontalLayout_3)


        self.verticalLayout_11.addWidget(self.groupBox_channels)

        self.groupBox_observations = QGroupBox(self.centralwidget)
        self.groupBox_observations.setObjectName(u"groupBox_observations")
        self.groupBox_observations.setFont(font1)
        self.verticalLayout_4 = QVBoxLayout(self.groupBox_observations)
        self.verticalLayout_4.setObjectName(u"verticalLayout_4")
        self.textEdit_observation = QTextEdit(self.groupBox_observations)
        self.textEdit_observation.setObjectName(u"textEdit_observation")
        self.textEdit_observation.setFont(font1)

        self.verticalLayout_4.addWidget(self.textEdit_observation)


        self.verticalLayout_11.addWidget(self.groupBox_observations)

        self.groupBox_roi = QGroupBox(self.centralwidget)
        self.groupBox_roi.setObjectName(u"groupBox_roi")
        self.groupBox_roi.setFont(font1)
        self.horizontalLayout_5 = QHBoxLayout(self.groupBox_roi)
        self.horizontalLayout_5.setObjectName(u"horizontalLayout_5")
        self.pushButton_roi = QPushButton(self.groupBox_roi)
        self.pushButton_roi.setObjectName(u"pushButton_roi")
        self.pushButton_roi.setFont(font1)

        self.horizontalLayout_5.addWidget(self.pushButton_roi)

        self.pushButton_apply_roi = QPushButton(self.groupBox_roi)
        self.pushButton_apply_roi.setObjectName(u"pushButton_apply_roi")
        self.pushButton_apply_roi.setFont(font1)

        self.horizontalLayout_5.addWidget(self.pushButton_apply_roi)


        self.verticalLayout_11.addWidget(self.groupBox_roi)

        self.horizontalLayout_8 = QHBoxLayout()
        self.horizontalLayout_8.setObjectName(u"horizontalLayout_8")
        self.groupBox = QGroupBox(self.centralwidget)
        self.groupBox.setObjectName(u"groupBox")
        self.groupBox.setFont(font1)
        self.horizontalLayout_9 = QHBoxLayout(self.groupBox)
        self.horizontalLayout_9.setObjectName(u"horizontalLayout_9")
        self.horizontalLayout_9.setContentsMargins(-1, -1, -1, 3)
        self.label_13 = QLabel(self.groupBox)
        self.label_13.setObjectName(u"label_13")

        self.horizontalLayout_9.addWidget(self.label_13)

        self.spinBox_zmin = QSpinBox(self.groupBox)
        self.spinBox_zmin.setObjectName(u"spinBox_zmin")

        self.horizontalLayout_9.addWidget(self.spinBox_zmin)

        self.horizontalSpacer_2 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_9.addItem(self.horizontalSpacer_2)

        self.label_14 = QLabel(self.groupBox)
        self.label_14.setObjectName(u"label_14")

        self.horizontalLayout_9.addWidget(self.label_14)

        self.spinBox_zmax = QSpinBox(self.groupBox)
        self.spinBox_zmax.setObjectName(u"spinBox_zmax")

        self.horizontalLayout_9.addWidget(self.spinBox_zmax)

        self.horizontalLayout_9.setStretch(0, 1)
        self.horizontalLayout_9.setStretch(1, 1)
        self.horizontalLayout_9.setStretch(2, 4)
        self.horizontalLayout_9.setStretch(3, 1)
        self.horizontalLayout_9.setStretch(4, 1)

        self.horizontalLayout_8.addWidget(self.groupBox)

        self.checkBox_show_distributions = QCheckBox(self.centralwidget)
        self.checkBox_show_distributions.setObjectName(u"checkBox_show_distributions")

        self.horizontalLayout_8.addWidget(self.checkBox_show_distributions)


        self.verticalLayout_11.addLayout(self.horizontalLayout_8)

        self.pushButton_procces = QPushButton(self.centralwidget)
        self.pushButton_procces.setObjectName(u"pushButton_procces")
        font3 = QFont()
        font3.setPointSize(12)
        font3.setBold(True)
        self.pushButton_procces.setFont(font3)

        self.verticalLayout_11.addWidget(self.pushButton_procces)

        self.verticalLayout_11.setStretch(0, 1)
        self.verticalLayout_11.setStretch(1, 2)
        self.verticalLayout_11.setStretch(2, 1)
        self.verticalLayout_11.setStretch(3, 1)
        self.verticalLayout_11.setStretch(4, 1)
        self.verticalLayout_11.setStretch(5, 1)
        self.verticalLayout_11.setStretch(6, 1)
        self.verticalLayout_11.setStretch(7, 1)
        self.verticalLayout_11.setStretch(8, 1)

        self.horizontalLayout_10.addLayout(self.verticalLayout_11)

        self.graphWidget = ImageView(self.centralwidget)
        self.graphWidget.setObjectName(u"graphWidget")
        self.graphWidget.setStyleSheet(u"background-color: rgb(0, 0, 0);")

        self.horizontalLayout_10.addWidget(self.graphWidget)

        self.horizontalLayout_10.setStretch(0, 1)
        self.horizontalLayout_10.setStretch(1, 2)
        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)

        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"MorphoScope", None))
        self.pushButton_loadImage.setText(QCoreApplication.translate("MainWindow", u"Load images", None))
        self.pushButton_clean.setText(QCoreApplication.translate("MainWindow", u"Clean", None))
        self.groupBox_properties.setTitle(QCoreApplication.translate("MainWindow", u"Image Properties", None))
        self.label_5.setText(QCoreApplication.translate("MainWindow", u"Image Size (X,Y,Z):", None))
        self.label_2.setText(QCoreApplication.translate("MainWindow", u"Voxel Size (X,Y,Z):", None))
        self.lineEdit_image_size_X.setText(QCoreApplication.translate("MainWindow", u"X", None))
        self.lineEdit_pixel_size_X.setText(QCoreApplication.translate("MainWindow", u"X", None))
        self.label_6.setText(QCoreApplication.translate("MainWindow", u"x", None))
        self.label_9.setText(QCoreApplication.translate("MainWindow", u"x", None))
        self.lineEdit_image_size_Y.setText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.lineEdit_pixel_size_Y.setText(QCoreApplication.translate("MainWindow", u"Y", None))
        self.label_7.setText(QCoreApplication.translate("MainWindow", u"x", None))
        self.label_10.setText(QCoreApplication.translate("MainWindow", u"x", None))
        self.lineEdit_image_size_Z.setText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.lineEdit_pixel_size_Z.setText(QCoreApplication.translate("MainWindow", u"Z", None))
        self.label_8.setText(QCoreApplication.translate("MainWindow", u"px", None))
        self.label_11.setText(QCoreApplication.translate("MainWindow", u"um", None))
        self.groupBox_visualizer.setTitle(QCoreApplication.translate("MainWindow", u"Visualizer", None))
        self.label_4.setText(QCoreApplication.translate("MainWindow", u"Channel:", None))
        self.label_12.setText(QCoreApplication.translate("MainWindow", u"Visualization Mode:", None))
        self.radioButton_plotZproject.setText(QCoreApplication.translate("MainWindow", u"Z project", None))
        self.radioButton_plotStack.setText(QCoreApplication.translate("MainWindow", u"Stack", None))
        self.groupBox_channels.setTitle(QCoreApplication.translate("MainWindow", u"Channels for quantification", None))
        self.label.setText(QCoreApplication.translate("MainWindow", u"Plasticity:     ", None))
        self.comboBox_filter_type_chP.setItemText(0, QCoreApplication.translate("MainWindow", u"Select Filter", None))
        self.comboBox_filter_type_chP.setItemText(1, QCoreApplication.translate("MainWindow", u"Threshold", None))
        self.comboBox_filter_type_chP.setItemText(2, QCoreApplication.translate("MainWindow", u"Gaussian Blur", None))
        self.comboBox_filter_type_chP.setItemText(3, QCoreApplication.translate("MainWindow", u"Median Filter", None))
        self.comboBox_filter_type_chP.setItemText(4, QCoreApplication.translate("MainWindow", u"Rolling Ball", None))

        self.pushButton_applyfilter_chP.setText(QCoreApplication.translate("MainWindow", u"Apply", None))
        self.pushButton_undofilter_chP.setText(QCoreApplication.translate("MainWindow", u"Undo", None))
        self.label_3.setText(QCoreApplication.translate("MainWindow", u"Fluoresence:", None))
        self.comboBox_filter_type_chF.setItemText(0, QCoreApplication.translate("MainWindow", u"Select Filter", None))
        self.comboBox_filter_type_chF.setItemText(1, QCoreApplication.translate("MainWindow", u"Threshold", None))
        self.comboBox_filter_type_chF.setItemText(2, QCoreApplication.translate("MainWindow", u"Gaussian Blur", None))
        self.comboBox_filter_type_chF.setItemText(3, QCoreApplication.translate("MainWindow", u"Median Filter", None))
        self.comboBox_filter_type_chF.setItemText(4, QCoreApplication.translate("MainWindow", u"Rolling Ball", None))

        self.pushButton_applyfilter_chF.setText(QCoreApplication.translate("MainWindow", u"Apply", None))
        self.pushButton_undofilter_chF.setText(QCoreApplication.translate("MainWindow", u"Undo", None))
        self.groupBox_observations.setTitle(QCoreApplication.translate("MainWindow", u"Observation", None))
        self.groupBox_roi.setTitle(QCoreApplication.translate("MainWindow", u"Region of Interest (ROI) Selection", None))
        self.pushButton_roi.setText(QCoreApplication.translate("MainWindow", u"Create ROI", None))
        self.pushButton_apply_roi.setText(QCoreApplication.translate("MainWindow", u"Apply ROI", None))
        self.groupBox.setTitle(QCoreApplication.translate("MainWindow", u"Slice Selection", None))
        self.label_13.setText(QCoreApplication.translate("MainWindow", u"From:", None))
        self.label_14.setText(QCoreApplication.translate("MainWindow", u"To:", None))
#if QT_CONFIG(tooltip)
        self.checkBox_show_distributions.setToolTip(QCoreApplication.translate("MainWindow", u"Display X, Y, Z distribution plots after processing", None))
#endif // QT_CONFIG(tooltip)
        self.checkBox_show_distributions.setText(QCoreApplication.translate("MainWindow", u"Show distribution plots", None))
        self.pushButton_procces.setText(QCoreApplication.translate("MainWindow", u"Process image and save", None))
    # retranslateUi


class MyMainWindow(QMainWindow):
    """
    Main application window for MorphoScope structural plasticity analysis.
    
    This class provides the graphical user interface and core functionality for:
    - Loading and visualizing 3D microscopy images
    - Applying preprocessing filters
    - Defining regions of interest (ROI)
    - Calculating structural plasticity metrics
    - Exporting results to CSV
    
    Attributes
    ----------
    ui : Ui_MainWindow
        Auto-generated UI from Qt Designer (PySide6)
    roi : pyqtgraph.PolyLineROI
        Current polygonal region of interest
    creating_roi : bool
        Flag indicating if ROI creation is active
    current_image_data : List[np.ndarray]
        Loaded image data (one array per channel)
    current_metadata_channel : List[str]
        Channel names/labels
    current_metadata_voxel_size : Tuple[float, float, float]
        Voxel dimensions (x, y, z) in micrometers
    current_metadata_dimension : Tuple[int, int, int]
        Image dimensions (width, height, depth)
    complexity_channel : np.ndarray
        Masked image for plasticity analysis
    fluor_channel : Optional[np.ndarray]
        Masked image for fluorescence analysis (if available)
    csv_file_path : str
        Path to output CSV file
    
    Notes
    -----
    - The UI must be generated from Qt Designer (.ui file) using pyside6-uic
    - Supports CZI, TIF, and LSM microscopy formats
    - Uses ImageProcessor for quantification (see image_processor.py)
    
    See Also
    --------
    image_processor.ImageProcessor : Core processing engine
    """
    
    def __init__(self):
        """
        Initialize the main window and connect UI elements to handlers.
        
        This method:
        1. Sets up the auto-generated UI
        2. Connects buttons to their respective functions
        3. Initializes instance variables
        """
        super().__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        
        # ================================================================
        # CONNECT UI ELEMENTS TO HANDLERS
        # ================================================================
        
        # Image loading and management
        self.ui.pushButton_loadImage.clicked.connect(self.load_images)
        self.ui.pushButton_clean.clicked.connect(self.clean_list)
        self.ui.listWidget_images.itemSelectionChanged.connect(self.on_image_selected)
        
        # ROI selection
        self.ui.pushButton_roi.clicked.connect(self.enable_polygonal_roi_creation)
        self.ui.pushButton_apply_roi.clicked.connect(self.apply_roi_mask)
        
        # Image display controls
        self.ui.comboBox_channel_selector.currentIndexChanged.connect(self.update_display)
        self.ui.radioButton_plotZproject.toggled.connect(self.update_display)
        self.ui.radioButton_plotStack.toggled.connect(self.update_display)
        
        # Filter controls
        self.ui.pushButton_applyfilter_chP.clicked.connect(
            lambda: self.apply_selected_filter('chP')
        )
        self.ui.pushButton_applyfilter_chF.clicked.connect(
            lambda: self.apply_selected_filter('chF')
        )
        self.ui.pushButton_undofilter_chP.clicked.connect(
            lambda: self.undo_filter('chP')
        )
        self.ui.pushButton_undofilter_chF.clicked.connect(
            lambda: self.undo_filter('chF')
        )
        
        # Processing
        self.ui.pushButton_procces.clicked.connect(self.process)
        
        # ================================================================
        # INITIALIZE INSTANCE VARIABLES
        # ================================================================
        
        # ROI management
        self.roi = None
        self.creating_roi = False
        self.temp_points = []
        self.polygon_points = []
        self.AArea = 0.0
        
        # Image data
        self.current_image_filepath = None
        self.current_image_data = None
        self.original_image_data = None  # Backup for undo filter
        self.current_metadata_channel = None
        self.current_metadata_dimension = None
        self.current_metadata_voxel_size = None
        
        # Processed channels
        self.complexity_channel = None
        self.fluor_channel = None
        
        # Output
        self.csv_file_path = None
        
        logger.info("MorphoScope main window initialized")
    
    
    # ====================================================================
    # IMAGE LOADING METHODS
    # ====================================================================
    
    def load_images(self):
        """
        Load microscopy images and set up CSV output file.
        
        This method:
        1. Opens file dialog to select image files
        2. Prompts for output CSV filename
        3. Creates CSV with headers if new
        4. Adds images to processing list
        
        Supported Formats
        -----------------
        - .czi : Zeiss CZI files
        - .tif : TIFF stacks
        - .lsm : Zeiss LSM files
        
        Notes
        -----
        - Duplicate files are automatically filtered out
        - CSV headers are written only for new files
        - All images share the same output CSV file
        """
        logger.info("Loading images...")
        
        # Open file selection dialog
        file_paths, _ = QFileDialog.getOpenFileNames(
            None, 
            "Select Images", 
            "", 
            "Image Files (*.tif *.czi *.lsm)"
        )
        
        if not file_paths:
            QMessageBox.information(
                None, 
                "No Files Selected", 
                "Please select one or more image files."
            )
            return
        
        # Get output filename from user
        output_file_name, ok = QInputDialog.getText(
            None, 
            "Output File Name", 
            "Enter the name of the output CSV file:"
        )
        
        if not ok or not output_file_name.strip():
            QMessageBox.warning(
                None, 
                "No Output File", 
                "You must provide a name for the output file."
            )
            return
        
        # Ensure .csv extension
        if not output_file_name.endswith(".csv"):
            output_file_name += ".csv"
        
        # Determine output folder (same as first image)
        output_folder = os.path.dirname(file_paths[0]) if file_paths else os.getcwd()
        self.csv_file_path = os.path.join(output_folder, output_file_name)
        
        file_exists = os.path.isfile(self.csv_file_path)
        
        # Create CSV and write headers if new file
        try:
            with open(self.csv_file_path, mode='a', newline='', encoding='utf-8') as csv_file:
                writer = csv.writer(csv_file)
                if not file_exists:
                    writer.writerow([
                        "Image filename",
                        "Spread x [pixel]",
                        "Spread y [pixel]",
                        "Spread z [pixel]",
                        "Spread x*y [pixel2]",
                        "Spread x*y*z [pixel3]",
                        "Spread x [um]",
                        "Spread y [um]",
                        "Spread z [um]",
                        "Spread x*y [um2]",
                        "Spread x*y*z [um3]",                        
                        "Axonal Volume",
                        "Fluorescence_px",
                        "Fluorescence_um",
                        "Observation"
                    ])
                    logger.info(f"Created new CSV file: {self.csv_file_path}")
        except Exception as e:
            logger.error(f"Failed to create CSV file: {e}")
            QMessageBox.critical(
                self, 
                "Error Creating File", 
                f"An error occurred while creating the file:\n{e}"
            )
            return
        
        # Filter out duplicate files
        existing_files = [
            self.ui.listWidget_images.item(i).text() 
            for i in range(self.ui.listWidget_images.count())
        ]
        new_files = [fp for fp in file_paths if fp not in existing_files]
        
        if not new_files:
            QMessageBox.information(
                self, 
                "No New Files", 
                "All selected files are already loaded."
            )
            return
        
        # Add new files to list
        self.ui.listWidget_images.addItems(new_files)
        
        logger.info(f"Loaded {len(new_files)} new images")
        QMessageBox.information(
            None, 
            "Images Loaded", 
            f"Loaded {len(new_files)} images. Ready to process."
        )
    
    
    def _load_czi(self, file_path: str):
        """
        Load a Zeiss CZI microscopy file.
        
        Parameters
        ----------
        file_path : str
            Path to .czi file
        
        Notes
        -----
        - Extracts metadata (voxel size, dimensions, channels)
        - Applies zero-padding for non-square images
        - Transposes data to (Z, X, Y) format
        
        Raises
        ------
        ValueError
            If required metadata keys are missing
        RuntimeError
            If file cannot be loaded
        """
        try:
            logger.info(f"Loading CZI file: {file_path}")
            
            with czi.open_czi(file_path) as czifile:
                # Extract metadata
                metadata = czifile.metadata
                
                width_px = int(metadata['ImageDocument']['Metadata']['Information']['Image']['SizeX'])
                height_px = int(metadata['ImageDocument']['Metadata']['Information']['Image']['SizeY'])
                zlim = int(metadata['ImageDocument']['Metadata']['Information']['Image']['SizeZ'])
                
                voxel_size = metadata['ImageDocument']['Metadata']['Scaling']['Items']['Distance']
                voxel_size_x_um = float(next(
                    item for item in voxel_size if item['@Id'] == 'X'
                )['Value']) * 1e6
                voxel_size_y_um = float(next(
                    item for item in voxel_size if item['@Id'] == 'Y'
                )['Value']) * 1e6
                voxel_size_z_um = float(next(
                    item for item in voxel_size if item['@Id'] == 'Z'
                )['Value']) * 1e6
                
                num_channels = int(
                    metadata['ImageDocument']['Metadata']['Information']['Image'].get('SizeC', 1)
                )
                channels = [f"Channel {i+1}" for i in range(num_channels)] if num_channels > 1 else ["Channel 1"]
                
                logger.debug(f"CZI metadata: {width_px}x{height_px}x{zlim}, "
                           f"{num_channels} channels, voxel: {voxel_size_x_um:.3f}x"
                           f"{voxel_size_y_um:.3f}x{voxel_size_z_um:.3f} µm")
                
                # Read data for each channel
                data = []
                for i in range(num_channels):
                    channel_data = [
                        czifile.read(plane={"C": i, "Z": z})[:, :, 0]
                        for z in range(zlim)
                    ]
                    
                    # Apply zero-padding for non-square images
                    if width_px != height_px:
                        max_dim = max(width_px, height_px)
                        channel_data = np.array(channel_data)
                        padded_image = np.zeros(
                            (zlim, max_dim, max_dim), 
                            dtype=channel_data.dtype
                        )
                        y_start = (max_dim - height_px) // 2
                        x_start = (max_dim - width_px) // 2
                        padded_image[:, y_start:y_start+height_px, x_start:x_start+width_px] = channel_data
                        channel_data = padded_image
                        
                        logger.debug(f"Applied zero-padding: {width_px}x{height_px} -> {max_dim}x{max_dim}")
                    
                    # Stack and transpose to (Z, X, Y)
                    data.append(np.stack(channel_data, axis=0).transpose(0, 2, 1))
                
                # Update dimensions if padded
                if width_px != height_px:
                    width_px = max_dim
                    height_px = max_dim
                
                # Store loaded data
                self.current_image_data = data
                self.current_metadata_channel = channels
                self.current_metadata_voxel_size = (voxel_size_x_um, voxel_size_y_um, voxel_size_z_um)
                self.current_metadata_dimension = (width_px, height_px, zlim)
                
                logger.info(f"✓ CZI file loaded successfully")
                
        except KeyError as ke:
            logger.error(f"Metadata key not found: {ke}")
            raise ValueError(f"Metadata key not found: {ke}")
        except Exception as e:
            logger.error(f"Unexpected error loading CZI: {e}", exc_info=True)
            raise RuntimeError(f"Unexpected error while loading CZI file: {e}")
    
    
    def _load_tif(self, file_path: str):
        """
        Load a TIFF microscopy file.
        
        Parameters
        ----------
        file_path : str
            Path to .tif file
        
        Notes
        -----
        - Attempts to read ImageJ metadata for voxel sizes
        - Falls back to default values if metadata unavailable
        - Handles 2D, 3D, and 4D TIFF formats
        - Transposes data to (Z, X, Y) format
        
        Warnings
        --------
        If metadata is missing, default values are used:
        - voxel_size_x = 0.1 µm
        - voxel_size_y = 0.1 µm  
        - voxel_size_z = 1.0 µm
        """
        try:
            logger.info(f"Loading TIF file: {file_path}")
            
            with tifffile.TiffFile(file_path) as tiff:
                # Try to read ImageJ metadata
                try:
                    imagej_metadata = tiff.imagej_metadata
                    if imagej_metadata is None:
                        raise ValueError("Metadata is not available")
                    
                    num_channels = imagej_metadata.get('channels', 1)
                    voxel_size_z = float(imagej_metadata.get('spacing', 1.0))
                    
                    logger.debug(f"ImageJ metadata: {num_channels} channels, z-spacing: {voxel_size_z} µm")
                    
                except (AttributeError, ValueError):
                    # Use default values if metadata unavailable
                    logger.warning("No ImageJ metadata found. Using default values.")
                    QMessageBox.warning(
                        None, 
                        "Metadata Warning", 
                        "No metadata found. Using default values."
                    )
                    num_channels = 1
                    voxel_size_z = 1.0
                
                # Assume isotropic XY voxels
                voxel_size_x = 0.1
                voxel_size_y = voxel_size_x
                
                # Read image data
                full_image = tiff.asarray()
                dims = full_image.ndim
                
                # Handle different dimensionalities
                if dims == 2:
                    # Single slice, single channel
                    image_data = full_image[np.newaxis, :, :]
                    
                elif dims == 3:
                    # Either (Z, Y, X) or (C, Y, X)
                    if num_channels > 1 and full_image.shape[0] == num_channels:
                        # Multi-channel, single Z
                        image_data = [full_image[c, :, :][np.newaxis, :, :] for c in range(num_channels)]
                    else:
                        # Single channel, multiple Z
                        image_data = full_image
                        
                elif dims == 4:
                    # Multi-channel, multiple Z: (Z, C, Y, X)
                    image_data = [full_image[:, c, :, :] for c in range(num_channels)]
                    
                else:
                    raise ValueError(f"Unexpected image dimensions: {dims}")
                
                # Transpose to (Z, X, Y)
                if isinstance(image_data, list):
                    image_data = [img.transpose(0, 2, 1) for img in image_data]
                else:
                    image_data = image_data.transpose(0, 2, 1)
                
                # Get dimensions
                if num_channels == 1:
                    height, width = image_data.shape[1:3]
                    zlim = image_data.shape[0]
                    data = [image_data]
                else:
                    height, width = image_data[0].shape[1:3]
                    zlim = image_data[0].shape[0]
                    data = image_data
                
                # Store loaded data
                self.current_image_data = data
                self.current_metadata_channel = [f"Channel {i+1}" for i in range(num_channels)]
                self.current_metadata_voxel_size = (voxel_size_x, voxel_size_y, voxel_size_z)
                self.current_metadata_dimension = (width, height, zlim)
                
                logger.info(f"✓ TIF file loaded successfully: {width}x{height}x{zlim}")
                QMessageBox.information(
                    None, 
                    "TIF Loaded", 
                    f"Successfully loaded {file_path}."
                )
                
        except Exception as e:
            logger.error(f"Failed to load TIF: {e}", exc_info=True)
            QMessageBox.critical(
                None, 
                "Error Loading TIF", 
                f"Failed to load TIF file:\n{e}"
            )
    
    
    def _load_lsm(self, file_path: str):
        """
        Load a Zeiss LSM microscopy file.
        
        Parameters
        ----------
        file_path : str
            Path to .lsm file
        
        Notes
        -----
        - LSM files are a special type of TIFF with Zeiss metadata
        - Extracts voxel sizes from LSM metadata
        - Transposes data to (Z, X, Y) format
        """
        try:
            logger.info(f"Loading LSM file: {file_path}")
            
            with tifffile.TiffFile(file_path) as tif:
                # Extract LSM metadata
                metadata = tif.lsm_metadata
                
                num_channels = metadata['DimensionChannels']
                zlim = metadata['DimensionZ']
                width = metadata['DimensionX']
                height = metadata['DimensionY']
                
                # Convert voxel sizes to micrometers
                voxel_size_x = metadata['VoxelSizeX'] * 1e6
                voxel_size_y = metadata['VoxelSizeY'] * 1e6
                voxel_size_z = metadata['VoxelSizeZ'] * 1e6
                
                logger.debug(f"LSM metadata: {width}x{height}x{zlim}, {num_channels} channels, "
                           f"voxel: {voxel_size_x:.3f}x{voxel_size_y:.3f}x{voxel_size_z:.3f} µm")
                
                # Read image data
                if num_channels == 1:
                    # Single channel
                    image_data = tif.asarray().transpose((0, 2, 1))  # (Z, X, Y)
                else:
                    # Multi-channel: (Z, C, Y, X)
                    full_image = tif.asarray(series=0)
                    image_data = [
                        full_image[:, c, :, :].transpose((0, 2, 1)) 
                        for c in range(num_channels)
                    ]
                
                # Store loaded data
                self.current_image_data = image_data if num_channels > 1 else [image_data]
                self.current_metadata_channel = [f"Channel {i+1}" for i in range(num_channels)]
                self.current_metadata_voxel_size = (voxel_size_x, voxel_size_y, voxel_size_z)
                self.current_metadata_dimension = (width, height, zlim)
                
                logger.info(f"✓ LSM file loaded successfully")
                QMessageBox.information(
                    None, 
                    "LSM Loaded", 
                    f"Successfully loaded {file_path}."
                )
                
        except Exception as e:
            logger.error(f"Failed to load LSM: {e}", exc_info=True)
            QMessageBox.critical(
                None, 
                "Error Loading LSM", 
                f"Failed to load LSM file:\n{e}"
            )
    
    
    # ====================================================================
    # UI MANAGEMENT METHODS
    # ====================================================================
    
    def clean_list(self):
        """
        Clear all loaded images and reset the application state.
        
        This method:
        - Clears the image list widget
        - Resets all image data variables
        - Removes any active ROI
        - Clears the display
        - Re-enables channel selectors
        """
        logger.info("Cleaning image list and resetting state")
        
        # Clear image list
        self.ui.listWidget_images.clear()
        
        # Reset data
        self.current_image_data = None
        self.original_image_data = None
        self.current_metadata_channel = None
        self.current_metadata_dimension = None
        self.current_metadata_voxel_size = None
        
        # Re-enable channel selectors
        self.ui.comboBox_plasticityChannel.setEnabled(True)
        self.ui.comboBox_fluoChannel.setEnabled(True)
        
        # Clear display
        self.ui.graphWidget.clear()
        
        # Cancel ROI creation if active
        if self.creating_roi:
            self.creating_roi = False
            self.temp_points = []
            if hasattr(self, "vertex_scatter"):
                self.ui.graphWidget.removeItem(self.vertex_scatter)
                self.vertex_scatter = None
            self.ui.pushButton_roi.setEnabled(True)
            self.ui.pushButton_roi.setStyleSheet("")
        
        # Remove existing ROI
        self.clear_polygonal_roi()
        if self.roi is not None:
            self.ui.graphWidget.removeItem(self.roi)
            self.roi = None
        
        logger.info("✓ Application state reset")
    
    
    def reset_image_data(self):
        """
        Reset current image data variables.
        
        Used when switching between images in the list.
        """
        self.current_image_filepath = None
        self.current_image_data = None
        self.original_image_data = None
        self.current_metadata_channel = None
        self.current_metadata_dimension = None
        self.current_metadata_voxel_size = None
        self.ui.graphWidget.clear()
    
    
    def on_image_selected(self):
        """
        Handle image selection from the list widget.
        
        This method:
        1. Resets previous image data
        2. Loads the selected image file
        3. Updates channel selectors
        4. Displays image metadata
        5. Shows initial visualization
        
        Notes
        -----
        - Maintains channel selection when switching between images of same format
        - Automatically detects file format and calls appropriate loader
        - Updates UI with image dimensions and voxel sizes
        """
        logger.info("Image selected from list")
        
        self.reset_image_data()
        
        # Get selected file
        selected_items = self.ui.listWidget_images.selectedItems()
        if not selected_items:
            return
        
        selected_file = selected_items[0].text()
        self.current_image_filepath = selected_file
        
        logger.info(f"Loading: {os.path.basename(selected_file)}")
        
        # Clear display and ROI
        self.ui.graphWidget.clear()
        
        # Cancel ROI creation if active
        if self.creating_roi:
            self.creating_roi = False
            self.temp_points = []
            if hasattr(self, "vertex_scatter"):
                self.ui.graphWidget.removeItem(self.vertex_scatter)
                self.vertex_scatter = None
            self.ui.pushButton_roi.setEnabled(True)
            self.ui.pushButton_roi.setStyleSheet("")
        
        self.clear_polygonal_roi()
        if self.roi is not None:
            self.ui.graphWidget.removeItem(self.roi)
            self.roi = None
        
        # Load image based on file extension
        try:
            if selected_file.endswith(".czi"):
                self._load_czi(selected_file)
            elif selected_file.endswith(".tif"):
                self._load_tif(selected_file)
            elif selected_file.endswith(".lsm"):
                self._load_lsm(selected_file)
            else:
                QMessageBox.warning(
                    None, 
                    "Unsupported Format", 
                    f"Unsupported file format: {selected_file}"
                )
                return
        except Exception as e:
            logger.error(f"Failed to load image: {e}", exc_info=True)
            QMessageBox.critical(
                None, 
                "Error Loading Image", 
                f"Error loading image:\n{e}"
            )
            return
        
        # Save previous channel selections
        prev_index_plasticity = self.ui.comboBox_plasticityChannel.currentIndex()
        prev_index_fluo = self.ui.comboBox_fluoChannel.currentIndex()
        prev_channel_count = self.ui.comboBox_plasticityChannel.count()
        
        # Update channel selectors
        self.ui.comboBox_plasticityChannel.setEnabled(True)
        self.ui.comboBox_fluoChannel.setEnabled(True)
        
        self.ui.comboBox_channel_selector.clear()
        self.ui.comboBox_channel_selector.addItems(self.current_metadata_channel)
        
        self.ui.comboBox_plasticityChannel.clear()
        self.ui.comboBox_plasticityChannel.addItems(self.current_metadata_channel)
        
        self.ui.comboBox_fluoChannel.clear()
        self.ui.comboBox_fluoChannel.addItem("No channel")
        self.ui.comboBox_fluoChannel.addItems(self.current_metadata_channel)
        
        # Restore previous selection if possible
        new_channel_count = len(self.current_metadata_channel)
        if prev_channel_count == new_channel_count:
            if 0 <= prev_index_plasticity < self.ui.comboBox_plasticityChannel.count():
                self.ui.comboBox_plasticityChannel.setCurrentIndex(prev_index_plasticity)
            if 0 <= prev_index_fluo + 1 < self.ui.comboBox_fluoChannel.count():
                self.ui.comboBox_fluoChannel.setCurrentIndex(prev_index_fluo)
        else:
            # Default selection for different channel count
            self.ui.comboBox_plasticityChannel.setCurrentIndex(0)
            self.ui.comboBox_fluoChannel.setCurrentIndex(0)
        
        # Update metadata display
        image_size_x, image_size_y, image_size_z = self.current_metadata_dimension
        voxel_size_x, voxel_size_y, voxel_size_z = self.current_metadata_voxel_size
        
        self.ui.lineEdit_image_size_X.setText(str(image_size_x))
        self.ui.lineEdit_image_size_Y.setText(str(image_size_y))
        self.ui.lineEdit_image_size_Z.setText(str(image_size_z))
        self.ui.lineEdit_pixel_size_X.setText(f"{voxel_size_x:.3f}")
        self.ui.lineEdit_pixel_size_Y.setText(f"{voxel_size_y:.3f}")
        self.ui.lineEdit_pixel_size_Z.setText(f"{voxel_size_z:.3f}")
        
        self.ui.spinBox_zmax.setValue(image_size_z - 1)
        
        # Display image
        if self.ui.radioButton_plotZproject.isChecked():
            # Show maximum Z projection
            max_projection = np.max(self.current_image_data[0], axis=0)
            self.ui.graphWidget.setImage(max_projection)
        elif self.ui.radioButton_plotStack.isChecked():
            # Show full stack
            self.ui.graphWidget.setImage(self.current_image_data[0])
        else:
            QMessageBox.warning(
                self, 
                "Display Mode", 
                "Please select either Z-projection or Stack view."
            )
        
        # Hide ROI button (not needed for this application)
        self.ui.graphWidget.ui.roiBtn.hide()
        
        logger.info(f"✓ Image loaded and displayed: {image_size_x}x{image_size_y}x{image_size_z}")
    
    
    def update_display(self):
        """
        Update the image display based on current channel and view mode selection.
        
        Called automatically when:
        - Channel selector changes
        - View mode radio button toggles (Z-projection vs Stack)
        
        View Modes
        ----------
        - Z-projection: Maximum intensity projection along Z axis
        - Stack: Full 3D stack (can scroll through slices)
        """
        # Check if image is loaded
        selected_items = self.ui.listWidget_images.selectedItems()
        if not selected_items:
            return
        
        selected_channel = self.ui.comboBox_channel_selector.currentIndex()
        
        if 0 <= selected_channel < len(self.current_image_data):
            channel_data = self.current_image_data[selected_channel]
            
            if self.ui.radioButton_plotZproject.isChecked():
                # Show maximum Z projection
                max_projection = np.max(channel_data, axis=0)
                self.ui.graphWidget.setImage(max_projection)
                logger.debug(f"Displaying Z-projection of channel {selected_channel}")
                
            elif self.ui.radioButton_plotStack.isChecked():
                # Show full stack
                self.ui.graphWidget.setImage(channel_data)
                logger.debug(f"Displaying full stack of channel {selected_channel}")
                
            else:
                QMessageBox.warning(
                    self, 
                    "Display Mode", 
                    "Please select either Z-projection or Stack view."
                )
    
    
    # ====================================================================
    # IMAGE FILTERING METHODS
    # ====================================================================
    
    def apply_selected_filter(self, channel: str):
        """
        Apply the selected filter to the specified channel.
        
        Parameters
        ----------
        channel : str
            Channel identifier: 'chP' (plasticity) or 'chF' (fluorescence)
        
        Available Filters
        -----------------
        - Threshold: Remove pixels below percentage of maximum intensity
        - Gaussian Blur: Smooth image with Gaussian kernel
        - Median Filter: Remove noise with median filter
        
        Notes
        -----
        - Original image is backed up before filtering (for undo)
        - User is prompted for filter parameters via dialog
        - Display is automatically updated after filtering
        """
        logger.info(f"Applying filter to channel: {channel}")
        
        # Get channel index and filter type
        if channel == 'chP':
            selected_filter = self.ui.comboBox_filter_type_chP.currentText()
            channel_index = self.ui.comboBox_plasticityChannel.currentIndex()
        elif channel == 'chF':
            selected_filter = self.ui.comboBox_filter_type_chF.currentText()
            channel_index = self.ui.comboBox_fluoChannel.currentIndex() - 1
        else:
            QMessageBox.warning(self, "Invalid Channel", "Invalid channel specified.")
            return
        
        # Validate channel index
        if not (0 <= channel_index < len(self.current_image_data)):
            QMessageBox.warning(
                self, 
                "Invalid Channel", 
                f"Invalid channel index: {channel_index}"
            )
            return
        
        # Backup original image (for undo)
        if self.original_image_data is None:
            self.original_image_data = {}
        if channel_index not in self.original_image_data:
            self.original_image_data[channel_index] = self.current_image_data[channel_index].copy()
            logger.debug(f"Backed up original image for channel {channel_index}")
        
        # Get image copy
        image = self.current_image_data[channel_index].copy()
        
        # Apply selected filter
        if selected_filter == "Threshold":
            threshold, ok = QInputDialog.getDouble(
                self, 
                "Threshold", 
                "Percentage of maximum (0-100):", 
                3.0, 0, 100, 1
            )
            if ok:
                max_val = np.max(image)
                image[image <= (threshold / 100.0) * max_val] = 0
                logger.info(f"Applied threshold: {threshold}% of max ({max_val:.2f})")
        
        elif selected_filter == "Gaussian Blur":
            sigma, ok = QInputDialog.getDouble(
                self, 
                "Gaussian Blur", 
                "Sigma:", 
                1.0, 0.1, 50.0, 1
            )
            if ok:
                image = gaussian(image, sigma=sigma, preserve_range=True)
                logger.info(f"Applied Gaussian blur: sigma={sigma}")
        
        elif selected_filter == "Median Filter":
            size, ok = QInputDialog.getInt(
                self, 
                "Median Filter", 
                "Kernel size (odd):", 
                3, 1, 99, 2
            )
            if ok:
                filtered_slices = [
                    median_filter(image[z], size=size) 
                    for z in range(image.shape[0])
                ]
                image = np.stack(filtered_slices, axis=0)
                logger.info(f"Applied median filter: size={size}")
        
        else:
            QMessageBox.warning(self, "Invalid Filter", "Please select a valid filter.")
            return
        
        # Save filtered image and update display
        self.current_image_data[channel_index] = image
        self.update_display()
        
        logger.info(f"✓ Filter applied successfully to channel {channel_index}")
    
    
    def undo_filter(self, channel: str):
        """
        Undo the last filter applied to the specified channel.
        
        Parameters
        ----------
        channel : str
            Channel identifier: 'chP' (plasticity) or 'chF' (fluorescence)
        
        Notes
        -----
        - Restores the original image from backup
        - Removes backup after restoring
        - Only one undo level is supported
        """
        logger.info(f"Undoing filter for channel: {channel}")
        
        # Get channel index
        if channel == 'chP':
            channel_index = self.ui.comboBox_plasticityChannel.currentIndex()
        elif channel == 'chF':
            channel_index = self.ui.comboBox_fluoChannel.currentIndex() - 1
        else:
            QMessageBox.warning(self, "Invalid Channel", "Invalid channel specified.")
            return
        
        # Check if backup exists
        if self.original_image_data is None or channel_index not in self.original_image_data:
            QMessageBox.information(
                self, 
                "No Backup", 
                "No original version saved for this channel."
            )
            return
        
        # Restore original image
        self.current_image_data[channel_index] = self.original_image_data[channel_index].copy()
        del self.original_image_data[channel_index]
        
        # Clean up backup dict if empty
        if not self.original_image_data:
            self.original_image_data = None
        
        self.update_display()
        
        logger.info(f"✓ Filter undone for channel {channel_index}")
    
    
    # ====================================================================
    # ROI SELECTION METHODS
    # ====================================================================
    
    def enable_polygonal_roi_creation(self):
        """
        Enable interactive polygonal ROI creation mode.
        
        In this mode:
        - Click to add vertices to polygon
        - Backspace/Delete/Space to remove last vertex
        - Escape to cancel ROI creation
        - Apply ROI button to finalize selection
        
        Notes
        -----
        - Vertices are displayed as red dots
        - Button turns green while in ROI creation mode
        - Mouse and keyboard events are captured
        """
        logger.info("Enabling polygonal ROI creation")
        
        self.creating_roi = True
        self.temp_points = []
        
        # Create scatter plot for vertices
        self.vertex_scatter = pyqtgraph.ScatterPlotItem(
            size=5, 
            pen=pyqtgraph.mkPen(None), 
            brush=pyqtgraph.mkBrush(255, 0, 0, 150)
        )
        self.ui.graphWidget.addItem(self.vertex_scatter)
        
        # Connect mouse click event
        self.ui.graphWidget.view.scene().sigMouseClicked.connect(self.add_vertex)
        
        # Override keyboard event handler
        self.ui.graphWidget.keyPressEvent = self.handle_key_press
        
        # Update button appearance
        self.ui.pushButton_roi.setStyleSheet("background-color: lightgreen;")
        self.ui.pushButton_roi.setEnabled(False)
        
        logger.info("✓ ROI creation mode active")
    
    
    def add_vertex(self, event):
        """
        Add a vertex to the polygon ROI on mouse click.
        
        Parameters
        ----------
        event : QGraphicsSceneMouseEvent
            Mouse click event
        
        Notes
        -----
        - Only active when creating_roi flag is True
        - Converts scene coordinates to image coordinates
        - Updates vertex visualization immediately
        """
        # Check if ROI creation is enabled
        if not getattr(self, 'creating_roi', False):
            return
        
        # Get click position in image coordinates
        pos = event.scenePos()
        img_pos = self.ui.graphWidget.view.mapSceneToView(pos)
        self.temp_points.append((img_pos.x(), img_pos.y()))
        
        # Update vertex visualization
        self.vertex_scatter.setData(
            [p[0] for p in self.temp_points], 
            [p[1] for p in self.temp_points]
        )
        
        logger.debug(f"Added vertex: ({img_pos.x():.1f}, {img_pos.y():.1f})")
    
    
    def handle_key_press(self, event):
        """
        Handle keyboard events during ROI creation.
        
        Parameters
        ----------
        event : QKeyEvent
            Keyboard event
        
        Key Bindings
        ------------
        - Backspace/Delete/Space/Left: Remove last vertex
        - Escape: Cancel ROI creation
        """
        key = event.key()
        
        # Remove last vertex
        if key in [Qt.Key_Backspace, Qt.Key_Delete, 
                   Qt.Key_Space, Qt.Key_Left]:
            if self.temp_points:
                removed_point = self.temp_points.pop()
                logger.debug(f"Removed vertex: {removed_point}")
                
                # Update visualization
                self.vertex_scatter.setData(
                    [p[0] for p in self.temp_points], 
                    [p[1] for p in self.temp_points]
                )
            else:
                logger.debug("No vertices to remove")
        
        # Cancel ROI creation
        elif key == Qt.Key_Escape:
            logger.info("ROI creation canceled")
            self.creating_roi = False
            self.temp_points = []
            
            if hasattr(self, "vertex_scatter"):
                self.ui.graphWidget.removeItem(self.vertex_scatter)
                self.vertex_scatter = None
            
            # Restore button
            self.ui.pushButton_roi.setEnabled(True)
            self.ui.pushButton_roi.setStyleSheet("")
    
    
    def clear_polygonal_roi(self):
        """
        Clear all temporary ROI data.
        
        Removes:
        - Vertex scatter plot
        - Temporary point list
        """
        if hasattr(self, 'vertex_scatter') and self.vertex_scatter is not None:
            self.ui.graphWidget.removeItem(self.vertex_scatter)
            self.vertex_scatter = None
        self.temp_points = []
    
    
    def apply_roi_mask(self):
        """
        Apply the defined ROI mask to the selected channels.
        
        This method:
        1. Finalizes the polygon ROI
        2. Calculates ROI area
        3. Creates binary mask
        4. Applies mask to plasticity channel
        5. Applies mask to fluorescence channel (if selected)
        6. Updates display
        7. Disables channel selectors (locked for processing)
        
        Notes
        -----
        - Minimum 3 vertices required for valid polygon
        - Mask is 2D, replicated across all Z slices
        - Masked regions are set to zero
        - Original images remain unchanged (copies are masked)
        
        Raises
        ------
        Warning
            If no ROI defined or less than 3 vertices
        """
        logger.info("Applying ROI mask...")
        
        # Validate ROI
        if self.roi is None:
            if not hasattr(self, 'temp_points') or len(self.temp_points) < 3:
                QMessageBox.warning(
                    self, 
                    "No ROI", 
                    "Please select an ROI or create one with at least 3 vertices."
                )
                return
            
            # Close polygon automatically
            self.temp_points.append(self.temp_points[0])
            self.roi = pyqtgraph.PolyLineROI(
                self.temp_points, 
                closed=True, 
                pen=pyqtgraph.mkPen('g', width=2)
            )
            self.ui.graphWidget.addItem(self.roi)
        
        # Deactivate ROI creation mode
        self.creating_roi = False
        self.ui.pushButton_roi.setEnabled(True)
        self.ui.pushButton_roi.setStyleSheet("")
        
        # Get ROI vertices
        roi_positions = self.roi.getLocalHandlePositions()
        vertices = [(pos.x(), pos.y()) for name, pos in roi_positions]
        
        # Calculate area and store vertices
        polygon = np.array(vertices + [vertices[0]])  # Close polygon
        self.AArea = Polygon(polygon).area
        self.polygon_points = vertices
        
        logger.info(f"ROI area: {self.AArea:.2f} pixels²")
        
        # Get selected channels
        plasticity_channel_index = self.ui.comboBox_plasticityChannel.currentIndex()
        self.ui.comboBox_plasticityChannel.setEnabled(False)  # Lock selection
        
        fluo_channel_index = self.ui.comboBox_fluoChannel.currentIndex()
        self.ui.comboBox_fluoChannel.setEnabled(False)  # Lock selection
        
        # Validate plasticity channel
        if 0 <= plasticity_channel_index < len(self.current_image_data):
            plasticity_image = self.current_image_data[plasticity_channel_index]
        else:
            QMessageBox.warning(
                self, 
                "Invalid Channel", 
                "Invalid plasticity channel selected."
            )
            return
        
        # Create binary mask from ROI
        slices, height, width = plasticity_image.shape
        
        # Create coordinate grid
        y, x = np.meshgrid(np.arange(width), np.arange(height), indexing='xy')
        points = np.vstack((x.ravel(), y.ravel())).T
        
        # Check which points are inside polygon
        path = MPLPath(polygon)
        mask = path.contains_points(points).reshape(height, width)
        
        # Expand mask to 3D (replicate for all Z slices)
        expanded_mask = np.stack([mask] * slices, axis=0)
        
        logger.info(f"ROI mask created: {np.sum(mask)} pixels, "
                   f"{slices} slices = {np.sum(expanded_mask)} total voxels")
        
        # Apply mask to plasticity channel
        masked_plasticity_image = np.copy(plasticity_image)
        masked_plasticity_image[~expanded_mask] = 0
        self.complexity_channel = masked_plasticity_image
        
        # Update display
        self.ui.graphWidget.setImage(masked_plasticity_image, autoLevels=True)
        
        # Apply mask to fluorescence channel (if selected)
        if fluo_channel_index == 0:
            # "No channel" selected
            fluo_image = None
        elif 0 < fluo_channel_index <= len(self.current_image_data):
            fluo_image = self.current_image_data[fluo_channel_index - 1]
        else:
            QMessageBox.warning(
                self, 
                "Invalid Channel", 
                "Invalid fluorescence channel selected."
            )
            return
        
        if fluo_image is not None:
            masked_fluo_image = np.copy(fluo_image)
            masked_fluo_image[~expanded_mask] = 0
            self.fluor_channel = masked_fluo_image
            logger.info("Fluorescence channel masked")
        else:
            self.fluor_channel = None
            logger.info("No fluorescence channel selected")
        
        logger.info(f"✓ ROI mask applied - Complexity channel: "
                   f"{np.count_nonzero(self.complexity_channel)} non-zero voxels, "
                   f"intensity: {np.sum(self.complexity_channel):.2f}")
    
    
    # ====================================================================
    # PROCESSING METHODS
    # ====================================================================
    
    def process(self):
        """
        Process the selected image to calculate structural plasticity metrics.
        
        This is the main processing pipeline that:
        1. Validates parameters and ROI
        2. Extracts Z-range subset
        3. Applies ROI mask
        4. Initializes ImageProcessor
        5. Calculates PCA rotation and spreads
        6. Calculates fluorescence (if second channel)
        7. Plots distributions (optional)
        8. Saves results to CSV
        9. Updates UI
        
        Workflow Steps
        --------------
        1. Parameter validation
        2. Progress dialog creation
        3. Image preparation (ROI masking, Z-subset)
        4. Processor initialization
        5. PCA and spread calculation
        6. Fluorescence calculation (optional)
        7. Distribution plotting (optional)
        8. Results export to CSV
        9. UI update and success message
        
        Notes
        -----
        - Shows progress dialog during processing
        - Can be canceled by user at any time
        - Errors are logged and displayed to user
        - Processed images are marked in the list (green=success, red=error)
        
        See Also
        --------
        ImageProcessor.process_image : Core processing algorithm
        _validate_processing_parameters : Parameter validation
        _prepare_masked_images : Image preparation
        _save_results_to_csv : Results export
        """
        try:
            logger.info("\n" + "="*70)
            logger.info("STARTING PROCESSING")
            logger.info("="*70)
            
            # ============================================================
            # STEP 1: VALIDATE PARAMETERS
            # ============================================================
            if not self._validate_processing_parameters():
                return
            
            # ============================================================
            # STEP 2: GET PARAMETERS
            # ============================================================
            try:
                voxel_size_x = float(self.ui.lineEdit_pixel_size_X.text())
                voxel_size_y = float(self.ui.lineEdit_pixel_size_Y.text())
                voxel_size_z = float(self.ui.lineEdit_pixel_size_Z.text())
                z_start = self.ui.spinBox_zmin.value()
                z_end = self.ui.spinBox_zmax.value()
            except ValueError as e:
                QMessageBox.critical(
                    self, 
                    "Invalid Parameters", 
                    f"Please check your numeric inputs:\n{str(e)}"
                )
                return
            
            logger.info(f"Parameters: voxel_x={voxel_size_x:.3f}µm, "
                       f"voxel_y={voxel_size_y:.3f}µm, voxel_z={voxel_size_z:.3f}µm, "
                       f"z_range=[{z_start}, {z_end}]")
            
            # Validate voxel sizes
            is_valid, error_msg = validate_parameters(
                voxel_size_x, voxel_size_y, voxel_size_z
            )
            if not is_valid:
                QMessageBox.warning(self, "Invalid Voxel Sizes", error_msg)
                return
            
            # ============================================================
            # STEP 3: CREATE PROGRESS DIALOG
            # ============================================================
            progress = QProgressDialog(
                "Processing image...", 
                "Cancel", 
                0, 100, 
                self
            )
            progress.setWindowModality(Qt.WindowModal)
            progress.setWindowTitle("Structural Plasticity Analysis")
            progress.setMinimumDuration(0)  # Show immediately
            progress.setValue(5)
            
            # ============================================================
            # STEP 4: PREPARE MASKED IMAGES
            # ============================================================
            progress.setLabelText("Applying ROI mask...")
            progress.setValue(15)
            
            masked_image_complexity, masked_image_fluor = self._prepare_masked_images(
                z_start, z_end
            )
            
            if masked_image_complexity is None:
                progress.close()
                return
            
            if progress.wasCanceled():
                logger.info("Processing canceled by user")
                progress.close()
                return
            
            progress.setValue(30)
            
            # ============================================================
            # STEP 5: INITIALIZE PROCESSOR
            # ============================================================
            progress.setLabelText("Initializing processor...")
            
            processor = ImageProcessor(
                voxel_size_x=voxel_size_x,
                voxel_size_y=voxel_size_y,
                voxel_size_z=voxel_size_z
            )
            
            progress.setValue(40)
            
            # ============================================================
            # STEP 6: PROCESS IMAGE (PCA + SPREADS)
            # ============================================================
            progress.setLabelText("Calculating PCA and spreads...")
            
            logger.info(f"Processing complexity channel. Shape: {masked_image_complexity.shape}")
            
            results = processor.process_image(
                image_3d=masked_image_complexity,
                mask_area_pixels=self.AArea
            )
            
            if progress.wasCanceled():
                logger.info("Processing canceled by user")
                progress.close()
                return
            
            progress.setValue(70)
            
            # ============================================================
            # STEP 7: CALCULATE FLUORESCENCE (if second channel exists)
            # ============================================================
            if masked_image_fluor is not None:
                progress.setLabelText("Calculating fluorescence...")
                logger.info("Calculating fluorescence for second channel")
                
                fluor_px, fluor_um = processor.calculate_fluorescence(
                    masked_image_fluor, 
                    self.AArea
                )
                
                results['fluorescence_px'] = fluor_px
                results['fluorescence_um'] = fluor_um
            
            progress.setValue(80)
            
            # ============================================================
            # STEP 8: PLOT DISTRIBUTIONS (if enabled)
            # ============================================================
            if self.ui.checkBox_show_distributions.isChecked():
                self._plot_distributions(
                    results['MMsum'], 
                    results['MMyy'], 
                    results['MMzz']
                )
            
            # ============================================================
            # STEP 9: SAVE RESULTS
            # ============================================================
            progress.setLabelText("Saving results...")
            progress.setValue(90)
            
            self._save_results_to_csv(results)
            
            # ============================================================
            # STEP 10: UPDATE UI
            # ============================================================
            self._update_ui_after_processing(success=True)
            
            progress.setValue(100)
            progress.close()
            
            # ============================================================
            # STEP 11: SHOW SUCCESS MESSAGE
            # ============================================================
            QMessageBox.information(
                self,
                "Processing Complete",
                f"Image processed successfully!\n\n"
                f"3D Spread: {results['spread_xyz_um']:.2f} µm³\n"
                f"Axonal Volume: {results['axonal_volume']:.2f}\n"
                f"Rotation Angle: {results['rotation_angle']:.2f}°"
            )
            
            logger.info("="*70)
            logger.info("PROCESSING COMPLETED SUCCESSFULLY")
            logger.info("="*70 + "\n")
            
        except Exception as e:
            logger.error(f"Processing error: {e}", exc_info=True)
            
            if 'progress' in locals():
                progress.close()
            
            QMessageBox.critical(
                self,
                "Processing Error",
                f"An error occurred during processing:\n\n{str(e)}\n\n"
                f"Check the log file (plasticity_analyzer.log) for details."
            )
            
            self._update_ui_after_processing(success=False)
    
    
    def _validate_processing_parameters(self) -> bool:
        """
        Validate all parameters before processing.
        
        Checks:
        - Image is loaded
        - ROI is defined (minimum 3 vertices)
        - ROI area is positive
        - Z range is valid
        - Voxel sizes are positive
        
        Returns
        -------
        bool
            True if all parameters are valid, False otherwise
        
        Notes
        -----
        Displays warning dialogs for validation failures
        """
        # Check if image is loaded
        if not hasattr(self, 'complexity_channel') or self.complexity_channel is None:
            QMessageBox.warning(self, "No Image", "Please load an image first.")
            return False
        
        # Check if ROI is defined
        if not hasattr(self, 'polygon_points') or len(self.polygon_points) < 3:
            QMessageBox.warning(
                self, 
                "No ROI", 
                "Please define a Region of Interest (ROI) first.\n"
                "Click 'Select ROI' and draw a polygon around the structure."
            )
            return False
        
        # Check if ROI area was calculated
        if not hasattr(self, 'AArea') or self.AArea <= 0:
            QMessageBox.warning(
                self, 
                "Invalid ROI", 
                "ROI area is zero. Please redraw the ROI."
            )
            return False
        
        # Validate Z range
        z_start = self.ui.spinBox_zmin.value()
        z_end = self.ui.spinBox_zmax.value()
        
        if z_start >= z_end:
            QMessageBox.warning(
                self, 
                "Invalid Z Range", 
                f"Z start ({z_start}) must be less than Z end ({z_end})."
            )
            return False
        
        if z_end > self.complexity_channel.shape[0]:
            QMessageBox.warning(
                self, 
                "Invalid Z Range",
                f"Z end ({z_end}) exceeds image depth "
                f"({self.complexity_channel.shape[0]})."
            )
            return False
        
        # Validate voxel sizes
        try:
            voxel_x = float(self.ui.lineEdit_pixel_size_X.text())
            voxel_y = float(self.ui.lineEdit_pixel_size_Y.text())
            voxel_z = float(self.ui.lineEdit_pixel_size_Z.text())
            
            if voxel_x <= 0 or voxel_y <= 0 or voxel_z <= 0:
                QMessageBox.warning(
                    self, 
                    "Invalid Voxel Size",
                    "All voxel sizes must be positive numbers."
                )
                return False
                
        except ValueError:
            QMessageBox.warning(
                self, 
                "Invalid Voxel Size",
                "Please enter valid numeric values for voxel sizes."
            )
            return False
        
        logger.info("✓ All parameters validated")
        return True
    
    
    def _prepare_masked_images(
        self, 
        z_start: int, 
        z_end: int
    ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """
        Prepare masked images for complexity and fluorescence channels.
        
        This method:
        1. Extracts Z-slice subset
        2. Returns already-masked images (mask was applied in apply_roi_mask)
        
        Parameters
        ----------
        z_start : int
            Starting Z slice (inclusive)
        z_end : int
            Ending Z slice (inclusive)
        
        Returns
        -------
        Tuple[Optional[np.ndarray], Optional[np.ndarray]]
            (masked_complexity, masked_fluor)
            - masked_complexity: Masked plasticity channel
            - masked_fluor: Masked fluorescence channel (None if not selected)
            Returns (None, None) if preparation fails
        
        Notes
        -----
        - Images are already masked (apply_roi_mask was called earlier)
        - Only Z-range extraction is performed here
        - Logs intensity statistics for verification
        """
        try:
            logger.info("Preparing masked images...")
            
            # Extract Z-range from complexity channel
            complexity_subset = self.complexity_channel[z_start:z_end+1, :, :]
            logger.info(f"Complexity channel shape: {complexity_subset.shape}")
            logger.info(f"Intensity: {np.sum(complexity_subset):.2f}, "
                       f"Non-zero: {np.count_nonzero(complexity_subset)}")
            
            # Extract Z-range from fluorescence channel (if exists)
            fluor_subset = None
            if self.fluor_channel is not None:
                fluor_subset = self.fluor_channel[z_start:z_end+1, :, :]
                logger.info(f"Fluorescence channel shape: {fluor_subset.shape}")
            
            return complexity_subset, fluor_subset
            
        except Exception as e:
            logger.error(f"Error preparing masked images: {e}", exc_info=True)
            QMessageBox.critical(
                self, 
                "Mask Error",
                f"Failed to prepare images:\n{str(e)}"
            )
            return None, None
    
    
    def _plot_distributions(
        self, 
        MMsum: np.ndarray, 
        MMyy: np.ndarray, 
        MMzz: np.ndarray
    ):
        """
        Plot X, Y, Z distributions for visual inspection.
        
        Creates a 3-panel figure showing:
        - Intensity distribution along X (horizontal)
        - Y-spread at each X position
        - Z-spread at each X position
        
        Parameters
        ----------
        MMsum : np.ndarray
            Intensity sum at each X position
        MMyy : np.ndarray
            Local Y-variance at each X position
        MMzz : np.ndarray
            Local Z-variance at each X position
        
        Notes
        -----
        - Plots are shown in a matplotlib window
        - User must close window to continue
        - Useful for quality control and troubleshooting
        """
        try:
            fig, axs = plt.subplots(1, 3, figsize=(15, 4))
            
            # X distribution (intensity)
            axs[0].plot(MMsum, linewidth=2, color='#2E86AB')
            axs[0].set_title(
                'Intensity Along X (Horizontal)', 
                fontsize=14, 
                fontweight='bold'
            )
            axs[0].set_xlabel('X Position (pixels)', fontsize=12)
            axs[0].set_ylabel('Intensity Sum', fontsize=12)
            axs[0].grid(True, alpha=0.3)
            axs[0].set_facecolor('#F8F9FA')
            
            # Y distribution (variance)
            axs[1].plot(MMyy, linewidth=2, color='#F18F01')
            axs[1].set_title('Y-Spread at Each X', fontsize=14, fontweight='bold')
            axs[1].set_xlabel('X Position (pixels)', fontsize=12)
            axs[1].set_ylabel('Variance (pixels²)', fontsize=12)
            axs[1].grid(True, alpha=0.3)
            axs[1].set_facecolor('#F8F9FA')
            
            # Z distribution (variance)
            axs[2].plot(MMzz, linewidth=2, color='#06A77D')
            axs[2].set_title('Z-Spread at Each X', fontsize=14, fontweight='bold')
            axs[2].set_xlabel('X Position (pixels)', fontsize=12)
            axs[2].set_ylabel('Variance (slices²)', fontsize=12)
            axs[2].grid(True, alpha=0.3)
            axs[2].set_facecolor('#F8F9FA')
            
            plt.suptitle('Distribution Analysis', fontsize=16, fontweight='bold', y=1.02)
            plt.tight_layout()
            plt.show()
            
            logger.info("Distribution plots displayed")
            
        except Exception as e:
            logger.warning(f"Failed to plot distributions: {e}")
    
    
    def _save_results_to_csv(self, results: dict):
        """
        Save processing results to CSV file.
        
        Appends a new row with:
        - Image filename
        - Spread metrics (pixels and µm)
        - Axonal volume
        - Fluorescence values
        - User observation/notes
        
        Parameters
        ----------
        results : dict
            Processing results from ImageProcessor
        
        Notes
        -----
        - CSV file path was set during load_images()
        - Observation text is taken from UI text edit widget
        - Results are appended to existing CSV (does not overwrite)
        """
        try:
            if not hasattr(self, 'csv_file_path') or not self.csv_file_path:
                logger.error("No CSV file path defined")
                QMessageBox.warning(
                    self, 
                    "No Output File",
                    "Output CSV file is not defined."
                )
                return
            
            observation = self.ui.textEdit_observation.toPlainText()
            image_name = os.path.basename(self.current_image_filepath)
            
            # Write to CSV
            with open(self.csv_file_path, mode='a', newline='', encoding='utf-8') as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow([
                    image_name,
                    results['spread_x_pixel'],
                    results['spread_y_pixel'],
                    results['spread_z_pixel'],
                    results['spread_xy_pixel'],
                    results['spread_xyz_pixel'],
                    results['spread_x_um'],
                    results['spread_y_um'],
                    results['spread_z_um'],
                    results['spread_xy_um'],
                    results['spread_xyz_um'],
                    results['axonal_volume'],
                    results['fluorescence_px'],
                    results['fluorescence_um'],
                    f'"{observation}"'
                ])
            
            logger.info(f"✓ Results saved to: {self.csv_file_path}")
            
        except Exception as e:
            logger.error(f"Error saving results: {e}", exc_info=True)
            QMessageBox.critical(
                self, 
                "Save Error",
                f"Failed to save results to CSV:\n{str(e)}"
            )
    
    
    def _update_ui_after_processing(self, success: bool = True):
        """
        Update UI elements after processing.
        
        Parameters
        ----------
        success : bool, optional
            Whether processing was successful (default: True)
        
        UI Updates
        ----------
        - Mark processed image in list:
            - Green background: Success
            - Red background: Error
        - Update status bar (if available)
        
        Notes
        -----
        Errors in UI update are logged but do not interrupt workflow
        """
        try:
            # Mark processed image in list
            selected_items = self.ui.listWidget_images.selectedItems()
            if selected_items:
                if success:
                    selected_items[0].setBackground(pyqtgraph.mkColor('lightgreen'))
                else:
                    selected_items[0].setBackground(pyqtgraph.mkColor("#FB889A"))
            
            # Update status bar if exists
            if hasattr(self, 'statusBar'):
                if success:
                    self.statusBar().showMessage("✓ Processing completed", 5000)
                else:
                    self.statusBar().showMessage("✗ Processing failed", 5000)
            
        except Exception as e:
            logger.warning(f"Failed to update UI: {e}")


# ========================================================================
# APPLICATION ENTRY POINT
# ========================================================================

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MyMainWindow()
    window.show()
    sys.exit(app.exec())