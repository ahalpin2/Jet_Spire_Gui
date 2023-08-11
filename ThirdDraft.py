
from IPython.core.display import display, HTML
#display(HTML("<style>.container { width:95% !important; }</style>"))
import importlib
import traceback
from scipy.interpolate import interp1d
from sys import path
import matplotlib.pyplot as plt
import sunpy.map
import sunpy.data.sample
import astropy.units as u
import numpy as np
import base64
import os
import cv2
import imageio as imageio
import contextlib
from astropy.coordinates import SkyCoord
from astropy.time import Time
from sunpy.net import Fido, attrs as a, vso
from sunpy.time import TimeRange
from ndcube import NDCube, NDCubeSequence, NDCollection 
from sunpy.coordinates.utils import GreatArc
from sunpy.map import Map
from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QPushButton, QWidget, QLabel, QFileDialog, QMessageBox, QComboBox
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QInputDialog, QLineEdit
from PyQt5.QtGui import QImage, QPixmap
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from PIL import Image
from matplotlib import patches
from matplotlib.cm import get_cmap
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1 import make_axes_locatable
from datetime import datetime
from sunpy.coordinates import Helioprojective
###algorythm files 
################################################################### 
## add in custom errors for a few common errors (directory, invalid region, dimensions is bad, time format, string format, plus more that come to you)
class HoverableLabel(QtWidgets.QLabel):
    mouseHoverEvent = QtCore.pyqtSignal(SkyCoord)
    def __init__(self, parent=None):
        super(HoverableLabel, self).__init__(parent)
class MainWindow(QtWidgets.QMainWindow):
    first_click = None
    def __init__(self):
        super(MainWindow, self).__init__()
        self.data_map = None
        self.fig, self.ax = plt.subplots()
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.label_image = HoverableLabel()
        self.setCentralWidget(self.label_image)
#setting the central window up
        self.setWindowTitle("My Window")
        self.resize(900, 600)
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        self.layout().setContentsMargins(0, 0, 0, 0)
# Create the buttons as instance attributes
        self.button1 = QPushButton("AIA Search", self)
        self.button2 = QPushButton("Image Creator", self)
        self.button3 = QPushButton("Movie Creator", self)
        self.button4 = QPushButton("Light Curve Creator", self)
        self.button5 = QPushButton("Select Input_Directory", self)
        self.button6 = QPushButton("Select Output_Directory", self)
        self.button7 = QPushButton("Velocity Calculator", self)
        self.button8 = QPushButton("Velocity Group", self)
#        self.button9 = QPushButton("EmToolKit", self)
        self.button10 = QPushButton("HMI search", self)
        self.button11 = QPushButton("EMT search", self)
        self.button7_load_image = QPushButton("Load Image", self)
        self.button7_load_image_2 = QPushButton("Load Image 2", self)       
# Set the positions of the buttons using move() method
        self.button1.setGeometry(10, 100, 180, 30)
        self.button2.setGeometry(10, 325, 140, 30)
        self.button3.setGeometry(10, 275, 140, 30)
        self.button4.setGeometry(10, 400, 180, 30)
        self.button5.setGeometry(10, 0, 180, 30)
        self.button6.setGeometry(10, 50, 180, 30)
        self.button7.setGeometry(10, 450, 180, 30)
        self.button8.setGeometry(10, 500, 180, 30)
#        self.button9.setGeometry(10, 550, 180, 30)
        self.button10.setGeometry(10, 200, 180, 30)
        self.button11.setGeometry(10, 150, 180, 30)
        self.button7_load_image.setGeometry(155, 275, 125, 30)
        self.button7_load_image_2.setGeometry(155, 325, 125, 30)
# Connect the clicked signal of each button to the corresponding function
        self.button1.clicked.connect(self.search)
        self.button2.clicked.connect(self.generate_images)
        self.button3.clicked.connect(self.create_gif)
        self.button4.clicked.connect(self.create_plot)
        self.button5.clicked.connect(self.select_input_directory)
        self.button6.clicked.connect(self.select_output_directory)
        self.button7.clicked.connect(self.Velocity_Calculator)
        self.button8.clicked.connect(self.Velocity_Group_Calculator)
#        self.button9.clicked.connect(self.emtoolkit)
        self.button10.clicked.connect(self.hmi_search)
        self.button11.clicked.connect(self.emt_search)
        self.button7_load_image.clicked.connect(self.load_image)
        self.button7_load_image_2.clicked.connect(self.load_image_2)
#constructing the line edits
        self.label1 = QLabel("Start Time:", self)
        self.label2 = QLabel("End Time:", self)
        self.label3 = QLabel("Wavelength:", self)
        self.label27 = QLabel("JSOC Email:", self)
        self.label4 = QLabel("Window Top X:", self)
        self.label7 = QLabel("Window Top Y:", self)
        self.label8 = QLabel("Window Bottom X:", self)
        self.label9 = QLabel("Window Bottom Y:", self)
        self.label10 = QLabel("Arc Start X :", self)
        self.label11 = QLabel("Arc Start Y:", self)
        self.label12 = QLabel("Arc End X:", self)
        self.label13 = QLabel("Arc End Y:", self)
        self.label5 = QLabel(self)
        self.label6 = QLabel(self)
        self.label16 = QLabel("Arc Start X 2:", self)
        self.label17 = QLabel("Arc Start Y 2:", self)
        self.label18 = QLabel("Time 1:", self)
        self.label19 = QLabel("Time 2:", self)
        self.label24 = QLabel("Velocity (km/s):", self)
        self.label20 = QLabel("Example: YYYY-MM-DD HH:MM:SS", self)
        self.label21 = QLabel("Example: YYYY-MM-DD HH:MM:SS", self)
        self.label22 = QLabel("Example: 960 ", self)
        self.label23 = QLabel("Range: (-1100, 1100) ", self)         
#geometry
        self.label1.setGeometry(210, 100, 230, 30)
        self.label2.setGeometry(210, 150, 230, 30)
        self.label3.setGeometry(210, 200, 130, 30)
        self.label27.setGeometry(410, 150, 100, 130)
        self.label4.setGeometry(285, 275, 130, 30)
        self.label7.setGeometry(285, 325, 130, 30)
        self.label8.setGeometry(480, 275, 130, 30)
        self.label9.setGeometry(480, 325, 130, 30)
        self.label10.setGeometry(210, 400, 130, 30)
        self.label11.setGeometry(210, 450, 130, 30)
        self.label12.setGeometry(210, 500, 130, 30)
        self.label13.setGeometry(400, 500, 130, 30)
        self.label5.setGeometry(200, 0, 250, 30)
        self.label6.setGeometry(200, 50, 400, 30)
        self.label16.setGeometry(400, 400, 95, 30)
        self.label17.setGeometry(400, 450, 95, 30)
        self.label18.setGeometry(585, 400, 75, 30)
        self.label19.setGeometry(585, 450, 75, 30)
        self.label24.setGeometry(585, 500, 105, 30)
        self.label20.setGeometry(460, 100, 250, 30)
        self.label21.setGeometry(460, 150, 250, 30)
        self.label22.setGeometry(690, 275, 250, 30)
        self.label23.setGeometry(690, 325, 250, 30)
#line edits
        self.line_edit1 = QLineEdit(self)
        self.line_edit2 = QLineEdit(self)
        self.line_edit3 = QLineEdit(self)
        self.line_edit4 = QLineEdit(self)
        self.line_edit5 = QLineEdit(self)
        self.line_edit6 = QLineEdit(self)
        self.line_edit7 = QLineEdit(self)
        self.line_edit8 = QLineEdit(self)
        self.line_edit9 = QLineEdit(self)
        self.line_edit10 = QLineEdit(self)
        self.line_edit11 = QLineEdit(self)
        self.line_edit12 = QLineEdit(self)
        self.line_edit13 = QLineEdit(self)
        self.line_edit14 = QLineEdit(self)
        self.line_edit15 = QLineEdit(self)
        self.line_edit17 = QLineEdit(self)
        self.line_edit18 = QLineEdit(self)
        self.line_edit14.setReadOnly(True)
        self.line_edit15.setReadOnly(True)
        self.line_edit18.setReadOnly(True)
#line edit geometry
        self.line_edit1.setGeometry(300, 100, 150, 30) #start time
        self.line_edit2.setGeometry(300, 150, 150, 30) #end time
        self.line_edit3.setGeometry(200, 200, 125, 30) #wavelength
        self.line_edit4.setGeometry(390, 275, 75, 30) #Image window
        self.line_edit5.setGeometry(390, 325, 75, 30) #Image window
        self.line_edit6.setGeometry(610, 275, 75, 30) #Image window
        self.line_edit7.setGeometry(610, 325, 75, 30) #Image window
        self.line_edit17.setGeometry(500, 200, 125, 30) #  HMI jsoc email
        self.line_edit8.setGeometry(300, 400, 75, 30) #jet start 1
        self.line_edit9.setGeometry(300, 450, 75, 30) #jet start 1
        self.line_edit10.setGeometry(300, 500, 75, 30) #jet end  
        self.line_edit11.setGeometry(500, 500, 75, 30) #jet end
        self.line_edit12.setGeometry(500, 400, 75, 30) # jet top 2
        self.line_edit13.setGeometry(500, 450, 75, 30) # jet top 2
        self.line_edit14.setGeometry(645, 400, 200, 30) #jet time 1
        self.line_edit15.setGeometry(645, 450, 200, 30) # jet time 2
        self.line_edit18.setGeometry(695, 500, 150, 30) # calculated velocity
#setting up the image display for the map
        self.label_image= QLabel(self)
#setting up a combo box for the search features
        self.combo_box2 = QComboBox(self)
        self.combo_box2.setGeometry(300,200,100,30)
        Wavelength_list = ["171","131","193","94","211","304","335",]
        self.combo_box2.addItems(Wavelength_list)
        self.combo_box2.setLineEdit(self.line_edit3)
#setting up input and out put directory stuff for all future buttons to use
        self.input_directory = ""
        self.output_directory = ""
###setting up the input and out folder selects for the following window
    def select_input_directory(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        directory = QFileDialog.getExistingDirectory(self, "Select Input Directory", options=options)
        if directory:
            self.input_directory = directory
            self.label5.setText("Input Directory: " + self.input_directory)
            self.label5.adjustSize()
    def select_output_directory(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        directory = QFileDialog.getExistingDirectory(self, "Select Output Directory", options=options)
        if directory:
            self.output_directory = directory
            self.label6.setText("Output Directory: " + self.output_directory)
            self.label6.adjustSize()
#search button script
    def search(self):
      try:
        start_time = self.line_edit1.text() 
        end_time = self.line_edit2.text()
        wavelength = float(self.line_edit3.text()) * u.angstrom
        instrument = 'AIA'
        search_results = Fido.search(a.Time(start_time, end_time), a.Instrument(instrument), a.Wavelength(wavelength))
        downloaded_files = Fido.fetch(search_results, path=self.output_directory, verify=False)
        print(downloaded_files.errors)  
      except Exception as e:
        series = a.jsoc.Series('hmi.synoptic_mr_polfil_720s')
        start_time
        error_message = str(e)
        QMessageBox.critical(None, "Error", "An error occurred: {}".format(error_message)) 
#filtered_searches for the EmToolKit
    def hmi_search(self):
      try:
        start_time = self.line_edit1.text() 
        end_time = self.line_edit2.text()  
        JSOC_Email = self.line_edit17.text()
        print(JSOC_Email)
        series = ('hmi.M_45s')        
        search_results = Fido.search(a.Time(start_time, end_time), a.jsoc.Series(series), a.jsoc.Notify(JSOC_Email))
        downloaded_files = Fido.fetch(search_results, path=self.output_directory, verify=False)
      except Exception as e:
        error_message = str(e)
        QMessageBox.critical(None, "Error", "An error occurred: {}".format(error_message))
#EMT_Search
    def emt_search(self):
      try:         
        date= self.line_edit1.text() 
# Commands for initial data download. Comment out once that's successful.
# VSO can sometimes be a bit flakey here, in my experience, may require multiple tries:
        paths = []
        passbands = np.array([94,131,171,193,211,335])*u.angstrom
        for band in passbands: 
            qry = Fido.search(a.Time(TimeRange(date,12*u.s)),a.Instrument('AIA'),a.Wavelength(band))[0,0]
            paths.append(Fido.fetch(qry,path=self.output_directory, verify=False))
      except Exception as e:
        error_message = str(e)
        QMessageBox.critical(None, "Error", "An error occurred: {}".format(error_message))
#Generated buttons script
    def generate_images(self): 
      try:       
        for filename in os.listdir(self.input_directory):
            f = os.path.join(self.input_directory, filename)
            if os.path.isfile(f):
                plt.figure()					
                data_map = sunpy.map.Map(f)
                bottom_left = SkyCoord(float(self.line_edit6.text())*u.arcsec, float(self.line_edit7.text())*u.arcsec, frame=data_map.coordinate_frame)
                top_right = SkyCoord(float(self.line_edit4.text())*u.arcsec, float(self.line_edit5.text())*u.arcsec, frame=data_map.coordinate_frame) 
                aia_map = data_map.submap(bottom_left, top_right=top_right)
                aia_map.plot()
                aia_map.draw_limb() 
                aia_map.draw_contours([10,20,30,40,50,60,70,80,90]*u.percent)
                plt.colorbar()
                output_filename = os.path.join(self.output_directory, filename.replace('.fits', '.png'))
                plt.savefig(output_filename)
                print("Generated Image:", output_filename)
      except Exception as e:
        error_message = str(e)
        QMessageBox.critical(None, "Error", "An error occurred: {}".format(error_message))
#generated gif button script
    def create_gif(self):
        try:
            jpg_dir = self.input_directory
            first_image = os.path.join(jpg_dir, sorted(os.listdir(jpg_dir))[0])
            first_imagee = cv2.imread(first_image)
            height, width, _ = first_imagee.shape
            video_name = "Sun.mp4"
            video_path = os.path.join(self.output_directory, video_name)
            mp4 = cv2.VideoWriter_fourcc(*'mp4v')
            writer = cv2.VideoWriter(video_path, mp4, 30, (width, height))
            for file_name in sorted(os.listdir(jpg_dir)):
                if file_name.endswith('.jpg'):
                    file_path = os.path.join(jpg_dir, file_name)
                    frame = cv2.imread(file_path)
                    writer.write(frame)
            writer.release()        
        except Exception as e:
            error_message = str(e)
            QMessageBox.critical(None, "Error", "An error occurred: {}".format(error_message))
#############################creating plot button script####################################
    def create_plot(self):
#      try:
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        output_file, _ = QFileDialog.getSaveFileName(self, "Save Plot", self.input_directory, "JPEG Files (*.jpg);;PNG Files (*.png)", options=options)
        plt.ion()
        df = self.input_directory
        Times = []  # x=np.linespace
        LC_Datas = []
        Angless = []  # a list of arrays
        Intensity = []
        for filename in os.listdir(df):
            f = os.path.join(df, filename)
            if os.path.isfile(f):		
                data_map = sunpy.map.Map(f)
                #sunpy 5.0 way of interpolating alongn the line
                line_coords = SkyCoord([round(float(self.line_edit8.text())),round(float(self.line_edit9.text()))],[round(float(self.line_edit10.text())),round(float(self.line_edit11.text()))], unit = (u.arcsec,u.arcsec), frame = data_map.coordinate_frame)
                intensity_coords = sunpy.map.pixelate_coord_path(data_map, line_coords)
                intensity  = sunpy.map.sample_at_coords(data_map, intensity_coords)
                angular_seperation = intensity_coords.separation(intensity_coords[0]).to(u.arcsec)
                Angless.append(angular_seperation)
                Intensity.append(intensity) 
                LC_Data = float(np.sum(intensity.value))
                LC_Datas.append(LC_Data)                
                hdul = fits.open(f)
                time = hdul[1].header['T_OBS']
                Times.append(time)
        NTimes = np.arange(len(os.listdir(df)))*12 #new time based on a numerical system
        NAngless = np.array(Angless) #converting the list of Angless to an array
        print("Length of NTimes:", len(NTimes))
        print("Length of Intensity:", len(Intensity))
#defining the fig and subplots        
        fig, ax = plt.subplots(1, 2, figsize=(20,15))	#creating a new plot object so I have so use ax	
        ax[0].plot(NTimes, LC_Datas)
        ax[0].set_xlabel(f'Time from' + Times[0] + '[s]')
        ax[0].set_ylabel(f'Intensity DNs')
        ax[0].set_title('Time vs Total Intensity')
 #       for NTimese, NAnglesse, Intensitye in zip(NTimes, NAngless, Intensity): #extra e on each idk
 #           t = np.linspace(Intensity[0], Intensity[-1]) *len(NAnglesse)
 #           ax[1].scatter([NTimese] * len(NAnglesse), NAnglesse, c=Intensitye, cmap='copper')
        ax[1].set_xticks(NTimes)
        intensity_map = np.array(Intensity).reshape(len(Times),-1)
        extent = [NTimes[0], NTimes[-1], 0, NAngless.shape[1]]
        im = ax[1].imshow(intensity_map, cmap='copper', extent= extent, aspect = 'auto')
        ax[1].set_xlabel(f'Time from' + Times[0] + '[s]')
        ax[1].set_ylabel(f'Distance along arc ')
        ax[1].locator_params(axis='x', nbins=10) 
        ax[1].set_title('Time vs Distance along Slit Intensity Map')	
#        scatter_plot = ax[1].scatter([], [], c=[], cmap='copper')
        cbar = fig.colorbar(im, ax=[1])
        cbar.set_label('Intensity')
        cbar.mappable.set_clim(np.min(Intensity), np.max(Intensity))				
        plt.savefig(output_file)
###adding in velocity calculation?
#      except Exception as e:
#        error_message = str(e)
#        QMessageBox.critical(None, "Error", "An error occurred: {}".format(error_message))
#############################################creating a mouse hover event hopefully
    def on_click(self, event, source):
        if event.button == QtCore.Qt.LeftButton:
            wcs = self.ax.wcs
            coord = SkyCoord(wcs.pixel_to_world(event.xdata, event.ydata), frame='helioprojective')
            if self.first_click is None:
                if source == 'load_image':
                     self.line_edit8.setText(str(coord.Tx)[:-6])
                     self.line_edit9.setText(str(coord.Ty)[:-6])
                elif source == 'load_image_2':
                     self.line_edit12.setText(str(coord.Tx)[:-6])
                     self.line_edit13.setText(str(coord.Ty)[:-6])
                self.first_click = (event.xdata, event.ydata)
            else:
                wcs = self.ax.wcs
                coord = SkyCoord(wcs.pixel_to_world(event.xdata, event.ydata), frame='helioprojective')
                print(coord)
                self.line_edit10.setText(str(coord.Tx)[:-6])
                self.line_edit11.setText(str(coord.Ty)[:-6])
                self.first_click = None
    def load_image(self):
      try:
            file_dialog = QFileDialog()
            image_path, _ = file_dialog.getOpenFileName(self, "Select FITS Image", "", "FITS Files (*.fits)")
            if image_path:
                self.data_map = Map(image_path)  # Set the data_map attribute
                hdul = fits.open(image_path)
                time = hdul[1].header['T_OBS']
                self.line_edit14.setText(str(time))
                hdul.close() 
                # Define the region of interest
                bottom_left = SkyCoord(float(self.line_edit6.text())*u.arcsec, float(self.line_edit7.text())*u.arcsec, frame=self.data_map.coordinate_frame)
                top_right = SkyCoord(float(self.line_edit4.text())*u.arcsec, float(self.line_edit5.text())*u.arcsec, frame=self.data_map.coordinate_frame) 
                aia_map = self.data_map.submap(bottom_left, top_right=top_right)
                wcs =  aia_map.wcs#
                self.ax = plt.subplot(projection=wcs)#2 part creating projection 
                aia_map.draw_limb(axes=self.ax) 
                self.ax.imshow(aia_map.data, origin='lower', aspect='auto', extent=(bottom_left.data.lon.value, top_right.data.lon.value, bottom_left.data.lat.value, top_right.data.lat.value))               
                plt.ion()
                plt.show()
                plt.gcf().canvas.mpl_connect('button_press_event', lambda event: self.on_click(event, "load_image"))     
      except Exception as e:
        error_message = str(e)
        QMessageBox.critical(None, "Error", "An error occurred: {}".format(error_message))
    def load_image_2(self):
      try:
            file_dialog = QFileDialog()
            image_path, _ = file_dialog.getOpenFileName(self, "Select FITS Image", "", "FITS Files (*.fits)")
            if image_path:
                self.data_map = Map(image_path)  # Set the data_map attribute
                hdul = fits.open(image_path)
                time = hdul[1].header['T_OBS']
                self.line_edit15.setText(str(time))
                hdul.close()                  
                # Define the region of interest
                bottom_left = SkyCoord(float(self.line_edit6.text())*u.arcsec, float(self.line_edit7.text())*u.arcsec, frame=self.data_map.coordinate_frame)
                top_right = SkyCoord(float(self.line_edit4.text())*u.arcsec, float(self.line_edit5.text())*u.arcsec, frame=self.data_map.coordinate_frame) 
                aia_map = self.data_map.submap(bottom_left, top_right=top_right)
                wcs =  aia_map.wcs#
                self.ax = plt.subplot(projection=wcs)#2 part creating projection 
                aia_map.plot(axes=self.ax)
                aia_map.draw_limb(axes=self.ax)
                self.ax.imshow(aia_map.data, origin='lower', aspect='auto', extent=(bottom_left.data.lon.value, top_right.data.lon.value, bottom_left.data.lat.value, top_right.data.lat.value))  
                plt.ion()
                plt.show() 
                plt.gcf().canvas.mpl_connect('button_press_event', lambda event: self.on_click(event, "load_image_2"))   
      except Exception as e:
        error_message = str(e)
        QMessageBox.critical(None, "Error", "An error occurred: {}".format(error_message))
    def Velocity_Calculator(self):
      try:
            x= float(self.line_edit8.text())
            y= float(self.line_edit9.text())
            x2= float(self.line_edit12.text())  #from the first load image sitch
            y2= float(self.line_edit13.text())  #from the 2nd load image sitch
            t= str(self.line_edit14.text())   #time from the fisrt image
            t2= str(self.line_edit15.text())
            time_format = "%Y-%m-%dT%H:%M:%S.%fZ"
            tt = datetime.strptime(t, time_format)#starting the process of taking the time format and finding the difference between the two
            tt2 = datetime.strptime(t2, time_format)                               
            dist = ((x2-x)**2 + (y2-y)**2)**.5 #setting up the distance  
            print(dist*u.arcsec)
            deltt = (tt2-tt).total_seconds()#datetime way of getting the difference between two times in seconds 
            print(deltt)
            arc2km = 725 #conversion constant
            velocity = ((dist * arc2km) / deltt)
            velocity_text = str(velocity)
            self.line_edit18.setText(velocity_text)
      except Exception as e:
        error_message = str(e)
        QMessageBox.critical(None, "Error", "An error occurred: {}".format(error_message))
    def Velocity_Group_Calculator(self):
        try: 
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            output_file, _ = QFileDialog.getSaveFileName(self, "Save Plot", self.input_directory, "JPEG Files (*.jpg);;PNG Files (*.png)", options=options)
            plt.ion()
            df = self.input_directory
            Times = []
            Angless = []  # a list of arrays
            Intensity = []
            first_data_map = None
            great_arc_path = None 
            for filename in os.listdir(df):
                f = os.path.join(df, filename)
                if os.path.isfile(f):
                    data_map = sunpy.map.Map(f) 
                    hdul = fits.open(f)
                    time = hdul[1].header['T_OBS']
                    Times.append(time)
                    #sunpy 5.0 way of interpolating alongn the line
                    line_coords = SkyCoord([round(float(self.line_edit8.text())),round(float(self.line_edit9.text()))],[round(float(self.line_edit10.text())),round(float(self.line_edit11.text()))], unit = (u.arcsec,u.arcsec), frame = data_map.coordinate_frame)
                    intensity_coords = sunpy.map.pixelate_coord_path(data_map, line_coords)
                    intensity  = sunpy.map.sample_at_coords(data_map, intensity_coords)
                    angular_seperation = intensity_coords.separation(intensity_coords[0]).to(u.arcsec)
                    Angless.append(angular_seperation)
                    Intensity.append(intensity)
                    IS = np.power(intensity, 0.3) 
                    if first_data_map is None: 
                        first_data_map = data_map 
                    if great_arc_path is None:
                        great_arc_path = intensity_coords      
            NTimes = np.arange(len(os.listdir(df)))*12 #new time based on a numerical system
            NAngless = np.array(Angless) 
            fig, axs = plt.subplots(3, 2, figsize = None)
            plt.tight_layout() 
            intensity_map = np.array(Intensity).reshape(len(Times),-1)
            extent = [NTimes[0], NTimes[-1], 0, NAngless.shape[1]]
            im1 = axs[0,0].imshow(intensity_map, cmap='copper', extent= extent, aspect = 'auto')
            axs[0,0].set_ylabel(f"Distance Along Slit [{angular_seperation.unit}]")
            axs[0,0].set_xlabel("Time" )
            axs[0,0].set_title('Spectrogram')
            axs[0,0].grid(True)
            cbar0 = fig.colorbar(im1, ax=axs[0,0])
            cbar0.set_label('Intensity')
            cbar0.mappable.set_clim(np.max(Intensity), np.min(Intensity))
#
            intensity_map2 = np.log(np.array(Intensity))     
            intensity_map2 = intensity_map2.reshape(len(Times),-1)
            im2 = axs[0,1].imshow(intensity_map2, cmap='copper', extent= extent, aspect = 'auto')
            axs[0,1].set_ylabel(f"Distance Along Slit [{angular_seperation.unit}]")
            axs[0,1].set_xlabel("Time")
            axs[0,1].set_title('Log Intensity Spectrogram')
            axs[0,1].grid(True)
            cbar1 = fig.colorbar(im2, ax=axs[0,1])
            cbar1.set_label('Intensity')
            cbar1.mappable.set_clim(np.log(np.max(Intensity)), np.log(np.min(Intensity)))              
#
            intensity_map3 = (np.array(Intensity))                    
            im3 = axs[1,0].imshow(intensity_map3, cmap='copper_r', extent= extent, aspect = 'auto')
            axs[1,0].set_ylabel(f"Distance Along Slit [{angular_seperation.unit}]")
            axs[1,0].set_xlabel("Time")
            axs[1,0].set_title('Intensity Inverted Spectrogram')
            axs[1,0].grid(True)
            cbar2 = fig.colorbar(im3, ax=axs[1,0])
            cbar2.set_label('Intensity')
            cbar2.mappable.set_clim(np.max(Intensity), np.min(Intensity))

            intensity_map4 = np.log(np.array(Intensity))     
            intensity_map4 = intensity_map4.reshape(len(Times),-1)
            im4 = axs[1,1].imshow(intensity_map4, cmap='copper_r', extent= extent, aspect = 'auto')             
            axs[1,1].set_ylabel(f"Distance Along Slit [{angular_seperation.unit}]")
            axs[1,1].set_xlabel("Time")
            axs[1,1].set_title('Intensity Inverted + Log Scaled Spectrogram')
            axs[1,1].grid(True)
            cbar2 = fig.colorbar(im4, ax=axs[1,1])
            cbar2.set_label('Intensity')
            cbar2.mappable.set_clim(np.max(Intensity), np.min(Intensity))

            ax = plt.subplot(3, 2, 5, projection=first_data_map)
            first_data_map.plot(axes= ax, clip_interval=(1, 99.99)*u.percent)
            ax.plot_coord(great_arc_path)
            ax.plot_coord(great_arc_path[0], marker="o", color="blue", label="start")
            ax.plot_coord(great_arc_path[1], marker="o", color="green", label="end")
            axs[2,1].axis('off')
            plt.show()            
#calculation time
            arc2km = 725  # 725km in 1 arcsec
            x2 = float(self.line_edit8.text())
            y2 = float(self.line_edit9.text())
            x1 = float(self.line_edit10.text())
            y1 = float(self.line_edit11.text())
            t = str(self.line_edit14.text())   # time from the first image
            t2 = str(self.line_edit15.text())
            time_format = "%Y-%m-%dT%H:%M:%S.%fZ"
            tt = datetime.strptime(t, time_format)  # starting the process of taking the time format and finding the difference between the two
            print("Parsed Successfully:", tt)
            tt2 = datetime.strptime(t2, time_format)                               
            dt_distance = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5 * arc2km
            dt_time = (tt2 - tt).total_seconds()
            velocity = (dt_distance / abs(dt_time)) * u.km / u.s
            velocity_text = str(velocity)
            self.line_edit18.setText(velocity_text)            
        except Exception as e:
         tb_str = traceback.format_tb(e.__traceback__)
         line_number = int(tb_str[-1].split(",")[1].split()[-1])
         error_message = str(e)
         error_with_line = f"An error occurred at line {line_number}: {error_message}"
         QMessageBox.critical(None, "Error", error_with_line)

if __name__ == '__main__':
    app = QApplication([])
    window = MainWindow()
    window.show()
    app.exec_()
