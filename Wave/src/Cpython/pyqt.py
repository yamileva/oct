import sys
import os
import pyqtgraph as pg
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
from pyqtgraph.Qt import QtWidgets, mkQApp
import pyqtgraph.exporters 
import imageio.v2 as imageio
import scipy.signal as ss
import numpy as np
import h5py

from PySide6.QtCore import Slot,QRectF
from PySide6.QtGui import QAction, QKeySequence, QShortcut, QFont, QIcon
from PySide6.QtWidgets import (QApplication, QHBoxLayout, QFileDialog,QDialog, QCheckBox, QGroupBox, QGridLayout, QSpinBox, QDoubleSpinBox,
							   QLabel, QMainWindow, QVBoxLayout, QPushButton, QMessageBox, QWidget)
from enum import Enum
from qt_material import apply_stylesheet

class Dimension(Enum):
	One,Two,Three = range(1,4)

class Database():
	def __init__(self,):
		self.dim = None
		self.wave = None
		self.xi = None
		self.phi = None
		self.scan = None
		self.pml_crd = None
		self.skin_crd = None
		self.h = None
		self.left_borders = None
		self.right_borders = None
		self.steps = None
		self.fname = None
		self.ptr_time = 0
		self.color = (45, 202, 240)
		self.imgLevels = None
		self.imgs = []
		self.size = 1000
		self.file_list = []
	def read_file(self,file_name, dim):
		import os
		self.dim = dim
		self.fname = file_name
		grid_settings = None
		if os.path.splitext(self.fname)[1] == ".h5":
			file = h5py.File(self.fname, 'r')
			self.wave = file['wave']
			p_file = h5py.File(os.path.splitext(self.fname)[0] + '_params.h5', 'r')
			grid_settings = p_file['grid_settings'][0]
			self.pml_crd = np.column_stack((p_file['pml_left'][0],p_file['pml_right'][0]))
			self.skin_crd = p_file['skin_crd'][0]
			self.layer_t,self.layer_x,self.layer_y = p_file['write_layer_N'][0]
			self.lamda, self.omega_min, self.omega_max, self.tau = p_file['light_param'][0]
			print("l: ",self.layer_t,self.layer_x,self.layer_y)
			# self.omega =  p_file['omega'][0]
			print(*self.skin_crd)
		else:
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("Please select the .h5 file format!")
			msg.exec()
			return
		if self.dim == Dimension.One:
			self.h = {'time': grid_settings[0], 'x':  grid_settings[1]}
			self.left_borders = {'time': grid_settings[2], 'x':  grid_settings[3]}
			self.right_borders = {'time': grid_settings[4], 'x':  grid_settings[5]}
			self.steps = {'time': (self.right_borders['time'] - self.left_borders['time'])/self.h['time'],
		 				 'x':  (self.right_borders['x'] - self.left_borders['x'])/self.h['x']}
		elif self.dim == Dimension.Two:
			self.h = {'time': grid_settings[0], 'x':  grid_settings[1],'y':  grid_settings[2]}
			self.left_borders = {'time': grid_settings[3], 'x':  grid_settings[4], 'y':  grid_settings[5]}
			self.right_borders = {'time': grid_settings[6], 'x':  grid_settings[7], 'y':  grid_settings[8]}
			self.steps = {'time': (self.right_borders['time'] - self.left_borders['time'])/self.h['time'],
		 				 'x':  (self.right_borders['x'] - self.left_borders['x'])/self.h['x'],
						 'y':  (self.right_borders['y'] - self.left_borders['y'])/self.h['y']}
		elif self.dim == Dimension.Three:
			self.h = {'time': grid_settings[0], 'x':  grid_settings[1],'y':  grid_settings[2], 'z':  grid_settings[3]}
			self.left_borders = {'time': grid_settings[4], 'x':  grid_settings[5], 'y':  grid_settings[6], 'z':  grid_settings[7]}
			self.right_borders = {'time': grid_settings[8], 'x':  grid_settings[9], 'y':  grid_settings[10], 'z':  grid_settings[11]}
	def plot(self, graph = None, in_curve=None):
		curve = in_curve
		if self.dim == Dimension.One:
			if curve == None:
				curve = graph.plot(pen={'color':self.color, 'width':2})
				for pmlx in self.pml_crd[0]:
					graph.addItem(pg.InfiniteLine(pos=[pmlx,0],movable=False, angle=90, pen={'color':(255,0,0),'width':1}, label='x={value:0.2f}',
								labelOpts={'position':0.1, 'color': (0,0,0), 'fill': (0,0,0,50), 'movable': False}))
				for sk_x in self.skin_crd:
					graph.addItem(pg.InfiniteLine(pos=[sk_x,0],movable=False, angle=90, pen={'color':(29, 233, 182),'width':1}))
			curve.setData(np.linspace(self.left_borders['x'],self.right_borders['x'],len(self.wave[0])),
						self.wave[self.ptr_time])
		elif self.dim == Dimension.Two:
			print(self.ptr_time * self.h['time'])
			data = self.wave[self.ptr_time].reshape(int((self.steps['x']/self.layer_x + 1)), int((self.steps['y']/self.layer_y + 1)))
			if curve == None:
				t_min = -0.01#-0.0003#np.min(self.wave)/4
				t_max = 0.01#0.0003#np.max(self.wave)/4
				self.imgLevels = (t_min, t_max)
				#print("time: ", self.ptr_time * self.h['time'],self.imgLevels)
				curve = pg.ImageItem(image = data, rect = QRectF( self.left_borders['x'], self.left_borders['y'], self.right_borders['x'], self.right_borders['y']),
							levels = self.imgLevels)
				graph.addItem(curve)
				graph.addColorBar(curve,colorMap='CET-D9', interactive=False)#,colorMap='CET-',orientation = 'h'
				for pmlx in self.pml_crd[0]:
					graph.addItem(pg.InfiniteLine(pos=[pmlx,0],movable=False, angle=90, pen={'color':(255,0,0),'width':1}, label='x={value:0.2f}',
							labelOpts={'position':0.1, 'color': (0,0,0), 'fill': (0,0,0,50), 'movable': False}))
				for pmly in self.pml_crd[1]:
					graph.addItem(pg.InfiniteLine(pos=[0,pmly],movable=False, angle=0, pen={'color':(255,0,0),'width':1}, label='y={value:0.2f}',
							labelOpts={'position':0.1, 'color': (0,0,0), 'fill': (0,0,0,50), 'movable': False}))
				for sk_x in self.skin_crd:
					graph.addItem(pg.InfiniteLine(pos=[sk_x,0],movable=False, angle=90, pen={'color':(29, 233, 182),'width':1}))
			else:
				curve.setImage(data, levels = self.imgLevels) #levels = self.imgLevels
		return curve

	def update_wave(self, ptr, graph, in_curve):
		self.ptr_time = ptr
		return self.plot(graph,in_curve)

	def scan_plot(self,x = None, y = None, graph = None, curve = None):
		def omega_f(t):
			return ((t/self.tau)*(self.omega_max - self.omega_min) + self.omega_min)
		delta_y = (35/1300) * 1e3
		radius = 0.2
		if self.dim == Dimension.One:
			if curve == None:
				curve = graph.plot(pen=self.color)
			self.scan = self.wave[:,int((x-self.left_borders['x'])/(self.h['x']*self.layer_x))]
			curve.setData(np.linspace(self.left_borders['time'],self.right_borders['time'],len(self.scan)),
		 				self.scan)	
		elif self.dim == Dimension.Two:
			self.scan = np.zeros(len(self.wave), dtype=np.float64)
			self.ref  = np.zeros(len(self.wave), dtype=np.float64)
			for t in range(len(self.wave)):
				data = self.wave[t].reshape(int((self.steps['x']/self.layer_x + 1)), int((self.steps['y']/self.layer_y + 1)))
				self.scan[t] = data[int((x-self.left_borders['x'])/(self.h['x']*self.layer_x))][int((y-self.left_borders['y'])/(self.h['y']*self.layer_y))]
				self.ref[t] = np.exp(-((self.pml_crd[1][0] + (self.pml_crd[1][1] - self.pml_crd[1][0])/2) - delta_y)**2/radius**2)*np.sin(omega_f(t * self.h['time']) * t * self.h['time'])
			if curve == None:
				curve = graph.plot(pen=self.color)
			curve.setData(np.linspace(self.left_borders['time'],self.right_borders['time'],len(self.scan)),self.scan)
		return curve
	def A_scan_plot(self, graph = None, curve = None):
		self.A_scan = 2*self.scan*self.ref
		from scipy.fft import fft, fftfreq
		yf = fft(self.A_scan)
		from scipy.signal import blackman
		w = blackman(len(self.A_scan))
		#w = np.kaiser(len(self.A_scan), 2)
		#w = np.hanning(len(self.A_scan) + 1)[:-1]
		ywf = fft(self.A_scan*w)
		from scipy.constants import speed_of_light
		#xf = np.arange(len(self.A_scan)) * speed_of_light / 2 / 1200e-9
		xf = fftfreq(len(self.scan), self.h['time'])
		if curve == None:
			curve = graph.plot(pen=self.color)
		curve.setData(xf, np.abs(ywf)) 
		return curve

	def Intensity_plot(self, x_range, graph = None, curve = None):
		if self.dim == Dimension.One:
			x_list = np.arange(x_range[0],x_range[1], 0.5)
			Intensity = np.zeros(len(x_list), dtype=np.float64)
			for i,x in enumerate(x_list):
				print("Intensity_plot(%): ", (i*100)/len(x_list))
				Intensity[i] = np.max(self.wave[:,int((x-self.left_borders['x'])/self.h['x'])])**2
			if curve == None:
					curve = graph.plot(pen=self.color)
			curve.setData(x_list,
						Intensity)
		elif self.dim == Dimension.Two:
			#print("Intensity_plot")
			#y = (30/1300)*1e3
			y = 30
			x_list = np.arange(x_range[0],x_range[1], 0.5)
			Intensity = np.zeros(len(x_list), dtype=np.float64)
			scan = np.zeros(len(self.wave), dtype=np.float64)
			for i,x in enumerate(x_list):
				print("Intensity_plot(%): ", (i*100)/len(x_list))
				for t in range(len(self.wave)):
					data = self.wave[t].reshape(int((self.steps['x']/self.layer_x + 1)), int((self.steps['y']/self.layer_y + 1)))
					scan[t] = data[int((x-self.left_borders['x'])/self.h['x'])][int((y-self.left_borders['y'])/self.h['y'])]
				Intensity[i] = np.max(scan)
			if curve == None:
					curve = graph.plot(pen=self.color)
			curve.setData(x_list,
						Intensity) 
		return curve
	def save_img(self,view, dir, step = 5):
		if self.ptr_time % step == 0:
			#ImageExporter
			exporter = pg.exporters.ImageExporter(view.scene())
			img_names = dir + str(self.ptr_time) + '.jpeg'
			exporter.export(img_names)
			self.imgs.append(img_names)
	def create_gif(self, dir):
		with imageio.get_writer(dir + 'wave' + str(self.dim.value) + '.gif', mode='I') as writer:
			for filename in self.imgs:
				#image = svg2rlg(filename)
				image = imageio.imread(filename)
				writer.append_data(image)
	def clear_image(self, dir):
		for file in os.scandir(dir):
			os.remove(file.path)

	def cpp_read(self, dim, cpp_datas):
		data, grid_settings, pml_crd, xi_data = cpp_datas
		self.dim = dim
		self.wave = data
		self.pml_crd = pml_crd
		self.xi = xi_data
		if self.dim == Dimension.One:
			self.h = {'time': grid_settings[0], 'x':  grid_settings[1]}
			self.left_borders = {'time': grid_settings[2], 'x':  grid_settings[3]}
			self.right_borders = {'time': grid_settings[4], 'x':  grid_settings[5]}
		elif self.dim == Dimension.Two:
			self.h = {'time': grid_settings[0], 'x':  grid_settings[1],'y':  grid_settings[2]}
			self.left_borders = {'time': grid_settings[3], 'x':  grid_settings[4], 'y':  grid_settings[5]}
			self.right_borders = {'time': grid_settings[6], 'x':  grid_settings[7], 'y':  grid_settings[8]}
		elif self.dim == Dimension.Three:
			self.h = {'time': grid_settings[0], 'x':  grid_settings[1],'y':  grid_settings[2], 'z':  grid_settings[3]}
			self.left_borders = {'time': grid_settings[4], 'x':  grid_settings[5], 'y':  grid_settings[6], 'z':  grid_settings[7]}
			self.right_borders = {'time': grid_settings[8], 'x':  grid_settings[9], 'y':  grid_settings[10], 'z':  grid_settings[11]}
	
class CustomDialog(QDialog):
	def __init__(self):
		super().__init__()

		self.setWindowTitle("Dim!")

		self.dim_list = [Dimension.One,Dimension.Two,Dimension.Three]
		self.select_dim = None

		self.layout = QVBoxLayout()
		message = QLabel("Please select the dimension of the graph.")
		message.setStyleSheet("font_family': 'Roboto")
		message.setStyleSheet("font-size: 20px")
		btn = QPushButton("I chose. I'm not lying!")
		btn.clicked.connect(self.close_dlg)

		self.check_groupbox = QGroupBox(self)
		self.check_layout = QHBoxLayout()
		self.check_dims = [QCheckBox("1d",self),QCheckBox("2d",self),QCheckBox("3d",self)]
		for check in self.check_dims:
			self.check_layout.addWidget(check)
			check.stateChanged.connect(self.on_statedChanged)

		self.check_groupbox.setLayout(self.check_layout)

		self.layout.addWidget(message)
		self.layout.addWidget(self.check_groupbox)
		self.layout.addWidget(btn)
		
		self.setLayout(self.layout)
	def close_dlg(self):
		for idx,item in enumerate(self.check_dims):
			if item.isChecked():
				self.select_dim = self.dim_list[idx]
		if (self.select_dim != None):
			self.close()
		
	def on_statedChanged(self,state):
		for item in self.check_dims:
			if item != state:
				item.setChecked(False)

class ApplicationWindow(QMainWindow):
	def __init__(self, parent=None):
		QMainWindow.__init__(self, parent)
		widget = QWidget()
		self.Glayout = QGridLayout()
		#
		self.img_dir = '../../img/'
		self.icon_dir = '../../icon/'
		self.gif_dir = '../../gif/'
		#time
		self.time = 0
		self.incr = 1
		self.running = False
		self.gif = True
		self.setWindowIcon(QIcon("../../icon/icon.png"))
		self.setWindowTitle("PML")
		#var
		self.fname = None
		self.dim = None
		self.database = Database()

		#graph
		self.first_graph = None			#wave
		self.scan_graph = None			#scan
		#self.A_scan_graph = None		#A scan
		self.Intensity_graph = None		#

		#curve
		self.first_curve = None			#wave
		self.scan_curve = None			#scan
		#self.A_scan_curve = None		#A scan
		self.Intensity_curve = None	#

		# Main menu bar
		self.menu = self.menuBar()
		self.menu_file = self.menu.addMenu("File")

		open_file = QAction("Open file", self, triggered=self.open_file_dialog)
		open_file.setShortcuts(QKeySequence('Ctrl+F'))
		exit = QAction("Exit", self, triggered=QApplication.quit)
	
		actionlist = [open_file,exit]
		for i in actionlist:
			self.menu_file.addAction(i)

		self.menu_about = self.menu.addMenu("&About")
		self.menu_about.addAction("Hot keys", self.about, shortcut='F1')

		#Hot keys
		self.shortcut_play = QShortcut(QKeySequence('Ctrl+P'), self)
		self.shortcut_play.activated.connect(self.run_wave)
		#print(self.shortcut_play.text())

		self.shortcut_stop = QShortcut(QKeySequence('Ctrl+S'), self)
		self.shortcut_stop.activated.connect(self.stop_wave)

		self.shortcut_restart = QShortcut(QKeySequence('Ctrl+R'), self)
		self.shortcut_restart.activated.connect(self.restart_wave)

		self.shortcut_faster = QShortcut(QKeySequence('Ctrl+M'), self)
		self.shortcut_faster.activated.connect(self.faster)

		self.shortcut_slowly = QShortcut(QKeySequence('Ctrl+L'), self)
		self.shortcut_slowly.activated.connect(self.slowly)


		self.shortcut_scan = QShortcut(QKeySequence('Ctrl+A'), self)
		self.shortcut_scan.activated.connect(self.open_close_subplot)


		self.shortcut_about = ['Ctrl+P - Start animate wave',
			 				'Ctrl+S - Stop animate wave',
							'Ctrl+R - Restart animate wave',
							'Ctrl+M - Faster animate wave',
							'Ctrl+L - Slowly animate wave',
							'Ctrl+A - Close/Open A-scan graph']
		# spinbox
		self.spinbox_x = QDoubleSpinBox()
		self.spinbox_x.setPrefix("X: ")
		self.spinbox_x.setSingleStep(0.5)
		self.spinbox_x.valueChanged.connect(self.value_change)  #####

		self.spinbox_y = QDoubleSpinBox()
		self.spinbox_y.setPrefix("Y: ")
		self.spinbox_y.setSingleStep(0.5)
		self.spinbox_y.valueChanged.connect(self.value_change)  #####

		self.spinbox = [self.spinbox_x, self.spinbox_y]
		# 
		tmp_n = 0
		for i,s in enumerate(self.spinbox):
			s.hide()
			self.Glayout.addWidget(s, i + 1, 0)
			tmp_n = i + 1

		self.btn = QPushButton("Show scan!")
		self.btn.clicked.connect(self.start_subplot)
		self.btn.hide()
		self.Glayout.addWidget(self.btn, tmp_n + 1, 0)
		# Central widget
		self.gr_wid = pg.GraphicsLayoutWidget(show=True)
		self.Glayout.addWidget(self.gr_wid, 0, 0)
		
		widget.setLayout(self.Glayout)
		self.setCentralWidget(widget)
	
		self.gr_wid.scene().sigMouseClicked.connect(self.mouse_clicked)	

	def open_file_dialog(self):
		self.restart_wave()
		self.gr_wid.clear()
		fname, _ = QFileDialog.getOpenFileName(self, "Выберите файл", '.', "All types (*.*)")
		self.fname = fname
		print(f"fname = {self.fname}")
		self.showDialog_dim()
		self.database.read_file(self.fname,self.dim)
		if self.database.wave != None:
			self.spinbox_x.setRange(self.database.left_borders['x'], self.database.right_borders['x'])
			self.spinbox_y.setRange(self.database.left_borders['y'], self.database.right_borders['y'])
			self.plot_all()
			self.scan_graph.hide()

	def main_plot(self):
		self.gr_wid.clear()
		self.first_graph = self.gr_wid.addPlot(title="Wave")
		self.first_curve = self.database.plot(self.first_graph)
		
	def subplot(self):
		self.spinbox_x.setValue(self.database.pml_crd[0][0])
		self.spinbox_y.setValue(self.database.pml_crd[1][0] + (self.database.pml_crd[1][1] - self.database.pml_crd[1][0])/2)
		self.gr_wid.nextRow()
		self.scan_graph = self.gr_wid.addPlot(title="Scan")
		self.scan_curve = self.database.scan_plot(self.spinbox_x.value(), self.spinbox_y.value(), self.scan_graph)
		##self.gr_wid.nextRow()
		#self.Intensity_graph = self.gr_wid.addPlot(title="Intensity")
		#self.Intensity_curve = self.database.Intensity_plot((self.database.pml_crd[0][0],self.database.pml_crd[0][1]), self.Intensity_graph)
		self.gr_wid.nextRow()
		self.A_scan_graph = self.gr_wid.addPlot(title="A-Scan")
		self.A_scan_graph.setLogMode(x=False, y=False)
		#self.A_scan_graph.addLegend(offset=(-100, 10), colCount = 2, labelTextSize ='12pt')
		self.A_scan_curve = self.database.A_scan_plot(self.A_scan_graph)
		# for text, spin in self.spins:
		# 	spin.show()
		# 	#label.hide()

	def cpp_plot(self,cpp_datas):
		if(len(cpp_datas) > 2):
			self.dim = cpp_datas[0]
			self.database.cpp_read(self, self.dim, cpp_datas)
		else:
			self.dim, self.fname = cpp_datas
			self.database.read_file(self.fname,self.dim)
		if self.database.wave != None:
			self.plot_all()
			self.scan_graph.hide()
			self.A_scan_graph.hide()
			#self.Intensity_graph.hide()

	def plot_all(self):
		self.main_plot()
		self.subplot()
		#	self.about()
	@Slot()
	def open_close_subplot(self):
		if self.scan_graph.isVisible():
			self.scan_graph.hide()
			self.A_scan_graph.hide()
			#self.Intensity_graph.hide()
			self.btn.hide()
			if self.dim == Dimension.One:
				self.spinbox_x.hide()
			elif self.dim == Dimension.Two:
				self.spinbox_x.hide()
				self.spinbox_y.hide()
		else:
			self.scan_graph.show()
			self.A_scan_graph.show()
			#self.Intensity_graph.show()
			self.btn.show()
			if self.dim == Dimension.One:
				self.spinbox_x.setRange(self.database.left_borders["x"], self.database.right_borders["x"])
				self.spinbox_x.show()
			elif self.dim == Dimension.Two:
				self.spinbox_x.setRange(self.database.left_borders["x"], self.database.right_borders["x"])
				self.spinbox_x.show()
				self.spinbox_y.setRange(self.database.left_borders["y"], self.database.right_borders["y"])
				self.spinbox_y.show()
		# if self.gr_wid.getItem(1,0) == None:
		# 	self.gr_wid.addItem(self.scan_graph, row = 1, col= 0)
		# else:
		# 	self.gr_wid.removeItem(self.scan_graph)

	def showDialog_dim(self):
		dlg = CustomDialog()
		dlg.exec()
		self.dim = dlg.select_dim
	@Slot()
	def run_wave(self):
		if self.time + self.incr < len(self.database.wave) - 1:
			self.first_curve.clear()
			#self.database.save_img(self.gr_wid,self.img_dir)
			self.database.update_wave(self.time,self.first_graph,self.first_curve)
			self.running = True
		while self.running == True and self.time + self.incr < len(self.database.wave) - 1:
			self.time += self.incr
			#self.database.save_img(self.gr_wid,self.img_dir)
			self.database.update_wave(self.time,self.first_graph,self.first_curve)
			QApplication.processEvents()
	@Slot()	
	def stop_wave(self):
		self.running = False
	@Slot()
	def restart_wave(self):
		self.stop_wave()
		#self.database.create_gif(self.gif_dir)
		#self.database.clear_image(self.img_dir)
		self.incr = 1
		self.time = 0
		self.database.update_wave(self.time,self.first_graph,self.first_curve)
	@Slot()
	def start_subplot(self):
		self.scan_curve.clear()
		self.scan_curve = self.database.scan_plot(self.spinbox_x.value(), self.spinbox_y.value(), self.scan_graph, self.scan_curve)
		self.A_scan_curve = self.database.A_scan_plot(self.A_scan_graph, self.A_scan_curve)
	@Slot()
	def faster(self):
		self.incr +=1
	@Slot()
	def slowly(self):
		if self.incr - 1 != 0:
			self.incr = self.incr - 1 
	@Slot()
	def about(self):
		# Hot keys names
		text = "Hot keys:\n"
		for keys in self.shortcut_about:
			text += keys + "\n"
		QMessageBox.about(self, "Информация", text)
	@Slot()
	def closeEvent(self, event):
		self.stop_wave()
		event.accept() # let the window close
	@Slot()
	def mouse_clicked(self,event):
		self.vb = self.first_graph.vb
		items = self.first_graph.scene().items(event.scenePos())
		mousePoint = self.vb.mapSceneToView(event._scenePos)
		if self.first_graph.sceneBoundingRect().contains(event._scenePos):
			mousePoint = self.vb.mapSceneToView(event._scenePos)
			index = int(mousePoint.x())
			#print(f'clicked plot X: {mousePoint.x()}, Y: {mousePoint.y()}')
	@Slot()
	def value_change(self):
		pass
	
def cpp_app(parametrs):
	
	if(len(parametrs) > 2):
		dim = parametrs[0]
		data = np.asarray(parametrs[1], dtype = np.float64)
		grid_settings = np.asarray(parametrs[2], dtype = np.float64)
		pml_crd = np.asarray(parametrs[3], dtype = np.float64)
		xi_data = parametrs[4]

		in_cpp_datas = (dim,data,grid_settings,pml_crd,xi_data)
	else:
		in_cpp_datas = parametrs
		
	if not QApplication.instance():
		app = QApplication(sys.argv)
	else:
		app = QApplication.instance()
	apply_stylesheet(app, theme='light_teal.xml', invert_secondary=True)
	w = ApplicationWindow()
	w.cpp_plot(in_cpp_datas)
	w.resize(1280, 720)
	w.show()
	app.exec()

## Start Qt event loop
if __name__ == '__main__':	
	# cpp_app((Dimension.One,"../../Data/temp/"))
	cpp_app((Dimension.Two,"../../Data/Save/hm2_28.h5"))
