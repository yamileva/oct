import sys
from PySide6.QtCore import Slot
from PySide6.QtWidgets import (QApplication, QMainWindow, QGridLayout, QTableWidgetItem, QWidget, QSpinBox,
							   QTableWidget, QHeaderView)
from qt_material import apply_stylesheet
import numpy as np
import bisect
from sklearn import preprocessing
from scipy.constants import speed_of_light, epsilon_0, pi,mu_0
from scipy.interpolate import splrep, BSpline
import json
import pyqtgraph as pg

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')


# TODO: перейти на частоту
class Dielectric_params(QMainWindow):
	def __init__(self, settings_file_path,lamda = 1000,degree_fit = 0, parent=None):
		QMainWindow.__init__(self, parent)
		widget = QWidget()

		self.permeability = {}
		self.conductivity = {}
		self.refraction = {}

		self.min_length = 800
		self.max_length = 1400

		self.settings_file = settings_file_path
		self.params = None

		self.xx = np.arange(0, 1, 0.001)

		self.chromophorf_data = {}
		self.chromophorf_func = {}

		self.absorpton = {}
		self.scattering = {}
		self.layers_name = ["SC", "LE", "PD", "UD", "RD"]

		self.Glayout = QGridLayout()
		self.win = pg.GraphicsLayoutWidget(show=True)
		self.win.setWindowTitle('Data')
		self.graph = self.win.addPlot(title="Chromophors")
		self.graph.addLegend(offset=(-100, 10), colCount = 2, labelTextSize ='12pt')
		self.graph.setLogMode(x=False, y=False)

		self.spinbox_s = QSpinBox()
		self.spinbox_s.setPrefix("S: ")
		self.spinbox_s.setSingleStep(1)
		self.spinbox_s.setRange(0, 20)
		self.spinbox_s.valueChanged.connect(self.value_change_interp)  #####

		# plot window goes on right side, spanning 3 rows
		self.Glayout.addWidget(self.win, 0, 0)
		self.Glayout.addWidget(self.spinbox_s, 0, 1)

		self.table = QTableWidget()
		header = self.table.horizontalHeader()
		header.setSectionResizeMode(QHeaderView.Stretch)
		self.table.setRowCount(3)
		self.table.setColumnCount(len(self.layers_name) + 1)

		self.spinbox = QSpinBox()
		self.spinbox.setPrefix("lamda: ")
		self.spinbox.setSingleStep(10)
		self.spinbox.setRange(self.min_length, self.max_length)
		self.spinbox.valueChanged.connect(self.value_change_calc)  ####

		self.Glayout.addWidget(self.table, 1, 0)
		self.Glayout.addWidget(self.spinbox, 1, 1)
		# self.setLayout(self.Glayout)
		widget.setLayout(self.Glayout)
		self.setCentralWidget(widget)

		self.load_data()
		self.interpolate_data(degree_fit)
		self.calc_dielectric_parameters(lamda)
		#self.set_table_data()

	@Slot()
	def value_change_interp(self):
		self.graph.clear()
		self.interpolate_data(self.spinbox_s.value())
		self.calc_dielectric_parameters(self.spinbox.value())
		self.set_table_data()

	@Slot()
	def value_change_calc(self):
		self.calc_dielectric_parameters(self.spinbox.value())
		self.set_table_data()

	def interpolate_data(self, k=0):
		for key, value in self.chromophorf_data.items():
			x = value[:, 0]
			y = value[:, 1]
			tck_s = splrep(x, y, s=k)
			self.chromophorf_func[key] = BSpline(*tck_s)
			self.graph.plot(*value.T, pen=None, symbol='o', symbolPen=self.colors[key], symbolSize=2.5,
							name="data_" + key)
			self.graph.plot(self.xx, BSpline(*tck_s)(self.xx), pen=self.colors[key], name="splrep_" + key)

	def parse(self):
		with open(self.settings_file) as user_file:
			self.params = json.load(user_file)

	def load_data(self):
		self.parse()
		self.chromophorf_data = {
			"HbO2": np.loadtxt(self.params["Path_HbO2"], dtype=np.float64),
			"Hb": np.loadtxt(self.params["Path_Hb"], dtype=np.float64),
			"Water": np.loadtxt(self.params["Path_Water"], dtype=np.float64),
			"Lipid": np.loadtxt(self.params["Path_Lipid"], dtype=np.float64),
		}
		self.colors = {
			"HbO2": (247, 0, 255),
			"Hb": (255, 0, 0),
			"Water": (0, 17, 255),
			"Lipid": (11, 255, 3)
		}
		self.prepare_data()

	def prepare_data(self):
		for key, value in self.chromophorf_data.items():
			# slice array to 800 - 1400 nm
			index = [bisect.bisect_right(value[:, 0], self.min_length),
					 bisect.bisect_right(value[:, 0], self.max_length)]
			self.chromophorf_data[key] = value[index[0]:index[1] + 1]

			# delete unique rows basic in lamda
			_, idx = np.unique(self.chromophorf_data[key][:, 0], return_index=True)
			self.chromophorf_data[key] = self.chromophorf_data[key][idx]

			# convert lamda to omega
			# chromophorf[key][:,0] = (c * 2 * pi)/chromophorf[key][:,0]

			# normalize
			scaler = preprocessing.MinMaxScaler()
			self.chromophorf_data[key][:, 0] = scaler.fit_transform(self.chromophorf_data[key][:, 0].reshape(-1, 1)).reshape(1, -1)
		self.chromophorf_func["Melanin"] = lambda x: (1.7 * (10**12)) * (x**(-3.48))
		self.chromophorf_func["Base"] = lambda x: (7.84 * (10**7)) * (x**(-3.225))

	# сюда приходит нормированное значение
	def calc_absorpton_layers(self, lamda_n):
		lamda = lamda_n * (self.max_length - self.min_length) + self.min_length
		S = self.params["FHb"] * self.params["FRBC"] * self.params["Ht"]

		# #print(self.chromophorf_func["Lipid"](lamda_n), self.params["fLipid"][0] )
		for index, layer in enumerate(self.layers_name):
			#print("-"*10)
			#print(layer)
			#print("refraction:",self.params["refraction"][index])
			#print("Base:", self.chromophorf_func["Base"](lamda))
			#print("fBlood:", self.params["fBlood"][index], "\tBlood: ",
				  #S * self.chromophorf_func["HbO2"](lamda_n) + (1 - S) * self.chromophorf_func["Hb"](lamda_n))
			#print("fWater: ", self.params["fWater"][index], "\tWater: ", self.chromophorf_func["Water"](lamda_n))
			#print("fMelanin: ", self.params["fMelanin"][index], "\tMelanin: ", self.chromophorf_func["Melanin"](lamda))
			#print("fLipid: ", self.params["fLipid"][index], "\tLipid: ", self.chromophorf_func["Lipid"](lamda_n))

			self.absorpton[layer] = self.chromophorf_func["Base"](lamda) \
			+ self.params["fBlood"][index] * (S * self.chromophorf_func["HbO2"](lamda_n) + (1 - S) * self.chromophorf_func["Hb"](lamda_n))\
			+ self.params["fWater"][index] * self.chromophorf_func["Water"](lamda_n)\
			+ self.params["fMelanin"][index] * self.chromophorf_func["Melanin"](lamda)\
			+ self.params["fLipid"][index] * self.chromophorf_func["Lipid"](lamda_n)

			#print(f"Absorption: {self.absorpton[layer]:.10f}")
			#print("-" * 10)
	def calc_scattering_layers(self, lamda):
		for index, layer in enumerate(self.layers_name):
			self.scattering[layer] = self.params["a"][index]*(self.params["f_ray"][index] * (lamda/500)**(-4)
					       + (1-self.params["f_ray"][index]) * (lamda/500)**self.params["b_mie"][index])

	# #print(self.absorpton)

	# сюда приходит ненормированное значение
	def calc_dielectric_parameters(self, lamda):
		lamda_n = (lamda - self.min_length) / (self.max_length - self.min_length)
		self.calc_absorpton_layers(lamda_n)
		self.calc_scattering_layers(lamda)
		print("-"*10)
		for index, layer in enumerate(self.layers_name):
			print("M-t:" + layer, self.absorpton[layer])
			self.refraction[layer] = [self.params["refraction"][index], (self.absorpton[layer] + self.scattering[layer]) * (lamda * (10**-7))/ (4 * pi)]
			#print("-" * 10)
			#print(f"Re refraction {layer}:",self.refraction[layer][0],f"\nIm refraction {layer}:",self.refraction[layer][-1])
			self.permeability[layer] = [self.refraction[layer][0] ** 2 - self.refraction[layer][-1] ** 2,
										2 * self.refraction[layer][0] * self.refraction[layer][-1]]
			self.conductivity[layer] = (speed_of_light/(lamda * (10**-9))) * (2*pi*epsilon_0) * self.permeability[layer][-1]
			#print(f"Conductivity {layer}:", self.conductivity[layer])
			#print("-" * 10)

		#print(epsilon_0)
		#print(mu_0)
	# #print(self.permeability)
	# #print(self.conductivity)

	def set_table_data(self):
		self.table.setHorizontalHeaderLabels(["Params"] + self.layers_name)
		self.table.setItem(0, 0, QTableWidgetItem("Permeability(Re)"))
		self.table.setItem(1, 0, QTableWidgetItem("Permeability(Im)"))
		self.table.setItem(2, 0, QTableWidgetItem("Conductivity"))
		for index, layer in enumerate(self.layers_name):
			self.table.setItem(0, index + 1, QTableWidgetItem(f"{self.permeability[layer][0]:.10f}"))
			self.table.setItem(1, index + 1, QTableWidgetItem(f"{self.permeability[layer][-1]:.10f}"))
			self.table.setItem(2, index + 1, QTableWidgetItem(f"{self.conductivity[layer]:.10f}"))
			##print(self.permeability[layer][0],end=",")
			##print(self.conductivity[layer],end=",")

def cpp_dielectric_params(lamda):
	if not QApplication.instance():
		app = QApplication(sys.argv)
	else:
		app = QApplication.instance()
	apply_stylesheet(app, theme='light_teal.xml', invert_secondary=True)
	w = Dielectric_params("src/Cpython/settings.json",lamda)
	
	tpl = (list(np.array(list(w.permeability.values()))[:,0]), list(w.conductivity.values()))
	print(tpl)
	return tpl

if __name__ == '__main__':
	if not QApplication.instance():
		app = QApplication(sys.argv)
	else:
		app = QApplication.instance()
	apply_stylesheet(app, theme='light_teal.xml', invert_secondary=True)
	w = Dielectric_params("src/Cpython/settings.json")
	w.resize(1280, 720)
	w.show()
	app.exec()

