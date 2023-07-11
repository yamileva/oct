import pyqtgraph.examples
pyqtgraph.examples.run()



import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore
from PySide6.QtCore import Slot,QRectF
from scipy.constants import speed_of_light, epsilon_0, pi,mu_0

app = pg.mkQApp("Plotting Example")
#mw = QtWidgets.QMainWindow()
#mw.resize(800,800)

win = pg.GraphicsLayoutWidget(show=True, title="Basic plotting examples")
win.resize(1000,600)
win.setWindowTitle('pyqtgraph example: Plotting')

lamda = 1400
period = (lamda * 1e-9) / speed_of_light
omega = 2 * pi * speed_of_light *  period / (lamda * 1e-9)

lx = 460
r_t = (2 * lx * 1e-6) /(speed_of_light * period)
tau = (5 * lamda * 1e-9)/(speed_of_light * period)
print("t Импульса: ",tau)
delta = tau/4
print("Delta : ",delta)
print("Omega: ",omega)
dt = 1/50
tList = np.arange(0,10,dt, dtype=np.float64)
x = np.zeros(len(tList), dtype=np.float64)
#print(tList)
def E(t):
    if t <= delta:
        return t/delta
    elif t > delta and t <= tau - delta:
        return 1
    elif t > tau - delta and t <= tau:
        return 1 - (t - tau + delta) * 1/delta
    elif t > tau:
        return 0
    
for i,t in enumerate(tList):
    x[i] =  np.sin(omega*t) * E(t)

p1 = win.addPlot(title="test")
curve = p1.plot(pen={'color':(45, 202, 240), 'width':2})
curve.setData(tList,x)


if __name__ == '__main__':
    pg.exec()


#import numpy as np

# crd = np.array([[2.,8.],[2.,4.]], dtype=np.float64)
# pml = np.array(np.meshgrid(crd[0], crd[1])).T.reshape(-1,2)

# print(pml)


# data = np.array([[1,3,5,5,8,1,2,3],[1,3,5,5,8,1,2,3]], dtype=np.float64)
# x_size = 2
# y_size = 4
# wave = data[0].reshape(x_size, y_size)
# print(wave)
# wave = data[0].reshape(x_size, y_size)
# wave = wave.reshape(x_size * y_size)
# print(wave)


# import numpy as np
# import pyqtgraph as pg
# from pyqtgraph.Qt import QtCore
# from PySide6.QtCore import Slot,QRectF
# import h5py
# pg.setConfigOption('background', 'w')
# pg.setConfigOption('foreground', 'k')

# file1 = h5py.File("D:/USATU/OCT/Wave/Data/wave2d_g.h5", 'r')
# file2 = h5py.File("D:/USATU/OCT/Wave/Data/wave2d.h5", 'r')

# wave_1 = file1['wave']
# wave_2 = file2['wave']

# grid_settings = file1['grid_settings'][0]
# pml_crd = np.column_stack((file1['pml_left'][0],file1['pml_right'][0]))
# skin_crd = file1['skin_crd'][0]

# h = {'time': grid_settings[0], 'x':  grid_settings[1],'y':  grid_settings[2]}
# left_borders = {'time': grid_settings[3], 'x':  grid_settings[4], 'y':  grid_settings[5]}
# right_borders = {'time': grid_settings[6], 'x':  grid_settings[7], 'y':  grid_settings[8]}
# steps = {'time': (right_borders['time'] - left_borders['time'])/h['time'],
#                 'x':  (right_borders['x'] - left_borders['x'])/h['x'],
#                 'y':  (right_borders['y'] - left_borders['y'])/h['y']}


# x = 4.5
# y = 2
# scan_1 = np.zeros(len(wave_1), dtype=np.float64)
# for i in range(len(wave_1)):
#     data = wave_1[i].reshape(int(steps['x'] + 1), int(steps['y'] + 1))
#     scan_1[i] = data[int((x-left_borders['x'])/h['x'])][int((y-left_borders['y'])/h['y'])]

# scan_2 = np.zeros(len(wave_2), dtype=np.float64)
# for i in range(len(wave_2)):
#     data = wave_2[i].reshape(int(steps['x'] + 1), int(steps['y'] + 1))
#     scan_2[i] = data[int((x-left_borders['x'])/h['x'])][int((y-left_borders['y'])/h['y'])]

# dif = scan_2 - scan_1

# app = pg.mkQApp("Plotting Example")
# #mw = QtWidgets.QMainWindow()
# #mw.resize(800,800)

# win = pg.GraphicsLayoutWidget(show=True, title="Basic plotting examples")
# win.resize(1000,600)
# win.setWindowTitle('pyqtgraph example: Plotting')

# t = np.linspace(left_borders['time'],right_borders['time'],len(dif))
# ##
# p1 = win.addPlot(title="scan_1 with Gauss")
# curve1 = p1.plot(pen={'color':(45, 202, 240), 'width':2})
# curve1.setData(t,scan_1)
# p1.addItem(curve1)
# win.nextRow()
# ##
# p2 = win.addPlot(title="scan_2")
# curve2 = p2.plot(pen={'color':(45, 202, 240), 'width':2})
# curve2.setData(t,scan_2)
# p2.addItem(curve2)
# win.nextRow()
# ##
# p3 = win.addPlot(title="scan_2 - scan_1")
# curve3 = p3.plot(pen={'color':(45, 202, 240), 'width':2})
# curve3.setData(t,dif)
# p3.addItem(curve3)
# win.nextRow()

# if __name__ == '__main__':
#     pg.exec()
