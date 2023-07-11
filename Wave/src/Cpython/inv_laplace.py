import sympy as sp
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore
from PySide6.QtCore import Slot, QRectF
from scipy.constants import speed_of_light, epsilon_0, pi, mu_0
from mpmath import *

app = pg.mkQApp("Plotting Example")
#mw = QtWidgets.QMainWindow()
#mw.resize(800,800)

win = pg.GraphicsLayoutWidget(show=True, title="Basic plotting examples")
win.resize(1000,600)
win.setWindowTitle('pyqtgraph example: Plotting')


# z, p = sp.symbols('z, p')
# t = sp.Symbol('t', positive=True)

lamda = 1000
period = (lamda * 1e-9) / speed_of_light

sigma1 = 0.0665302480707485
sigma2 = 0.3012867941819113
eps1 = 2.3408999999983005
eps2 = 1.795599999954565

a1 = (mu_0 * sigma1 * (lamda * 1e-9) ** 2) / period
a2 = (mu_0 * sigma2 * (lamda * 1e-9) ** 2) / period
b1 = (mu_0 * epsilon_0 * eps1 * (lamda * 1e-9) ** 2) / period ** 2
b2 = (mu_0 * epsilon_0 * eps2 * (lamda * 1e-9) ** 2) / period ** 2

# sqrt(p ** 2 * b1 + p * a1) = sqrt(p ** 2 * b1 + p * a1)
# sqrt(p ** 2 * b2 + p * a2) = sqrt(p ** 2 * b2 + p * a2)

sigma = 0.3
delta_z = 2

m = 5
L1 = 0
L2 = 10


U1 = lambda p,z :(-exp(sqrt(p ** 2 * b1 + p * a1) ** 2 * sigma ** 2 / 4 + (z - 2 * L1 + delta_z) * sqrt(p ** 2 * b1 + p * a1)) * (b1 * p + a1) * erf((-sqrt(p ** 2 * b1 + p * a1) * sigma ** 2 + 2 * L1 - 2 * delta_z) / sigma / 2) + exp(sqrt(p ** 2 * b1 + p * a1) * (sqrt(p ** 2 * b1 + p * a1) * sigma ** 2 - 4 * delta_z + 4 * z) / 4) * (b1 * p + a1) * erf((-sqrt(p ** 2 * b1 + p * a1) * sigma ** 2 + 2 * delta_z - 2 * z) / sigma / 2) + exp(sqrt(p ** 2 * b1 + p * a1) * (sqrt(p ** 2 * b1 + p * a1) * sigma ** 2 - 4 * delta_z + 4 * z) / 4) * (b1 * p + a1) * erf((sqrt(p ** 2 * b1 + p * a1) * sigma ** 2 + 2 * L1 - 2 * delta_z) / sigma / 2) - exp(-sqrt(p ** 2 * b1 + p * a1) * (-sqrt(p ** 2 * b1 + p * a1) * sigma ** 2 / 4 + z - delta_z)) * (b1 * p + a1) * erf((sqrt(p ** 2 * b1 + p * a1) * sigma ** 2 + 2 * delta_z - 2 * z) / sigma / 2) + (erf((-sqrt(p ** 2 * b1 + p * a1) * sigma ** 2 + 2 * L1 - 2 * delta_z) / sigma / 2) * (b1 * p + a1) * (sqrt(p ** 2 * b2 + p * a2) - sqrt(p ** 2 * b1 + p * a1)) * exp(sqrt(p ** 2 * b1 + p * a1) ** 2 * sigma ** 2 / 4 + (m - 2 * L1 + delta_z) * sqrt(p ** 2 * b1 + p * a1) + sqrt(p ** 2 * b2 + p * a2) * (-2 * L2 + m)) - (b1 * p + a1) * (erf((-sqrt(p ** 2 * b1 + p * a1) * sigma ** 2 + 2 * delta_z - 2 * m) / sigma / 2) + erf((sqrt(p ** 2 * b1 + p * a1) * sigma ** 2 + 2 * L1 - 2 * delta_z) / sigma / 2)) * (sqrt(p ** 2 * b2 + p * a2) - sqrt(p ** 2 * b1 + p * a1)) * exp(sqrt(p ** 2 * b1 + p * a1) ** 2 * sigma ** 2 / 4 + (m - delta_z) * sqrt(p ** 2 * b1 + p * a1) + sqrt(p ** 2 * b2 + p * a2) * (-2 * L2 + m)) + erf((sqrt(p ** 2 * b1 + p * a1) * sigma ** 2 + 2 * delta_z - 2 * m) / sigma / 2) * (b1 * p + a1) * (sqrt(p ** 2 * b1 + p * a1) + sqrt(p ** 2 * b2 + p * a2)) * exp(sqrt(p ** 2 * b1 + p * a1) ** 2 * sigma ** 2 / 4 + (-m + delta_z) * sqrt(p ** 2 * b1 + p * a1) + sqrt(p ** 2 * b2 + p * a2) * (-2 * L2 + m)) + erf((-sqrt(p ** 2 * b1 + p * a1) * sigma ** 2 + 2 * L1 - 2 * delta_z) / sigma / 2) * (b1 * p + a1) * (sqrt(p ** 2 * b1 + p * a1) + sqrt(p ** 2 * b2 + p * a2)) * exp(sqrt(p ** 2 * b1 + p * a1) ** 2 * sigma ** 2 / 4 + (m - 2 * L1 + delta_z) * sqrt(p ** 2 * b1 + p * a1) - sqrt(p ** 2 * b2 + p * a2) * m) - (b1 * p + a1) * (erf((-sqrt(p ** 2 * b1 + p * a1) * sigma ** 2 + 2 * delta_z - 2 * m) / sigma / 2) + erf((sqrt(p ** 2 * b1 + p * a1) * sigma ** 2 + 2 * L1 - 2 * delta_z) / sigma / 2)) * (sqrt(p ** 2 * b1 + p * a1) + sqrt(p ** 2 * b2 + p * a2)) * exp(sqrt(p ** 2 * b1 + p * a1) ** 2 * sigma ** 2 / 4 + (m - delta_z) * sqrt(p ** 2 * b1 + p * a1) - sqrt(p ** 2 * b2 + p * a2) * m) + erf((sqrt(p ** 2 * b1 + p * a1) * sigma ** 2 + 2 * delta_z - 2 * m) / sigma / 2) * (b1 * p + a1) * (sqrt(p ** 2 * b2 + p * a2) - sqrt(p ** 2 * b1 + p * a1)) * exp(sqrt(p ** 2 * b1 + p * a1) ** 2 * sigma ** 2 / 4 + (-m + delta_z) * sqrt(p ** 2 * b1 + p * a1) - sqrt(p ** 2 * b2 + p * a2) * m) + 2 * sqrt(p ** 2 * b1 + p * a1) * (p * b2 + a2) * (exp(sqrt(p ** 2 * b2 + p * a2) * (sqrt(p ** 2 * b2 + p * a2) * sigma ** 2 - 4 * delta_z) / 4) * erf((-sqrt(p ** 2 * b2 + p * a2) * sigma ** 2 + 2 * delta_z - 2 * m) / sigma / 2) + exp(sqrt(p ** 2 * b2 + p * a2) * (sqrt(p ** 2 * b2 + p * a2) * sigma ** 2 - 4 * delta_z) / 4) * erf((sqrt(p ** 2 * b2 + p * a2) * sigma ** 2 + 2 * L2 - 2 * delta_z) / sigma / 2) - exp(sqrt(p ** 2 * b2 + p * a2) ** 2 * sigma ** 2 / 4 + (-2 * L2 + delta_z) * sqrt(p ** 2 * b2 + p * a2)) * erf((sqrt(p ** 2 * b2 + p * a2) * sigma ** 2 + 2 * delta_z - 2 * m) / sigma / 2) - exp(sqrt(p ** 2 * b2 + p * a2) ** 2 * sigma ** 2 / 4 + (-2 * L2 + delta_z) * sqrt(p ** 2 * b2 + p * a2)) * erf((-sqrt(p ** 2 * b2 + p * a2) * sigma ** 2 + 2 * L2 - 2 * delta_z) / sigma / 2))) / ((-sqrt(p ** 2 * b2 + p * a2) + sqrt(p ** 2 * b1 + p * a1)) * exp(sqrt(p ** 2 * b1 + p * a1) * (-2 * L1 + m) + sqrt(p ** 2 * b2 + p * a2) * (-2 * L2 + m)) + (-sqrt(p ** 2 * b2 + p * a2) - sqrt(p ** 2 * b1 + p * a1)) * exp(sqrt(p ** 2 * b1 + p * a1) * (-2 * L1 + m) - sqrt(p ** 2 * b2 + p * a2) * m) + (sqrt(p ** 2 * b1 + p * a1) + sqrt(p ** 2 * b2 + p * a2)) * exp(sqrt(p ** 2 * b2 + p * a2) * (-2 * L2 + m) - sqrt(p ** 2 * b1 + p * a1) * m) + exp(-m * (sqrt(p ** 2 * b1 + p * a1) + sqrt(p ** 2 * b2 + p * a2))) * (sqrt(p ** 2 * b2 + p * a2) - sqrt(p ** 2 * b1 + p * a1))) * (exp(-sqrt(p ** 2 * b1 + p * a1) * z) - exp(sqrt(p ** 2 * b1 + p * a1) * (-2 * L1 + z)))) / sqrt(p ** 2 * b1 + p * a1) / 4


y = lambda p: U1(p,2)
#res = ilt.invert(y,t,7,1)
print(invertlaplace(y,10,method='talbot'))
# func =  sp.sin(z*t)
# F = sp.laplace_transform(func, t, p)
# print(F)
#f = sp.inverse_laplace_transform(U1, p, t)
dx = 0.5
dt = 0.5
x_l = np.arange(L1, m + 1, dx, dtype=np.float64)
t_l = np.arange(1, 15, dt, dtype=np.float64)
u = np.zeros((len(t_l),len(x_l)), dtype=np.float64)
u_tmp = np.zeros(len(x_l), dtype=np.float64)
#print(x_l)
#print(t_l)

# for i,t in enumerate(t_l):
#     for j,x in enumerate(x_l):
#         y = lambda p: U1(p,x)
#         u_tmp[j] = invertlaplace(y,t,method='talbot')
#     u[i] = u_tmp
#     print(t)

# np.save('E:/Wave/src/Cpython/test2.npy', u)

u = np.load('E:/Wave/src/Cpython/test2.npy')
p1 = win.addPlot(title="test")
curve = p1.plot(pen={'color':(45, 202, 240), 'width':2})
ptr = 0
#print(u[ptr])
curve.setData(x_l, u[-1])

# def update():
#     global curve, u, ptr, x_l
#     curve.setData(x_l,u[ptr%10])
#     ptr += 1
# timer = QtCore.QTimer()
# timer.timeout.connect(update)
# timer.start(0)

pg.exec()