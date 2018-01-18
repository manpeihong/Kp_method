import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA

ks = np.arange(0, 0.05, 0.001)
k2 = 0.00001
k = np.sqrt(k2) * 3.7

a1 = 1
a2 = 0
a3 = 0
a = np.sqrt(a1 * a1 + a2 * a2 + a3 * a3)
energies = []
kx = k * a1 / a
ky = k * a2 / a
kz = k * a3 / a

L = -32
M = -5.3
N = -32.4
delta = 0.29

H11 = L * kx * kx + M * (ky * ky + kz * kz)
H12 = N * kx * ky
H13 = N * kx * kz
H21 = N * kx * ky
H22 = L * ky * ky + M * (kx * kx + kz * kz)
H23 = N * ky * kz
H31 = N * kx * kz
H32 = N * ky * kz
H33 = L * kz * kz + M * (kx * kx + ky * ky)

error = 1
Ek = 0
Eks = []
errors = []

for Ek in np.arange(0, 0.45, 0.001):
    _H11 = H11 + k * k + Ek
    _H22 = H22 + k * k + Ek
    _H33 = H33 + k * k + Ek

    total = _H11 * _H22 * _H33 + 2 * H12 * H23 * H13 - _H11 * H23 * H23 - _H22 * H13 * H13 - _H33 * H12 * H12 \
          - delta / 3 * (_H11 * _H22 + _H11 * _H33 + _H22 * _H33 - H12 * H12 - H13 * H13 - H23 * H23)
    Eks.append(Ek)
    errors.append(total)


print(errors)

plt.plot(Eks, errors)
plt.xlabel('Ek')
plt.ylabel('Error')
plt.title('About as simple as it gets, folks')
plt.grid(True)
# plt.savefig("test.png")
# plt.show()

a = np.matrix([[1, 0], [0, 1]])
b = np.matrix([[0, 1], [1, 0]])

Jx = np.matrix([[0, np.sqrt(3)/2, 0, 0], [np.sqrt(3)/2, 0, 1, 0], [0, 1, 0, np.sqrt(3)/2], [0, 0, np.sqrt(3)/2, 0]])
Jy = np.matrix([[0, -np.sqrt(3)/2 * 1j, 0, 0], [np.sqrt(3)/2 * 1j, 0, -1j, 0], [0, 1j, 0, -np.sqrt(3)/2 * 1j], [0, 0, np.sqrt(3)/2 * 1j, 0]])
Jz = np.matrix([[1.5, 0, 0, 0], [0, 0.5, 0, 0], [0, 0, -0.5, 0], [0, 0, 0, -1.5]])
J = Jx + Jy + Jz
J2 = np.dot(Jx, Jx) + np.dot(Jy, Jy) + np.dot(Jz, Jz)

H = [[a, b], [b, a]]

H2 = np.matrix([[1, 0, 0, 1], [0, 1, 1, 0], [0, 1, 1, 0], [1, 0, 0, 1]])
print(np.dot(Jx, Jy))

w, v = LA.eig(H)

print("*"*60)


hbar = 6.582e-16
k = 0.00e10
m = 9.11e-31
EP0 = 25.7
EP1 = 5.9
EQ = 13.5
P02 = m * EP0 / 3
P12 = m * EP1 / 3
Q2 = m * EQ / 3
P0 = np.sqrt(P02)
P1 = np.sqrt(P12)
Q = np.sqrt(Q2)
Deltac = 0.171
Delta = 0.341

Ec = 3.14
Eg = 1.519
Joues = 1.60218e-19

bracket = 1 / (Ec + Eg - Deltac) - 1 / (Ec + Eg)

result = hbar * hbar * k * k / np.sqrt(2) / m / m * Q * P1 * bracket * 1000 * Joues

HHfactor = hbar * hbar * k * k / 2 / m * Joues

hps11 = Eg + HHfactor
hps12 = 1j * P0 * hbar * k / m * np.sqrt(Joues)
hps13 = 0
hps14 = 1j * P0 * hbar * k / m / np.sqrt(2) * np.sqrt(Joues)
hps21 = - hps12
hps22 = HHfactor
hps23 = hps24 = hps31 = hps32 = 0
hps33 = HHfactor
hps34 = 0
hps41 = - hps14
hps42 = hps43 = 0
hps44 = - delta + HHfactor

hpsmatrix = np.array([[hps11, hps12, hps13, hps14], [hps21, hps22, hps23, hps24],
                      [hps31, hps32, hps33, hps34], [hps41, hps42, hps43, hps44]])

deltahps11 = HHfactor * P12 / m * (2 / (Eg - Ec - Eg) + 1 / (Eg - Ec - Eg + Deltac))
deltahps12 = 0
deltahps13 = HHfactor * 1j * P1 * Q / m * (np.sqrt(2) / (Eg - Ec - Eg) - np.sqrt(2) / (Eg - Ec - Eg + Deltac))
deltahps14 = 0
deltahps21 = 0
deltahps22 = HHfactor * Q2 / m / (0 - Ec - Eg)
deltahps23 = 0
deltahps24 = HHfactor * - Q2 / m / (- Ec - Eg) * np.sqrt(2)
deltahps31 = HHfactor * -1j * P1 * Q / m * (np.sqrt(2) / (0 - Ec - Eg) - np.sqrt(2) / (0 - Ec - Eg + Deltac))
deltahps32 = 0
deltahps33 = HHfactor * Q2 / m / ((0 - Ec - Eg) + 2 / (0 - Ec - Eg + Deltac))
deltahps34 = 0
deltahps41 = 0
deltahps42 = HHfactor * -np.sqrt(2) * Q2 / m / (-delta - Ec - Eg)
deltahps43 = 0
deltahps44 = HHfactor * Q2 / m / (-delta - Ec - Eg)

deltahpsmatrix = np.array([[deltahps11, deltahps12, deltahps13, deltahps14], [deltahps21, deltahps22, deltahps23, deltahps24],
                      [deltahps31, deltahps32, deltahps33, deltahps34], [deltahps41, deltahps42, deltahps43, deltahps44]])

w, v = LA.eig(hpsmatrix + deltahpsmatrix)

print(w)

