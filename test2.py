import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA

ks = np.arange(0, 0.2, 0.001)
ks = ks * 1e10

C_band = []
LH_band = []
HH_band = []
SO_band = []
C_SO_band = []
C_H_band = []
C_L_band = []

hbar = 6.582e-16
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

for k in ks:
    HHfactor = hbar * hbar * k * k / 2 / m * Joues
    Hfactor = hbar * k / m * np.sqrt(Joues)

    hps11 = Eg + HHfactor
    hps12 = 1j * P0 * Hfactor
    hps13 = 0
    hps14 = 1j * P0 / np.sqrt(2) * Hfactor
    hps21 = - hps12
    hps22 = HHfactor
    hps23 = hps24 = hps31 = hps32 = 0
    hps33 = HHfactor
    hps34 = 0
    hps41 = - hps14
    hps42 = hps43 = 0
    hps44 = - Delta + HHfactor

    hps_matrix = np.array([[hps11, hps12, hps13, hps14], [hps21, hps22, hps23, hps24],
                          [hps31, hps32, hps33, hps34], [hps41, hps42, hps43, hps44]])

    hp11 = hp22 = hp44 = hp55 = Ec + Eg
    hp33 = hp66 = Ec + Eg - Deltac

    Hp = np.array([[hp11, 0, 0, 0, 0, 0], [0, hp22, 0, 0, 0, 0], [0, 0, hp33, 0, 0, 0],
                         [0, 0, 0, hp44, 0, 0], [0, 0, 0, 0, hp55, 0], [0, 0, 0, 0, 0, hp66]])

    hpsp11 = 1j * P1 * Hfactor
    hpsp12 = 0
    hpsp13 = 1j * P1 / np.sqrt(2) * Hfactor
    hpsp21 = 0
    hpsp22 = - Q / np.sqrt(2) * Hfactor
    hpsp23 = 0
    hpsp31 = Q / np.sqrt(2) * Hfactor
    hpsp32 = 0
    hpsp33 = - Q * Hfactor
    hpsp41 = 0
    hpsp42 = Q * Hfactor
    hpsp43 = 0

    hpsp_matrix = np.array([[hpsp11, hpsp12, hpsp13], [hpsp21, hpsp22, hpsp23],
                            [hpsp31, hpsp32, hpsp33], [hpsp41, hpsp42, hpsp43]])

    hpsp_matrix_conjT = hpsp_matrix.conj().T

    zero_matrix44 = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    hps_temp1 = np.concatenate((hps_matrix, zero_matrix44), axis=1)
    hps_temp2 = np.concatenate((zero_matrix44, hps_matrix), axis=1)
    Hps = np.concatenate((hps_temp1, hps_temp2), axis=0)

    zero_matrix43 = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
    hpsp_temp1 = np.concatenate((hpsp_matrix, zero_matrix43), axis=1)
    hpsp_temp2 = np.concatenate((zero_matrix43, hpsp_matrix), axis=1)
    Hpsp = np.concatenate((hpsp_temp1, hpsp_temp2), axis=0)

    Hpsp_conjT = Hpsp.conj().T

    H_temp1 = np.concatenate((Hps, Hpsp), axis=1)
    H_temp2 = np.concatenate((Hpsp_conjT, Hp), axis=1)
    H = np.concatenate((H_temp1, H_temp2), axis=0)

    w, v = LA.eig(H)

    w = sorted(w)

    SO_band.append(w[0].real)
    LH_band.append(w[2].real)
    HH_band.append(w[4].real)
    C_band.append(w[6].real)
    C_SO_band.append(w[8].real)
    C_H_band.append(w[10].real)
    C_L_band.append(w[12].real)

plt.plot(ks, SO_band)
plt.plot(ks, LH_band)
plt.plot(ks, HH_band)
plt.plot(ks, C_band)
plt.plot(ks, C_SO_band)
plt.plot(ks, C_H_band)
plt.plot(ks, C_L_band)

plt.xlabel('k')
plt.ylabel('E')
plt.title('14band K.p')
plt.grid(True)
# plt.savefig("test.png")
plt.show()
