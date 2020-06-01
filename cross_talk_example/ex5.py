#!/usr/bin/python3

from cross_talk import *
import math

def func1(time):
    return 10*(1-math.e**(-1e6 * time))

def func2(time):
    return 10* math.e**(-1e6 * (time-1e-5))

Rs = 50
Rne = 50
Rfe = 50
Rl = 50
Rsh = 0.0898
Lsh = 2.48e-6
LGR = 1.98e-6
cm = 20.3e-12  # F/m
c1 = 40.6e-12 - cm # F/m
circ = xtalk_circuit(Rs, Rl, Rne, Rfe)
circ.L = 0.2
circ.set_shield_both_end_parameters(Rsh, Lsh, LGR)
circ.plot_bode_Vne_sbe()
