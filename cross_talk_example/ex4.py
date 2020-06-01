#!/usr/bin/python3

from cross_talk import *
import math

def func1(time):
    return 10*(1-math.e**(-1e6 * time))

def func2(time):
    return 10* math.e**(-1e6 * (time-1e-5))

Rs = 50
Rne = 500
Rfe = 200
Rl = 100
L =  0.2 
cm = 20.3e-12  # F/m
c1 = 40.6e-12 - cm # F/m
c2 = 29.7e-12 - cm # F/m
circ = xtalk_circuit(Rs, Rl, Rne, Rfe)
circ.lm = 0.69e-6 #H/m
circ.l2 = 1.1e-6
circ.l1 = 1.4e-6
circ.set_c_parasitic(cm,c1,c2)
circ.L = L

circ.print_data_freq_elaboration(100e6)

circ.set_Vs_time_signal([func1,func2],[10e-6,20e-6])
circ.print_data_time_elaboration()
