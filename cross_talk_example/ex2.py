#!/usr/bin/python3

from cross_talk import *

# xtalk of unshielded cable
Rs = 0
Rne = 50
Rfe = 50
Rl = 50
L =  5 
cm = 18.7e-12  # F/m
c1 = 37.4e-12 - cm # F/m
c2 = 24.98e-12 - cm # F/m
circ = xtalk_circuit(Rs, Rl, Rne, Rfe)
circ.set_c_parasitic(cm,c1,c2)
circ.L = 4
# calculus of Vne and Vfe from previous data
circ.print_data_freq_elaboration(1e6)
