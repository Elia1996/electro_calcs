# electro_calcs
Install module using install.sh file, this file download plotting repository 
and set PYTHONPATH variable in ~/.bashrc in order to be able to use
this module everywere in terminal
## cross_talk.py script
		The main implementation in this script is the following class:

### xtalk_circuit
		This class consider the following circuit as base for all formulas and calculations:

									Generator Conductor
				┌----/\/\/\------------------------>--------------------------┐
				|      Rs                         Ig(z,f)       +             |
				|                     Receptor Conductor      Vg(z,f)         |
				|          ┌------------>----------------------------┐        |
		Vs(f) / + \        |  +      Ir(z,f)   +                  +  |        |
			  \_-_/        \                                         /        \
				|      Rne /  Vne            Vr(z,f)             Vfe \ Rfe    / Rl
				|          \                                         /        \
				|          |  -                -      Ig+Ir       -  |        |
				└----------┴---------------------------->------------┴--------┘
									Refecence conductor
							|-----------------------------------------|
												L
							|-----------------------------------------|------> z
							z=0                                       z=L

		In order to define all variable these functions are provided
		- set_c_parasitic(cm, c1, c2)
		- set_l_parasitic(lm, l1, l2)
		- set_shield_both_end_parameters( Rsh, Lsh, LGR)

		
		Starting from this basic circuit we have two possible modification,
		1) adding a shield and the case in which the reference has a 
		considerably resistence R0, so the main three functions are:

		- xtalk() -> calculate capacitance and inductance parameter for main circuit
			- Vne_Vfe_xtalk(f_list) -> compute Vne/Vs and Vfe/Vs using xtalk data
			- Vne_Vfe_xtalk(f_list) -> compute Vne/Vs and Vfe/Vs using xtalk data in dB
			- Vne_Vfe_xtalk_valid_set -> Return a list of xtalk value ranging in valid frequency
			- print_data_freq_elaboration -> UEFULLL!!! print a complete sreenshot of circuit defined
		- cw_xtalk() -> calulate capacitance and inductance parameter for circuit with R0
		- xtalk_shielded_both_end() -> calulate capacitance and inductance parameter 
			for both ended shilded receptor conductor
			- Vne_Vfe_xtalk_dB_shielded_both_end(f_list) -> find Vne/Vs and Vfe/Vs for shilded cable
			- Vne_Vfe_xtalk_valid_set_shielded_both_end() -> Return a list of xtalk value 
				ranging in valid frequency

		There is also some function useful to plot the corrisponding cross talk
		created by a signal in time, so passing in time domain we use derivative.
		The first function to call is:

		- set_Vs_time_signal(self, func_list, validity_interval_list, step=1000)
				this fuction receive a Vs defined in many interval and described
				by the function list. So an example of this could be:
						def func1(time):
							return 10*(1-math.e**(-1e6 * time))

						def func2(time):
							return 10* math.e**(-1e6 * (time-1e-5))
						circ = xtalk_circuit(Rs, Rl, Rne, Rfe) # definition of circuit
						circ.set_c_parasitic(cm,c1,c2)
						circ.L = L
				----->	circ.set_Vs_time_signal([func1,func2],[10e-6,20e-6])
						circ.print_data_time_elaboration() 
		
		After we could use simply:
		- print_data_time_elaboration()
		That print all plot: Vs, d_Vs, d_Vne, d_Vfe, in  order to see the 
		signal added to receptor conductor


### other function
	Two important function to transform l1,l2,lm to c1,c2,cm and viceversa are:
	- ind_to_cap_pul(lm, l1, l2)
	- cap_to_ind_pul(cm, c1,c2)
	While a function to calculate Inductance capacitance for two wires above grounds
	- lc_twag(**kwargs)
