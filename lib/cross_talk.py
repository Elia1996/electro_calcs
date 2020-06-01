#!/usr/bin/env python3

from scipy.constants import mu_0
import scipy.constants as scik
import math
from engineering_notation import EngNumber
import sys
import numpy
from  plotter import *


def normalized_list(n_lista):
    mean = 0
    i=0
    for n in n_lista:
        if n!=0:
            mean += math.log10(math.fabs(n))
        i+=1
    mean=int(mean/i)
    while mean%3 != 0:
        mean+=mean/math.fabs(mean)
    return {"order":mean, "lista":[x/(10**mean) for x in n_lista ]}

def getSettedPlot():
    plot_style = { 'style' : 'Solarize_Light2', 
				'marker' : '.', 
				'linewidth' : 0.5,
				'color' : '#00a1dc'
			}
    objp=plotter(40,False,**plot_style)
    objp.title_font['size']=18
    objp.title_font['fontweight']='bold'
    objp.ticksize(12)
    objp.x_font['color']='#667B83'
    objp.x_font['size']='12'
    objp.y_font['size']='12'
    objp.y_font['color']='#667B83'
    objp.x_font['weight']='bold'
    objp.x_font['fontweight']='bold'
    objp.y_font['weight']='bold'
    objp.y_font['fontweight']='bold'
    return objp

def easyPlot(title, x, xlab, y, ylab):
    objp  = getSettedPlot()
    objp.plotXY(xlab, x, ylab, y, title)

def easyPlotSemilogx(title, x, xlab, y, ylab):
    objp  = getSettedPlot()
    objp.semilogxXY(xlab, x, ylab, y, title)


def iselement(name, **kwargs):
    if name in kwargs:
        return kwargs[name]
    else:
        sys.exit()

def islistelement(name_list, **kwargs):
    lista = []
    for i in name_list:
        lista.append(iselement(i, **kwargs))
    return lista

def ind_to_cap_pul(lm, l1, l2):
    """Calculate cm c1 and c2 from pul inductance
     _             _                         _        _
    | c1+cm    -cm  |            1          |  l2  -lm |
    |               | = --------------------| 
    |_ -cm    c2+cm_|   c^2 * (l1*l2 - lm^2)|_-lm   l1_|

    """
    den = (scik.c**2*(l1*l2 - lm**2))
    cm = lm / den
    c1 = ( l2 / den ) - cm
    c2 = ( l1 / den ) - cm

    return [cm, c1, c2]

def cap_to_ind_pul(cm, c1,c2):
    """Calculate l1, l2 and lm from pul capacitante
    according to the following formula:
     _             _                         _        _
    | c1+cm    -cm  |            1          |  l2  -lm |
    |               | = --------------------| 
    |_ -cm    c2+cm_|   c^2 * (l1*l2 - lm^2)|_-lm   l1_|
    But we should calculate inverse of both matrix, so we remember that
     _     _ -1          _      _
    | a   b |       1   |  d  -b | 
    |       |   = ----- |        |
    |_c   d_|     ab-bc |_-c   a_|

    And so 
     _             _                                                    _        _
    | c2+cm     cm  |           1                 c^2 * (l1*l2 - lm^2) |  l1   lm |
    |               | ----------------------   = --------------------- |          | 
    |_  cm    c1+cm_| (c1-cm)(c2-cm) - cm^2          (l1*l2 - lm^2)    |_ lm   l2_|
            |                                                                |
            |                                                                |
            |         -1             1                                       |
             ---->   C   --------------------------- = L  <------------------      
                         c^2 [(c1-cm)(c2-cm) - cm^2]
                                    |
                                    |        -1   1
                                     --->   C    --- = L
                                                  K 
    
    Son now we know how to calculate L from C matrix
    """
    K = scik.c**2  * ((c1+cm)*(c2+cm) - cm**2 )
    l1 = (c2+cm) / K
    l2 = (c1+cm) / K
    lm = cm / K 

    return [lm, l1, l2]


def lc_twag(**kwargs):
    [s, hr, rwr, hg, rwg] = islistelement(["s", "hr", "rwr", "hg", "rwg"], **kwargs)
    """Inductance capacitance two wires above grounds
    Electronic schema is the following:
                   s
            |--------------|
          .---.          .---.
      _  / rwr \        / rwg \ _   rw(r,g) are radius
      |  \     /        \     / |           of wire r,g
      |   '---'          '---'  |
    h1|                         |h2
      |                         |
    -------------------------------
    
    two_wires_above_grounds calculate:
    [mutual_inductance, inductance_wr_pul, inductance_wg_pul]
    where:
    mutual_inductance:  mutual inductance between the two wire
    inductance_wr_pul:  inductance per unit length of wire reference
    inductance_wg_pul:  inductance per unit length of wire g
    mutual_capacitance: mutual capacitance between the two wire
    capacitance_wr_pul: capacitance per unit length of wire reference
    capacitance_wg_pul: capacitance per unit length of wire g
    
    Ipotesys:
        NO insulation
        s sufficently higher then rwr and rwg
    Example:
        ./cross_talk.py F=twag S=0.02 HR=0.02 RWR=0.001 HG=0.02 RWG=0.001
    """
    inductance_wr_pul = (mu_0/(2*math.pi))*math.log(2*hr/rwr,math.e)
    inductance_wg_pul = (mu_0/(2*math.pi))*math.log(2*hg/rwg,math.e)
    mutual_inductance = (mu_0/(4*math.pi))*math.log(1+4*hr*hg/s**2)
    
    cap_list = ind_to_cap_pul(mutual_inductance, inductance_wg_pul, inductance_wg_pul)

    return [mutual_inductance, inductance_wr_pul, inductance_wg_pul] + cap_list

class xtalk_circuit:
    def __init__(self, Rs, Rl, Rne, Rfe, R0=0):
        self.Rs = Rs
        self.Rl = Rl
        self.Rne = Rne
        self.Rfe = Rfe
        self.R0 = R0
        self.parasitic_set=0
        self.L = 0 # length of wires
        [ self.lm, self.l1, self.l2, self.cm, self.c1, self.c2 ] = [0,0,0,0,0,0]

    #### Following function are different way to generate c1,c2,cm,l1,l2,lm
    def set_geometry(self, geometry, **kwargs):
        if geometry == "twag": # two wire above ground
           # s, hr, rwr, hg, rwg
           [self.lm, self.l1, self.l2, self.cm, self.c1, self.c2] = lc_twag(**kwargs)       
        else:
            print("Geometry not implemented")
            sys.exit()
        self.parasitic_set=1
        return 0

    def set_c_parasitic(self, cm, c1, c2 ):
        self.parasitic_set=1
        self.cm = cm
        self.c1 = c1
        self.c2 = c2
        if self.lm == 0 and self.l1 ==0 and self.l2 ==0:
            [self.lm, self.l1, self.l2] = cap_to_ind_pul(cm,c1,c2)

    def set_l_parasitic(self, lm, l1, l2):
        self.parasitic_set=1
        self.lm = lm
        self.l1 = l1
        self.l2 = l2
        if self.cm == 0 and self.c1 ==0 and self.c2 ==0:
            [self.cm, self.c1, self.c2] = ind_to_cap_pul(lm,l1,l2)
    ### Citcuit modify setting
    def set_shield_both_end_parameters(self, Rsh, Lsh, LGR):
        """
        Rsh is the radius of the shield
        Lsh is the total inductance of shield
        LGR inductance from generator wire to receptor wire
        """
        self.Rsh = Rsh
        self.Lsh = Lsh 
        self.LGR = LGR

    
    ### Time domain set
    def set_Vs_time_signal(self, func_list, validity_interval_list, step=1000):
        """
        func_list is a list of function int the form:
            func(time) 
        that return value of signal in Volt for a given instant of time

        validity_interval_list is the list of validity of previous function,
        if for example we have:
                  _
                 |  10^2 * e^(-10^6 *t )                0s <= t <= 10us
        Vne(t) = |
                 |_ -10^-2 * e^(-10^6(t-10^-5))       10us<= t<= 20us
        we will have:
            validity_interval_list = [10e-6, 20e-6]

        step are the number of point that will be sampled in period
        """
        self.ts =  validity_interval_list[-1]/step # sample time
        self.time = [ t*validity_interval_list[-1]/step for t in range(0,step+1) ] # list of time
        dizionario = normalized_list(self.time)
        order = dizionario['order']
        self.time_norm = dizionario['lista']
        self.time_order= order
        t=0
        i = 0 
        self.Vs = []
        for func in func_list:
            while self.time[t] < validity_interval_list[i]:
                self.Vs.append(func(self.time[t])) 
                t+=1
            i+=1
        self.Vs.append(func(self.time[t]))

    ### Plotting
    def time_plotting(self, title, signal_list, ylabelname):
        """Normalize y axes and plot
        """
        dizionario = normalized_list(signal_list)
        order = dizionario['order']
        norm_sig = dizionario['lista']
        ylabel = str(ylabelname+" [ "+str(EngNumber(10**dizionario['order']))+ "V]")
        easyPlot(title,self.time_norm, "time ["+str(EngNumber(10**self.time_order))+"s]", norm_sig, ylabel)

    def plot_Vs(self):
        self.time_plotting("Vs", self.Vs, "Voltage")
    
    def plot_d_Vs(self):
        try:
            self.d_Vs
        except AttributeError:
            self.Vs_derivative()
        self.time_plotting("Derivata di Vs",self.d_Vs, "Voltage")
    
    def plot_d_Vne(self):
        try:
            self.Vne_intime
        except AttributeError:
            self.Vne_Vfe_xtalk_time()
        self.time_plotting("Vne nel dominio del tempo", self.Vne_intime, "Voltage")
    
    def plot_d_Vfe(self):
        try:
            self.Vfe_intime
        except AttributeError:
            self.Vne_Vfe_xtalk_time()
        self.time_plotting("Vfe nel dominio del tempo", self.Vfe_intime, "Voltage")

    def plot_bode_Vne_sbe(self):
        """ Plot shielded both end"""
        diz = self.Vne_Vfe_xtalk_valid_set_shielded_both_end()
        easyPlotSemilogx("Bode of Vne, shilded circuit, Fsh="+str(EngNumber(self.Fsh))+"Hz",diz['f'], "freq [Hz]" ,diz['xtalk']['Vne'], "Voltage dB")


    ### time domain
    def Vs_derivative(self):
        """ Find the derivative of vs 
        """
        try:
            self.d_Vs=[]
            flag=0
            for t,vs in zip(self.time, self.Vs):
                if flag==0:
                    t_old = t
                    vs_old = vs
                    flag=1
                elif flag==1:
                    self.d_Vs.append( (vs-vs_old)/(t-t_old))
                    self.d_Vs.append( (vs-vs_old)/(t-t_old))
                    flag=2
                else:
                    self.d_Vs.append( (vs-vs_old)/(t-t_old))
                    t_old = t
                    vs_old = vs

        except NameError:
            print("Error in Vs_derivative, call get_Vs_time_signal before")
            sys.exit()

    def Vne_Vfe_xtalk_time(self):
        """Finf Vne and Vfe value in time domain and save it in:
        self.Vne_intime
        self.Vfe_intime
        """
        [ Mne_ind, Mne_cap,  Mfe_ind, Mfe_cap, dominating_effect] = self.xtalk()
        self.Vs_derivative()
        self.Vne_intime = [v*(Mne_ind+Mne_cap) for v in self.d_Vs]
        self.Vfe_intime = [v*(Mfe_ind+Mfe_cap) for v in self.d_Vs]

    ### Find xtalk shilded cable
    def xtalk_shielded_both_end(self):
        [ Rs, Rl, Rne, Rfe, Rsh, Lsh, LGR] = [self.Rs, self.Rl, self.Rne, self.Rfe, self.Rsh, self.Lsh, self.LGR]
        """
        In this case the receptor conductor have a shield connected to ground at both end
        so the capacitive effect is 0!!!!!
        remain only inductive effect so will be returned following list:
        [ Fsh, Mne_ind_LF,  Mfe_ind_LF, Mne_ind_HF, Mfe_ind_HF ]
        Fsh is the transition frequency 
        Mne_ind_LF,  Mfe_ind_LF  low frequency = Vne_ind/( jw * Vs ) f<Fsh  DEPEND from FREQUENCY   
        Mne_ind_HF,  Mfe_ind_HF  high frequency = Vne_ind/Vs   f>=Fsh  CONSTANT
        """
        try:
            self.Rsh
        except AttributeError:
            print("Error, call set_shield_both_end_parameters")
            sys.exit()
        Fsh = Rsh/(2*math.pi*Lsh)
        Mne_ind_LF = (Rne*LGR)/((Rne+Rfe)*(Rs+Rl))
        Mne_ind_HF = (Rne*LGR*Rsh)/((Rne+Rfe)*(Rs+Rl)*Lsh)
        Mfe_ind_LF = -(Rfe*LGR)/((Rne+Rfe)*(Rs+Rl))
        Mfe_ind_HF = -(Rfe*LGR*Rsh)/((Rne+Rfe)*(Rs+Rl)*Lsh)
        return [ Fsh, Mne_ind_LF,  Mfe_ind_LF, Mne_ind_HF, Mfe_ind_HF ]
    
    def Vne_Vfe_xtalk_dB_shielded_both_end(self, f_list):
        """ Finf Near-End and Far-End cross talk using setted constant
        for f < Fsh     Vne/Vs = jw*(Mne_ind_LF)  
                        Vfe/Vs = jw*(Mfe_ind_LF)
        fir f>= Fsh     Vne/Vs = Mne_ind_HF
                        Vfe/Vs = Mnfe_ind_HF
        """
        [ Fsh, Mne_ind_LF,  Mfe_ind_LF, Mne_ind_HF, Mfe_ind_HF ] = self.xtalk_shielded_both_end()
        Vne = []
        Vfe = []
        self.Fsh = Fsh
        for f in f_list:
            if f < Fsh:
                Vne.append(   20*math.log10( math.fabs( 2*math.pi* f * Mne_ind_LF  ) )   )
                Vfe.append(   20*math.log10( math.fabs( 2*math.pi* f * Mfe_ind_LF ) )   )
            else:
                Vne.append(   20*math.log10( math.fabs( Mne_ind_HF ) )  )
                Vfe.append(   20*math.log10( math.fabs( Mfe_ind_HF ) )  )
        
        return { "Vne":Vne, "Vfe":Vfe }
    
    def Vne_Vfe_xtalk_valid_set_shielded_both_end(self, fmin=10, fstep=100):
        """Return a list of xtalk value ranging in valid frequency
        """
        f_list = list(numpy.logspace(math.log10(fmin), math.log10(int(self.freq_validity())), fstep))
        value = self.Vne_Vfe_xtalk_dB_shielded_both_end(f_list)
        return { "f": f_list, "xtalk": value }


    ### Find xtalk  non shilded

    def xtalk(self):
        if self.parasitic_set==0:
            print("Error, call set_geometry before call xtal!!")
            sys.exit()

        [ Rs, Rl, Rne, Rfe, Lm, Cm, L] = [self.Rs, self.Rl, self.Rne, self.Rfe, self.lm, self.cm, self.L]
        """
        Returned list:
            [Mne_ind   Mne_cap  Mfe_ind  Mfe_cap  dominating_effect]

        This function take as reference the following circuit:

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
                        |-----------------------------------------|---------------------> z
                        z=0                                       z=L
            
        The resulting equivalent circuit is the following:
                        
                    jw Lm IGcd  
                        ----
            ┌----------|+  -|----┬-----------------------┐
            |   +       ----    _|_                  +   |
            \                  / ^ \ jw Cm VGdc          /
         Rne/  Vne             \_|_/                Vfe  \ Rfe
            \                    |                       /
            |   -                |                   -   |
            └--------------------┴-----------------------┘
        So starting from these equation:ù
        
                Rne                       Rne*Rfe
        Vne = ----------- jw Lm IGdc   +  ----------- jmw Cm VGdc
            Rne + Rfe                   Rne + Rfe
            
                    Rfe                       Rne*Rfe
        Vfe = - ----------- lw Lm IGdc   +  ----------- lw Cm VGdc
                Rne + Rfe                   Rne + Rfe 

            ( inductive coupling )      ( capacitive coupling )

        Using equation of IGdc = Vs/(Rl+Rs) and VGdc = Vs*Rl/(Rs+Rl) we find:

        Vne = Vs * jw * [ Mne_ind + Mne_cap ]

        Vfe = Vs * jw * [ Mfe_ind + Mfe_cap ]

        Where:
                    Rne        Lm                             Rne*Rfe      Rl*Cm
        Mne_ind = ─────────── ─────────              Mne_cap = ─────────── ─────────
                    Rne + Rfe   Rs + Rl                         Rne + Rfe    Rs + Rl 
    

                        Rfe       Lm                             Rne*Rfe     Rl*Cm
        Mfe_ind = - ─────────── ───────              Mfe_cap = ─────────── ──────── = Mne_cap
                    Rne + Rfe  Rs + Rl                         Rne + Rfe   Rs + Rl 
        
        This function calculate previous coefficient, but olso find if we have a dominating capacitive
        or indictive effect
                                                            if
        Inductive Effect is dominating if M_ind > M_cap    --->  (Rfe*Rl)/(Lm/Cm) < 1
                                                                (Rne*Rl)/(Lm/Cm) < 1

        Capacitive Effect is dominating if M_cap > M_ind   so previous condition are false

        """             
        Mne_ind = (Rne*Lm*L)/((Rne+Rfe)*(Rs+Rl))
        Mne_cap = (Rne*Rfe*Rl*Cm*L)/((Rne+Rfe)*(Rs+Rl))
        Mfe_ind = - Mne_ind*Rfe/Rne
        Mfe_cap = Mne_cap

        if (Mfe_cap > Mfe_ind and Mne_cap > Mfe_cap):
            return [Mne_ind, Mne_cap,  Mfe_ind, Mfe_cap, "capacitive"]
        elif (Mfe_cap < Mfe_ind and Mne_cap < Mfe_cap):
            return [Mne_ind, Mne_cap,  Mfe_ind, Mfe_cap, "inductive"]
        else:
            return [Mne_ind, Mne_cap,  Mfe_ind, Mfe_cap, "similar"]

    def cw_xtalk(self):
        [ Rne, Rfe, Rs, Rl, R0] = [self.Rne, self.Rfe, self.Rs, self.Rl, self.R0]
        """We in this case consider circuit with R0 as resistance of reference conductor
                                Generator Conductor
            ┌----/\/\/\------------------------>--------------------------┐
            |      Rs                         Ig(z,f)       +             |
            |                     Receptor Conductor      Vg(z,f)         |
           _|          ┌------------>----------------------------┐        |
    Vs(f) / + \        |  +      Ir(z,f)   +                  +  |        |
          \_-_/        \                                         /        \
            |      Rne /  Vne            Vr(z,f)             Vfe \ Rfe    / Rl
            |          \                                         /        \
            |          |  -        R0       -      Ig+Ir      -  |        |
            └----------┴---------/\/\/\------------->------------┴--------┘
                                -   V0   +     Refecence conductor
        V0 = R0*Ig = Vs * R0/(Rs+Rl)
    
        So:
                        Rne * R0                                  Rfe * R0
        Mne_cw = ─────────────────────           Mfe_cw = - ──────────────────────
                (Rne + Rfe)*(Rs + Rl)                      (Rne +  Rfe)*(Rs + Rl)

        where:
        
        Vne = Vs * Mne_cw
        Vfe = Vs * Mfe_cw

        """
        Mne_cw = (Rne * R0)/((Rne+Rfe)*(Rs+Rl))
        Mfe_cw = (Rfe * R0)/((Rne+Rfe)*(Rs+Rl))

        return [Mne_cw, Mfe_cw]

    def freq_validity(self):
        """Return mazimum frequency for which xtalk analysis is valid
        L < lambda/10 
        lambda = c/f
        so:
        fmax << c/L
        """
        if self.L==0:
            return 0
        return scik.c/(10*self.L)

    def Vne_Vfe_xtalk(self, f_list):
        """ Finf Near-End and Far-End cross talk using setted constant
        Vne/Vs = jw*(Mne_ind + Mne_cap) + Mne_cw
        Vfe/Vs = jw*(Mfe_ind + Mfe_cap) + Mfe_cw
        """
        [ Mne_ind, Mne_cap,  Mfe_ind, Mfe_cap, dominating_effect] = self.xtalk()
        [ Mne_cw, Mfe_cw ] = self.cw_xtalk()
        Vne = []
        Vfe = []
        for f in f_list:
            Vne.append(  2*math.pi* f * (Mne_ind+Mne_cap) + Mne_cw    )
            Vfe.append(  2*math.pi* f * (Mfe_ind+Mfe_cap) + Mfe_cw    )
     
        return { "Vne":Vne, "Vfe":Vfe }
    
    def Vne_Vfe_xtalk_dB(self, f_list):
        """ Finf Near-End and Far-End cross talk using setted constant
        Vne/Vs = jw*(Mne_ind + Mne_cap) + Mne_cw
        Vfe/Vs = jw*(Mfe_ind + Mfe_cap) + Mfe_cw
        """
        [ Mne_ind, Mne_cap,  Mfe_ind, Mfe_cap, dominating_effect] = self.xtalk()
        [ Mne_cw, Mfe_cw ] = self.cw_xtalk()
        Vne = []
        Vfe = []
        for f in f_list:
            Vne.append(   20*math.log10( math.fabs( 2*math.pi* f * (Mne_ind+Mne_cap) + Mne_cw ) )   )
            Vfe.append(   20*math.log10( math.fabs( 2*math.pi* f * (Mfe_ind+Mfe_cap) + Mfe_cw ) )   )
        
        return { "Vne":Vne, "Vfe":Vfe }

    def Vne_Vfe_xtalk_valid_set(self, fmin=10, fstep=100):
        """Return a list of xtalk value ranging in valid frequency
        """
        f_list =[f for f in range(fmin, int(self.freq_validity()), fstep) ]
        value = Vne_Vfe_xtalk(f_list)
        return { "f": f_list, "xtalk": value }


    def print_data_freq_elaboration(self, xtalk_freq):
        print("Not shilded circuit!!")
        print("lm: {0}, l1: {1}, l2:{2}".format(EngNumber(self.lm), EngNumber(self.l1), EngNumber(self.l2)))
        print("cm: {0}, c1: {1}, c2:{2}".format(EngNumber(self.cm), EngNumber(self.c1), EngNumber(self.c2)))
        [Mne_ind, Mne_cap,  Mfe_ind, Mfe_cap, dominating_effect] = self.xtalk()
        print("Mne_ind: {0}         Mne_cap:{1}  \nMfe_ind:{2}       Mfe_cap:{3}\nDominating_effect:{4}".\
                format(EngNumber(Mne_ind), EngNumber(Mne_cap),  EngNumber(Mfe_ind),EngNumber(Mfe_cap), dominating_effect))
        print("frequence validity: {0}".format(EngNumber(self.freq_validity())))
        xtalk = self.Vne_Vfe_xtalk([xtalk_freq])
        print("Linear Vne and Vfe values")
        print("     Vne at {1}Hz: {0}".format(EngNumber(xtalk["Vne"][0]), EngNumber(xtalk_freq)))
        print("     Vfe at {1}Hz: {0}".format(EngNumber(xtalk["Vfe"][0]), EngNumber(xtalk_freq)))
        xtalk = self.Vne_Vfe_xtalk_dB([xtalk_freq])
        print("Vne and Vfe values n dB")
        print("     Vne at {1}Hz: {0}dB".format(EngNumber(xtalk["Vne"][0]), EngNumber(xtalk_freq)))
        print("     Vfe at {1}Hz: {0}dB".format(EngNumber(xtalk["Vfe"][0]), EngNumber(xtalk_freq)))
    
    def print_data_time_elaboration(self):
        try:
            self.time
        except AttributeError:
            print("Error, Call set_Vs_time_signal before")
            sys.exit()
        print("Remember that time axes have order : {0}".format(EngNumber(10**self.time_order)))
        self.plot_Vs()
        self.plot_d_Vs()
        self.plot_d_Vne()
        self.plot_d_Vfe()

#from pwn import *
#
#if args['F']:
#    if args.F == "twag": # two wire above ground
#        if args['S'] and args['HR'] and args['RWR'] and args['HG'] and args['RWG']:
#            listas = two_wires_above_grounds(float(args.S), 
#                                            float(args.HR), 
#                                            float(args.RWR),
#                                            float(args.HG), 
#                                            float(args.RWG))
#            print("lm: {0}H/m\nlr: {1}H/m\nlg: {2}H/m\n"
#                    .format(EngNumber(listas[0]),EngNumber(listas[1]),EngNumber(listas[2])))
#            print("cm: {0}F/m\ncr: {1}F/m\ncg: {2}F/m\n"
#                    .format(EngNumber(listas[3]),EngNumber(listas[4]),EngNumber(listas[4])))
#
#
