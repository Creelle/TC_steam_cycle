## Project LMECA2150-Thermal cycle
#Template for the Steam turbine
#
# Author: Paolo Thiran & Gauthier Limpens
# Version: 2020
#
# Students can modify this script.
# However, the arguments for the ST defined in ST_arguments.py CANNOT be modified



import numpy as np;
import ST_arguments as ST_arg;
import matplotlib.pyplot as plt
from pyXSteam.XSteam import XSteam # see documentation here: https://pypi.org/project/pyXSteam/

# Add your packages here:



def psychrometrics(Tdb,absolute_humidity):
    """
        psychrometrics gives the dew points temperature [°C] for given inputs:
         - dry bulb temperature (Tdb [°C])
         - absolute humidity ( absolute_humidity [kg_water/kg_dry_air]).
        Equations and data based on ASHRAE 2013 fundamentals.
        Example: psychrometrics(30,0.01) gives 14.07°C
    """
    P_atm = 101.325;#kPa
    Pw=absolute_humidity*P_atm/(0.621945+absolute_humidity); #partial pressure of water wapor
    alpha=np.log(Pw);
    Tdp=6.54 + 14.526*alpha+0.7389*(alpha**2)+0.09486*(alpha**3)+0.4569*(Pw**0.1984); # valid for Tdp between 0 C and 93 C
    return Tdp;

def ST(ST_inputs):
    """
     ST Steam power plants modelisation
     ST computes the thermodynamics states for a Steam
     power plant (combustion, exchanger, cycle) turbine based on several
     inputs (imposed in ST_arguments) and based on a given electricity production P_e.
     It returns the main results. It can as well plots graphs if input
     argument DISPLAY = true (<=> DISPLAY=1)

     INPUTS (some inputs can be dependent on others => only one of these 2 can
             be activated). Refer to Fig 2.33 from reference book (in english)
     P_E = electrical power output target [kW]
     nsout     [-] : Number of feed-heating
     reheat    [-] : Number of reheating
     T_max     [°C] : Maximum steam temperature
     T_cond_out[°C] : Condenser cold outlet temperature
     p3_hp     [bar] : Maximum pressure
     drumFlag  [-] : if =1 then drum if =0 then no drum.
     eta_mec   [-] : mechanic efficiency of shafts bearings
     comb is a class containing combustion data :
           -comb.Tmax     [°C] : maximum combustion temperature
           -comb.Lambda   [-] : air excess
           -comb.x        [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
           -comb.y        [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
     T_exhaust [°C] : Temperature of exhaust gas out of the chimney
     p4       [bar] : High pressure after last reheating
     x6        [-] : Vapor ratio [gaseous/liquid] (in french : titre)
     T_0       [°C] : Reference temperature
     T_ext     [°C] : External temperature (atmospheric)
     TpinchSub [°C] : Temperature pinch at the subcooler
     TpinchEx  [°C] : Temperature pinch at a heat exchanger
     TpinchCond[°C] : Temperature pinch at condenser
     TpinchHR  [°C] : Temperature pinch at the heat recovery system
     Tdrum     [°C] : minimal drum temperature
     eta_SiC    [-] : Internal pump efficiency
     eta_SiT    [-] : Isotrenpic efficiency for Turbine. It can be a vector of 2 values :
                                eta_SiT(1)=eta_SiT_HP,eta_SiT(2)=eta_SiT_others
     DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then do not plot.
     """
    steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS); # m/kg/sec/°C/bar/W
    arg_in = ST_inputs;

    ## Check input arguments
    # ======================
    # if value == -1 <=> value not defined. Thus, it will be defined here:
    Pe = arg_in.Pe;
    if Pe ==-1.:
        Pe = 250e3;#250 kWe
    nsout = arg_in.nsout;
    if nsout ==-1.:
        nsout = 1;#15°C
    reheat = arg_in.reheat
    if reheat == -1:
        reheat = 1.;# Number of reheating
    T_max = arg_in.T_max;
    if T_max == -1.:
        T_max = 520 #°C
    T_cond_out = arg_in.T_cond_out;
    if T_cond_out == -1.:
        T_cond_out = 30#°C
    p3_hp = arg_in.p3_hp;
    if p3_hp == -1.:
        p3_hp =40 #bar
    eta_mec = arg_in.eta_mec;
    if eta_mec == -1.:
        eta_mec = 0.95 #[-]

    comb = arg_in.combustion
    if comb.Tmax == -1:
        comb.Tmax = 1400#°C
    if comb.Lambda == -1:
        comb.Lambda = 2;
    if comb.x == -1:
        comb.x = 0;
    if comb.y == -1:
        comb.y = 4;

    T_exhaust = arg_in.T_exhaust;
    if T_exhaust == -1.:
        T_exhaust = 600 #°C
    p4 = arg_in.p4;
    if p4 == -1.:
        p4 = 0.0503 #bar
    x6 = arg_in.x6;
    if x6 == -1.:
        x6 = 0.91 #[-]
    T0 = arg_in.T_0;
    if T0 == -1.:
        T0 = 15#°C
    T_ext = arg_in.T_ext;
    if T_ext ==-1.:
        T_ext = 15;#15°C

    #Bon il reste a faire TpinchSub,TpinchEx, TpinchCond,TpinchHR,Tdrum

    eta_SiC = arg_in.eta_SiC;
    if eta_SiC == -1.:
        eta_SiC = 0.89
    eta_SiT = arg_in.eta_SiT;
    if eta_SiT == -1.:
        eta_SiT = 0.89

    # Put all temperatures in Kelvin
    T_max,T_cond_out, T_exhaust, T0, T_ext = T_max+273.15,T_cond_out+273.15, T_exhaust+273.15, T0+273.15, T_ext+273.15 #K
    ## cycle definition
    # =================
    # Your job
    # YOUR WORK HERE
    # 1 basic rankine (4 points)
    # 2 Add reheating
    # 3 Add bleedings (Heat exchangers are a hard part)
    # 4 Add drum (Iterative process required)
    # 3bis  (in //), implement a combustion for CHyOx
    # 3tris (in //), implement heat recovery (see component 'Ra' in Fig 2.2 of reference book (english version)).
    #        To help you, answer the following questions (approximately):
    #        With several bleedings, what is the water temperature entering the boiler (T_2)?
    #        What is the highest combustion temperature allowed? Why? (T_max)
    #       Without a heat recovery, how much energy will I loose?
    #           For this point, if we assume that cp is constant => loss = cp * (T_2 - T_ext).
    #                           the useful energy was => useful = cp * (T_3 - T_2).
    #                           are the losses negligible compare to the useful energy?
    #        As a conclusion, what can we do to decrease this loss?

    """
    1) Pump
    """
    T1= T_cond_out
    p1 = steamTable.psat_t(T1-273.15)
    h1= steamTable.hL_p(p1)
    s1= steamTable.sL_p(p1)
    x1 = 0
    v1 = steamTable.vL_p(p1)



    p2=p3_hp#bar
    #s2=s1
    # h2s=steamTable.h_ps(p2,s1)
    # h2 = (h2s-h1)/eta_SiC+h1
    h2=v1*(p2-p1)*10**2/eta_SiC+h1#kJ/kg
    T2= steamTable.t_ph(p2,h2)
    s2 = steamTable.s_ph(p2,h2)
    x2= None # eau non saturée


    """
    2) Boiler
    """
    # s'occuper ici de la combustion

    T3=T_max #K
    p3= p3_hp #bar
    h3=steamTable.h_pt(p3,T3-273.15)
    s3=steamTable.s_pt(p3,T3-273.15)
    x3 = None # vapeur surchauffée

    Q1 = h3-h2 #kJ/kg_v

    """
    3) Turbine
    """


    p4=p4
    h4s = steamTable.h_ps(p4,s3)
    h4= h3-(h3-h4s)*eta_SiT
    x4 = x6
    s4 = steamTable.s_ph(p4,h4)
    T4 = steamTable.t_ph(p4,h4)#°C

    """
    4) Condenser
    """
    Q2 = h4-h1 # kJ/kg_v

    """
    5) Mechanical work:
    """
    Wm_t = h3-h4
    Wm_c = h2-h1
    Wm = Wm_t-Wm_c

    """
    6) Cycle efficiency and massflows
    """
    eta_cyclen = Wm/Q1
    eta_toten = 0
    mv = Pe/(Wm*eta_mec) #kg_v/s

    eta_gen = 0

    """
    7) Computation of enegy losses :
    """
    P_cond = Q2*mv#kW
    P_mec = Pe-Wm*mv
    P_gen = 0
    """
    8) Exergy efficiencies
    """
    eta_cyclex =0
    eta_totex =0
    eta_gex=0
    eta_combex = 0
    eta_chemex = 0
    eta_condex = 0
    eta_transex =0
    """
    Last) Define output arguments
    """
    outputs = ST_arg.ST_outputs();
    outputs.eta[0:]= [eta_cyclen,eta_toten,eta_cyclex,eta_totex,eta_gen,eta_gex,eta_combex,eta_chemex,eta_condex,eta_transex]
    outputs.daten[0:]=[P_gen, P_mec, P_cond]
    outputs.dat[0:]= [[T1-273.15,T2-273.15,T3-273.15,T4-273.15],[p1,p2,p3,p4],[h1,h2,h3,h4],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    return outputs;

# Example:
steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS); # m/kg/sec/°C/bar/W
print (steamTable.h_pt(50,300));
print (steamTable.s_pt(50,300));
print (steamTable.h_ps(0.05,steamTable.s_pt(50,300)));


ST_inputs = ST_arg.ST_inputs();
results = ST(ST_inputs);
print(results.dat[:])
