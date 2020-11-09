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
    arg_in = ST_inputs;

    ## Check input arguments
    # ======================
    # if value == -1 <=> value not defined. Thus, it will be defined here:
    Pe = arg_in.Pe;
    if Pe ==-1.:
        Pe = 250e3;#250 kWe
    T_ext = arg_in.T_ext;
    if T_ext ==-1.:
        T_ext = 288.15;#15°C
    nsout = arg_in.nsout;
    if nsout ==-1.:
        nsout = 3;#15°C



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



    ## define output arguments
    # ======================
    outputs = ST_arg.ST_outputs();
    outputs.eta[1] = 0.35;
    outputs.Xmassflow = np.array(nsout);

    ## Generate graph to export:
    # ==========================
    # random data:
    t = np.arange(0.0, 2.0, 0.01)
    s1 = np.sin(2*np.pi*t)
    s2 = np.sin(4*np.pi*t)

    # My 1st figure:
    fig1 = plt.figure(1)
    plt.plot(t, s1)

    # My 2nd figure:
    fig2 = plt.figure(2)
    plt.plot(t, s2)

    # ...
    outputs.fig = [fig1, fig2]

    if (arg_in.DISPLAY == 1):
        plt.show()

    return outputs;

# Example:
steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS); # m/kg/sec/°C/bar/W
print (steamTable.h_pt(50,300));
print (steamTable.s_pt(50,300));
print (steamTable.h_ps(0.05,steamTable.s_pt(50,300)));


ST_inputs = ST_arg.ST_inputs();
ST_inputs.DISPLAY = 1;
results = ST(ST_inputs);
results.fig[1].show();
