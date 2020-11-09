## Project LMECA2150-Thermal cycle
# Template for the Steam turbine arguments
#
# Author: Paolo Thiran & Gauthier Limpens
# Version: 2020
#
# !!! This script CANNOT be modified by the students.

import numpy;

class ST_inputs:
    """
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
    
    class Combustion:
        def __init__(self):
            self.Tmax   =-1.;#   [°C] : maximum combustion temperature
            self.Lambda =-1.;#   [-] : air excess
            self.x      =-1.;#   [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
            self.y      =-1.;#   [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
    
    def __init__(self, comb = Combustion()):
        self.Pe = 250.0e3;# 250e3 [kW]
        self.nsout = -1.;# Number of feed-heating 
        self.reheat = -1.;# Number of reheating
        self.T_max = -1.;# Maximum steam temperature
        self.T_cond_out = -1.;# Condenser cold outlet temperature
        self.p3_hp  = -1.;# Maximum pressure
        self.drumFlag = -1.;# if =1 then drum if =0 then no drum. 
        self.eta_mec = -1.;# mechanic efficiency of shafts bearings
        self.combustion = comb;        
        self.T_exhaust = -1.;# Temperature of exhaust gas out of the chimney
        self.p4       = -1.;# High pressure after last reheating
        self.x6       = -1.;# Vapor ratio [gaseous/liquid] (in french : titre)
        self.T_0      = -1.;# Reference temperature
        self.T_ext    = -1.;# External temperature (atmospheric)
        self.TpinchSub = -1.;# Temperature pinch at the subcooler
        self.TpinchEx  = -1.;# Temperature pinch at a heat exchanger
        self.TpinchCond = -1.;# Temperature pinch at condenser 
        self.TpinchHR   = -1.;# Temperature pinch at the heat recovery system
        self.Tdrum      = -1.;# minimal drum temperature
        self.eta_SiC    = -1.;# Internal pump efficiency
        self.eta_SiT    = -1.;# Isotrenpic efficiency for Turbines.
        self.DISPLAY = -1;
    
class ST_outputs:
    """
     OUPUTS : 
     ETA is a vector with :
       -eta(1) : eta_cyclen, cycle energy efficiency
       -eta(2) : eta_toten, overall energy efficiency
       -eta(3) : eta_cyclex, cycle exegy efficiency
       -eta(4) : eta_totex, overall exergy efficiency
       -eta(5) : eta_gen, Steam generator energy efficiency
       -eta(6) : eta_gex, Steam generator exergy efficiency
       -eta(7) : eta_combex, Combustion exergy efficiency
       -eta(8) : eta_chemex, Chimney exergy efficiency (losses)
       -eta(9) : eta_condex, Condenser exergy efficiency
       -eta(10): eta_transex, Steam bleeding heat exchanger overall exergy efficiency
       FYI : eta(i) \in [0;1] [-]
     Xmassflow is a vector with each feedheating massflow [kg/s] (with respect to figure 
               2.33, page 91 "Thermal Power Plants" English version).
               Xmassflow(1) = mass flow at 6_1 etc...
     DATEN is a vector with : 
       -daten(1) : perte_gen [kW]
       -daten(2) : perte_mec [kW]
       -daten(3) : perte_cond [kW]
     DATEX is a vector with :
       -datex(1) : perte_mec    [kW]
       -datex(2) : perte_totex  [kW]
       -datex(3) : perte_rotex  [kW]
       -datex(4) : perte_combex [kW]
       -datex(5) : perte_condex [kW]
       -datex(6) : perte_chemex [kW]
       -datex(7) : perte_transex[kW]
     DAT is a matrix containing :
     dat = {T_1       , T_2       , ...       , T_6_I,     T_6_II, ... ;  [°C]
            p_1       , p_2       , ...       , p_6_I,     p_6_II, ... ;  [bar]
            h_1       , h_2       , ...       , h_6_I,     h_6_II, ... ;  [kJ/kg]
            s_1       , s_2       , ...       , s_6_I,     s_6_II, ... ;  [kJ/kg/K]
            e_1       , e_2       , ...       , e_6_I,     e_6_II, ... ;  [kJ/kg]
            x_1       , x_2       , ...       , x_6_I,     x_6_II, ... ;   };[-]
     MASSFLOW is a vector containing : 
       -massflow(1) = m_a, air massflow [kg/s]
       -massflow(2) = m_v, water massflow at 2 [kg/s]
       -massflow(3) = m_c, combustible massflow [kg/s] 
       -massflow(4) = m_f, exhaust gas massflow [kg/s]
     
     COMBUSTION is a structure with :
       -combustion.LHV    : the Lower Heat Value of the fuel [kJ/kg]
       -combustion.e_c    : the combustible exergie         [kJ/kg]
       -combustion.Lambda : the air excess                   [-]
       -combustion.Cp_g   : heat capacity of exhaust gas at T_max [kJ/kg/K]
       -combustion.fum  : is a vector of the exhaust gas composition :
           -fum(1) = m_O2f  : massflow of O2 in exhaust gas [kg/s]
           -fum(2) = m_N2f  : massflow of N2 in exhaust gas [kg/s]
           -fum(3) = m_CO2f : massflow of CO2 in exhaust gas [kg/s]
           -fum(4) = m_H2Of : massflow of H2O in exhaust gas [kg/s] 
    
     HEAT_RECOVERED is a structure containing:
       - heat_recovered.T_hot_in   : inlet  temperature of exhaust gases after exchangers with cycle [°C]
       - heat_recovered.T_hot_out  : outlet temperature of exhaust gases after exchangers with cycle [°C]
       - heat_recovered.T_cold_in  : inlet  temperature of air before heat recovery [°C]
       - heat_recovered.T_cold_out : outlet temperature of air after heat recovery[°C]
       - heat_recovered.T_dew : dew point temperature at the outlet of the chimney
     According Fig 2.2 in reference book (english version):
            T_hot_in = T_b
            T_hot_out = T_c
            T_cold_in = T_atmospheric 
            T_cold_out = T_a (assuming state 'a' before mixture with fuel)
    
    
     FIG is a vector of all the figure you plot. Before each figure, define a
     figure environment such as:  
      "FIG(1) = figure;
      plot(x,y1);
      [...]
       FIG(2) = figure;
      plot(x,y2);
      [...]"
      Your vector FIG will contain all the figure plot during the run of this
      code (whatever the size of FIG).
    """
    class Combustion:
        def __init__(self):
            self.LHV =-1.; # [kJ/kg]
            self.e_c =-1.; # [kJ/kg]
            self.Lambda =-1.; #   [-]
            self.Cp_g =-1.; #  [kJ/kg/K]
            self.fum = numpy.zeros(4);
            
    class Heat_recovery:
        def __init__(self):
            self.T_hot_in   = -1.; # [kJ/kg]
            self.T_hot_out  = -1.; # [kJ/kg]
            self.T_cold_in  = -1.; #   [-]
            self.T_cold_out = -1.; #  [kJ/kg/K]
            self.T_dew      = -1.;
            
    def __init__(self, comb = Combustion()):
        self.eta = numpy.zeros(10); # [-] \in [0,1]
        self.Xmassflow = list();# must be defined by the user # [kg/s]
        self.daten = numpy.zeros(3);# [kW]
        self.datex = numpy.zeros(7);# [kW]
        self.dat = list();#numpy.zeros((5, 4));
        self.massflow = numpy.zeros(4); # must be defined by the user
        self.combustion = comb;
        self.fig = list();
