import numpy as np


class boiler_input:
    """
     combustion(lambda,x_O2a,x_N2a) computes the combustion of methane with air given the excess air (lambda).
     It returns the molar fraction of the different

     INPUTS
     lambda_comb : excess air coefficient
     x_O2a : molar fraction of oxygen concentration in air
     x_N2a : molar fraction of nitrogen concentration in air
     T_in  : gas temperature at the inlet
     h_in  : enthalpy of the gas at the inlet
     LHV  : Fuel 'Low Heating Value'. CH4 here [kJ/kg_CH4]
     HHV : Fuel 'High Heating Value'. CH4 here [kJ/kg_CH4]
     inversion : boolean in order to know if it has to compute T_out (False) or lambda (True)
     ftype :  the fuel type
    """
    def __init__(self, Lambda = 2,#excess air
                     x_O2a = 0.21,# molar fraction
                     x_N2a = 0.79,# molar fraction
                     T_in = 15,#°C
                     inversion =False,
                     T_out = 1400, #°Cpour trouver lambda en fonction de T_out mettre true
                     HHV = 55695, # [kJ/kg_CH4],
                     LHV =50150,
                     T_ext = 15,#°C
                     # T_exhaust = 150,#°C
                     TpinchHR = 30,#°C
                     T_exhaust = 50,#°C
                     x = 0,
                     y =4,
                     Tdb = 30, # °C
                     absolute_humidity = 0.01,
                     Q=3361*31):# [kW]
        # combustion
        self.Lambda = Lambda;
        self.x_O2a = x_O2a;
        self.x_N2a = x_N2a ;
        self.T_in = T_in;
        self.inversion = inversion;
        self.T_out = T_out;
        self.HHV = HHV;
        self.LHV = LHV;
        self.x = x
        self.y = y
        self.Tdb = Tdb
        self.T_exhaust=T_exhaust
        self.absolute_humidity = absolute_humidity

        #exchanger
        self.Q = Q
        self.T_ext = T_ext
        # self.T_exhaust = T_exhaust
        self.TpinchHR = TpinchHR

class boiler_output:
    """
     combustion(lambda,x_O2a,x_N2a) computes the combustion of methane with air given the excess air (lambda).
     It returns the mass fraction of the different

     OUTPUTS
     m_O2f : mass fraction of oxygen in exhaust gases [-]
     m_N2f : mass fraction of nitrogen in exhaust gases [-]
     m_CO2f : mass fraction of carbon dioxyde in exhaust gases [-]
     m_H2Of : mass fraction of water steam in exhaust gases [-]
     T_out  : outlet gas temperature [K]

    """
    def __init__(self,m_O2f = 0.1,# molar fraction
                     m_N2f = 0.4,# molar fraction
                     m_CO2f = 0.2,#°C
                     m_H2Of = 0.3,# enthalpy
                     air_massflow = 1, #kg/s
                     LHV = 50000,#Low heating value [kJ/kg_f]
                     Lambda=2,
                     T_out = 1200,#°C
                     e_c = 1., #kJ/kg_ch4
                     eta_combex = 1.,
                     Cp_g = 1,# [kJ/kg/K]
                     T_hot_in = 150,#T_b
                     T_hot_out= 100,#T_c
                     T_cold_in = 15,#T_atmospherique
                     T_cold_out= 50,# T_a before combustion
                     T_dew = 2,# °C
                     T_exhaust = 150, # °C
                     P_chimney = 1, #kW
                     L_comb = 1,
                     L_HR =1,
                     L_exhaust = 1,
                     e_boiler_in = 1,
                     e_boiler_out = 1,
                     eta_gen = 1,
                     eta_chemex = 1,
                     eta_transex_HR=0.5):

        self.m_O2f = m_O2f;
        self.m_N2f = m_N2f ;
        self.m_CO2f = m_CO2f;
        self.m_H2Of = m_H2Of;


        self.Lambda = Lambda;
        self.T_out = T_out#°C
        self.e_c=e_c;
        self.Cp_g = Cp_g;#[kJ/kg/K]
        self.LHV = LHV;

        self.boiler_massflow = np.zeros(4);# air, vapor/water , combutible, exhaust
        self.T_hot_in = T_hot_in;
        self.T_hot_out = T_hot_out;
        self.T_cold_in = T_cold_in;
        self.T_cold_out = T_cold_out;
        self.T_dew = T_dew;
        self.T_exhaust = T_exhaust;

        self.P_chimney = P_chimney
        self.L_comb = L_comb
        self.L_HR = L_HR
        self.L_exhaust = L_exhaust

        self.eta_gen = eta_gen
        self.eta_transex_HR = eta_transex_HR;
        self.eta_combex = eta_combex;
        self.eta_chemex = eta_chemex;

        self.e_boiler_in = e_boiler_in
        self.e_boiler_out = e_boiler_out

# class exchanger_input:
#     """
#     this class computes the necessary inputs for a simple heat exchanger
#     INPUTS :
#     U : transmission coefficient
#     Mflow_air_in : the flow of air (cold) coming in [kg/s]
#     Mflow_f_in : the flow of hot gasses coming in
#     T_air_in: cold source temperature [°C]
#     T_f_in : hot source temperature [°C]
#     comb_lambda : excess air coefficient
#     courant :  contre-courant = -1 ; co-courant = 1
#     """
#     def __init__(self, U = 3,
#                      Mflow_air_in = 45,
#                      Mflow_f_in = 50,
#                      T_air_in = 400, #°C
#                      T_f_in = 700, # °C
#                      comb_lambda = 2,
#                      courant = -1):
#         self.U = U;
#         self.Mflow_air_in = Mflow_air_in
#         self.Mflow_f_in = Mflow_f_in
#         self.T_air_in = T_air_in
#         self.T_f_in = T_f_in
#         self.comb_lambda = comb_lambda
#         self.courant = courant;
#
# class exchanger_output:
#     """
#     T_air_out: cold source tempearature at the outlet
#     T_f_out: hot source temperature at the outlet
#     eta_transex : exergetic heat exchanger efficiency
#     Surf : surface of the heat exchanger [m**2]
#     Q: exchanged heat betweent hot and cold source [kJ]
#     """
#     def __init__(self,T_air_out = 450, # [°C]
#                      T_f_out = 550, # [°C]
#                      eta_transex = 0.5,
#                      Surf = 50, # surface d'échange [m**2]
#                      Q = 1000): #[kJ]
#         self.T_air_out = T_air_out
#         self.T_f_out = T_f_out
#         self.eta_transex = eta_transex
#         self.Surf = Surf
#         self.Q = Q
