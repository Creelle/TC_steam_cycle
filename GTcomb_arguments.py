
class comb_input:
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
     k_cc: pressure loss in the combustion chamber
     r :  compression ratio in the compressor
     inversion : boolean in order to know if it has to compute T_out (False) or lambda (True)
    """
    def __init__(self, lambda_comb = 2,#excess air
                     x_O2a = 0.21,# molar fraction
                     x_N2a = 0.79,# molar fraction
                     T_in = 600,#°C
                     T_in_comb = 15,# °C
                     h_in = 650,# enthalpy[kJ/kg_air]
                     HHV = 55695, # [kJ/kg_CH4]
                     inversion =False,
                     T_out = 1000, #°Cpour trouver lambda en fonction de T_out mettre true
                     k_cc = 1.0,
                     r = 1,
                     LHV =50150):# [kJ/kg_CH4]

        self.lambda_comb = lambda_comb;
        self.x_O2a = x_O2a;
        self.x_N2a = x_N2a ;
        self.T_in = T_in;
        self.T_in_comb = T_in_comb;
        self.h_in = h_in;
        self.LHV = LHV;
        self.HHV = HHV;
        self.inversion = inversion
        self.k_cc = k_cc
        self.r = r
        self.T_out = T_out

class comb_output:
    """
     combustion(lambda,x_O2a,x_N2a) computes the combustion of methane with air given the excess air (lambda).
     It returns the mass fraction of the different

     OUTPUTS
     R_f   : ideal gas constant (R*) for specific gas (R/Mm_f) [kJ/kg/K]
     m_O2f : mass fraction of oxygen in exhaust gases [-]
     m_N2f : mass fraction of nitrogen in exhaust gases [-]
     m_CO2f : mass fraction of carbon dioxyde in exhaust gases [-]
     m_H2Of : mass fraction of water steam in exhaust gases [-]
     T_out  : outlet gas temperature [K]

    """
    def __init__(self, R_f = -1.,#50 MW
                     m_O2f = 0.1,# molar fraction
                     m_N2f = 0.4,# molar fraction
                     m_CO2f = 0.2,#°C
                     m_H2Of = 0.3,# enthalpy
                     Mm_af = -1.,
                     lambda_comb=1,
                     ma1 =1,
                     e_c = 1., #kJ/kg_ch4
                     eta_combex = 1.,
                     T_out = -1.):# °C
        self.R_f = R_f;
        self.m_O2f = m_O2f;
        self.m_N2f = m_N2f ;
        self.m_CO2f = m_CO2f;
        self.m_H2Of = m_H2Of;
        self.T_out = T_out;
        self.Mm_af = Mm_af;
        self.lambda_comb = lambda_comb
        self.ma1 =ma1
        self.e_c=e_c
        self.eta_combex = eta_combex


class exchanger_input:
    """
    this class computes the necessary inputs for a simple heat exchanger
    INPUTS :
    U : transmission coefficient
    Mflow_air_in : the flow of air (cold) coming in [kg/s]
    Mflow_f_in : the flow of hot gasses coming in
    T_air_in: cold source temperature [°C]
    T_f_in : hot source temperature [°C]
    comb_lambda : excess air coefficient
    courant :  contre-courant = -1 ; co-courant = 1
    """
    def __init__(self, U = 3,
                     Mflow_air_in = 45,
                     Mflow_f_in = 50,
                     T_air_in = 400, #°C
                     T_f_in = 700, # °C
                     comb_lambda = 2,
                     courant = -1):
        self.U = U;
        self.Mflow_air_in = Mflow_air_in
        self.Mflow_f_in = Mflow_f_in
        self.T_air_in = T_air_in
        self.T_f_in = T_f_in
        self.comb_lambda = comb_lambda
        self.courant = courant;

class exchanger_output:
    """
    T_air_out: cold source tempearature at the outlet
    T_f_out: hot source temperature at the outlet
    eta_transex : exergetic heat exchanger efficiency
    Surf : surface of the heat exchanger [m**2]
    Q: exchanged heat betweent hot and cold source [kJ]
    """
    def __init__(self,T_air_out = 450, # [°C]
                     T_f_out = 550, # [°C]
                     eta_transex = 0.5,
                     Surf = 50, # surface d'échange [m**2]
                     Q = 1000): #[kJ]
        self.T_air_out = T_air_out
        self.T_f_out = T_f_out
        self.eta_transex = eta_transex
        self.Surf = Surf
        self.Q = Q
