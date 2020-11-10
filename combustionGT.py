import numpy as np;
import GTcomb_arguments as GT_arg;
import useful
"""
Convention :  All temperatures in input are given in °C but are immediatly changed into K
inside of the function. output temperatures are also in °C
"""

def combustionGT(comb_input):
    """
     combustion GT handles the combustion chamber inside the gas turbine
     It is called in the main function GT.

     INPUTS
     comb_input see in GTcomb_arguments.py
     Important: combustionGT can find T_out fromb the excess air coefficient lambda
     or the opposite. This has to be specified in the inputs via the flag inversion
     If inversion == True , then lambda is found from T_out

     OUTPUTS
     R_f   : ideal gas constant (R*) for specific gas (R/Mm_f) [kJ/kg/K]
     m_O2f : mass fraction of oxygen in exhaust gases [-]
     m_N2f : mass fraction of nitrogen in exhaust gases [-]
     m_CO2f : mass fraction of carbon dioxyde in exhaust gases [-]
     m_H2Of : mass fraction of water steam in exhaust gases [-]
     T_out  : outlet gas temperature [K]
    """
    lambda_comb = comb_input.lambda_comb
    x_O2a = comb_input.x_O2a
    x_N2a = comb_input.x_N2a
    coeff = x_N2a/x_O2a
    T_in  = comb_input.T_in +273.15 #[K]
    T_in_comb = comb_input.T_in_comb + 273.15  #[K]
    h_in  = comb_input.h_in #kJ/kg_air
    LHV   =comb_input.LHV #kJ/kg_ch4]
    HHV = comb_input.HHV #kJ/kg_CH4
    T0 =  288.15 #[K]
    inversion = comb_input.inversion #boolean set to false
    T_out = comb_input.T_out +273.15 #[K]
    kcc = comb_input.k_cc
    r= comb_input.r
    Mm_CH4 = 0.016; Mm_O2 = 0.032; Mm_N2 = 0.028; Mm_H2O = 0.018; Mm_CO2 = 0.044
    CH4=useful.CH4

    molar_mass = np.array([0.028,0.044,0.018,0.032]) #kg/mol N2- CO2 - H2O - O2
    Mm_a = x_O2a * Mm_O2 + x_N2a * Mm_N2 # [kg/mol_air]
    ma1 =  Mm_a/Mm_CH4 * 2/x_O2a # kg_air/kg_CH4 = proportion of air  vs combustible
    mass_conc0 = np.array([x_N2a,0,0,x_O2a])*molar_mass/Mm_a

    # At the exit of the combustion:
    coeff_stochio = np.array([2*lambda_comb*coeff,1,2,2*(lambda_comb-1)]) # N2- CO2 - H2O - O2
    total_n = sum(coeff_stochio) # total moles
    molar_conc = coeff_stochio/total_n # concentration of elements
    Mm_af = sum(molar_conc*molar_mass) #kg/mol_air
    mass_conc = molar_conc*molar_mass/Mm_af #[-] kg_co2/kg_tot


    #  find lambda or T_out if the flag is up by iterations
    iter = 1
    error = 1
    dt = 0.1

    h_f0 =  useful.janaf_integrate_air(useful.cp_air,mass_conc,Mm_af,T0-15,T0,dt)
    hc= useful.janaf_integrate(useful.cpCH4,273.15,T_in_comb,0.001)/Mm_CH4
    ha = useful.janaf_integrate_air(useful.cp_air,mass_conc0,Mm_a,T0-15,T_in,0.0001) #attention useful.cp_air [J/kg_air] #basically h_in=ha

    if (inversion == False):
        T_out = 1273.15 #[K] first estimation
        while iter < 50 and error > 0.01 :
            cp_f = useful.cp_mean_air(useful.cp_air,mass_conc,Mm_af,T0,T_out,dt)
            T_out_final = (T0 + ((1000*LHV + hc + lambda_comb*ma1*ha)/((lambda_comb*ma1+1)*cp_f)) - h_f0/(cp_f))
            iter = iter + 1
            error = abs(T_out_final - T_out)
            T_out = T_out_final
            # print("Nombre d'itérations : ",iter)
            # print("T_out : ",T_out,"K")

    if (inversion == True):
        lambda_comb = 2 # first estimation
        while iter <50 and error > 0.01 :
            coeff_stochio = np.array([2*lambda_comb*coeff,1,2,2*(lambda_comb-1)]) # N2- CO2 - H2O - O2
            total_n = sum(coeff_stochio) # total moles
            molar_conc = coeff_stochio/total_n # concentration of elements mol_co2/mol_tot
            Mm_af = sum(molar_conc*molar_mass)
            mass_conc = molar_conc*molar_mass/Mm_af

            cp_f = useful.cp_mean_air(useful.cp_air,mass_conc,Mm_af,T0,T_out,dt)
            h_f0 = useful.janaf_integrate_air(useful.cp_air,mass_conc,Mm_af,T0-15,T0,dt)

            lambda_comb_final = (cp_f*(T_out-T0) + h_f0 - LHV*1000 - hc)/(ma1*(ha - h_f0 + cp_f*(T0-T_out)))
            iter = iter + 1
            error = abs(lambda_comb_final-lambda_comb)
            lambda_comb = lambda_comb_final

    # exergie and eta_combex ==> see formula page 28
    e_c = HHV+ 15 * (useful.cp_mean(useful.cpCH4,273.15,T0,dt)/0.016 + mass_conc[3]/0.032*useful.cp_mean(useful.cpO2,273.15,T0,dt) - mass_conc[1]/0.044*useful.cp_mean(useful.cpCO2,273.15,T0,dt) - mass_conc[2]/0.018*useful.cp_mean(useful.cpH2O,273.15,T0,dt))/1000 - T0*(CH4.S(273.15)/0.016+useful.cp_mean(useful.cpCH4,273.15,T0,dt)/0.016*np.log(T0/273.15))/1000 # kJ/kg_ch4

    #influence of preheating
    e_a = useful.cp_mean_air(useful.cp_air,mass_conc0,Mm_a,T0,T_in,dt)*(T_in-T0) - T0*(useful.janaf_integrate_air(useful.cp_air_T,mass_conc0,Mm_a,T0,T_in,dt)) #attention useful.cp_air [J/kg_air/K] => e_a = J/kg_air
    e_cr = useful.cp_mean(useful.cpCH4,T0,T_in,dt)*(T_in-T0)/Mm_CH4 - (useful.janaf_integrate(useful.cpCH4_T,T0,T_in,dt)*T0/Mm_CH4) #J/kg_CH4
    e_r = e_a*((lambda_comb*ma1)/(lambda_comb*ma1+1)) + e_cr*(1/(lambda_comb*ma1+1))

    Rf=8.31/Mm_af # [J/kg/K]
    delta_sf = useful.janaf_integrate_air(useful.cp_air_T,mass_conc,Mm_af,T0,T_out,dt)#- Rf*np.log(kcc*r)
    e_f = useful.janaf_integrate_air(useful.cp_air,mass_conc,Mm_af,T0,T_out,dt) - T0*delta_sf #ici j ai changé
    eta_combex = (e_f-e_r)*(lambda_comb*ma1+1)/(e_c*1000)

    # final outputs
    outputs = GT_arg.comb_output();
    outputs.Mm_af = Mm_af
    outputs.lambda_comb = lambda_comb
    outputs.ma1 = ma1
    outputs.T_out = T_out - 273.15 #°C
    outputs.R_f = Rf
    outputs.m_N2f,outputs.m_CO2f,outputs.m_H2Of,outputs.m_O2f = mass_conc  #[-]
    outputs.eta_combex =eta_combex
    outputs.e_c = e_c
    return outputs;
