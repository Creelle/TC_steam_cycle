import numpy as np
import ST_arguments as ST_arg
import STboiler_arguments as STboiler_arg
import ST_useful as useful
from thermochem import janaf
db = janaf.Janafdb();

#ici contre courant
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

"""
Le but ici sera de calculer ce qui se passe dans le boiler,
Dans un premier temps, on considere qu 'il n y a pas de prechauffage de l'air
donc pas de heat recovery
A la page 43, cela veut dire que le point a vaut le point b
L'entrée est donnée par Q qui est la chaleur nécessaire a echangée ,
Les temperatures de l 'air exterieur , T2 et T3 sont donnée '
Pour la temperature dans la chambre de combustion , on nous donne soit Lambda ou T_max
"""
def boiler(STboiler_input):
    arg_in = STboiler_input

    #Chimneu
    Tdb = arg_in.Tdb
    absolute_humidity = arg_in.absolute_humidity
    T_dew =psychrometrics(Tdb,absolute_humidity)+273.15#K

    T_exhaust = arg_in.T_exhaust+273.15

    #Combustion
    Lambda = arg_in.Lambda
    xO2a=arg_in.x_O2a
    xN2a = arg_in.x_N2a
    T_in = arg_in.T_in +273.15 #K
    inversion = arg_in.inversion
    T_out = arg_in.T_out + 273.15#K

    x = arg_in.x
    y =  arg_in.y
    HHV = arg_in.HHV
    LHV = arg_in.LHV

    #Exchanger
    Q= arg_in.Q #kW
    T_ext =arg_in.T_ext+273.15#K
    TpinchHR = arg_in.TpinchHR

    molar_mass_f = 0.018
    coeff = xN2a/xO2a
    T0 = 288.15 #K

    Mm_f = 0.012 + y*0.001 + x*0.016 ; Mm_O2 = 0.032; Mm_N2 = 0.028; Mm_H2O = 0.018; Mm_CO2 = 0.044

    mc = 0.012/Mm_f ; mh = y*0.001/Mm_f ; mo = (x*0.016)/Mm_f

    LHV = (38.2*mc + 84.9*(mh - mo/8))*1000 #[kJ/kg_fuel] --> formule modified Dulong's formula.
    HHV = LHV
    if y!=0:
        HHV = LHV + (y/2)*40.752/Mm_f #[kJ/kg_fuel] --> we're adding the latent heat of water.

    """
    1) Calcul de la combustion sans préchauffage :
    """

    #at the entering of the combustion
    molar_mass = np.array([Mm_N2,Mm_CO2,Mm_H2O,Mm_O2]) #kg/mol N2- CO2 - H2O - O2
    Mm_a = xO2a * Mm_O2 + xN2a * Mm_N2 # [kg/mol_air]
    ma1 =  Mm_a/Mm_f * (1 + (y-2*x)/4)/xO2a # kg_air/kg_CH4 = proportion of air  vs combustible
    mass_conc0 = np.array([xN2a,0,0,xO2a])*molar_mass/Mm_a

    # At the exit of the combustion:
    coeff_stochio = np.array([(1 + (y-2*x)/4)*Lambda*coeff,1,y/2,(1 + (y-2*x)/4)*(Lambda-1)]) # N2- CO2 - H2O - O2
    total_n = sum(coeff_stochio) # total moles
    molar_conc = coeff_stochio/total_n # concentration of elements
    Mm_af = sum(molar_conc*molar_mass) #kg/mol_air
    mass_conc = molar_conc*molar_mass/Mm_af #[-] kg_co2/kg_tot

    #  find lambda (inversion ==True) or T_out(inversion == False) by iterations

    T_in = T_ext+60 #K # préchauffage # premiere estimation

    #T_in = T_ext #sans préchauffage
    error2 = 1
    iter2 = 1

    while error2>0.1 and iter2<100:
        dt = 0.1
        iter = 1
        error = 1

        h_f0 =  useful.janaf_integrate_air(useful.cp_air,mass_conc,Mm_af,T0-15,T0,dt)
        #on neglige hc
        ha = useful.janaf_integrate_air(useful.cp_air,mass_conc0,Mm_a,T0-15,T_in,0.01) #attention useful.cp_air [J/kg_air] #basically h_in=ha

        if (inversion == False):
            T_out = 1273.15 #[K] first estimation
            while iter < 50 and error > 0.01 :
                cp_f = useful.cp_mean_air(useful.cp_air,mass_conc,Mm_af,T0,T_out,dt)
                T_out_final = (T0 + ((1000*LHV + Lambda*ma1*ha)/((Lambda*ma1+1)*cp_f)) - h_f0/(cp_f))
                iter = iter + 1
                error = abs(T_out_final - T_out)
                T_out = T_out_final
                # print("Nombre d'itérations : ",iter)
                # print("T_out : ",T_out,"K")

        if (inversion == True):
            Lambda = 2 # first estimation
            while iter <50 and error > 0.01 :
                coeff_stochio = np.array([(1 + (y-2*x)/4)*Lambda*coeff,1,y/2,(1 + (y-2*x)/4)*(Lambda-1)]) # N2- CO2 - H2O - O2
                total_n = sum(coeff_stochio) # total moles
                molar_conc = coeff_stochio/total_n # concentration of elements
                Mm_af = sum(molar_conc*molar_mass) #kg/mol_air
                mass_conc = molar_conc*molar_mass/Mm_af #[-] kg_co2/kg_tot

                cp_f = useful.cp_mean_air(useful.cp_air,mass_conc,Mm_af,T0,T_out,dt)
                h_f0 = useful.janaf_integrate_air(useful.cp_air,mass_conc,Mm_af,T0-15,T0,dt)

                Lambda_final = (cp_f*(T_out-T0) + h_f0 - LHV*1000)/(ma1*(ha - h_f0 + cp_f*(T0-T_out)))
                iter = iter + 1
                error = abs(Lambda_final-Lambda)
                Lambda = Lambda_final


        """
        2) Determine massflows
        """
        massflow_f = Q*1000/useful.janaf_integrate_air(useful.cp_air,mass_conc,Mm_af,T_in+TpinchHR,T_out,dt) #kg_f/kg_vapor
        massflow_a = massflow_f/(1+1/(Lambda*ma1))
        massflow_c = massflow_f-massflow_a #kg/s
        #print('massflow_c',massflow_c)

        """
        3) From the massflow, calculate a new T_in
        """

        Cp_f = useful.cp_mean_air(useful.cp_air,mass_conc,Mm_af,T_exhaust,T_in+TpinchHR,dt)
        Cp_a = useful.cp_mean_air(useful.cp_air,mass_conc0,Mm_a,T_ext,T_in,dt)
        T_in_new = (massflow_f*Cp_f*T_exhaust-massflow_f*Cp_f*TpinchHR-massflow_a*Cp_a*T_ext)/(massflow_f*Cp_f-massflow_a*Cp_a)
        error2 = abs(T_in-T_in_new)
        iter2+=1
        if iter2 ==100:
            print("Max iteration occured, there may have been no convergence")
        T_in = T_in_new


    T_hot_in = T_in+TpinchHR

    """
    4)  calcul des puissances
    """
    hair_ext = useful.janaf_integrate_air(useful.cp_air,mass_conc0,Mm_a,T0,T_ext,dt)/1000 #kJ/kg
    #print("puissance de l \'air entrant",hair_in*massflow_a)
    hair_out =  useful.janaf_integrate_air(useful.cp_air,mass_conc,Mm_af,T0,T_exhaust,dt)/1000 #kJ/kg
    # print("puissance de l \'air sortant",hair_out*massflow_f)
    # print("comparison",Q+hair_out*massflow_f-hair_in*massflow_a,massflow_c*LHV)
    """
    5) Exergie and losses in the boiler
    """
    CH4=useful.CH4
    ec = HHV- T0*(CH4.S(273.15)/0.016/1000 ) #kJ/kg_c

    s_air_ext = useful.janaf_integrate_air(useful.cp_air_T,mass_conc0,Mm_a,T0,T_ext,dt)
    hair_ext = useful.janaf_integrate_air(useful.cp_air,mass_conc0,Mm_a,T0,T_ext,dt)/1000
    e_air_ext = hair_ext-T0*s_air_ext/1000#kJ/kg_air

    s_air_in = useful.janaf_integrate_air(useful.cp_air_T,mass_conc0,Mm_a,T0,T_in,dt)
    hair_in = useful.janaf_integrate_air(useful.cp_air,mass_conc0,Mm_a,T0,T_in,dt)/1000
    e_air_in = hair_in-T0*s_air_in/1000#kJ/kg_air

    #boiler in
    s_f = useful.janaf_integrate_air(useful.cp_air_T,mass_conc,Mm_af,T0,T_out,dt)
    hf = useful.janaf_integrate_air(useful.cp_air,mass_conc,Mm_af,T0,T_out,dt)/1000
    e_f = hf-T0*s_f/1000#kJ/kg_f

    #boiler out
    s_f_out = useful.janaf_integrate_air(useful.cp_air_T,mass_conc,Mm_af,T0,T_hot_in,dt)
    hf_out = useful.janaf_integrate_air(useful.cp_air,mass_conc,Mm_af,T0,T_hot_in,dt)/1000
    e_f_out = hf_out-T0*s_f_out/1000#kJ/kg_f

    #exhaust
    s_f_exhaust = useful.janaf_integrate_air(useful.cp_air_T,mass_conc,Mm_af,T0,T_exhaust,dt)
    hf_exhaust = useful.janaf_integrate_air(useful.cp_air,mass_conc,Mm_af,T0,T_exhaust,dt)/1000
    e_f_exhaust = hf_exhaust-T0*s_f_exhaust/1000#kJ/kg_f

    L_comb = massflow_c*ec+massflow_a*e_air_in-massflow_f*e_f #kW
    L_HR = massflow_a*(e_air_ext-e_air_in)+massflow_f*(e_f_out-e_f_exhaust)
    L_exhaust = massflow_f*e_f_exhaust

    """
    6) Calculation of exergetic efficiencies
    """
    eta_combex = (massflow_f*e_f-massflow_a*e_air_in)/(massflow_c*ec)
    eta_chemex = (e_f-e_f_out)*massflow_f/(massflow_f*e_f-massflow_a*e_air_in)
    eta_transex_HR = -(e_air_in-e_air_ext)/(e_f_exhaust-e_f_out)

    #les pertes dans le boiler seront comptés dans ST.py
    """
    Last) Define the outputs :
    """
    outputs = STboiler_arg.boiler_output()
    outputs.m_N2f,outputs.m_CO2f,outputs.m_H2Of,outputs.m_O2f = mass_conc  #[-]
    outputs.Lambda = Lambda
    outputs.T_out = T_out -273.15 #°C
    outputs.Cp_g = useful.cp_air(T_out,mass_conc,Mm_af)/1000 # [kJ/kg/K]
    outputs.LHV = LHV

    outputs.T_hot_out = T_exhaust-273.15 #°C
    outputs.T_cold_in = T_ext-273.15 #°C
    outputs.T_cold_out = T_in-273.15#°C
    outputs.T_hot_in = T_hot_in-273.15 #°C
    outputs.T_dew = T_dew - 273.15 #°C
    outputs.T_exhaust = T_exhaust-273.15 #°C

    outputs.e_c =ec #kJ/kg_c
    outputs.P_chimney = hair_out*massflow_f#kW
    outputs.eta_gen = Q/(LHV*massflow_c+massflow_a*hair_ext)
    outputs.eta_combex = eta_combex
    outputs.eta_transex_HR =eta_transex_HR
    outputs.eta_chemex =eta_chemex

    outputs.L_comb=L_comb
    outputs.L_HR =L_HR
    outputs.L_exhaust = L_exhaust

    outputs.e_boiler_in = e_f #kJ/kg_f
    outputs.e_boiler_out = e_f_out #kJ/kg_f

    outputs.boiler_massflow[0:]= [massflow_a,0,massflow_c,massflow_f]
    print('boiler temperature',T_out,T_hot_in)


    return outputs

# results = boiler(STboiler_arg.boiler_input(inversion=True))
# print(results.Lambda)
# results2 = boiler(STboiler_arg.boiler_input(inversion=False,Lambda = results.Lambda))
# print(results2.T_out)
