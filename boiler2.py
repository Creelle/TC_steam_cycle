import numpy as np
import ST_arguments as ST_arg
import STboiler_arguments as STboiler_arg
import useful
from thermochem import janaf
db = janaf.Janafdb();

#ici contre courant


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

    #Combustion
    Lambda = arg_in.Lambda
    xO2a=arg_in.x_O2a
    xN2a = arg_in.x_N2a
    T_in = arg_in.T_in +273.15 #K
    inversion = arg_in.inversion
    T_out = arg_in.T_out + 273.15#K
    ftype = arg_in.ftype
    x = arg_in.x
    y =  arg_in.y
    HHV = arg_in.HHV
    LHV = arg_in.LHV

    #Exchanger
    Q= arg_in.Q
    T_ext =arg_in.T_ext+273.15#K
    T_exhaust = arg_in.T_exhaust+273.15#K
    TpinchHR = arg_in.TpinchHR

    # molar_mass_f = 0.018
    coeff = xN2a/xO2a
    T0 = 288.15 #K

    # if x==0 and y!=1 :
    #     ftype = 'C'+ 'H' + str(y)
    # if x==0 and y==1 :
    #     ftype = 'C'+ 'H'
    # if y==0 and x!=1 :
    #     ftype = 'C' + 'O' + str(x)
    # if y==0 and x==1 :
    #     ftype = 'C' + 'O'
    # if x!=0 and y!=0 :
    #     ftype = 'C'+ 'H' + str(y) + 'O' + str(x)

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
    molar_mass = np.array([0.028,0.044,0.018,0.032]) #kg/mol N2- CO2 - H2O - O2
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

    T_in = T_ext+60 #K # préchauffage
    print(T_in)
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
            #Lambda = 2 # first estimation
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
        2) Heat recovery add a preheater for the air before entering the combustion
        chamber
        """
        # on connait T_ext, TpinchHR et T_exhaust
        # T_out ou lambda nous est donné par la combustion
        # Tout va dependre T_cold_out donc il faudra boucler

        #T_hot_in =T_in +TpinchHR#K premiere estimation



        """
        3) Determine massflows
        """
        massflow_f = Q*1000/useful.janaf_integrate_air(useful.cp_air,mass_conc,Mm_af,T_in+TpinchHR,T_out,dt) #kg_f/kg_vapor
        massflow_a = massflow_f/(1+1/(Lambda*ma1))
        massflow_c = massflow_f-massflow_a #kg/s
        print('massflow_c',massflow_c)

        Cp_f = useful.cp_mean_air(useful.cp_air,mass_conc,Mm_af,T_exhaust,T_in+TpinchHR,dt)
        Cp_a = useful.cp_mean_air(useful.cp_air,mass_conc0,Mm_a,T_ext,T_in,dt)
        T_in_new = (massflow_f*Cp_f*T_exhaust-massflow_f*Cp_f*TpinchHR-massflow_a*Cp_a*T_ext)/(massflow_f*Cp_f-massflow_a*Cp_a)
        #print(T_in_new,T_ext,T_exhaust,T_hot_in)
        error2 = abs(T_in-T_in_new)
        iter2+=1
        if iter2 ==100:
            print("Max iteration occured, there may have been no convergence")
        T_in = T_in_new


    T_hot_in = T_in+TpinchHR
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
    outputs.T_hot_in = T_hot_in-273.15

    outputs.boiler_massflow[0:]= [massflow_a,0,massflow_c,massflow_f]

    #pour l instant il manque e_c , eta_combex,
    return outputs

results = boiler(STboiler_arg.boiler_input(inversion=True))
print(results.boiler_massflow)
