import numpy as np;

import GTcomb_arguments as GT_comb_arg;
import useful

Mm_CH4 = 0.016; Mm_O2 = 0.032; Mm_N2 = 0.028; Mm_H2O = 0.018; Mm_CO2 = 0.044

def heatexchanger(exchanger_input,T_air_out):
    """
    When using this method, we know the mass flow of the two gas that are exchanging heat. For know we are only using this function
    with air and the flue gases.
    We also know the air temperature (T_air_in) and the flue gas temperature (T_f_in) coming in the exchanger.
    The user has to indicate the temperature of the air at the exit of the exchanger : T_air_out. Once heated,
    this air will go in the combustion chamber.

    The purpose of the function is to compute the temperature of the flue gases out.
    It will also compute the exchanger area S (given a Heat transfer coefficient U in the argument of exchanger input) using Hausbrand
    equation and the exergy efficiency of the exchanger : eta_transex.

    For the future, we can improve this function in three ways :
        1) for now, heatexchanger only takes air and flue gases as an input. We could make the function accept more than one type of compound
        2) a better estimation of U and S could be done but we have not the solution for this problem now.
        3) We could use NTU for defining the temperatures in and out.
    WARNING : the function has temperature limits : when those limits are exceeded an error message will be given the command
    consol.
    """
    U = exchanger_input.U
    Mflow_air_in = exchanger_input.Mflow_air_in
    Mflow_f_in = exchanger_input.Mflow_f_in
    T_air_in = exchanger_input.T_air_in +273.15 #[K]
    T_f_in = exchanger_input.T_f_in +273.15 #[K]
    lambda_comb = exchanger_input.comb_lambda

    x_O2a = 0.21
    x_N2a = 0.79
    coeff = x_N2a/x_O2a
    T0 =  288.15 #[K]

    # at the air input

    Mm_a = x_O2a * Mm_O2 + x_N2a * Mm_N2 # [kg/mol_air]
    mass_conc1=np.array([x_N2a*Mm_N2/Mm_a,0,0,x_O2a*Mm_O2/Mm_a])
    ma1 =  Mm_a/Mm_CH4 * 2/x_O2a # kg_air/kg_CH4 = proportion d air entrant vs combustible

    # at the flue gases input
    coeff_stochio = np.array([2*lambda_comb*coeff,1,2,2*(lambda_comb-1)]) # N2- CO2 - H2O - O2
    total_n = sum(coeff_stochio) # nombre de moles total
    molar_conc = coeff_stochio/total_n # concentration des elements mol_co2/mol_t
    molar_mass = np.array([0.028,0.044,0.018,0.032]) #kg/mol N2- CO2 - H2O - O2
    Mm_af = sum(molar_conc*molar_mass) #somme ponderé des masse molaire
    mass_conc2 = molar_conc*molar_mass/Mm_af #[-] kg_co2/kg_tot


    dt = 0.01
    iter = 1
    error = 20

    #Heat
    T_air_out = T_air_out +273.15 # K
    Q= Mflow_air_in * useful.janaf_integrate_air(useful.cp_air,mass_conc1,Mm_a,T_air_in,T_air_out,dt)
    print("Q : ",Q,'[J]')

    #first estimation of T_f_out before the iterations
    #this esyimation is made with cp constants (without integration)
    cp_f = useful.cp_air(T_f_in,mass_conc2,Mm_af) #J/kg_fumée
    T_f_out = (-Q/(Mflow_f_in*cp_f)) + T_f_in
    T_f_out_secours = T_f_out

    #We used a third degree polynomial to interpolate the cp values of each component and integrate those values on a temperature interval.
    #Once the polynomial computed, we can iterate the calcul and find a more accurate T_f_out.

    while iter<50 and error>0.01 :
        iter = iter+1
        heat_const_N2 = useful.cp_Iconstants('N2',T_f_out,T_f_in)
        heat_const_CO2 = useful.cp_Iconstants('CO2',T_f_out,T_f_in)
        heat_const_H2O = useful.cp_Iconstants('H2O',T_f_out,T_f_in)
        heat_const_O2 = useful.cp_Iconstants('O2',T_f_out,T_f_in)
        int_N2 = (1/2)*heat_const_N2[2][2]*(T_f_in**2-T_f_out**2) + (1/3)*heat_const_N2[2][1]*(T_f_in**3-T_f_out**3) + (1/4)*heat_const_N2[2][0]*(T_f_in**4-T_f_out**4)
        int_CO2 = (1/2)*heat_const_CO2[2][2]*(T_f_in**2-T_f_out**2) + (1/3)*heat_const_CO2[2][1]*(T_f_in**3-T_f_out**3) + (1/4)*heat_const_CO2[2][0]*(T_f_in**4-T_f_out**4)
        int_H2O = (1/2)*heat_const_H2O[2][2]*(T_f_in**2-T_f_out**2) + (1/3)*heat_const_H2O[2][1]*(T_f_in**3-T_f_out**3) + (1/4)*heat_const_H2O[2][0]*(T_f_in**4-T_f_out**4)
        int_O2 = (1/2)*heat_const_O2[2][2]*(T_f_in**2-T_f_out**2) + (1/3)*heat_const_O2[2][1]*(T_f_in**3-T_f_out**3) + (1/4)*heat_const_O2[2][0]*(T_f_in**4-T_f_out**4)

        hfout = np.array([int_N2,int_CO2,int_H2O,int_O2])
        hf = np.dot(hfout,mass_conc2)/Mm_af #J/kg_fumée

        A_variables = np.array([heat_const_N2[2][3],heat_const_CO2[2][3],heat_const_H2O[2][3],heat_const_O2[2][3]])
        A_var = np.dot(A_variables,mass_conc2)/Mm_af

        T_f_out_final = ((hf*Mflow_f_in - Q)/(A_var*Mflow_f_in))+T_f_in
        error = abs(T_f_out_final-T_f_out)
        T_f_out = T_f_out_final

        if iter==50 :
            print("Le heat exchanger ne converge peut-être pas")
            T_f_out = T_f_out_secours
    print("Nombres d'itération : ",iter)
    print("T_f_out : ",T_f_out-273.15,'C')

    cp_f = useful.cp_air(T_f_out,mass_conc2,Mm_af)
    print('second in',cp_f)

    #Hausbrand equation depending on two situations
    if Mflow_air_in*useful.cp_air(T_air_in,mass_conc1,Mm_a)>Mflow_f_in*cp_f :
        Deltag_T = T_f_in - T_air_out
        Deltap_T = T_f_out - T_air_in
        DeltaM_T = (Deltag_T - Deltap_T)/np.log(Deltag_T/Deltap_T)
        S = Q/(U*DeltaM_T)
    else :
        Deltag_T = T_f_out - T_air_in
        Deltap_T = T_f_in - T_air_out
        DeltaM_T = (Deltag_T - Deltap_T)/np.log(Deltag_T/Deltap_T)
        S = Q/(U*DeltaM_T)

    print("S : ",S)

    """
    MEAN TEMPERATURE CALCULATIONS
    """

    deltas_f = useful.janaf_integrate_air(useful.cp_air_T,mass_conc2,Mm_af,T_f_out,T_f_in,dt)
    deltah_fout = useful.janaf_integrate_air(useful.cp_air,mass_conc2,Mm_af,273.15,T_f_out,dt)
    deltah_fin = useful.janaf_integrate_air(useful.cp_air,mass_conc2,Mm_af,273.15,T_f_in,dt)
    Tfm = (-deltah_fout+deltah_fin)/deltas_f #[K] mean temperature des fumées

    ###################################

    Tam = (useful.janaf_integrate_air(useful.cp_air,mass_conc1,Mm_a,273.15,T_air_out,dt)-useful.janaf_integrate_air(useful.cp_air,mass_conc1,Mm_a,273.15,T_air_in,dt))/ useful.janaf_integrate_air(useful.cp_air_T,mass_conc1,Mm_a,T_air_in,T_air_out,dt)

    """
    Without dissipative terms, exergetic efficiency of the heat transfer is :
    """

    eta_transex = ((Tam - 288.15)*Tfm)/(Tam*(Tfm-288.15))
    print("eta_transex :",eta_transex)

    outputs = GT_comb_arg.exchanger_output();
    outputs.T_air_out = T_air_out-273.15 #[°C]
    outputs.T_f_out = T_f_out -273.15 #°C
    outputs.Q = Q/1000#[kJ]
    outputs.eta_transex = eta_transex
    outputs.Surf = S
    return outputs
