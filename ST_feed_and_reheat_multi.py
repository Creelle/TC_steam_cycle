## Project LMECA2150-Thermal cycle
#Template for the Steam turbine
#
# Author: Paolo Thiran & Gauthier Limpens
# Version: 2020
#
# Students can modify this script.
# However, the arguments for the ST defined in ST_arguments.py CANNOT be modified

import numpy as np;
from scipy.optimize import fsolve
import ST_arguments as ST_arg;
import STboiler_arguments as STboiler_arg;
from boiler import boiler
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
     figures a faire TS, HS, piechart energy, piechart exergy
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
        nsout = 0;#15°C
    reheat = arg_in.reheat
    if reheat == -1:
        reheat = 0;# Number of reheating
    T_max = arg_in.T_max;
    if T_max == -1.:
        T_max = 520 #°C
    T_cond_out = arg_in.T_cond_out;
    if T_cond_out == -1.:
        T_cond_out = 25#°C
    p3_hp = arg_in.p3_hp;
    if p3_hp == -1.:
        p3_hp =50#bar
    eta_mec = arg_in.eta_mec;
    if eta_mec == -1.:
        eta_mec = 0.95 #[-]

    comb = arg_in.combustion

    inversion = True  # on doit calculer Lambda a partir de Tmax
    comb.Tmax = 273.15+1400;

    if comb.x == -1:
        comb.x = 0;
    if comb.y == -1:
        comb.y = 4;

    T_exhaust = arg_in.T_exhaust;
    if T_exhaust == -1.:
        T_exhaust = 150 #°C
    p4 = arg_in.p4;
    if p4 == -1.:
        p4 = 30 #bar
    x6 = arg_in.x6;
    if x6 == -1.:
        x6 = 0.91 #[-]
    T0 = arg_in.T_0;
    if T0 == -1.:
        T0 = 15#°C
    T_ext = arg_in.T_ext;
    if T_ext ==-1.:
        T_ext = 15;#15°C
    Tdrum = arg_in.Tdrum;
    if Tdrum ==-1.:
        Tdrum = 15;#15°C

    #Bon il reste a faire TpinchSub,TpinchEx, TpinchCond,TpinchHR
    # ce sont des differences de temperatures donc on ne change pas

    TpinchSub = arg_in.TpinchSub;
    if TpinchSub ==-1.:
        TpinchSub = 15;#delta K
    TpinchEx = arg_in.TpinchEx;
    if TpinchEx ==-1.:
        TpinchEx = 15;#delta K
    TpinchCond = arg_in.TpinchCond;
    if TpinchCond ==-1.:
        TpinchCond = 5;#delta K
    TpinchHR = arg_in.TpinchHR;
    if TpinchHR ==-1.:
        TpinchHR = 100;#delta K

    eta_SiC = arg_in.eta_SiC;
    if eta_SiC == -1.:
        eta_SiC = 0.88
    eta_SiT = arg_in.eta_SiT;
    if eta_SiT == -1.:
        eta_SiT = 0.88

    # Put all temperatures in Kelvin
    T_max,T_cond_out, T_exhaust, T0, T_ext, T_drum= T_max+273.15,T_cond_out+273.15, T_exhaust+273.15, T0+273.15, T_ext+273.15, Tdrum +273.15 #K
    if nsout ==0:
        results = np.zeros((6,6+3*(nsout)+2*reheat))#7,10,1,2,3,[reheat],6
    else :
        results = np.zeros((6,6+3*(nsout)+2+2*reheat))#results 7,10,1,2,3,[reheat :41,51,...,4last,pre ], 6, 61,62,...,6n, 81,82,...,8n,101,102,...,10n(1),91

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

    Q1 = 0 #chaleur fournie au boiler
    Q2 = 0 #chaleur échangée au condenseur
    Q_boiler_exergie = 0 #bilan exergie au boiler
    Wm_t = 0 #kJ/kg_v
    Wm_tmax = 0 #kJ/kg_v
    L_turbine_mv =0 #il faudra faire * mv plus tard


    """
    1) Initial
    """

    T7= T_cond_out+TpinchCond
    p7 = steamTable.psat_t(T7-273.75)
    h7= steamTable.hL_p(p7)
    s7= steamTable.sL_p(p7)
    x7 = 0
    v7 = steamTable.vL_p(p7)
    e7 = h7-T0*s7#kJ/kg
    results[:,0]=T7-273.75,p7,h7,s7,x7,e7

    """
    2) Extraction pump
    """

    p10=p7*500
    h10=v7*(p10-p7)*10**2/eta_SiC+h7#kJ/kg
    T10= steamTable.t_ph(p10,h10)+273.15#K
    s10 = steamTable.s_ph(p10,h10)
    x10= None # eau non saturée
    e10 = h10-T0*s10#kJ/kg
    results[:,1]=T10-273.105,p10,h10,s10,x10,e10

    """
    3) Boiler
    """

    T3=T_max #K
    p3= p3_hp #bar
    h3=steamTable.h_pt(p3,T3-273.15)
    s3=steamTable.s_pt(p3,T3-273.15)
    x3 = None # vapeur surchauffée
    e3 = h3-T0*s3#kJ/kg
    results[:,4]=T3-273.35,p3,h3,s3,x3,e3

    """
    4) Reheating
    """

    """
    Here we begin the code vectorization because we want to add the number of reheating that we want
    """
    for i in range (0,reheat):

            h_pre = results[2][4+2*i]
            s_pre = results[3][4+2*i]
            e_pre = results[5][4+2*i]
            p4i = p3 - (p3-p4)*(i+1)/reheat #donné dans les arguments
            s4is = s_pre #détente isentropique
            h4is = steamTable.h_ps(p4i,s4is) #attention on est ni à saturation, ni sous la courbe! => à checker
            h4i = h_pre-(h_pre-h4is)*eta_SiT
            x4i = None #check que ça doit pas passer sous la courbe?
            T4i = steamTable.t_ph(p4i,h4i)+273.15#K
            s4i = steamTable.s_ph(p4i,h4i)
            e4i = h4i-T0*s4i #kJ/kg#kJ/kg
            results[:,4+2*i+1]=T4i-273.15,p4i,h4i,s4i,x4i,e4i

            Wm_t += h_pre-h4i
            Wm_tmax += e_pre-e4i
            L_turbine_mv += T0*(s4i-s_pre)

            p5i = p4i # on considère que la combustion se fait sans perte de charge
            T5i = T3
            h5i = steamTable.h_pt(p5i,T5i-273.15)
            s5i = steamTable.s_pt(p5i,T5i-273.15)
            x5i = None
            e5i = h5i-T0*s5i
            results[:,4+2*i+2]=T5i-273.15,p5i,h5i,s5i,x5i,e5i

            Q1 += (h5i-h4i)
            Q_boiler_exergie += e5i-e4i

    """
    5) Turbine : after the reheat
    """
    s_pre = results[3][4+2*reheat]
    h_pre = results[2][4+2*reheat]
    e_pre = results[5][4+2*reheat]

    p6=p7
    h6s = steamTable.h_ps(p7,s_pre)

    h62= h_pre-(h_pre-h6s)*eta_SiT
    x6=x6
    h6 = x6*steamTable.hV_p(p6)+(1-x6)*steamTable.hL_p(p6)
    s6 = steamTable.s_ph(p6,h6)
    s62 = x6*steamTable.sV_p(p6)+(1-x6)*steamTable.sL_p(p6)
    T6 = steamTable.t_ph(p6,h6)+273.15#K
    e6 = h6-T0*s6 #kJ/kg#kJ/kg
    results[:,5+2*reheat]=T6-273.15,p6,h6,s6,x6,e6

    """
    6) FWH
    """
    if nsout!= 0:

        #trouver l 'etat 2'
        T2_prime = steamTable.tsat_p(p3)+273.15
        h2_prime = steamTable.hL_p(p3)
        s2_prime =steamTable.sL_p(p3)

        #calcul de l'etat 6i et 8i
        deltah_bleedings = (h_pre-h6)/(nsout+1)

        h6i = np.zeros(nsout)
        for i in range(nsout):
            h6i[i] = h6+(i+1)*deltah_bleedings
        h6is = (h6i-h_pre)/eta_SiT+h_pre

        p8i = np.zeros(nsout)
        T8i = np.zeros(nsout)
        h8i = np.zeros(nsout)
        s8i = np.zeros(nsout)

        T6i = np.zeros(nsout)
        s6i = np.zeros(nsout)


        for i in range(nsout):
            p8i[i] = steamTable.p_hs(h6is[i],s_pre) #its a mistake but okay
            T8i[i] = steamTable.tsat_p(p8i[i])+273.15
            h8i[i] = steamTable.hL_t(T8i[i]-273.15)
            s8i[i] = steamTable.sL_t(T8i[i]-273.15)

            T6i[i] = steamTable.t_ph(p8i[i],h6i[i])+273.15
            s6i[i] = steamTable.s_ph(p8i[i],h6i[i])
        p6i = p8i
        x6i = None
        x8i = 0
        e8i = h8i-T0*s8i
        e6i = h6i - T0*s6i

        #calcul de l état 10i
        T10i =T8i-TpinchEx #T102,T103,...,T1
        p10i = p10*np.ones(nsout)
        h10i =np.zeros(nsout)
        s10i = np.zeros(nsout)
        for i in range(nsout):
            h10i[i] = steamTable.h_pt(p10,T10i[i]-273.15)
            s10i[i] = steamTable.s_pt(p10,T10i[i]-273.15)
        x10i = None
        e10i = h10i-T0*s10i

        #calcul de l 'état 1 '
        T1,p1,h1,s1,x1,e1 = T10i[-1],p10i[-1],h10i[-1],s10i[-1],x10i,e10i[-1]
        v1 =  steamTable.vL_p(p1)

        def function_FHW_T81(x):

            T101 = x[0]
            Xis = x[1:] #X1,X2,X3

            """
            This function computes the energy balance at the exchanger condenser and
            the exchanger subcooler for one feed heating
            INPUTS :  x (T81) estimated temperature at the exit of the condenser
            y (T101) estimated temperature between the two exchangers for the cold fluid
            """

            h101 = steamTable.h_pt(p10,T101-273.15)

            T10ia,h10ia = np.append(T101,T10i),np.append(h101,h10i) # (h101,h102,...h10n,h1) longueur n+1

            p81 = p8i[0]
            T9 = T10+TpinchSub #K
            h9 = steamTable.h_pt(p81,T9-273.15)


            if nsout == 1:
                X1 = Xis
                Fone = X1*(h8i[0]-h9)-(1+X1)*(h101-h10)
                Flast = X1*(h6i[-1]-h8i[-1])-(1+X1)*(h10i[-1]-h101)
                Functions = Fone
                Functions = np.append(Fone,Flast)
            else:

                sumXi = sum(Xis)
                Fone = sumXi*(h8i[0]-h9)-(1+sumXi)*(h101-h10)
                Fa = np.zeros(nsout-1)
                for i in range(len(Fa)):
                    a = i+1
                    sumXi_aplus_n = sum(Xis[a:])
                    Fa[i] = sumXi_aplus_n*(h8i[a]-h8i[a-1])+Xis[a-1]*(h6i[a-1]-h8i[a-1])-(1+sumXi)*(h10ia[a]-h10ia[a-1])# c est chaud pour les notations
                Flast = Xis[-1]*(h6i[-1]-h8i[-1])-(1+sumXi)*(h10i[-1]-h10i[-2])
                Functions = np.append(Fone,Fa)
                Functions = np.append(Functions,Flast)
            return Functions

        #premiere estimation de T81
        def initial(nsout):

            T101 = T10+15
            initial = 0.1*np.ones(nsout)
            initial = np.append(T101,initial)
            return initial

        # print(function_FHW_T81(initial(nsout)))

        Solutions= fsolve(function_FHW_T81,initial(nsout))
        T101 = Solutions[0]
        X_bleedings = Solutions[1:]

        if (X_bleedings[0]<0):
            print('opposite flow in the subcooler T TpinchSub different from T pinch Ex')
            print('modify the pressure in the extraction pump p10')
        if(T10i[-1]>steamTable.tsat_p(p10)+273.15):
            print('vapor formed before the activation pump, increase pressure at alimentation pump p10')

        # print('Solutions',T101,X_bleedings,sum(X_bleedings))
        # print(T10i)
        # print(steamTable.tsat_p(p10)+273.15)

        """
        8) FWH states
        """

        #state 91
        T91 = T10+TpinchSub #K
        h91 = steamTable.h_pt(p8i[0],T91-273.15)
        s91 = steamTable.s_pt(p8i[0],T91-273.15)
        p91 = p8i[0]
        x91 = None
        e91 = h91-T0*s91

        #state 91_postvanne
        h91_postvanne = h91
        x91_postvanne = steamTable.x_ph(p7,h91)

        #state 101
        h101 = steamTable.h_pt(p10,T101-273.15)
        s101 = steamTable.s_pt(p10,T101-273.15)
        p101 = p10
        x101 = None
        e101 = h101-T0*s101


        # Condesor :
        Q2 +=  sum(X_bleedings)/(1+sum(X_bleedings))*(h91-h7)

        #turbine work :
        sumXbleedings = sum(X_bleedings)
        Wm_t += 1/(1+sumXbleedings)*(h_pre-h6)
        Wm_tmax += 1/(1+sumXbleedings)*(e_pre-e6)
        L_turbine_mv += 1/(1+sumXbleedings)*T0*(s6-s_pre)
        for i in range(nsout):
            Wm_t += X_bleedings[i]/(1+sumXbleedings)*(h_pre-h6i[i])
            Wm_tmax += X_bleedings[i]/(1+sumXbleedings)*(e_pre-e6i[i])
            L_turbine_mv += X_bleedings[i]/(1+sumXbleedings)*T0*(s6i[i]-s_pre)


        #formater le resultat rappel : results 7,10,1,2,3,[reheat :41,51,...,4last,pre ], 6, 61,62,...,6n, 81,82,...,8n,101,102,...,10n(1),91
        results[:,6:6+nsout]=T6i-273.15,p6i,h6i,s6i,-1*np.ones(nsout),e6i

        results[:,6+nsout:6+2*nsout]=T8i-273.15,p8i,h8i,s8i,x8i*np.ones(nsout),e8i

        results[:,6+2*nsout]=T101-273.15,p101,h101,s101,x101,e101

        results[:,6+2*nsout+1:6+3*nsout+1] =T10i-273.15,p10i,h10i,s10i,-1*np.ones(nsout),e10i

        results[:,6+3*nsout+1]=T91-273.15,p91,h91,s91,x91,e91


    else:
        Wm_t += (h_pre-h6)
        Wm_tmax += (e_pre-e6)
        L_turbine_mv += T0*(s6-s_pre)
        T1,p1,h1,s1,x1,v1,e1 = T10,p10,h10,s10,x10,v7,e10
        h91 = h7
        e91 = e7
    try:
        sumXbleedings = sum(X_bleedings)
    except:
        print('nsout = 0')
        sumXbleedings = 0

    """
    7) Alimentation pump
    """

    results[:,2] = T1,p1,h1,s1,x1,e1

    p2=p3_hp#bar
    h2=v1*(p2-p1)*10**2/eta_SiC+h1#kJ/kg
    T2= steamTable.t_ph(p2,h2)+273.15#K
    s2 = steamTable.s_ph(p2,h2)
    x2= None # eau non saturée
    e2 = h2-T0*s2#kJ/kg
    results[:,3]=T2-273.25,p2,h2,s2,x2,e2

    Q1 += h3-h2
    Q_boiler_exergie += e3-e2

    """
    9) Condenser
    """

    Q2 += 1/(1+sumXbleedings)*(h6-h7)


     # kJ/kg_v
    massflow_condenser_coeff = Q2/(steamTable.CpL_t(T_cond_out-273.15)*(T_cond_out-T_ext))

    h_cond_water_in = steamTable.h_pt(1.01325,T_ext-273.15)
    h_cond_water_out = steamTable.h_pt(1.01325,T_cond_out-273.15)
    s_cond_water_in = steamTable.s_pt(1.01325,T_ext-273.15)
    s_cond_water_out = steamTable.s_pt(1.01325,T_cond_out-273.15)
    e_cond_water_in = h_cond_water_in-T0*s_cond_water_in
    e_cond_water_out = h_cond_water_out-T0*s_cond_water_out

    """
    10) Mechanical work:
    """

    Wm_c = h2-h1+h10-h7 #kJ/kg_v #*(1+X)/(1+X)
    Wm_cmax = e2-e1+e10-e7
    Wm = Wm_t-Wm_c
    Wm_max = Wm_tmax-Wm_cmax

    """
    11) massflows
    """
    mv = Pe/(Wm*eta_mec) #kg_v/s
    Q_boiler = mv*Q1 #kW

    try:
        m_prin = 1/(1+sumXbleedings)*mv
        m_bleedings = X_bleedings/(1+sumXbleedings)*mv
    except:
        m_prin = mv
        m_bleedings = np.zeros(1)

    """
    12) Boiler combustion and heat recovery
    """
    boiler_inputs = STboiler_arg.boiler_input(inversion=inversion, Lambda = comb.Lambda, T_out = comb.Tmax-273.15,
                                            T_exhaust =T_exhaust-273.15,TpinchHR = TpinchHR,T_ext = T_ext-273.15,
                                            Q = Q_boiler)
    boiler_outputs = boiler(boiler_inputs)
    ma,dummy,mc,mf = boiler_outputs.boiler_massflow[0:]

    """
    13) Energetic efficiencies
    """
    eta_cyclen = Wm/Q1
    eta_gen = boiler_outputs.eta_gen
    eta_toten = eta_mec*eta_gen*eta_cyclen

    """
    14) Computation of energy losses :
    """
    P_cond = Q2*mv#kW
    Pf_mec = Wm*mv-Pe
    P_boiler = Q_boiler
    P_chimney = boiler_outputs.P_chimney
    #P_prim = P_boiler+P_chimney #-P_air_in mais =0 = +/- LHV*mc
    P_prim = boiler_outputs.LHV*mc
    print("energie chequ up",P_prim,P_chimney+Pf_mec+P_cond+Pe)
    # print("P_chimney",P_chimney,"Q_boiler",Q_boiler,"P_cond",P_cond,"Pf_mec",Pf_mec,"Pe",Pe)
    """
    15) Exergy efficiencies
    """
    ec = boiler_outputs.e_c
    e_boiler_in = boiler_outputs.e_boiler_in
    e_boiler_out = boiler_outputs.e_boiler_out#kJ/kg_f aka e_ech

    eta_cyclex =Wm/Q_boiler_exergie
    eta_rotex = Wm/Wm_max

    eta_combex = boiler_outputs.eta_combex
    eta_chemex = boiler_outputs.eta_chemex

    eta_transex = -mv*(e10-e1)/(np.dot(m_bleedings,e6i)-sum(m_bleedings)*e91)#efficiency at the bleed exchangers

    eta_gex=mv*Q_boiler_exergie/(mc*ec)

    eta_totex = eta_cyclex*eta_gex*eta_mec

    eta_condex = 0

    """
    16) Computation of exergy losses
    """


    L_comb = boiler_outputs.L_comb #kW
    L_HR = boiler_outputs.L_HR
    L_exhaust = boiler_outputs.L_exhaust

    L_boiler = mv*(-Q_boiler_exergie)+mf*(e_boiler_in-e_boiler_out) #à checekr

    L_turbine = L_turbine_mv*mv
    L_pump = T0*(s2-s1)*mv+T0*(s10-s7)*mv#kW

    if nsout == 0 :
        L_cond=mv*(e6-e7)+mv*massflow_condenser_coeff*(e_cond_water_in-e_cond_water_out)#vanne included
        L_exchanger_soutex = 0
    else :
        L_cond = m_prin*(e6-e7)+sum(m_bleedings)*(e91-e7)+mv*massflow_condenser_coeff*(e_cond_water_in-e_cond_water_out)#vanne included
        L_exchanger_soutex = m_prin*(e10-e1)+np.dot(m_bleedings,e6i)-sum(m_bleedings)*e91
    print('L_cond',L_cond)
    print('L_exchanger_soutex',L_exchanger_soutex)

    #Il manque ce qui sort de la tour de refroidissement

    ec = boiler_outputs.e_c
    print('exegie chequ up',ec*mc,Pe+Pf_mec+L_turbine+L_pump+L_boiler+L_cond+L_exhaust+L_HR+L_comb+L_exchanger_soutex)

    L_rotex = L_pump+L_turbine
    L_totex = mc*ec-Pe


    """
    Last) Define output arguments
    """
    outputs = ST_arg.ST_outputs();
    outputs.eta[0:]= [eta_cyclen,eta_toten,eta_cyclex,eta_totex,eta_gen,eta_gex,eta_combex,eta_chemex,eta_condex,eta_transex]
    outputs.daten[0:]=[P_chimney, Pf_mec, P_cond]
    outputs.datex[0:]=[Pf_mec,L_totex,L_rotex,L_comb,L_cond,L_exhaust,L_exchanger_soutex]
    outputs.dat= results
    outputs.massflow = boiler_outputs.boiler_massflow #[ma,0,mc,mf]
    outputs.massflow[1] = mv
    outputs.Xmassflow = m_bleedings

    #combustion
    outputs.combustion.LHV = boiler_outputs.LHV #[kJ/kg_f]
    outputs.combustion.e_c = ec
    outputs.combustion.Lambda = boiler_outputs.Lambda
    outputs.combustion.Cp_g = boiler_outputs.Cp_g # [kJ/kg/K]
    outputs.combustion.fum =np.array([boiler_outputs.m_O2f,boiler_outputs.m_N2f,boiler_outputs.m_CO2f,boiler_outputs.m_H2Of])*mf

    #heat recovery
    outputs.HR.T_hot_in = boiler_outputs.T_hot_in
    outputs.HR.T_hot_out = boiler_outputs.T_hot_out
    outputs.HR.T_cold_in = boiler_outputs.T_cold_in
    outputs.HR.T_cold_out = boiler_outputs.T_cold_out
    outputs.HR.T_dew = boiler_outputs.T_dew #still gave to define

    """
    18) Pie charts and cycle graphs
    """
    # pie chart of the energie flux in the cycle
    fig,ax =  plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
    data = [Pe,Pf_mec,P_cond,P_chimney]
    labels = ['Useful power {v} [MW]'.format(v=round(Pe/1000)),'Mechanical losses {v} [MW]'.format(v=round(Pf_mec/1000)),'Condensor losses {v} [MW]'.format(v=round(P_cond/1000)),'Chimney losses {v} [MW]'.format(v=round(P_chimney/1000))]

    ax.pie(data,labels = labels,autopct='%1.2f%%',startangle = 90)
    ax.set_title("Primary energetic flux "+ str(round(P_prim/1000)) + "[MW]")

    # pie chart of the exergie flux in the cycle
    fig2,ax =  plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
    data = [Pe,Pf_mec,L_turbine+L_pump,L_cond,L_exchanger_soutex,L_boiler,L_HR,L_comb,L_exhaust]
    labels = ['Useful power {v} [MW]'.format(v=round(Pe/1000)),'Mechanical losses {v} [MW]'.format(v=round(Pf_mec/1000)),' \n \n Turbine and \n pump losses {v} [MW]'.format(v=round((L_turbine+L_pump)/1000)),
              '\n \n Condenser losses {v} [MW]'.format(v=round(L_cond/1000)),'\n \n Bleed heating losses {v} [MW]'.format(v=round(L_exchanger_soutex/1000)),
              'Boiler losses {v} [MW]'.format(v=round(L_boiler/1000)),'Heat recovery losses {v} [MW]'.format(v=round(L_HR/1000)),
              'Combustion losses {v} [MW]'.format(v=round(L_comb/1000)),'Chimney losses {v} [MW]'.format(v=round(L_exhaust/1000))]
    plt.savefig('figures/energie_pie.png')

    ax.pie(data,labels = labels,autopct="%1.2f%%",startangle = 90)
    ax.set_title("Primary exergetic flux "+ str(round(ec*mc/1000)) + "[MW]")
    plt.savefig('figures/exergie_pie.png')

    fig = [fig,fig2]
    outputs.fig = fig
    if (ST_inputs.DISPLAY == 1):
        plt.show()

    return outputs;


ST_inputs = ST_arg.ST_inputs();
ST_inputs.Pe = 35.0e3 #[kW]
ST_inputs.DISPLAY = 1
ST_inputs.nsout = 4
ST_inputs.reheat = 3
answers = ST(ST_inputs);
