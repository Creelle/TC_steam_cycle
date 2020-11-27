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
        nsout = 1;#15°C
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
        p3_hp =40 #bar
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
        p4 = 20 #bar
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
        TpinchSub = 10;#delta K
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
    print(nsout)
    if (nsout == 0):
        results = np.zeros((6,4+2*reheat)) # 6 initial states, states 1,2,3,6 if no reheat then 7 = 1 states 4 and 5 is for the reheating

    elif(nsout == 1):
        results = np.zeros((6,10+2*reheat)) #7,10,1,2,3,[reheat],6,61,81,91,101 # est ce qu on met 101 quelque part?
    else:
        results = np.zeros((6,10+2*(nsout-1)+2*reheat))
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

    p10=p3_hp/4#bar choix
    #s2=s1
    # h2s=steamTable.h_ps(p2,s1)
    # h2 = (h2s-h1)/eta_SiC+h1
    h10=v7*(p10-p7)*10**2/eta_SiC+h7#kJ/kg
    T10= steamTable.t_ph(p10,h10)+273.15#K
    s10 = steamTable.s_ph(p10,h10)
    x10= None # eau non saturée
    e10 = h10-T0*s10#kJ/kg
    results[:,1]=T10-1073.105,p10,h10,s10,x10,e10

    """
    3) Boiler
    """
    # s'occuper ici de la combustion

    T3=T_max #K
    p3= p3_hp #bar
    h3=steamTable.h_pt(p3,T3-273.15)
    s3=steamTable.s_pt(p3,T3-273.15)
    x3 = None # vapeur surchauffée
    e3 = h3-T0*s3#kJ/kg
    results[:,4]=T3-273.35,p3,h3,s3,x3,e3

    """
    4) Turbine
    """

    p6=p7
    h6s = steamTable.h_ps(p7,s3)
    h62= h3-(h3-h6s)*eta_SiT
    x6=x6
    h6 = x6*steamTable.hV_p(p6)+(1-x6)*steamTable.hL_p(p6)
    s6 = steamTable.s_ph(p6,h6)
    s62 = x6*steamTable.sV_p(p6)+(1-x6)*steamTable.sL_p(p6)
    T6 = steamTable.t_ph(p6,h6)+273.15#K
    e6 = h6-T0*s6 #kJ/kg#kJ/kg
    results[:,5+2*reheat]=T6-273.15,p6,h6,s6,x6,e6

    """
    4) FWH
    """
    #trouver l 'etat 2'
    T2_prime = steamTable.tsat_p(p3)+273.15
    h2_prime = steamTable.hL_p(p3)
    s2_prime =steamTable.sL_p(p3)

    #premiere estimation de T81
    T81 = 0.5*(T7+T2_prime)#K
    T101 = 0.5*(T10+T7)

    def function_FHW_T81(x):

        y=x[1]
        x=x[0]
        """
        This function computes the energy balance at the exchanger condenser and
        the exchanger subcooler for one feed heating
        INPUTS :  x (T81) estimated temperature at the exit of the condenser
        y (T101) estimated temperature between the two exchangers for the cold fluid
        """

        h81 = steamTable.hL_t(x-273.15) #kJ/kg_v
        p81 = steamTable.psat_t(x-273.15) #bar

        #find state 61
        p61=p81 #isobare
        h61s = steamTable.h_ps(p81,s3)
        h61 = -eta_SiT*(h3-h61s)+h3

        T1 = T81-TpinchEx#K
        h1=steamTable.h_pt(p10,T1-273.15)

        h101 = steamTable.h_pt(p10,y-273.15)

        T91 = T10+TpinchSub #K
        h91 = steamTable.h_pt(p81,T91-273.15)

        #enthalpie ratio
        X= (h1-h10)/(h61-h1+h10-h91)
        f1 = X*(h61-h81)-(1+X)*(h1-h101)
        f2 = X*(h81-h91)-(1+X)*(h101-h10)

        return np.array([f1,f2])

    T81,T101 = fsolve(function_FHW_T81,np.array([T81,T101]))

    """
    5) Alimentation pump
    """
    T1= T81-TpinchEx#K
    p1 = p10
    h1= steamTable.h_pt(p1,T1-273.15)
    s1= steamTable.s_pt(p1,T1-273.15)
    x1 = None
    v1 = steamTable.v_pt(p1,T1-273.15)
    e1 = h1-T0*s1#kJ/kg
    results[:,2]=T1-273.15,p1,h1,s1,x1,e1


    p2=p3_hp#bar
    #s2=s1
    # h2s=steamTable.h_ps(p2,s1)
    # h2 = (h2s-h1)/eta_SiC+h1
    h2=v1*(p2-p1)*10**2/eta_SiC+h1#kJ/kg
    T2= steamTable.t_ph(p2,h2)+273.15#K
    s2 = steamTable.s_ph(p2,h2)
    x2= None # eau non saturée
    e2 = h2-T0*s2#kJ/kg
    results[:,3]=T2-273.25,p2,h2,s2,x2,e2

    """
    6) FWH states
    """

    # state 81
    h81 = steamTable.hL_t(T81-273.15)
    p81 = steamTable.psat_t(T81-273.15)#bar
    s81 = steamTable.sL_t(T81-273.15) #kJ/kg
    x81=0
    e81 = h81-T0*s81
    i=0 # i =range(nsout)
    results[:,7]=T81-273.15,p81,h81,s81,x81,e81

    #state 6I
    p61=p81 #isobare
    h61s = steamTable.h_ps(p61,s3)
    h61 = -eta_SiT*(h3-h61s)+h3
    T61 = steamTable.t_ph(p61,h61)+273.15
    s61 = steamTable.s_ph(p61,h61)
    x61 = None
    e61 = h61-T0*s61
    results[:,6]=T61-273.35,p61,h61,s61,x61,e61

    #state 91
    T91 = T10+TpinchSub #K
    h91 = steamTable.h_pt(p81,T91-273.15)
    s91 = steamTable.s_pt(p81,T91-273.15)
    p91 = p81
    x91 = None
    e91 = h91-T0*s91
    results[:,8]=T91-273.35,p91,h91,s91,x91,e91

    #state 91_postvanne
    h91_postvanne = h91
    x91_postvanne = steamTable.x_ph(p7,h91)

    #state 101
    h101 = steamTable.h_pt(p10,T101-273.15)
    s101 = steamTable.s_pt(p10,T101-273.15)
    p101 = p10
    x101 = None
    e101 = h101-T0*s101
    results[:,9]=T101-273.35,p101,h101,s101,x101,e101

    #bleed proportion
    X= (h1-h10)/(h61-h1+h10-h91)

    Q1 = (h3-h2) #kJ/kg_v*(1+X)/(1+X)

    """
    5) Condenser
    """
    Q2 = 1/(1+X)*(h6-h7) + X/(1+X)*(h91-h7)
     # kJ/kg_v
    massflow_condenser_coeff = Q2/(steamTable.CpL_t(T_cond_out-273.15)*(T_cond_out-T_ext))

    h_cond_water_in = steamTable.h_pt(1.01325,T_ext-273.15)
    h_cond_water_out = steamTable.h_pt(1.01325,T_cond_out-273.15)
    s_cond_water_in = steamTable.s_pt(1.01325,T_ext-273.15)
    s_cond_water_out = steamTable.s_pt(1.01325,T_cond_out-273.15)
    e_cond_water_in = h_cond_water_in-T0*s_cond_water_in
    e_cond_water_out = h_cond_water_out-T0*s_cond_water_out

    """
    6) Mechanical work:
    """
    Wm_t = 1/(1+X)*(h3-h6) + X/(1+X)*(h3-h61) #kJ/kg_v
    Wm_tmax = 1/(1+X)*e3-e6 + X/(1+X)*(e3-e61) #kJ/kg_v
    Wm_c = h2-h1+h10-h7 #kJ/kg_v #*(1+X)/(1+X)
    Wm_cmax = e2-e1+e10-e7
    Wm = Wm_t-Wm_c
    Wm_max = Wm_tmax-Wm_cmax

    """
    7) massflows
    """
    mv = Pe/(Wm*eta_mec) #kg_v/s
    Q_boiler = mv*Q1 #kW

    """
    8) Boiler combustion and heat recovery
    """
    boiler_inputs = STboiler_arg.boiler_input(inversion=inversion, Lambda = comb.Lambda, T_out = comb.Tmax-273.15,
                                            T_exhaust =T_exhaust-273.15,TpinchHR = TpinchHR,T_ext = T_ext-273.15,
                                            Q = Q_boiler)
    boiler_outputs = boiler(boiler_inputs)
    ma,dummy,mc,mf = boiler_outputs.boiler_massflow[0:]

    """
    9) Energetic efficiencies
    """
    eta_cyclen = Wm/Q1
    eta_gen = boiler_outputs.eta_gen
    eta_toten = eta_mec*eta_gen*eta_cyclen

    """
    10) Computation of energy losses :
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
    11) Exergy efficiencies
    """
    ec = boiler_outputs.e_c
    e_boiler_in = boiler_outputs.e_boiler_in
    e_boiler_out = boiler_outputs.e_boiler_out#kJ/kg_f aka e_ech

    eta_cyclex =Wm/(e3-e2)
    eta_rotex = Wm/Wm_max

    eta_combex = boiler_outputs.eta_combex
    eta_chemex = boiler_outputs.eta_chemex
    eta_transex =mv*(e3-e2)/(mf*(e_boiler_in-e_boiler_out)) # echange au boiler


    eta_gex=mv*(e3-e2)/(mc*ec)
    # eta_gex2= eta_combex*eta_chemex*eta_transex
    # print(eta_gex,eta_gex2)
    eta_totex = eta_cyclex*eta_gex*eta_mec
    #eta_totex2 = Pe/(mc*ec)
    #print(eta_totex,eta_totex2)
    #print("eta_combex",eta_combex,"eta_chemex",eta_chemex,'eta_transex',eta_transex,"eta_gex",eta_gex,"eta_rotex",eta_rotex,'eta_cyclex',eta_cyclex,"eta_totex",eta_totex,'eta_gen',eta_gen)
    eta_condex = 0

    """
    12) Computation of exergy losses
    """
    L_comb = boiler_outputs.L_comb #kW
    L_HR = boiler_outputs.L_HR
    L_exhaust = boiler_outputs.L_exhaust

    L_boiler = mv*(e2-e3)+mf*(e_boiler_in-e_boiler_out)

    m_prin = 1/(1+X)*mv
    m_bleed = X/(1+X)*mv

    L_turbine = T0*(s6-s3)*m_prin + T0*(s6-s61)*m_bleed
    L_pump = T0*(s2-s1)*mv+T0*(s10-s7)*mv#kW

    L_cond=m_prin*(e6-e7)+m_bleed*(e91-e7)+mv*massflow_condenser_coeff*(e_cond_water_in-e_cond_water_out)#vanne included

    L_exchanger_soutex = m_bleed*(e61-e91)+mv*(e10-e1)

    #Il manque ce qui sort de la tour de refroidissement

    ec = boiler_outputs.e_c
    print('exegie chequ up',ec*mc,Pe+Pf_mec+L_turbine+L_pump+L_boiler+L_cond+L_exhaust+L_HR+L_comb+L_exchanger_soutex)

    """
    Last) Define output arguments
    """
    outputs = ST_arg.ST_outputs();
    outputs.eta[0:]= [eta_cyclen,eta_toten,eta_cyclex,eta_totex,eta_gen,eta_gex,eta_combex,eta_chemex,eta_condex,eta_transex]
    outputs.daten[0:]=[P_chimney, Pf_mec, P_cond]
    outputs.datex[0:]=[Pf_mec,0,0,L_comb,L_cond,L_exhaust,L_exchanger_soutex]
    outputs.dat= results
    outputs.massflow = boiler_outputs.boiler_massflow #[ma,0,mc,mf]
    outputs.massflow[1] = mv
    outputs.Xmassflow = np.array([m_bleed])

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
    13) Pie charts and cycle graphs
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
              'Condenser losses {v} [MW]'.format(v=round(L_cond/1000)),'Feed heating losses {v} [MW]'.format(v=round(L_exchanger_soutex/1000)),
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

# Example:
# steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS); # m/kg/sec/°C/bar/W
# print (steamTable.h_pt(50,300));
# print (steamTable.s_pt(50,300));
# print (steamTable.h_ps(0.05,steamTable.s_pt(50,300)));


ST_inputs = ST_arg.ST_inputs();
ST_inputs.Pe = 35.0e3 #[kW]
ST_inputs.DISPLAY = 1
answers = ST(ST_inputs);
print(answers.massflow,answers.Xmassflow)
