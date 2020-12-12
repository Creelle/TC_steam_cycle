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
from pyXSteam.XSteam import XSteam # see documentation here: https://pypi.org/project/pyXSteam
plt.rcParams.update({'font.size': 20})


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
    fig3,ax3 = plt.subplots()

    ## Check input arguments
    # ======================
    # if value == -1 <=> value not defined. Thus, it will be defined here:
    Pe = arg_in.Pe;
    if Pe ==-1.:
        Pe = 250e3;#250 kWe
    nsout = arg_in.nsout;
    if nsout ==-1.:
        nsout = 4;#15°C
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
    if p3_hp == -1.:
        p3_hp =100#bar
    p4 = arg_in.p4;
    if p4 == -1.:
        p4 = 50 #bar
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
        Tdrum = 120;#°C

    #Bon il reste a faire TpinchSub,TpinchEx, TpinchCond,TpinchHR
    # ce sont des differences de temperatures donc on ne change pas

    TpinchSub = arg_in.TpinchSub;
    if TpinchSub ==-1.:
        TpinchSub = 10;#delta K
    TpinchEx = arg_in.TpinchEx;
    if TpinchEx ==-1.:
        TpinchEx = 10;#delta K
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

    if (nsout<4):
        print("not possible , nsout should be >= 4")
        return
    # Put all temperatures in Kelvin
    T_max,T_cond_out, T_exhaust, T0, T_ext, Tdrum= T_max+273.15,T_cond_out+273.15, T_exhaust+273.15, T0+273.15, T_ext+273.15, Tdrum +273.15 #K
    results = np.zeros((6,8+3*(nsout)+2+2*reheat+1))#results 7,10,71,11,1,2,3,[reheat :41,51,...,4last,pre ], 6, 61,62,,...,6nsout,6_drum (= 6_IP),...,61_IP,62_IP,...,81,82,...,8n,101,102,...,10n(1),91

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
    p7 = steamTable.psat_t(T7-273.15)
    h7= steamTable.hL_p(p7)
    s7= steamTable.sL_p(p7)
    x7 = 0
    v7 = steamTable.vL_p(p7)
    e7 = h7-T0*s7#kJ/kg
    results[:,0]=T7-273.75,p7,h7,s7,x7,e7

    """
    3) Boiler
    """

    T3=T_max #K
    p3= p3_hp #bar
    h3=steamTable.h_pt(p3,T3-273.15)
    s3=steamTable.s_pt(p3,T3-273.15)
    x3 = None # vapeur surchauffée
    e3 = h3-T0*s3#kJ/kg
    results[:,6]=T3-273.15,p3,h3,s3,x3,e3

    """
    4) Reheating
    """

    """
    Here we begin the code vectorization because we want to add the number of reheating that we want
    """
    for i in range (0,reheat):
            p_pre = results[1][6+2*i]
            h_pre = results[2][6+2*i]
            s_pre = results[3][6+2*i]
            e_pre = results[5][6+2*i]
            p4i = p3 - (p3-p4)*(i+1)/reheat #donné dans les arguments
            s4is = s_pre #détente isentropique
            h4is = steamTable.h_ps(p4i,s4is) #attention on est ni à saturation, ni sous la courbe! => à checker
            h4i = h_pre-(h_pre-h4is)*eta_SiT
            x4i = None #check que ça doit pas passer sous la courbe?
            T4i = steamTable.t_ph(p4i,h4i)+273.15#K
            s4i = steamTable.s_ph(p4i,h4i)
            e4i = h4i-T0*s4i #kJ/kg#kJ/kg
            results[:,6+2*i+1]=T4i-273.15,p4i,h4i,s4i,x4i,e4i

            Wm_t += h_pre-h4i
            Wm_tmax += e_pre-e4i
            L_turbine_mv += T0*(s4i-s_pre)

            p5i = p4i # on considère que la combustion se fait sans perte de charge
            T5i = T3
            h5i = steamTable.h_pt(p5i,T5i-273.15)
            s5i = steamTable.s_pt(p5i,T5i-273.15)
            x5i = None
            e5i = h5i-T0*s5i
            results[:,6+2*i+2]=T5i-273.15,p5i,h5i,s5i,x5i,e5i

            Q1 += (h5i-h4i)
            Q_boiler_exergie += e5i-e4i

            #graph
            T4i5i = np.linspace(T4i-273.15,T5i-273.15,100)
            S4i5i = np.zeros(len(T4i5i))
            for i in range(0,len(T4i5i)):
                S4i5i[i] = steamTable.s_pt(p4i,T4i5i[i])
            ax3.plot(S4i5i,T4i5i,'g')

            ppre4i = np.linspace(p_pre,p4i,100)
            Sprep = np.linspace(s_pre,s_pre,100)
            Tpre4i = np.zeros(len(ppre4i))
            Spre4i = np.zeros(len(ppre4i))
            hpre4is = np.zeros(len(ppre4i))
            hpre4i = np.zeros(len(ppre4i))
            for i in range(0,len(ppre4i)):
                hpre4is[i] = steamTable.h_ps(ppre4i[i],s_pre)
                hpre4i[i] = h_pre-(h_pre-hpre4is[i])*eta_SiT
                Tpre4i[i] = steamTable.t_ph(ppre4i[i],hpre4i[i])
                Spre4i[i] = steamTable.s_ph(ppre4i[i],hpre4i[i])
            ax3.plot(Spre4i,Tpre4i,'g')#,Sprep,Tpre4i,'b--')

    """
    5) Turbine : after the reheat : Intermediate pressure
    """
    T_pre = results[0][6+2*reheat]
    p_pre = results[1][6+2*reheat]
    s_pre = results[3][6+2*reheat]
    h_pre = results[2][6+2*reheat]
    e_pre = results[5][6+2*reheat]

    T6_IP= Tdrum
    p6_IP = 2 #first estimation

    def Function_p6_IP(x):
        h6s_IP = steamTable.h_ps(x,s_pre)
        h6_IP= h_pre-(h_pre-h6s_IP)*eta_SiT

        T6_IP = steamTable.t_ph(x,h6_IP)+273.15
        return Tdrum-T6_IP
    p6_IP = fsolve(Function_p6_IP,p6_IP)[0]
    print('p6_IP',p6_IP)
    h6s_IP = steamTable.h_ps(p6_IP,s_pre)
    h6_IP= h_pre-(h_pre-h6s_IP)*eta_SiT
    s6_IP = steamTable.s_ph(p6_IP,h6_IP)
    x6_IP = steamTable.x_ph(p6_IP,h6_IP)
    e6_IP = h6_IP-T0*s6_IP

    #Wm_t +=
    #Mm_tmax +=
    #L_turbine_mv +=

    """
    5.B) Low pressure turbine
    """
    p6=p7
    h6s = steamTable.h_ps(p7,s6_IP)

    h62= h_pre-(h6_IP-h6s)*eta_SiT
    x6=x6
    h6 = x6*steamTable.hV_p(p6)+(1-x6)*steamTable.hL_p(p6)
    s6 = steamTable.s_ph(p6,h6)
    s62 = x6*steamTable.sV_p(p6)+(1-x6)*steamTable.sL_p(p6)
    T6 = steamTable.t_ph(p6,h6)+273.15#K
    e6 = h6-T0*s6 #kJ/kg#kJ/kg
    results[:,7+2*reheat]=T6-273.15,p6,h6,s6,x6,e6

    """
    2) Extraction pump
    """

    p10=p6_IP
    h10=v7*(p10-p7)*10**2/eta_SiC+h7#kJ/kg
    T10= steamTable.t_ph(p10,h10)+273.15#K
    s10 = steamTable.s_ph(p10,h10)
    x10= None # eau non saturée
    e10 = h10-T0*s10#kJ/kg
    results[:,1]=T10-273.15,p10,h10,s10,x10,e10


    """
    3) Drum pump
    """

    T71 = steamTable.tsat_p(p6_IP)+273.15
    h71 = steamTable.hL_p(p6_IP)
    p71 = p6_IP
    s71 = steamTable.sL_p(p6_IP)
    v71 = steamTable.vL_p(p6_IP)
    x71 = 0
    e71 = h71-T0*s71
    results[:,2] = T71-273.15,p71,h71,s71,x71,e71


    """
    6) FWH - and drum function
    """
    if nsout!= 0:
        if nsout%2 == 0:
            nsout_IP = int(nsout/2)
            nsout = int(nsout/2)
        else :
            nsout_IP = int((nsout-1)/2)
            nsout = int(nsout-nsout_IP)
        print('nsout,nsout_IP',nsout,nsout_IP)
        #nsout_IP=0
        #calcul de l'etat 6i et 8i
        deltah_bleedings = (h6_IP-h6)/(nsout+1)
        deltah_bleedings_IP = (h_pre-h6_IP)/(nsout_IP+1)


        h6i = np.zeros(nsout)
        h6i_IP = np.zeros(nsout_IP)
        for i in range(nsout):

            h6i[i] = h6+(i+1)*deltah_bleedings
        for i in range(nsout_IP):
            h6i_IP[i] = h6_IP+(i+1)*deltah_bleedings_IP
        # h6i_IP = h6i[nsout+2:]
        # h6_IPbis = h6i[nsout+1]
        # h6i =h6i[0:nsout]
        #print(h6_IPbis,h6_IP,"here")
        #print(h6i_IP,h6_IP,h6i)

        h6is = (h6i-h6_IP)/eta_SiT+h6_IP
        h6is_IP = (h6i_IP-h_pre)/eta_SiT+h_pre

        p8i = np.zeros(nsout)
        T8i = np.zeros(nsout)
        h8i = np.zeros(nsout)
        s8i = np.zeros(nsout)

        p8i_IP = np.zeros(nsout_IP)
        T8i_IP = np.zeros(nsout_IP)
        h8i_IP = np.zeros(nsout_IP)
        s8i_IP = np.zeros(nsout_IP)

        T6i = np.zeros(nsout)
        s6i = np.zeros(nsout)
        x6i = np.zeros(nsout)

        T6i_IP = np.zeros(nsout_IP)
        s6i_IP = np.zeros(nsout_IP)
        x6i_IP = np.zeros(nsout_IP)

        for i in range(nsout):
            p8i[i] = steamTable.p_hs(h6is[i],s6_IP) #its a mistake but okay
            T8i[i] = steamTable.tsat_p(p8i[i])+273.15
            h8i[i] = steamTable.hL_t(T8i[i]-273.15)
            s8i[i] = steamTable.sL_t(T8i[i]-273.15)

            T6i[i] = steamTable.t_ph(p8i[i],h6i[i])+273.15
            s6i[i] = steamTable.s_ph(p8i[i],h6i[i])
            x6i[i] = steamTable.x_ph(p8i[i],h6i[i])
        for i in range(nsout_IP):
            p8i_IP[i] = steamTable.p_hs(h6is_IP[i],s_pre) #its a mistake but okay
            T8i_IP[i] = steamTable.tsat_p(p8i_IP[i])+273.15
            h8i_IP[i] = steamTable.hL_t(T8i_IP[i]-273.15)
            s8i_IP[i] = steamTable.sL_t(T8i_IP[i]-273.15)

            T6i_IP[i] = steamTable.t_ph(p8i_IP[i],h6i_IP[i])+273.15
            s6i_IP[i] = steamTable.s_ph(p8i_IP[i],h6i_IP[i])
            x6i_IP[i] = steamTable.x_ph(p8i_IP[i],h6i_IP[i])

        p6i = p8i
        x8i = 0
        e8i = h8i-T0*s8i
        e6i = h6i - T0*s6i

        p6i_IP = p8i_IP
        x8i_IP = 0
        e8i_IP = h8i_IP-T0*s8i_IP
        e6i_IP = h6i_IP - T0*s6i_IP

        #graph
        for i in range(nsout):
            if x6i[i] >= 1:
                #je calcule les états intermédiaire entre T6i et T6i à saturation
                T6isat = steamTable.tsat_p(p8i[i])+1
                S6isat = steamTable.sV_t(T6isat)
                T6i_6isat = np.linspace(T6i[i]-273.15,T6isat)
                S6i_6isat = np.zeros(len(T6i_6isat))
                for j in range(0,len(T6i_6isat)):
                    S6i_6isat[j] = steamTable.s_pt(p8i[i],T6i_6isat[j])
                ax3.plot(S6i_6isat,T6i_6isat,'k')

                T6isat_8i = np.linspace(T6isat,T8i[i]-273.15,100)
                S6isat_8i = np.linspace(S6isat,s8i[i],100)
                ax3.plot(S6isat_8i,T6isat_8i,'k')
            if x6i[i] < 1 :
                T6isat_8i = np.linspace(T6i[i]-273.15,T8i[i]-273.15,100)
                S6isat_8i = np.linspace(s6i[i],s8i[i],100)
                ax3.plot(S6isat_8i,T6isat_8i,'k')
        #vannes
        for i in range(1,nsout):
            P88sat = np.linspace(p8i[i],p8i[i-1],100)
            T88sat = np.zeros(100)
            S88sat = np.zeros(100)
            for j in range(len(S88sat)):
                S88sat[j] = steamTable.s_ph(P88sat[j],h8i[i])
                T88sat[j] = steamTable.t_ph(P88sat[j],h8i[i])

            ax3.plot(S88sat,T88sat,'-m')

        #graph
        for i in range(nsout_IP):
            if x6i_IP[i] >= 1:
                #je calcule les états intermédiaire entre T6i et T6i à saturation
                T6isat = steamTable.tsat_p(p8i_IP[i])+1
                S6isat = steamTable.sV_t(T6isat)
                T6i_6isat = np.linspace(T6i_IP[i]-273.15,T6isat)
                S6i_6isat = np.zeros(len(T6i_6isat))
                for j in range(0,len(T6i_6isat)):
                    S6i_6isat[j] = steamTable.s_pt(p8i_IP[i],T6i_6isat[j])
                ax3.plot(S6i_6isat,T6i_6isat,'k')

                T6isat_8i = np.linspace(T6isat,T8i_IP[i]-273.15,100)
                S6isat_8i = np.linspace(S6isat,s8i_IP[i],100)
                ax3.plot(S6isat_8i,T6isat_8i,'k')
            if x6i_IP[i] < 1 :
                T6isat_8i = np.linspace(T6i_IP[i]-273.15,T8i_IP[i]-273.15,100)
                S6isat_8i = np.linspace(s6i_IP[i],s8i_IP[i],100)
                ax3.plot(S6isat_8i,T6isat_8i,'k')
        #vannes
        for i in range(1,nsout_IP):
            P88sat = np.linspace(p8i_IP[i],p8i_IP[i-1],100)
            T88sat = np.zeros(100)
            S88sat = np.zeros(100)
            for j in range(len(S88sat)):
                S88sat[j] = steamTable.s_ph(P88sat[j],h8i_IP[i])
                T88sat[j] = steamTable.t_ph(P88sat[j],h8i_IP[i])

            ax3.plot(S88sat,T88sat,'-m')

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

        T10i_IP =T8i_IP-TpinchEx #T102,T103,...,T1

        """
        Drum pump
        """
        p11 =steamTable.psat_t(T10i_IP[-1]-273.15)+1
        h11=v71*(p11-p71)*10**2/eta_SiC+h71#kJ/kg
        T11= steamTable.t_ph(p11,h11)+273.15#K
        s11 = steamTable.s_ph(p11,h11)
        x11= None # eau non saturée
        e11 = h11-T0*s11#kJ/kg
        results[:,3] = T11-273.15,p11,h11,s11,x11,e11

        p10i_IP = p11*np.ones(nsout_IP)
        h10i_IP =np.zeros(nsout_IP)
        s10i_IP = np.zeros(nsout_IP)
        for i in range(nsout_IP):
            h10i_IP[i] = steamTable.h_pt(p11,T10i_IP[i]-273.15)
            s10i_IP[i] = steamTable.s_pt(p11,T10i_IP[i]-273.15)
        x10i_IP = None
        e10i_IP = h10i_IP-T0*s10i_IP


        #calcul du dernier etat du set 1 avant le drum
        T10n,p10n,h10n,s10n,x10n,e10n = T10i[-1],p10i[-1],h10i[-1],s10i[-1],x10i,e10i[-1]
        #calcu de dernier du set 2 (avant la pompe d 'alimentation')
        T1,p1,h1,s1,x1,e1 = T10i_IP[-1],p10i_IP[-1],h10i_IP[-1],s10i_IP[-1],x10i_IP,e10i_IP[-1]
        v1 = steamTable.vL_p(p1)


        def function_FHW_LP(x):

            T101 = x[0]
            Xis = x[1:nsout+1] #X1,X2,X3

            Xdrum = x[nsout+1]
            Xis_IP = x[nsout+2:]


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

            h10ia_IP,T10ia_IP = np.append(h11,h10i_IP),np.append(T11,T10i_IP)
            # print('here',T10, T10ia,T10ia_IP)
            # print(steamTable.tsat_p(p10)+273.15,steamTable.tsat_p(p11)+273.15)
            #print(T2,T3)
            print("evolution",T10-273.15,T101-273.15,T10i-273.15,T71-273.15,T11-273.15,T10i_IP-273.15)


            #set1
            sumXi = sum(Xis)
            Fone = sumXi*(h8i[0]-h9)-(1+sumXi)*(h101-h10)
            Fa = np.zeros(nsout-1)


            for i in range(len(Fa)):
                a = i+1
                sumXi_aplus_n = sum(Xis[a:])
                Fa[i] = sumXi_aplus_n*(h8i[a]-h8i[a-1])+Xis[a-1]*(h6i[a-1]-h8i[a-1])-(1+sumXi)*(h10ia[a]-h10ia[a-1])# c est chaud pour les notations

            Flast = Xis[-1]*(h6i[-1]-h8i[-1])-(1+sumXi)*(h10i[-1]-h10i[-2])

            #set2
            if nsout_IP !=0:
                sumXi2 = sum(Xis_IP)
                Flast2 = (1+Xdrum+sumXi+sumXi2)*(h10i_IP[-1]-h10i_IP[-2])-Xis_IP[-1]*(h6i_IP[-1]-h8i_IP[-1])
                Fa_IP = np.zeros(nsout_IP-1)

                for i in range(len(Fa_IP)):

                    a = i+1
                    sumXi_IP_aplus_n = sum(Xis_IP[a:])
                    Fa_IP[i] = Xis_IP[a-1]*(h6i_IP[a-1]-h8i_IP[a-1])+sumXi_IP_aplus_n*(h8i_IP[a]-h8i_IP[a-1])-(1+Xdrum+sumXi+sumXi2)*(h10ia_IP[a]-h10ia_IP[a-1])
                    #Drum

                Fdrum = Xdrum*(h6_IP)+(1+sumXi)*h10i[-1]+sumXi2*h8i_IP[0]-(1+Xdrum+sumXi+sumXi2)*h71
                # print(h6_IP,h10i[-1],h8i_IP[0],h71,'here')
                # print(h6,h6i,h6_IP,h6i_IP,h_pre)
                Functions = np.append(Fone,Fa)
                Functions = np.append(Functions,Flast)
                Functions = np.append(Functions,Fa_IP)
                Functions = np.append(Functions,Flast2)
                Functions = np.append(Functions,Fdrum)

            else :

                Fdrum =  Xdrum*(h6_IP)+(1+sumXi)*h10i[-1]-(1+Xdrum+sumXi)*h71
                Functions = np.append(Fone,Fa)
                Functions = np.append(Functions,Flast)
                Functions = np.append(Functions,Fdrum)
                #print(Functions)
            return Functions
        print('au drum', "T6_IP",T6_IP-273.15,"T10n",T10i[-1]-273.15,"T8i_IP[0]",T8i_IP[0]-273.15,'T71',T71-273.15)
        #premiere estimation de T81
        print('here',steamTable.tsat_p(5.1))

        def initial(nsout,nsout_IP):

            T101 = T10+15
            initial_set1 = 0.1*np.ones(nsout)
            Xdrum = 0.1
            initial_set2 = 0.1*np.ones(nsout_IP)
            initial = np.append(T101,initial_set1)
            initial=np.append(initial,Xdrum)
            initial = np.append(initial,initial_set2)
            return initial
        #print('here',function_FHW_LP(initial(2,2)))


        Solutions= fsolve(function_FHW_LP,initial(nsout,nsout_IP))
        T101 = Solutions[0]

        X_bleedings = Solutions[1:]
        X_bleedings_LP = Solutions[1:nsout+1]
        Xdrum =Solutions[nsout+1]
        X_bleedings_IP = Solutions[nsout+2:]

        sum_X_LP = sum(X_bleedings_LP)
        sum_X_IP = sum(X_bleedings_IP)

        if (X_bleedings[0]<0):
            print('opposite flow in the subcooler T TpinchSub different from T pinch Ex')
            print('modify the pressure in the extraction pump p10 or lower TpinchEx or TpinchSub')
        if(T10i[-1]>steamTable.tsat_p(p10)+273.15):
            print('vapor formed before the activation pump, increase pressure at alimentation pump p10')

        print('Solutions',T101,X_bleedings,sum(X_bleedings))
        print(T10i)
        print(steamTable.tsat_p(p10)+273.15)

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

        #échange de chaleur isobare entre 8 et 9 + vanne

        T89 = np.linspace(T91,T8i[0],100)
        S89 = np.zeros(100)
        for i in range(100):
            S89[i] = steamTable.s_pt(p8i[0],T89[i])
        ax3.plot(S89,T89,'-m')

        #vanne
        P97 = np.linspace(p7,p91,100)
        T97 = np.zeros(100)
        S97 = np.zeros(100)
        for i in range(len(S97)):
            S97[i] = steamTable.s_ph(P97[i],h91)
            T97[i] = steamTable.t_ph(P97[i],h91)

        ax3.plot(S97,T97,'-m')


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
        Q2 +=  (sum_X_LP)/(1+sum_X_IP+Xdrum+sum_X_LP)*(h91-h7)

        #turbine work  :
        sumXbleedings = sum(X_bleedings)
        Wm_t += 1/(1+sumXbleedings)*(h6_IP-h6)
        Wm_t += (1+sum_X_LP+Xdrum)/(1+sumXbleedings)*(h_pre-h6_IP)
        Wm_tmax += 1/(1+sumXbleedings)*(e6_IP-e6)
        Wm_tmax += (1+sum_X_LP+Xdrum)/(1+sumXbleedings)*(e_pre-e6)
        L_turbine_mv += 1/(1+sumXbleedings)*T0*(s6-s6_IP)
        L_turbine_mv += (1+sum_X_LP+Xdrum)/(1+sumXbleedings)*T0*(s6_IP-s_pre)

        #LP
        for i in range(nsout):
            Wm_t += X_bleedings_LP[i]/(1+sumXbleedings)*(h6_IP-h6i[i])
            Wm_tmax += X_bleedings_LP[i]/(1+sumXbleedings)*(e6_IP-e6i[i])
            L_turbine_mv += X_bleedings_LP[i]/(1+sumXbleedings)*T0*(s6i[i]-s6_IP)
        #IP
        for i in range(nsout_IP):
            Wm_t += X_bleedings_IP[i]/(1+sumXbleedings)*(h_pre-h6i_IP[i])
            Wm_tmax += X_bleedings_IP[i]/(1+sumXbleedings)*(e_pre-e6i_IP[i])
            L_turbine_mv += X_bleedings_IP[i]/(1+sumXbleedings)*T0*(s6i_IP[i]-s_pre)

        print("done",nsout)
        print(T6i_IP-273.15,p6_IP,h6i_IP,s6i_IP,x6i_IP,e6i_IP)
        #formater le resultat rappel : results 7,10,71,11,1,2,3,[reheat :41,51,...,4last,pre ], 6, 61,62,...,6n, 81,82,...,8n,101,102,...,10n(1),91
        results[:,8+2*reheat:8+nsout+2*reheat]=T6i-273.15,p6i,h6i,s6i,x6i,e6i
        results[:,8+nsout+2*reheat] = T6_IP-273.15,p6_IP,h6_IP,s6_IP,1,e6_IP
        results[:,8+nsout+2*reheat+1:8+nsout+2*reheat+1+nsout_IP]=T6i_IP-273.15,p6i_IP,h6i_IP,s6i_IP,x6i_IP,e6i_IP

        results[:,8+nsout+2*reheat+1+nsout_IP:8+2*nsout+2*reheat+1+nsout_IP]=T8i-273.15,p8i,h8i,s8i,x8i*np.ones(nsout),e8i
        results[:,8+2*nsout+2*reheat+1+nsout_IP:8+2*nsout+2*reheat+1+2*nsout_IP]=T8i_IP-273.15,p8i_IP,h8i_IP,s8i_IP,x8i_IP*np.ones(nsout_IP),e8i_IP

        results[:,8+2*nsout+2*reheat+1+2*nsout_IP]=T101-273.15,p101,h101,s101,x101,e101

        results[:,8+2*nsout+2*reheat+2+2*nsout_IP:8+3*nsout+2*reheat+2+2*nsout_IP] =T10i-273.15,p10i,h10i,s10i,np.zeros(nsout),e10i
        results[:,8+3*nsout+2*reheat+2+2*nsout_IP:8+3*nsout+2*reheat+2+3*nsout_IP] =T10i_IP-273.15,p10i_IP,h10i_IP,s10i_IP,np.zeros(nsout_IP),e10i_IP

        results[:,8+3*nsout+2*reheat+2+3*nsout_IP]=T91-273.15,p91,h91,s91,x91,e91


    else:
        Wm_t += (h_pre-h6)
        Wm_tmax += (e_pre-e6)
        L_turbine_mv += T0*(s6-s_pre)
        T1,p1,h1,s1,x1,v1,e1 = T10,p10,h10,s10,x10,v7,e10
        h91 = h7
        e91 = e7

        p11 = 5*p71
        h11=v71*(p11-p71)*10**2/eta_SiC+h71#kJ/kg
        T11= steamTable.t_ph(p11,h11)+273.15#K
        s11 = steamTable.s_ph(p11,h11)
        x11= None # eau non saturée
        e11 = h11-T0*s11#kJ/kg

    LP_flowratio = 1+sum_X_LP
    flowratio = 1+sum_X_LP+sum_X_IP+Xdrum

    """
    7) Alimentation pump
    """

    results[:,4] = T1,p1,h1,s1,x1,e1

    p2=p3_hp#bar
    h2=v1*(p2-p1)*10**2/eta_SiC+h1#kJ/kg
    T2= steamTable.t_ph(p2,h2)+273.15#K
    s2 = steamTable.s_ph(p2,h2)
    x2= None # eau non saturée
    e2 = h2-T0*s2#kJ/kg
    results[:,5]=T2-273.25,p2,h2,s2,x2,e2

    Q1 += h3-h2
    Q_boiler_exergie += e3-e2
    print('Q1',Q1)

    """
    9) Condenser
    """

    Q2 += 1/(flowratio)*(h6-h7)

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

    Wm_c = (h2-h1)+(h11-h71)+LP_flowratio/flowratio*(h10-h7) #kJ/kg_v #*(1+X)/(1+X)
    Wm_cmax = e2-e1+(e11-e71)+LP_flowratio/flowratio*(e10-e7)
    Wm = Wm_t-Wm_c
    Wm_max = Wm_tmax-Wm_cmax

    """
    11) massflows
    """
    mv = Pe/(Wm*eta_mec) #kg_v/s
    Q_boiler = mv*Q1 #kW

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
    eta_transex =mv*Q_boiler_exergie/(mf*(e_boiler_in-e_boiler_out)) # echange au boiler

    eta_gex=mv*Q_boiler_exergie/(mc*ec)
    # eta_gex2= eta_combex*eta_chemex*eta_transex
    # print(eta_gex,eta_gex2)
    eta_totex = eta_cyclex*eta_gex*eta_mec
    #eta_totex2 = Pe/(mc*ec)
    #print(eta_totex,eta_totex2)
    #print("eta_combex",eta_combex,"eta_chemex",eta_chemex,'eta_transex',eta_transex,"eta_gex",eta_gex,"eta_rotex",eta_rotex,'eta_cyclex',eta_cyclex,"eta_totex",eta_totex,'eta_gen',eta_gen)
    eta_condex = (e_cond_water_in-e_cond_water_out)*massflow_condenser_coeff/(e7-e6)

    """
    16) Computation of exergy losses
    """


    L_comb = boiler_outputs.L_comb #kW
    L_HR = boiler_outputs.L_HR
    L_exhaust = boiler_outputs.L_exhaust

    L_boiler = mv*(-Q_boiler_exergie)+mf*(e_boiler_in-e_boiler_out) #à checekr

    m_prim_LP = 1/(1+sumXbleedings)*mv
    m_bleedings_LP = X_bleedings_LP/(1+sumXbleedings)*mv
    m_drum = Xdrum/(1+sumXbleedings)*mv
    m_bleedings_IP = X_bleedings_IP/(1+sumXbleedings)*mv

    L_turbine = L_turbine_mv*mv
    L_pump = T0*(s2-s1)*mv+T0*(s11-s71)*mv+T0*(s10-s7)*(1+sum(m_bleedings_LP))/(1+sumXbleedings)*mv#kW



    L_cond = m_prim_LP*(e6-e7)+sum(m_bleedings_LP)*(e91-e7)+mv*massflow_condenser_coeff*(e_cond_water_in-e_cond_water_out)#vanne included
    L_exchanger_soutex = m_prim_LP*e10+np.dot(m_bleedings_LP,e6i)+np.dot(m_bleedings_IP,e6i_IP)+m_drum*e6_IP-(mv*e1+sum(m_bleedings_LP)*e91)  #prend en compte les pertes au drum
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
    outputs.Xmassflow = m_bleedings_LP,m_drum,m_bleedings_IP

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
              '\n \n \n Condenser losses {v} [MW]'.format(v=round(L_cond/1000)),'\n \n \n \n Bleed heating losses {v} [MW]'.format(v=round(L_exchanger_soutex/1000)),
              'Boiler losses {v} [MW]'.format(v=round(L_boiler/1000)),'Heat recovery losses {v} [MW]'.format(v=round(L_HR/1000)),
              'Combustion losses {v} [MW]'.format(v=round(L_comb/1000)),'Chimney losses {v} [MW]'.format(v=round(L_exhaust/1000))]
    plt.savefig('figures/energie_pie.png')

    ax.pie(data,labels = labels,autopct="%1.2f%%",startangle = 90)
    ax.set_title("Primary exergetic flux "+ str(round(ec*mc/1000)) + "[MW]")
    plt.savefig('figures/exergie_pie.png')

    """
    19) Cycle graphs
    """

    #tracer la cloche

    T= np.linspace(5,373,1000)
    S_L=np.zeros(len(T))
    S_V = np.zeros(len(T))
    for i in range(len(T)):
        S_L[i]=steamTable.sL_t(T[i])
        S_V[i]=steamTable.sV_t(T[i])
    ax3.plot(S_L,T,'-r')
    ax3.plot(S_V,T,'-r')
    ax3.plot([S_L[-1],S_V[-1]],[T[-1],T[-1]],'-r')



    T2p = steamTable.tsat_p(p2)  ; S2p = steamTable.sL_p(p2) #°C,kJ/kg
    T2pp =  T2p ; S2pp = steamTable.sV_p(p2)
    T2p2pp = T2p*np.ones(100)
    S2p2pp = np.linspace(S2p,S2pp,100)
    ax3.plot(S2p2pp,T2p2pp,'g')

    S22p = np.linspace(s2,S2p,100)
    T22p = np.zeros(len(S22p))
    for i in range(0,len(S22p)):
        T22p[i] = steamTable.t_ps(p2,S22p[i])
    ax3.plot(S22p,T22p,'g')

    T2pp3 = np.linspace(T2pp,T3-273.15,100)+1
    S2pp3 = np.zeros(len(T2pp3))
    for i in range(0,len(T2pp3)):
        S2pp3[i] = steamTable.s_pt(p3,T2pp3[i])
    ax3.plot(S2pp3,T2pp3,'-g')

    p36 = np.linspace(p_pre,p6,100)
    S3p = np.linspace(s_pre,s_pre,100)
    T36 = np.zeros(len(p36))
    S36 = np.zeros(len(p36))
    h36s = np.zeros(len(p36))
    h36 = np.zeros(len(p36))
    for i in range(0,len(p36)):
        h36s[i] = steamTable.h_ps(p36[i],s_pre)
        h36[i] = h_pre-(h_pre-h36s[i])*eta_SiT
        T36[i] = steamTable.t_ph(p36[i],h36[i])
        S36[i] = steamTable.s_ph(p36[i],h36[i])
    ax3.plot(S36,T36,'g')#,S3p,T36,'b--')

    T67 = np.linspace(T36[-1],T7-273.15,100)
    S67 = np.linspace(S36[-1],s7,100)
    ax3.plot(S67,T67,'-b')

    #PE pump
    ax3.plot([s7,s10],[T7-273.15,T10-273.15])

    T110 =np.linspace(T10-273.15,T71-273.15,100)
    S110 = np.zeros(100)
    for i in range(len(T110)):
        S110[i]=steamTable.s_pt(p10,T110[i])
    ax3.plot(S110,T110,'-b')

    #PA pump
    ax3.plot([s1,s2],[T1-273.15,T2-273.15])

    #drum pump
    ax3.plot([s71,s11],[T71-273.15,T11-273.15],'-g',)

    #set 2 exchangers
    T110 =np.linspace(T11-273.15,T1-273.15,100)
    S110 = np.zeros(100)
    for i in range(len(T110)):
        S110[i]=steamTable.s_pt(p11,T110[i])
    ax3.plot(S110,T110,'-c')

    #what happens at the drum
    if x6_IP >= 1:
        #je calcule les états intermédiaire entre T6i et T6i à saturation
        T6isat = steamTable.tsat_p(p71)+1
        S6isat = steamTable.sV_t(T6isat)
        T6i_6isat = np.linspace(T6_IP-273.15,T6isat)
        S6i_6isat = np.zeros(len(T6i_6isat))
        for j in range(0,len(T6i_6isat)):
            S6i_6isat[j] = steamTable.s_pt(p71,T6i_6isat[j])
        ax3.plot(S6i_6isat,T6i_6isat,'tab:orange',label='Drum')

        T6isat_8i = np.linspace(T6isat,T71-273.15,100)
        S6isat_8i = np.linspace(S6isat,s71,100)
        ax3.plot(S6isat_8i,T6isat_8i,'tab:orange')
    if x6_IP< 1 :
        T6isat_8i = np.linspace(T6_IP-273.15,T71-273.15,100)
        S6isat_8i = np.linspace(s6_IP,s71,100)
        ax3.plot(S6isat_8i,T6isat_8i,'tab:orange',label='Drum')

    P88sat = np.linspace(p6_IP,p8i_IP[0],100)
    T88sat = np.zeros(100)
    S88sat = np.zeros(100)
    for j in range(len(S88sat)):
        S88sat[j] = steamTable.s_ph(P88sat[j],h8i_IP[0])
        T88sat[j] = steamTable.t_ph(P88sat[j],h8i_IP[0])

    ax3.plot(S88sat,T88sat,'-m')


    #ax3.plot(S22p,T22p,'g',S2p2pp,T2p2pp,'g',S2pp3,T2pp3,'g',S36,T36,'g',S3p,T36,'b--',S67,T67,'g')
    ax3.set_xlabel('Entropy [kJ/kg/K]')
    ax3.set_ylabel('Temperature [°C]')
    ax3.grid(True)
    #ax3.set_title('T S graph of the steam turbine cycle')
    ax3.legend()

    fig = [fig,fig2]
    outputs.fig = fig
    if (ST_inputs.DISPLAY == 1):
        plt.show()

    return outputs;


# ST_inputs = ST_arg.ST_inputs();
# ST_inputs.Pe = 250.0e3 #[kW]
# ST_inputs.DISPLAY = 1
# ST_inputs.nsout = 5
# ST_inputs.reheat = 5
# ST_inputs.p3_hp=100
# ST_inputs.p4 = 30
#
# answers = ST(ST_inputs);
# print(answers.Xmassflow)
# print(answers.massflow)
