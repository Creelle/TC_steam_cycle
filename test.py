from thermochem import janaf
db = janaf.Janafdb();
import numpy as np;
import matplotlib.pyplot as plt;

import GT_arguments as GT_arg;
import combustionGT as comb;

O2 = db.getphasedata('O2','g');
N2 = db.getphasedata('N2','g');
CO2 = db.getphasedata('CO2',phase ='g');
H2O = db.getphasedata('H2O',phase ='g');
Mm_O2 = 0.032;#kg/mol
Mm_N2 = 0.028;#kg/mol
conc_O2 = 0.21;# 21% in molar
conc_N2 = 0.79;# 79% in molar

#hello

def air_mixture(T):#kJ/kg/K
    Mm_a = conc_O2 * Mm_O2 + conc_N2 * Mm_N2;
    m_O2 = (conc_O2*Mm_O2)/Mm_a;# mass proportion of O2
    m_N2 = (conc_N2*Mm_N2)/Mm_a;
    cp_a = m_O2 * O2.cp(T) + N2.cp(T) * m_N2;#J/mol/K
    Cp = cp_a/Mm_a/1000;#kJ/kg/K
    R = 8.31/Mm_a/1000
    gamma = Cp/(Cp-R)
    return Cp,gamma;

def cp_air(T,conc_mass,Mm_a):
    cps = np.array([N2.cp(T),CO2.cp(T),H2O.cp(T),O2.cp(T)])
    molar_mass = np.array([0.028,0.044,0.018,0.032])
    cp_air = np.dot(conc_mass,cps);#J/mol/K
    return cp_air/Mm_a #J/kg


#fonction qui donne l enthalpie (kJ/kg), T temperature concentration massique mass : array (N2 , CO2, H20,O2)
#molar mass
def air_enthalpy(T,conc_mass,Mm_a): #==> a chequer si c est pas diviser par Mm_a ou divisé par molar_mass
    enthalpies = np.array([N2.hef(T),CO2.hef(T),H2O.hef(T),O2.hef(T)])
    molar_mass = np.array([0.028,0.044,0.018,0.032])
    h_air = sum(conc_mass*enthalpies);#kJ/mol
    return h_air/Mm_a #kJ/kg

def air_entropy(T,conc_mass,Mm_a):
    entropies = np.array([N2.S(T),CO2.S(T),H2O.S(T),O2.S(T)])
    S_air = sum(conc_mass*entropies);#J/mol/K
    return S_air/Mm_a #kJ/kg
def cp_air_T(T,conc_mass,Mm_a):#J/kg/K

    return cp_air(T,conc_mass,Mm_a)/T;

# def exergy_air(T,conc_mass,Mm_a):
#     T0=288.15
#     #molar_mass = np.array([0.028,0.044,0.018,0.032])
#     enthalpies = np.array([N2.hef(T),CO2.hef(T),H2O.hef(T),O2.hef(T)])
#     entropies = np.array([N2.S(T),CO2.S(T),H2O.S(T),O2.S(T)])
#     enthalpies0 = np.array([N2.hef(T0),CO2.hef(T0),H2O.hef(T0),O2.hef(T0)])
#     entropies0 = np.array([N2.S(T0),CO2.S(T0),H2O.S(T0),O2.S(T0)])
#     exergies = (enthalpies-enthalpies0)*1000-T0*(entropies-entropies0) #J/mol
#     e_air = sum(conc_mass*exergies)/1000/Mm_a #kJ/kg
#     return e_air #kJ/kg

def janaf_integrate_air(f,conc_mass,Mm_a,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values,conc_mass,Mm_a)*dt)
def cp_mean_air(f,conc_mass,Mm_a,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values,conc_mass,Mm_a)/len(values)) #  cp_mean [J/kg/K]
def janaf_integrate(f,T1,T2,dt): #==> pour calculer enthalpie
    values = np.arange(T1,T2,dt)
    return sum(f(values)*dt) # int(cp)dt [J/mol/K]]


def GT_simple(GT_input):
    """
     GT Gas turbine modelisation
     GT(P_e,options,display) compute the thermodynamics states for a Gas
     turbine based on several inputs (given in OPTION) and based on a given
     electricity production P_e. It returns the main results. It can as well
     plots graphs if input argument DISPLAY = true (<=> DISPLAY=1)

     INPUTS (some inputs can be dependent on others => only one of these 2 can
             be activated) Refer to Fig 3.1 from reference book (in english)
     P_E = electrical power output target [kW]
     OPTIONS is a structure containing :
       -options.T_ext [°C] : External temperature
       -options.r     [-] : Comperssion ratio
                            chamber
       -options.T_3   [°C] : Temperature after combustion (before turbine)
       -option.eta_PiC[-] : Intern polytropic efficiency (Rendement
                            polytropique interne) for compression
       -option.eta_PiT[-] : Intern polytropic efficiency (Rendement
                            polytropique interne) for expansion
        a faire:
        merge combustionGT avec simpleGT
        exergie pour simpleGT
        optimiser le taux de compression jouant sur le taux de compression et le lambda
        preheating ==> modelisation d un echangeur
        (humidity chequ)
        pychart
        des graphes T s et pv des etats dans la turbine
        (faire plusieurs etages de compression et analyse au niveau exergetique et energetique pour avoir
        si ca change quelque chose)

        a faire dans l immediat
        changer la formule de l entropy_air
        janaf integrate air
        formule d exergetique
        rendements exergetique
    """
    arg_in = GT_input;

    ## Check input arguments
    # ======================
    Pe = arg_in.Pe;
    if Pe ==-1.:
        Pe = 50e3;#50MWe
    T_ext = arg_in.T_ext;
    if T_ext ==-1.:
        T_ext = 288.15;#15°C
    r = arg_in.r;
    if r ==-1.:
        r = 10;#compression ratio = 10;
    T3 =arg_in.T3;
    eta_pic = arg_in.eta_PiC;
    if eta_pic ==-1.:
        eta_pic = 0.9;#max temperature = 1050°C
    eta_pit = arg_in.eta_PiT;
    if eta_pit ==-1.:
        eta_pit = 0.9;#max temperature = 1050°C
    T0=arg_in.T_0
    k_mec = arg_in.k_mec
    kcc = arg_in.k_cc

    """
    ## preliminary data (air) ==> find gamma
    # ======================
    # cp air at 15°C (298K): [kJ/mol/K]
    """
    Cp_a,gamma= air_mixture(T0)
    Mm_a = Mm_a = conc_O2 * Mm_O2 + conc_N2 * Mm_N2;
    conc_mass1=np.array([conc_N2*Mm_N2/Mm_a,0,0,conc_O2*Mm_O2/Mm_a])
    Ra = 8.31/Mm_a
    #coeff polytroique (compresseur :m>gamma , turbine : gamma >m) en premiere estimation
    m_t = (-eta_pit*(gamma-1)/gamma+1)**(-1)
    m_c = (1-1/eta_pic*(gamma-1)/gamma)**(-1)
    #on va recalculer m_c et m_t en utilisant la definition du livre page 118 (3.20) (3.25 en faisant des iterations)


    """
    1) compressor

    """
    T1=T_ext # a changer lors du preaheating
    p1 = 1.0 #bar
    h1 = air_enthalpy(T1,conc_mass1,Mm_a)+air_enthalpy(T0+10,conc_mass1,Mm_a)- air_enthalpy(T0,conc_mass1,Mm_a) #car la ref est pris a 25°c et non 25°C
    s1 = air_entropy(T1,conc_mass1,Mm_a)-air_entropy(T0,conc_mass1,Mm_a) #car T0 est ma reference
    s12 = janaf_integrate_air(cp_air_T,conc_mass1,Mm_a,T0-15,T1,0.001)
    print("s1",s1,s12)
    e1 = h1-T0*s1/1000 #kJ/kg_in


    p2 = r*p1

    # calcul de T2 par iteration
    T2 = T1*(r)**((m_c-1)/m_c) #premiere estimation
    iter = 1
    error = 1
    dt = 0.1

    while iter < 50 and error >0.01 :#pg119
        exposant_c=1/eta_pic*(Ra)/cp_mean_air(cp_air,conc_mass1,Mm_a,T1,T2,0.01)  # na-1/na :  formule du livre (3.22)
        T2_new = T1*r**(exposant_c)
        iter=iter+1
        error = abs(T2_new-T2)
        T2=T2_new

    s2 = air_entropy(T2,conc_mass1,Mm_a)-air_entropy(T0,conc_mass1,Mm_a)-Ra*np.log(r)
    #s22 = janaf_integrate_air(cp_air_T,conc_mass1,Mm_a,T0,T2,0.001)
    h2 = air_enthalpy(T2,conc_mass1,Mm_a)+ air_enthalpy(T0+10,conc_mass1,Mm_a)- air_enthalpy(T0,conc_mass1,Mm_a)
    e2 = h2-T0*s2/1000

    deltah_c = h2-h1 #kJ/kg
    deltah_c2 = janaf_integrate_air(cp_air,conc_mass1,Mm_a,T1,T2,0.001)/1000 #kJ/kg
    print('enthalpy comparaison',deltah_c,deltah_c2)
    print(T2-T1)

    deltas_c1 = s2-s1
    deltas_c2 = janaf_integrate_air(cp_air_T,conc_mass1,Mm_a,T1,T2,0.001)
    print('entropy comparaison',deltas_c1,deltas_c2)

    delta_ex_c = e2-e1
    delta_ex_c2 = deltah_c - T0 * (deltas_c1)/1000
    print('exergie comparaison 1-2', delta_ex_c,delta_ex_c2)

    """
     2 ) combustion
    """

    p3 = p2*kcc

    comb_outputs = comb.combustionGT(GT_arg.comb_input(h_in=h2,T_in = T2,inversion=True,T_out=T3 ))
    T3=comb_outputs.T_out
    lambda_comb = comb_outputs.lambda_comb
    ma1 = comb_outputs.ma1
    Mm_af = comb_outputs.Mm_af
    Rf = comb_outputs.R_f
    conc_mass2 = np.array([comb_outputs.m_N2f,comb_outputs.m_CO2f,comb_outputs.m_H2Of,comb_outputs.m_O2f])
    print('lambda_comb',lambda_comb)
    h3 = air_enthalpy(T3,conc_mass2,Mm_af) +air_enthalpy(T0+10,conc_mass2,Mm_af)- air_enthalpy(T0,conc_mass2,Mm_af)#kJ/kg_f
    # h32 = cp_mean_air(cp_air,conc_mass2,Mm_af,T0,T3,0.001)*(T3-T0)
    # h33 = janaf_integrate_air(cp_air,conc_mass2,Mm_af,T0,T3,0.001)
    # print('h3',h3,h32,h33)
    massflow_coefficient = 1+1/(ma1*lambda_comb) #kg_fu/kg_air
    #print('h3-h2',massflow_coefficient*h3-h2, massflow_coefficient*janaf_integrate_air(cp_air,conc_mass2,Mm_af,T2,T3,0.001))
    s3 = air_entropy(T3,conc_mass2,Mm_af)-air_entropy(T0,conc_mass2,Mm_af)-Rf*np.log(kcc*r) #J/K/kg_f
    e3 = h3-T0*s3/1000 #kJ/kg_f
    delta_exer_comb = massflow_coefficient*e3-e2 #kJ/kg_air
    print('exergie 2-3',delta_exer_comb)
    """
    3)  detente
    """
    p4 = p3/(r*kcc)


    # calcul de T4 par iteration
    T4 = T3*(1/(r*kcc))**((m_t-1)/m_t)
    iter = 1
    error = 1
    dt = 0.1

    while iter < 50 and error >0.01 :#pg119
        exposant_t=eta_pit*(Rf)/cp_mean_air(cp_air,conc_mass2,Mm_af,T4,T3,0.01)  # na-1/na :  formule du livre (3.25)
        T4_new = T3*(1/(kcc*r))**(exposant_t)
        iter=iter+1
        error = abs(T4_new-T4)
        T4=T4_new


    h4 = air_enthalpy(T4,conc_mass2,Mm_af) +air_enthalpy(T0+10,conc_mass2,Mm_af)- air_enthalpy(T0,conc_mass2,Mm_af)# kJ/kg_f # pour fixer la ref a 15°C
    # h42 = janaf_integrate_air(cp_air,conc_mass2,Mm_af,T0,T4,0.001)
    # print('h4',h4,h42)
    deltah_t = h4-h3 #<0# kJ/kg_f
    s4 = air_entropy(T4,conc_mass2,Mm_af)-air_entropy(T0,conc_mass2,Mm_af) # kJ/kg_f
    e4 = h4-T0*s4/1000# kJ/kg_f
    delta_exer_t = e4-e3# kJ/kg_f
    deltas_t = s4-s3# kJ/kg_f
    deltas_t2 = -janaf_integrate_air(cp_air_T,conc_mass2,Mm_af,T4,T3,0.001)# kJ/kg_f
    print('exergie 3-4',delta_exer_t, deltah_t-T0*deltas_t/1000,deltah_t-T0*deltas_t2/1000)
    """
    4) travail moteur
    """
    #travail moteur
    Wm = -(deltah_c+(1+1/(lambda_comb*ma1))*deltah_t) #kJ/kg_in
    #apport calorifique
    Q_comb = massflow_coefficient*h3-h2 #kJ/kg_in
    #formula chequ of the book
    # 3.27#Wm2 = cp_mean_air(cp_air,conc_mass1,Mm_a,T1,T2,0.001)*T1*(massflow_coefficient*(cp_mean_air(cp_air,conc_mass2,Mm_af,T4,T3,0.001)/cp_mean_air(cp_air,conc_mass1,Mm_a,T1,T2,0.001))*T3/T1*(1-(1/(kcc*r))**exposant_t)-(r**exposant_c-1))
    Q_comb2 =  cp_mean_air(cp_air,conc_mass1,Mm_a,T1,T2,0.001)*T1*(massflow_coefficient*cp_mean_air(cp_air,conc_mass2,Mm_af,T2,T3,0.001)/cp_mean_air(cp_air,conc_mass1,Mm_a,T1,T2,0.001)*T3/T1-r**exposant_c)
    print('test',Q_comb,Q_comb2)
    """
    5) rendements cyclen et mass flux
    """
    ##====================
    # calcul des rendements
    eta_cyclen  =Wm/Q_comb
    eta_mec =1-k_mec* (massflow_coefficient*abs(deltah_t)+deltah_c)/(massflow_coefficient*abs(deltah_t)-deltah_c) # Pe/Pm = 1-k_mec*(Pmt+Pmc)/(Pmt-Pmc)
    print('eta_cyclen',eta_cyclen, 'eta_mec',eta_mec)
    #massflow calcul # on neglige m  flow combustion
    mf_in = Pe/(Wm*eta_mec)#kg/s
    mf_out = mf_in*massflow_coefficient
    mf_c = mf_out-mf_in
    """
    5) calcul des flux de puissances
    """
    """
    6) calcul de n_mec, n_toten , n_gen
    """
    """
    7) calcul des pertes
    """
    """
    8) calcul des rendements exergetique
    """
    """
    last) define output arguments
    """
    outputs = GT_arg.GT_outputs();
    outputs.eta[0] = eta_cyclen;
    #outputs.eta[1] = eta_toten;
    outputs.dat[0:]= [[T1,T2,T3,T4],[p1,p2,p3,p4],[h1,h2,h3,h4],[s1,s2,s3,s4],[e1,e2,e3,e4]]
    outputs.massflow[0:] = [mf_in,mf_c,mf_out]
    outputs.combustion.fum[0:]=np.array([comb_outputs.m_O2f,comb_outputs.m_N2f,comb_outputs.m_CO2f,comb_outputs.m_H2Of])*mf_out

    return outputs;
"""
attention, la temperature de reference dans janaf n est pas 288.15 mais 298.15
"""

GT_simple_outputs = GT_simple(GT_arg.GT_input(Pe = 230e3,T_ext=288.15,r=18.,T3 = 1673.15));
print(GT_simple_outputs.dat)
print(GT_simple_outputs.massflow)
