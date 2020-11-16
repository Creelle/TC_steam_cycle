from thermochem import janaf
db = janaf.Janafdb();
import numpy as np;

"""
This file contains all the useful function for GT,GT2, exchanger and combustionGT.
It is using the data contained in Janafdb from thermochem.
"""
# import all the phase data from each elements thar are used in the gas turbine.
O2 = db.getphasedata('O2','g');
N2 = db.getphasedata('N2','g');
CO2 = db.getphasedata('CO2',phase ='g');
H2O = db.getphasedata('H2O',phase ='g');
CH4=db.getphasedata('CH4',phase ='g');
CO = db.getphasedata('CO',phase ='g');
Mm_O2 = 0.032;#kg/mol
Mm_N2 = 0.028;#kg/mol
conc_O2 = 0.21;# 21% in molar
conc_N2 = 0.79;# 79% in molar

"""
air_mixture(T) takes temperature T in Kelvin f.
It returns the global Cp, gamma and R for atmospheric air at different temperatures.
"""
def air_mixture(T):#kJ/kg/K
    Mm_a = conc_O2 * Mm_O2 + conc_N2 * Mm_N2;
    m_O2 = (conc_O2*Mm_O2)/Mm_a;# mass proportion of O2
    m_N2 = (conc_N2*Mm_N2)/Mm_a;
    cp_a = m_O2 * O2.cp(T) + N2.cp(T) * m_N2;#J/mol/K
    Cp = cp_a/Mm_a/1000;#kJ/kg/K
    R = 8.31/Mm_a/1000
    gamma = Cp/(Cp-R)
    return Cp,gamma,R;

"""
cp_air takes as input temperatures in kelvin, the vector conc_mass wich contains
the different mass concentrations of each compound of the gas and Mm_a wich is the molar mass
of the gas.
It returns the the cp of the gas divided by the molar mass at a given temperature.
"""
def cp_air(T,conc_mass,Mm_a):
    cps = np.array([N2.cp(T),CO2.cp(T),H2O.cp(T),O2.cp(T)])
    cp_air = np.dot(conc_mass,cps);#J/mol/K
    return cp_air/Mm_a #J/kg/K

"""
The functions bellow return the cp and the cp/T of each compounds involved in the Gas
turbine.
"""
def cpCH4(T):
    return CH4.cp(T)
def cpCH4_T(T):
    return CH4.cp(T)*(1/T)
def cpO2(T):
    return O2.cp(T)
def cpCO2(T):
    return CO2.cp(T)
def cpN2(T):
    return N2.cp(T)
def cpH2O(T):
    return H2O.cp(T)

"""
cp_mean computes the mean cp for a given compound f and an interval of temperatures T1,T2.
dt caracterises the precision of the function.
"""
def cp_mean(f,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values)/len(values)) #  cp_mean [J/mol/K]

"""
cp_air_T makes the same thing as cp_air but it divides the value by the temperature
"""
def cp_air_T(T,conc_mass,Mm_a):#J/kg/K
    return cp_air(T,conc_mass,Mm_a)/T;

"""
janaf_integrate_air makes the integration of the cp of a gas for a given temperature interval
[T1,T2].
"""
def janaf_integrate_air(f,conc_mass,Mm_a,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values,conc_mass,Mm_a)*dt)

"""
cp_mean_air makes the same thing as cp_mean but it is not for a single compound.
It has to be used for a gas.
"""
def cp_mean_air(f,conc_mass,Mm_a,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values,conc_mass,Mm_a)/len(values)) #  cp_mean [J/kg/K]

"""
janaf_integrate_air makes the integration of the cp of a compound for a given temperature interval
[T1,T2].
"""
def janaf_integrate(f,T1,T2,dt):
    values = np.arange(T1,T2,dt)
    return sum(f(values)*dt) # int(cp)dt [J/mol/K]]

"""
cp_Iconstants takes a compound and a temperature interval.
It interpolates a lot of cp of the compound into this interval by a third degree polynomial.
Once all the coefficients of the polynomial have been found, the integration is calculated.
It returns I the integration ,cp_mean (not used) and the coefficient of integration.
"""
def cp_Iconstants(M,T_0,T_1):
    # donne la valeur de cp en J/mol.
    #cp = A + BT + CT^2 + DT^3 (cp évolue en T^3 (sauf au début => à débattre dans le rapport : faire un calcul d'erreur))
    n=1000 # n est arbitraire : j'ai pris n points pour avoir le plus de précision possible
    #db.getphasedata est différente pour l'oxygène que pour les autres gazs donc je dois mettre un if.
    if M=='O2' or M=='N2':
        chim = db.getphasedata(M,'g');
        T = np.linspace(T_0,T_1,n) #je divise en espaces réguliers l intervalle de températures
        cp_values = np.zeros(n)
        for i in range (0,n):
            cp_values[i] = chim.cp(T[i])  #je prends toutes les valeurs de cp correspondant aux températures
        heat_const = np.polyfit(T,cp_values,3) #je trouve les constantes A,B,C,D de de l'interpolation polynomiale grâce à polyfit
        cp_mean =  sum(cp_values)/len(cp_values)
        I = heat_const[3]*(T_1-T_0) + (1/2)*heat_const[2]*(T_1**2-T_0**2) + (1/3)*heat_const[1]*(T_1**3-T_0**3) + (1/4)*heat_const[0]*(T_1**4-T_0**4)
        #j'intègre cp sur l'intervalle de température
    else:
        chim = db.getphasedata(M,phase ='g');
        T = np.linspace(T_0,T_1,n)
        cp_values = np.zeros(n)
        for i in range (0,n):
            cp_values[i] = chim.cp(T[i])
        heat_const = np.polyfit(T,cp_values,3)
        cp_mean =  sum(cp_values)/len(cp_values)
        I = heat_const[3]*(T_1-T_0) + (1/2)*heat_const[2]*(T_1**2-T_0**2) + (1/3)*heat_const[1]*(T_1**3-T_0**3) + (1/4)*heat_const[0]*(T_1**4-T_0**4)
        #je retourne l'intégrale
    return I,cp_mean,heat_const
# T_null = 0.1
# T_reference = 298.15
# dt = 0.1
# LHV = (CO2.DeltaH(273.15)+2*H2O.DeltaH(273.15)-CH4.DeltaH(273.15)-2*O2.DeltaH(273.15))/0.016
# HHV = LHV-2*40752/0.016
# print(HHV/1000)
# LHV = (CO2.DeltaH(273.15)-0.5*O2.DeltaH(273.15)-CO.DeltaH(273.15))/0.028
# print(LHV/1000)
