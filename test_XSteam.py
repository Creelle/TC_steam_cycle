# Install XSteam:
# pip install pyXSteam

from pyXSteam.XSteam import XSteam
steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
# t   Temperature (°C)
# p   Pressure (bar)
# h   Enthalpy (kJ/kg)
# s   Specific entropy (kJ/(kg °C))
# x   Vapor fraction

## Example :
# get T° from other data
#Saturation temperature at 1.014 bar : 100°C
print(steamTable.tsat_p(30))
#Temperature at pressure = 300 bar and enthalpy of 3000 kJ/kg/K : 100°C
print(steamTable.t_ph(300,3000))

print("gamma",steamTable.tsat_p(1.01417),steamTable.CpV_p(1.01417)/steamTable.CvV_p(1.01417))
print(steamTable.CpL_t(15))
# get T° from other data
#Saturation pressure at 100°C : 1.014 bar
print(steamTable.psat_t(100))

#les inputs et outputs de XSteam sont en °C
# ...
