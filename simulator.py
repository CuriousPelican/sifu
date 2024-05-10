# Afin de simplifier la comprehension, les notations respectent au mieux celles introduites dans les travaux de Dean Wheeler

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

############################## Paramètres ##############################

# Remarque: Seuls les retours & entrées (cette section) utilisateur peuvent utiliser des unités hors SI pour simplifier la lecture.
# Les entrées sont ensuites converties en unités SI et toutes les variables et le programme fonctionnent avec des unités SI.

# Simulation

# Bouteille (lorsque l'on parle de la la bouteille, est a comprendre la chambre de pression & la tuyère)
Vb = 1 # volume (en L) de la bouteille
Rb = 40 # rayon (en mm)  de la bouteille
Rt = 4.5 # rayon (en mm) de la tuyere 
z_ret = 110 # hauteur (en mm) du rétrécissement de la bouteille (ie debut tuyère - début corps droit)
mb = 100 # masse a vide (en g) de la bouteille
Cd = 0.5 # coef de trainée (sans unité)

# Remplissage & pression
P = 7 # pression de remplissage de la bouteille (en bar)
Vl = 0.5 # volume (en L) du liquide de propulsion (eau) - < volume_bouteille

# Environnement
T = 25 # temperature (en °C)
P_atm = 1.01325 # pression atmosphérique (en bar)



############################## Conversion SI & Constantes ##############################

# Conversions SI
Vb /= 1000
Rb /= 1000
Rt /= 1000
z_ret /= 1000
mb /= 1000
P *= 100000
Vl /= 1000
T += 273.15
P_atm *= 100000

# Constantes
rho_l = 1000 # 999.972 kg/m^3 pour l'eau
R = 8.31 # 8.3144621 J/K/mol - constante gazs parfaits
M_air = 0.029 # 28.965 g/mol - masse molaire air (en kg/mol)
rho_air = 1.292*273.15/T # 1.292 kg/m^3 pour T=0°C



############################## Bouteille & son volume ##############################

# Referentiel Bouteille - z est la distance axiale par rapport à la sotrie de la tuyère (flux supposé unidirectionnel)

# Renvoie le rayon a z donné - basé sur puis arctan entre -param_tan et param_tan et 
def R(z):
  param_tan = 1.8 # meilleurs parametres trouvés pour correspondre a mes bouteilles
  z_goulot = 0.26*z_ret
  if z>=z_ret:
    return Rb
  elif z<=z_goulot:
    return Rt
  else:
    z_atan = param_tan*(2*(z-z_goulot)/(z_ret-z_goulot)-1)
    return Rt + (Rb-Rt)/2 * (1+np.arctan(z_atan)/np.arctan(param_tan))

# Renvoie la section a z donné
def A(z):
  return np.pi*R(z)**2


V_ret = integrate.quad(A, 0, z_ret)[0] # volume (en L) de z=0 a z=z_ret
print(f"V_ret = {round(V_ret*1000, 3)} L (volume de z=0 à z=z_ret={z_ret*1000} mm)")

z_max = (Vb-V_ret)/(np.pi*Rb**2) + z_ret # hauteur (en m) de la bouteille
print(f"Donc z_max={round(z_max*100, 1)}cm (hauteur totale bouteille)")

# Test courbe rétrécissement
# plt.plot([R(z) for z in np.linspace(0,z_ret,100)])
# plt.show()

H = "hum" # H(t) hauteur interface air-liquide
u = "hum" # u(z,t) vitesse axiale fluide/bouteille
rho = "hum" # rho(z,t) masse volumique



############################## 4 Later ##############################

# Referentiel Terrestre - y est la distance verticale entre le sol et la fusée (vol supposé purement vertical)
y = [] # y(t) altitude fusée
v = [] # v(t)
a = [] # a(t)

# Renvoie la norme de la force de trainée
def F_drag(Cd,rho,A,v):
  return 1/2*Cd*rho*A*v**2

# Renvoie la norme de la force de poussée
def F_thrust():
  return "hum"