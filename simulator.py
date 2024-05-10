# Afin de simplifier la comprehension, les notations respectent au mieux celles introduites dans les travaux de Dean Wheeler

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

############################## Paramètres & Coonstantes ##############################

# Simulation


# Bouteille
volume_bouteille = 1 # volume (en L) de la bouteille
rayon_bouteille = 40 # rayon (en mm)  de la bouteille
rayon_tuyere = 4.5 # rayon (en mm) de la tuyere 
hauteur_retrecissement = 80 # hauteur (en mm) du retrecissement de la bouteille (ie debut bouteille - début corps)
masse_vide = 100 # masse a vide (en g) de la bouteille
Cd = 0.5 # coef de trainée (sans unité)

# Remplissage & pression
P = 7 # pression de remplissage de la bouteille (en bar)
volume_liquide = 0.5 # volume (en L) du liquide de propulsion (eau) - < volume_bouteille

# Environnement
T = 25 # temperature (en °C)
P_atm = 1.01325 # pression atmosphérique (en bar)

# Conversions
volume_bouteille /= 1000
rayon_bouteille /= 1000
rayon_tuyere /= 1000
hauteur_retrecissement /= 1000
masse_vide /= 1000
P *= 100000
volume_liquide /= 1000
T += 273.15
P_atm *= 100000

# Constantes
rho_liquide = 1000 # 999.972 kg/m^3 pour l'eau
R = 8.31 # 8.3144621 J/K/mol - constante gazs parfaits
M_air = 0.029 # 28.965 g/mol - masse molaire air (en kg/mol)
rho_air = 1.292*273.15/T # 1.292 kg/m^3 pour T=0°C

############################## Bouteille & son volume ##############################
# Referentiel Bouteille - z est la distance axiale par rapport au fond de la tuyère (flux supposé unidirectionnel)

# Renvoie le rayon a z donné - basé sur arctan entre -2 et 2
def R(z):
  param_tan = 3
  if z>=hauteur_retrecissement:
    return rayon_bouteille
  else:
    z_atan = param_tan*(2*z/hauteur_retrecissement-1)
    return rayon_tuyere + (rayon_bouteille-rayon_tuyere)/2 * (1+np.arctan(z_atan)/np.arctan(param_tan))

def R2(z):
  hauteur_goulot = 0.2*hauteur_retrecissement
  if z>=hauteur_retrecissement:
    return rayon_bouteille
  elif z<=hauteur_goulot:
    return rayon_tuyere
  else:
    z_atan = 4*(z-hauteur_goulot)/(hauteur_retrecissement-hauteur_goulot)-2
    return rayon_tuyere + (rayon_bouteille-rayon_tuyere)/2 * (1+np.arctan(z_atan)/np.arctan(2))

# Renvoie la section a z donné
def A(z):
  return 2*np.pi*R(z)**2

V_retrecissement = integrate.quad(A, 0, hauteur_retrecissement)
print(hauteur_retrecissement)
print(V_retrecissement)

plt.plot([R(z) for z in np.linspace(0,2*hauteur_retrecissement,100)])
plt.show()

z_max = "hum" # hauteur bouteille
H = [] # H(t) hauteur interface air-liquide
R = [] # R(z) rayon de la bouteille, fonction de z
u = "hum" # u(z,t) vitesse axiale fluide/bouteille
rho = "hum" # rho(z,t) masse volumique
A = [np.pi()*R[i]**2 for i in range(len(R))] # A(z) Section de la bouteille





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