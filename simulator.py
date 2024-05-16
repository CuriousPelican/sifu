# Afin de simplifier la comprehension, les notations respectent au mieux celles introduites dans les travaux de Dean Wheeler

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import time

############################## Paramètres ##############################

# Remarque: Seuls les retours & entrées (cette section) utilisateur peuvent utiliser des unités hors SI pour simplifier la lecture.
# Les entrées sont ensuites converties en unités SI et toutes les variables et le programme fonctionnent avec des unités SI.

# Simulation
dz = 0.01 # pas (en mm) utilisé pour l'integration du volume de la bouteille, defaut 10µm

# Bouteille (lorsque l'on parle de la la bouteille, est a comprendre la chambre de pression & la tuyère)
Vb = 1 # volume (en L) de la bouteille
Rb = 40 # rayon (en mm)  de la bouteille
Rt = 4.5 # rayon (en mm) de la tuyere 
z_ret = 110 # hauteur (en mm) du rétrécissement de la bouteille (ie debut tuyère - début corps droit)
mb = 100 # masse a vide (en g) de la bouteille
Cd = 0.5 # coef de trainée (sans unité)

# Remplissage & pression
P = 7 # pression de remplissage de la bouteille (en bar)
Vw0 = 0.5 # volume (en L) initial du liquide de propulsion (eau) - doit être < volume_bouteille

# Lanceur
R_TO = 4 # rayon (en mm) extérieur du tube de lancement
R_TI = 3 # rayon (en mm) intérieur du tube de lancement
L_T = 200 # longueur (en mm) du tube de lancement
Vl = 0.5 # volume (en L) du lanceur

# Environnement
T = 25 # temperature (en °C)
P_atm = 1.01325 # pression atmosphérique (en bar)



############################## Conversion SI & Constantes ##############################

# Conversions SI
dz /= 1000
Vb /= 1000
Rb /= 1000
Rt /= 1000
z_ret /= 1000
mb /= 1000
P *= 100000
Vw0 /= 1000
T += 273.15
P_atm *= 100000
R_TO /= 1000
R_TI /= 1000
L_T /= 1000
Vl /= 1000

# Constantes
rho_w = 1000 # 999.972 kg/m^3 pour l'eau
R = 8.31 # 8.3144621 J/K/mol - constante gazs parfaits
M_air = 0.029 # 28.965 g/mol - masse molaire air (en kg/mol)
rho_air = 1.292*273.15/T # 1.292 kg/m^3 pour T=0°C



############################## R(z), A(z) ##############################

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

# Test courbe rétrécissement
# plt.plot([R(z) for z in np.linspace(0,z_ret,100)])
# plt.show()

############################## DISCRETISATION H(Vw,e), rho(z), m_tot(H), Vw(H) ##############################

# Methode simpson & discretisation z
# Renvoie la liste de z, A(z) & V(z) en partant de 0 discretisés, ainsi que V_ret & z_max
def init_bouteille():
  # Calculs avant rétrécissement
  z_liste1 = np.arange(0,z_ret+dz,dz)
  A_liste1 = np.array([A(z) for z in z_liste1])
  V_liste1 = integrate.cumulative_simpson(A_liste1,dx=dz) # liste des volumes, V[i] = volume entre 0 et z_liste1[i]

  # V_ret & z_max
  V_ret = V_liste1[-1] # volume (en m^3) de z=0 a z=z_ret
  print(f"V_ret = {round(V_ret*1000, 3)} L (volume de z=0 à z=z_ret={z_ret*1000} mm)")
  z_max = (Vb-V_ret)/(np.pi*Rb**2) + z_ret # hauteur (en m) de la bouteille
  print(f"Donc z_max={round(z_max*100, 1)}cm (hauteur totale bouteille)")

  # Calculs après rétrécissement
  z_liste2 = np.arange(z_ret+dz,z_max+dz,dz)
  A_liste2 = np.array([A(z) for z in z_liste2])
  V_liste2 = [np.pi*Rb**2*(z-z_ret) for z in z_liste2]

  z_liste = np.concatenate((z_liste1,z_liste2))
  A_liste = np.concatenate((A_liste1,A_liste2))
  V_liste = np.concatenate((V_liste1,V_liste2))

  return (z_liste,A_liste,V_liste,V_ret,z_max)

(z_liste,A_liste,V_liste,V_ret,z_max) = init_bouteille()

############################## Lanceur & (y,v,a) ##############################

A_TO = np.pi*R_TO**2 # section ext tube lancement
A_TI= np.pi*R_TI**2 # section int tube lancement

# Referentiel Terrestre - y est la distance verticale entre le sol et la fusée (vol supposé purement vertical)
y = 0 # y(t) altitude fusée
v = 0 # v(t)
a = 0 # a(t)



############################## 4 Later ##############################

u = "hum" # u(z,t) vitesse axiale fluide/bouteille

# Renvoie la norme de la force de trainée
def F_drag(Cd,rho,A,v):
  return 1/2*Cd*rho*A*v**2

# Renvoie la norme de la force de poussée
def F_thrust():
  return "hum"



############################## UNUSED ##############################

# UNUSED - FUNCTION OBJECT H(Vw,e), rho(z), m_tot(H), Vw(H)
"""
# Methode fonctionelle avec scipy

V_ret = integrate.quad(A, 0, z_ret)[0] # volume (en m^3) de z=0 a z=z_ret
print(f"V_ret = {round(V_ret*1000, 3)} L (volume de z=0 à z=z_ret={z_ret*1000} mm)")
z_max = (Vb-V_ret)/(np.pi*Rb**2) + z_ret # hauteur (en m) de la bouteille
print(f"Donc z_max={round(z_max*100, 1)}cm (hauteur totale bouteille)")

# Renvoie H, hauteur de liquide avec une précision e
# pour une volume de liquide V donné par recherche par dichotomie sur V

def H(Vw,e):
  def f(z):
    return integrate.quad(A, 0, z)[0]-Vw
  a,b = 0,z_max
  while (f(b)-f(a))>e:
    m = (a+b)/2
    if f(m)*f(b)>0:
      b = m
    else:
      a = m
  return m,f(m)
H0 = H(Vw0,10e-4)[0] # H(t) hauteur interface air-liquide

# Renvoie rho à z & H donnés
def rho(z,H):
  if z>=H:
    return rho_air
  else:
    return rho_w

# Masse linéique axiale rho(z,H)*A(z)
def dm(z,H):
  return rho(z,H)*A(z)

# Renvoie la masse totale de la fusée à H donné
def m_tot(H):
  return mb + integrate.quad(dm, 0, z_max, args=(H))[0]

# Renvoie le volume de liquide dans la fusée à H donné
def Vw(H):
  return integrate.quad(dm, 0, H, args=(H))[0]
"""

# UNUSED - Intégrale du rayon entre a et b (objectif était le volume mais expression analytique trop complexe)
def R_int(a,b):
  param_tan = 1.8 # meilleurs parametres trouvés pour correspondre a mes bouteilles
  z_goulot = 0.26*z_ret
  # Primitive de R
  def G(z):
    a = -(z_ret-z_goulot)*np.log((param_tan**2+1)*(z_ret**2+z_goulot**2)+2*z_ret*z_goulot*(param_tan**2-1)+4*param_tan**2*z*(z-z_goulot-z_ret))
    b = 2*param_tan*(z_ret+z_goulot-2*z)*np.arctan(param_tan*(z_ret+z_goulot-2*z)/(z_ret-z_goulot))
    c = (Rb+Rt)*z/2
    return (Rb-Rt)*(a+b) + c
  return(G(b)-G(a))