import numpy as np
import math as math
import matplotlib.pyplot as plt
from scipy.integrate import odeint

tau     =   10    # valeur taux compression [-]
D       =   0.08    # valeur alesage   [m]
C       =   0.09    # valeur course [m]
L       =   0.18    # valeur longueur bielle@ [m]
mpiston =   0.8   # valeur masse piston@ [kg]
mbielle =   0.9   # valeur masse bielle  [kg]
Q       =   2800000  # valeur chaleur emise par fuel par kg de melange admis@ [J/kg_inlet gas]

R = C / 2 # Longueur de manivelle
V_critique = (( math.pi * D**2) / 4 ) * (2 * R) # OK
masse_air = 0.029 * V_critique/(8.31451*303.15) # OK
beta = L / R # rapport entre la longueur de la bielle et la longueur de la manivelle
V_min = (1 / ( tau - 1)) * V_critique
V_max = (tau / (tau -1)) * V_critique

###############################
#    Runge-Kutta (pour p)     #
###############################

def integrator(Xstart,Ustart,Xend,h,f):
    """
    h le pas entre chaque X
    f la fonction à intégrer
    """
    imax = int((Xend-Xstart)/h)
    X = Xstart + np.arange(imax+1)*h
    U = np.zeros((imax+1,3)); U[0,:] = Ustart
    for i in range(imax):
        K1 = f(U[i,:]       )
        K2 = f(U[i,:]+K1*h/2)
        K3 = f(U[i,:]+K2*h/2)
        K4 = f(U[i,:]+K3*h  )
        U[i+1,:] = U[i,:] + h*(K1+2*K2+2*K3+K4)/6     
    return X,U

def vol (theta):
    
    """""""""""""""""
    Calcul du V_output en fonction de thêta
    
    """""""""""""""""
    V_output = np.zeros(theta.size)

    for i in range (len(theta)):
        V_output[i] = V_critique/2 * (1 - math.cos(theta[i]*math.pi/180) + beta - np.sqrt(beta**2 - math.sin(theta[i]*math.pi/180)**2)) + V_min
        
    """""""""""""""""
    Graphe de V_output en fct de theta
    
    """""""""""""""""
    plt.plot(theta, V_output/V_max)
    
    # Axes
    fs_text = 16 # Taille du texte
    plt.xlabel("$Thêta$ [degré]", fontsize=fs_text)
    plt.ylabel("$V_{output}$ [m^3]", fontsize=fs_text)
    
    # Titre
    plt.title("Volume/V_max en fonction de l'angle de vilebrequin", fontsize=fs_text)
    plt.show()

    return V_output

def q(theta):
    """""""""""""""""
    Calcul du Q_output en fonction de thêta 

    """""""""""""""""
    Q_output = np.zeros(theta.size)
    Q_max = float('-inf')
    Q_min = float('inf')
    for i in range (len(theta)):
        if (theta[i] < thetaC or theta [i] > thetaC + deltaThetaC):
            Q_output[i] = 0
        
        else:
            fraction = (theta[i] + thetaC)/deltaThetaC # Pas besoin de mettre en radiant car fraction.
            Q_output[i] = (Q *s* masse_air)/ 2 * (1 - math.cos(math.pi * fraction))
            if (Q_output[i] > Q_max) :
                Q_max = Q_output[i] 
                print(i)
            if (Q_output[i] < Q_min) :
                Q_min = Q_output[i]
    
    
    """""""""""""""""
    Graphe de Q_output en fct de theta
    
    """""""""""""""""
    plt.plot(theta, Q_output)
    
    # Axes
    fs_text = 16 # Taille du texte
    plt.xlabel("$Thêta$ [degré]", fontsize=fs_text)
    plt.ylabel("$Q_{output}$ [J]", fontsize=fs_text)
    
    # Titre
    plt.title("Q_max = " + str(Q_max) + "    Q_min = " + str(Q_min))
    plt.suptitle("Chaleur en fonction de l'angle de vilebrequin", fontsize=fs_text)
    plt.show()

    return Q_output



def pression(theta, s, thetaC, deltaThetaC) : 
    """""""""""""""""
    Calcul du p_output en fonction de thêta 

    """""""""""""""""
    def model(p,theta):

        theta_rad = math.radians(theta)
        racine = np.sqrt(beta*beta - (math.sin(theta_rad)*math.sin(theta_rad)))
        V_formule = V_critique/2 * (1 - math.cos(theta_rad) + beta - racine) + V_min # Je peux pas prendre vokl(theta) var renvoie un array et j'ai besoin d'un float.
        fraction = (theta-thetaC)/deltaThetaC
        
        dVdtheta = V_critique/2*(math.sin(theta_rad) + (math.sin(theta_rad)*math.cos(theta_rad))/(racine))
        
        dQdtheta = (math.pi*Q*masse_air*s)/(2*deltaThetaC) * (math.sin(math.pi*fraction))
        
        if (thetaC <= theta < thetaC + deltaThetaC) : # Si dans combustion
            dpdt = - 1.3 *p/V_formule * dVdtheta + (1.3 - 1)/V_formule * dQdtheta
        else :
            dpdt = - 1.3 *p/V_formule * dVdtheta
        
        return dpdt
 
    p_output = odeint(model,s,theta)
    plt.plot(theta,p_output)
    plt.title("Pression")
    plt.show()
     

def F_tete(theta, p_output, w) : 
    """""""""""""""""
    Calcul du F_tete_output en fonction de thêta 

    """""""""""""""""
    F_tete_output = np.zeros(theta.size)
    for k in range(len(theta)) : 
        F_tete_output[k] = -(math.pi*D**2)/4 *p_output[k] - (mpiston + mbielle)*R*w**2*math.cos(math.radians(theta[k]))
     
            
    plt.plot(theta, F_tete_output)
    plt.title("F_tete_output")
    plt.show()

    return F_tete_output


def F_pied(theta, p_output, w) : 
    """""""""""""""""
    Calcul du F_pied_output en fonction de thêta 

    """""""""""""""""
    F_pied_output = np.zeros(theta.size)
    for j in range(len(theta)) : 
        F_pied_output[j] = (math.pi*D**2)/4 *p_output[j] - mpiston*R*w**2*math.cos(math.radians(theta[j]))
        
    plt.plot(theta, F_pied_output)
    plt.title("F_pied_output")
    plt.show()

    return F_pied_output


def t() : 
    """""""""""""""""
    Calcul du t 

    """""""""""""""""
    pass
        

def myfunc(rpm, s, theta, thetaC, deltaThetaC):
    """ 
    dimBielle dimensionnement d'une bielle
    dimBielle(rpm, s, theta, thetaC, deltaThetaC) calcules les données thermodynamiques
    et les forces d'un système Piston-Bielle-Vilebrequin afin de dimensionner la section
    d'une bielle.

    INPUTS :
    rpm : vitesse angulaire du moteur [rotation per minute]
    s : surcharge du moteur [-]
    theta : angle auxquels renvoyer les données [°]
    thetaC : angle d'allumage [°]
    deltaThetaC : durée de la combustion (en angle) [°] 
    theta : angle auxquels renvoyer les données [°]
    (Les angles sont donnés entre 0 et 720°)

    OUTPUTS :
    t : section de la bielle [m]
    V(theta) : Volume de la chambre de combustion en fonction de theta [m3]
    Q(theta) : chaleur dégagée par la combustion en fonction de theta [J]
    F_pied(theta) : [N]
    F_tete(theta) : [N]
    F_inertie(theta) : [N]
    (une force est positive si dirigée vers le haut).
    """

    w = rpm/60*2*math.pi # Transformation de revolution/minute en rad/s.    

    V_output = vol (theta)
    Q_output = q(theta)
    p_output = pression(theta, s, thetaC, deltaThetaC)
    F_tete_output = F_tete(theta, p_output, w)
    F_pied_output = F_pied(theta, p_output, w)
    t = t()
    
    #V_output = np.ones_like(theta)
    #Q_output      = 2*np.ones_like(theta);
    #F_pied_output = 3*np.ones_like(theta);
    #F_tete_output = 4*np.ones_like(theta);
    #p_output      = 5*np.ones_like(theta);
    #t             = 6;

    return ( F_pied_output, F_tete_output, p_output, t  )



rpm = 2000
s = 1.5
Pin = s*1E5
thetaC = 25
deltaThetaC = 50
theta = np.linspace(-180., 180., 100)

#print(myfunc(2555,1.9,np.linspace(-1*180,180,10000), 26,43))
print(myfunc(rpm, s, theta, thetaC, deltaThetaC))

