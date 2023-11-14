import os, socket, math, sys
from devsim.python_packages.simple_dd import *
from devsim import *
sys.path.append('/home/ccc/devsim/QS_Micro_packages/')
from QSmodel_create import *
from QSsimple_physics import *


ece_name="Electron_Equation"
hce_name="Hole_____Equation"

from scipy.special import erfinv
from scipy.special import erfcinv 

#####
##### Constants
#####
### Gaussian Fermi Integral
def gfi(zeta,s):
    S = s**2
    sqrt2 = math.sqrt(2)
    H = sqrt2 / s * erfcinv(math.exp(-0.5 * S))
    if zeta < -S:
        sqrt2_pi = sqrt2/math.sqrt(math.pi)
        K = 2 * (1. - H / s * sqrt2_pi * math.exp(0.5 * S * (1. - H**2)))
        value = math.exp(0.5 * S + zeta) / (math.exp(K*(zeta+S)) + 1)
    else:
        value = 0.5 * math.erfc(-zeta / (s*sqrt2) * H)
    return value


### differential of Gaussian Fermi Integral
def difgfi(zeta,s):
    S = s**2
    sqrt2 = math.sqrt(2)
    H = sqrt2 / s * erfcinv(math.exp(-0.5 * S))
    if zeta<-S:
        sqrt2_pi = sqrt2/math.sqrt(math.pi)
        K = 2 * (1. - H / s * sqrt2_pi * math.exp(0.5 * S * (1. - H**2)))
        den_inv = 1. / (math.exp(K * (S + zeta)) + 1.)
        dvalue = math.exp(0.5 * S + zeta) * den_inv * (1. - K*math.exp(K * (S+zeta)) * den_inv)
    else:
        dvalue = H / s / math.sqrt(2*math.pi) * math.exp(-0.5 * (H*zeta)**2/S)
    return dvalue

### inverse function for Gaussian Fermi Integral
def gfiinv(g,s):
    # using the Newton method
    # The initial guess
    # perhaps the degenerate approximation
    # or provided from the last call
    # improves speed
    # g=np.float128(g)
    sqrt2 = np.float128(sqrt(2))
    S = np.float128(s**2)
    H = np.float128(sqrt2 / s * erfcinv(math.exp(-0.5 * S)))
    x = np.float128(s * sqrt(2) * math.erfinv(g*2-1) / H)

    if g >= 0.5 * erfc(S / (s*sqrt2) * H):
        return x 
    #x = s * sqrt(2) * erfcinv(g*2)
    #x = 0.5*erfinv(g*2-1)
    # x = 0.0
    print("x= %s"%x)
    rerr = sys.float_info.max
    i = 0
    
    while rerr > 1e-12 and i < 20:
        # calculate the error in our solution
        f = gfi(x, s) - g
        fp = difgfi(x, s)
        upd = -f / fp
        x += upd
        # place any bounds on x here
        # calculate relative error
        rerr = abs(upd)/(abs(x) + 1e-12)
        print("DEBUG %d %1.15e %1.15e" % (i, x, rerr))
        i += 1
    return x


def SetSchottkyContactParamaters2(device, region, contact, ThermVelocities_e=1e11, ThermVelocities_h=1e11):
  set_parameter(device=device, region=region, name="ThermVelocities_e",  value=ThermVelocities_e)
  set_parameter(device=device, region=region, name="ThermVelocities_h",  value=ThermVelocities_h)

def SetSchottkyContactParamaters(device, region, contact, ThermVelocities_h=1e11, ThermVelocities_e=1e11, Bandtype="Silicon"):
  #### Bandtype: Silicon or Gaussian
  # ThermVelocities 1e7 cm/s, for  organic 8.64e6 cm/s
  POffSet = GetContactPotentilOffset(device, contact)
  V_t = get_parameter(device=device, region=region, name="V_t")
  if Bandtype == "Silicon":
    n_i = get_parameter(device=device, region=region, name="n_i")
    V_t = get_parameter(device=device, region=region, name="V_t")
    Equi_Electrons= n_i * math.exp( POffSet/V_t)
    Equi_Holes    = n_i * math.exp(-POffSet/V_t)  # ThermVelocities_e=1.20e-6-/Equi_Electrons
    # ThermVelocities_h=1.20e-6/Equi_Holes
  elif Bandtype == "Gaussian":
    Lumo_BDensity = get_parameter(device=device, region=region, name = "GaussianLumo_BDensity")    
    Lumo_Width    = get_parameter(device=device, region=region, name = "GaussianLumo_Width")
    Lumo_Center   = get_parameter(device=device, region=region, name = "GaussianLumo_Center")
    Equi_Electrons= Lumo_BDensity*gfi((POffSet-Lumo_Center)/V_t, Lumo_Width/V_t)

    Homo_BDensity = get_parameter(device=device, region=region, name = "GaussianHomo_BDensity")
    Homo_Width    = get_parameter(device=device, region=region, name = "GaussianHomo_Width")
    Homo_Center   = get_parameter(device=device, region=region, name = "GaussianHomo_Center")
    Equi_Holes    = Homo_BDensity*gfi((Homo_Center-POffSet)/V_t, Homo_Width/V_t)
    print(Equi_Electrons, Equi_Holes)
    # input()
  set_parameter(device=device, region=region, name="ThermVelocities_e",  value = ThermVelocities_e)
  set_parameter(device=device, region=region, name="ThermVelocities_h",  value = ThermVelocities_h)
  set_parameter(device=device, region=region, name="Equi_Electrons",     value = Equi_Electrons)
  set_parameter(device=device, region=region, name="Equi_Holes",         value = Equi_Holes)

def CreateSchottkyDriftDiffusionAtContact(device, contact, element_contact=False): 
  electrons_schottky_equation ="ElectronCharge*ThermVelocities_e*(Electrons-Equi_Electrons) * ContactSurfaceArea/NodeVolume"
  holes_schottky_equation     ="ElectronCharge*ThermVelocities_h*(Holes-Equi_Holes) * ContactSurfaceArea/NodeVolume"
  ### The denominator in these are used to counteract the integration of "NodeVolume"
  # electrons_schottky_equation="pow(T,2)*ThermVelocities_e*(Electrons-Equi_Electrons)*ContactSurfaceArea/NodeVolume"
  # holes_schottky_equation="pow(T,2)*ThermVelocities_h*(Holes-Equi_Holes)*ContactSurfaceArea/NodeVolume"
  electrons_schottky_name="{0}_electrons_schottky".format(contact)
  holes_schottky_name="{0}_holes_schottky".format(contact)
  CreateContactNodeModel(device, contact, electrons_schottky_name, electrons_schottky_equation)
  CreateContactNodeModelDerivative(device, contact, electrons_schottky_name, electrons_schottky_equation, "Electrons")
  CreateContactNodeModel(device, contact, holes_schottky_name, holes_schottky_equation)
  CreateContactNodeModelDerivative(device, contact, holes_schottky_name, holes_schottky_equation, "Holes")
  if element_contact:
    contact_equation(device=device, contact=contact, name=ece_name, #variable_name="Electrons",
                node_model=electrons_schottky_name,  element_model="ElementElectronCurrent", element_current_model="ElementElectronCurrent")
    contact_equation(device=device, contact=contact, name=hce_name, #variable_name="Holes",
                node_model=holes_schottky_name,  element_model="ElementHoleCurrent", element_current_model="ElementHoleCurrent")

  else:
    contact_equation(device=device, contact=contact, name=ece_name, #variable_name="Electrons",
                  node_model=electrons_schottky_name, edge_model="ElectronCurrent", edge_current_model="ElectronCurrent")
    contact_equation(device=device, contact=contact, name=hce_name, #variable_name="Holes",
                  node_model=holes_schottky_name, edge_model="HoleCurrent", edge_current_model="HoleCurrent")

def CreateSchottkyElementElectronsContactEquation(device, contact):
  electrons_schottky_name="{0}_electrons_schottky".format(contact)
  contact_equation(device=device, contact=contact, name=ece_name, #variable_name="Electrons",
                  node_model=electrons_schottky_name,  element_model="ElementElectronCurrent", element_current_model="ElementElectronCurrent")

def CreateSchottkyElementHolesContactEquation(device, contact):
  holes_schottky_name="{0}_holes_schottky".format(contact)
  contact_equation(device=device, contact=contact, name=hce_name, #variable_name="Holes",
                  node_model=holes_schottky_name,  element_model="ElementHoleCurrent", element_current_model="ElementHoleCurrent")

def CreateSchottkyPotentilContact(device, contact, PotentialOffset):
  PotentialOffset=get_parameter(device=device, name = "%s_PotentialOffset"%contact)

def SetSchottkEquiCharge(device, region, contact):
  InitialContact_electrons="Initial%s_electrons"%contact
  InitialContact_holes="Initial%s_holes"%contact
  set_node_values(device=device, region=region, name=InitialContact_electrons, init_from="Electrons")
  set_node_values(device=device, region=region, name=InitialContact_holes, init_from="Holes")
