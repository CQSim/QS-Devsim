import os, socket, tarfile, traceback, math, sys
from devsim.python_packages.simple_dd import *
from devsim import *
import pandas as pd
sys.path.append('/home/ccc/devsim/QS_Micro_packages/')
from QSmodel_create import *
hostname = socket.gethostname()


def InDisplayHost():
  HostList=("lenovo.localdomain","C-Twist","ccc-OptiPlex-3040","Turbo-King","TurboKing","ccc-Inspiron-660s")
  return hostname in HostList

def InNoneCompressHost():
  HostList=("pvd","C-Twist")
  return hostname in HostList

#TODO: make this a class so that paramters can be changed
contactcharge_node="contactcharge_node"
PotentialEdgeFlux="PotentialEdgeFlux"
ece_name="Electron_Equation"
hce_name="Hole_____Equation"
celec_model = "(n_i + 0.5*abs(NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"
chole_model = "(n_i + 0.5*abs(-NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"



def SetExtendedPrecision():
  set_parameter(name = "extended_solver",   value=True)
  set_parameter(name = "extended_model",    value=True)
  set_parameter(name = "extended_equation", value=True)
  set_parameter(name = "threads_available", value=4   )
  set_parameter(name = "threads_task_size", value=2048)

def SetCentiMeterBasicParameters():
  q     = 1.6021766e-19 # Coul
  Kb    = 1.3806503e-23 # J/K
  eps0  = 8.8541878e-14 # F/cm
  set_parameter(name = "ElectronCharge",    value=q   )
  set_parameter(name = "Kb",                value=Kb  )
  set_parameter(name = "eps0",              value=eps0)

def SetMicroMeterBasicParameters():
  q     = 1.6021766e-19 # Coul
  Kb    = 1.3806503e-23 # J/K
  eps0  = 8.8541878e-18 # F/um
  set_parameter(name = "ElectronCharge",    value=q   )
  set_parameter(name = "Kb",                value=Kb  )
  set_parameter(name = "eps0",              value=eps0)

def SetOxideParameters(device, region, rPermittivity=1.0):
  Permittivity = rPermittivity * get_parameter(name="eps0")
  set_parameter(device=device, region=region, name="Permittivity",   value=Permittivity  )


def SetDopingParameters(device, region, bulk_doping=0, body_doping=0, drain_doping=0, source_doping=0):
  set_parameter(device=device, region=region, name="bulk_doping",    value=bulk_doping   )
  set_parameter(device=device, region=region, name="drain_doping",   value=drain_doping  )
  set_parameter(device=device, region=region, name="source_doping",  value=source_doping )
  set_parameter(device=device, region=region, name="body_doping",    value=body_doping   )

def SetMicroMeterSiliconParameters(device, region, 
                          T             = 300   , #K
                          rPermittivity = 11.1  ,
                          n_i           = 1.0e-2, # /um^3
                          nMobility     = 1e8   , # um^2/sV
                          pMobility     = 1e8   , # um^2/sV
                          nLifeTime     = 1e-5  , # s
                          pLifeTime     = 1e-5  , # s
                          ):
  '''
    Sets physical parameters assuming constants
  '''
  #### TODO: make T a free parameter and T dependent parameters as models
  Permittivity = rPermittivity * get_parameter(name="eps0")
  ElectronCharge = get_parameter(name = "ElectronCharge")
  Kb  = get_parameter(name="Kb")
  V_t = Kb * T/ElectronCharge
  WriteSiliconParameters(device, region, Permittivity, n_i,n_i,n_i, T, V_t, nMobility, pMobility, nLifeTime, pLifeTime)

def SetCentiMeterSiliconParameters(device, region, 
                          T             = 300   , #K
                          rPermittivity = 11.1  ,
                          n_i           = 1.0e10, # /cm^3
                          nMobility     = 1     , # cm^2/sV
                          pMobility     = 1     , # cm^2/sV
                          nLifeTime     = 1e-5  , # s
                          pLifeTime     = 1e-5  , # s
                          ):
  '''
    Sets physical parameters assuming constants
  '''
  #### TODO: make T a free parameter and T dependent parameters as models
  Permittivity = rPermittivity * get_parameter(name="eps0")
  ElectronCharge = get_parameter(name = "ElectronCharge")
  Kb  = get_parameter(name="Kb")
  V_t = Kb * T/ElectronCharge
  WriteSiliconParameters(device, region,Permittivity, n_i,n_i,n_i, T, V_t, nMobility, pMobility, nLifeTime, pLifeTime)


def WriteSiliconParameters(device, region,Permittivity, n_i, n1, p1, T, V_t, nMobility, pMobility, nLifeTime, pLifeTime):
  set_parameter(device=device, region=region, name="Permittivity",   value=Permittivity  )
  set_parameter(device=device, region=region, name="T",              value=T    )
  set_parameter(device=device, region=region, name="n_i",            value=n_i  )
  set_parameter(device=device, region=region, name="n1",             value=n1)
  set_parameter(device=device, region=region, name="p1",             value=p1)

  # set_parameter(device=device, region=region, name="kT",             value=kT)
  set_parameter(device=device, region=region, name="V_t",            value=V_t)
  set_parameter(device=device, region=region, name="mu_n_eff",       value=nMobility)
  set_parameter(device=device, region=region, name="mu_p_eff",       value=pMobility)

  #default SRH parameters
  set_parameter(device=device, region=region, name="taun",           value=nLifeTime)
  set_parameter(device=device, region=region, name="taup",           value=pLifeTime)




#####
##### Constants
#####


def CreatePotentialAndFlux(device, region):
  '''
    Create electric field model in oxide
    Creates Potential solution variable if not available
  '''
  CreateSolution(device, region, "Potential")
  
  efield="(Potential@n0 - Potential@n1)*EdgeInverseLength"
  CreateEdgeModel(device, region, "ElectricField", efield)
  CreateEdgeModelDerivatives(device, region, "ElectricField", efield, "Potential")

  dfield="Permittivity*ElectricField"
  CreateEdgeModel(device, region, "PotentialEdgeFlux", dfield)
  CreateEdgeModelDerivatives(device, region, "PotentialEdgeFlux", dfield, "Potential")



def CreateSiliconPotentialOnly(device, region):
  '''
    Creates the physical models for a Silicon region
  '''
  IntrinsicElectrons  = "n_i*exp(Potential/V_t)"
  IntrinsicHoles      = "n_i^2/IntrinsicElectrons"
  IntrinsicCharge     = "kahan3(IntrinsicHoles, -IntrinsicElectrons, NetDoping)"
  PotentialIntrinsicCharge = "-ElectronCharge * IntrinsicCharge"

  # require NetDoping
  for i in (
       ("IntrinsicElectrons",       IntrinsicElectrons      ),
       ("IntrinsicHoles",           IntrinsicHoles          ),
       ("IntrinsicCharge",          IntrinsicCharge),
       ("PotentialIntrinsicCharge", PotentialIntrinsicCharge)):
    n = i[0]
    e = i[1]
    CreateNodeModel(device, region, n, e)
    CreateNodeModelDerivative(device, region, n, e, "Potential")

  equation(device=device, region=region, name="PotentialEquation", variable_name="Potential",
           node_model="PotentialIntrinsicCharge", edge_model="PotentialEdgeFlux", variable_update="log_damp")


def CreateSiliconPotentialOnlyContact(device, region, contact, is_circuit=False, offset=False):
  '''
    Creates the potential equation at the contact
    if is_circuit is true, than use node given by GetContactBiasName
  '''
  # Means of determining contact charge
  # Same for all contacts
  if not InNodeModelList(device, region, "contactcharge_node"):
    CreateNodeModel(device, region, "contactcharge_node", "ElectronCharge*IntrinsicCharge")
  #### TODO: This is the same as D-Field

  set_parameter(device=device, name=GetContactBiasName(contact), value=0.0)

  if offset:
    contact_model = "Potential -{0} + ifelse(NetDoping > 0, -V_t*log({1}/n_i), \
    V_t*log({2}/n_i))+{3}_PotentialOffset".format(GetContactBiasName(contact), celec_model, chole_model, contact)
  else:
    contact_model = "Potential -{0} + ifelse(NetDoping > 0, -V_t*log({1}/n_i), \
    V_t*log({2}/n_i))".format(GetContactBiasName(contact), celec_model, chole_model)
  
  contact_model_name = GetContactNodeModelName(contact)
  CreateContactNodeModel(device, contact, contact_model_name, contact_model)
  # Simplify it too complicated
  CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_model_name,"Potential"), "1")
  if is_circuit:
    CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_model_name,GetContactBiasName(contact)), "-1")

  if is_circuit:
    contact_equation(device=device, contact=contact, name="PotentialEquation", variable_name="Potential",
                   node_model=contact_model_name, edge_model="",
                   node_charge_model="contactcharge_node", edge_charge_model="PotentialEdgeFlux",
                   node_current_model="", edge_current_model="", circuit_node=GetContactBiasName(contact))
  else:
    contact_equation(device=device, contact=contact, name="PotentialEquation", variable_name="Potential",
                   node_model=contact_model_name, edge_model="",
                   node_charge_model="contactcharge_node", edge_charge_model="PotentialEdgeFlux",
                   node_current_model="", edge_current_model="")

def CreateSRH(device, region):
  USRH="(Electrons*Holes - n_i^2)/(taup*(Electrons + n1) + taun*(Holes + p1))"
  Gn = "-ElectronCharge * USRH"
  Gp = "+ElectronCharge * USRH"
  CreateNodeModel(device, region, "USRH", USRH)
  CreateNodeModel(device, region, "ElectronGeneration", Gn)
  CreateNodeModel(device, region, "HoleGeneration", Gp)
  for i in ("Electrons", "Holes"):
    CreateNodeModelDerivative(device, region, "USRH", USRH, i)
    CreateNodeModelDerivative(device, region, "ElectronGeneration", Gn, i)
    CreateNodeModelDerivative(device, region, "HoleGeneration", Gp, i)

def CreateECE(device, region, mu_n):
  CreateElectronCurrent(device, region, mu_n)
  NCharge = "ElectronCharge * Electrons"
  CreateNodeModel(device, region, "NCharge", NCharge)
  CreateNodeModelDerivative(device, region, "NCharge", NCharge, "Electrons")
  if InNodeModelList(device, region,"ElectronGeneration"):
    equation(device=device, region=region, name=ece_name, variable_name="Electrons",
      time_node_model = "NCharge", node_model="ElectronGeneration",
      edge_model="ElectronCurrent", variable_update="positive")
  else:
    equation(device=device, region=region, name=ece_name, variable_name="Electrons",
      time_node_model = "NCharge",
      edge_model="ElectronCurrent", variable_update="positive")

def CreateHCE(device, region, mu_p):
  CreateHoleCurrent(device, region, mu_p)
  PCharge = "-ElectronCharge * Holes"
  CreateNodeModel(device, region, "PCharge", PCharge)
  CreateNodeModelDerivative(device, region, "PCharge", PCharge, "Holes")
  if InNodeModelList(device, region,"HoleGeneration"):
    equation(device=device, region=region, name=hce_name, variable_name="Holes",
              time_node_model = "PCharge", node_model="HoleGeneration",
              edge_model="HoleCurrent", variable_update="positive")
  else:
    equation(device=device, region=region, name=hce_name, variable_name="Holes",
              time_node_model = "PCharge",
              edge_model="HoleCurrent", variable_update="positive")


#Create Net Charge Node Model
def CreatePE(device, region):
  pne = "-ElectronCharge*kahan3(Holes, -Electrons, NetDoping)"
  CreateNodeModel(device, region, "PotentialNodeCharge", pne)
  CreateNodeModelDerivative(device, region, "PotentialNodeCharge", pne, "Electrons")
  CreateNodeModelDerivative(device, region, "PotentialNodeCharge", pne, "Holes")

  equation(device=device, region=region, name="PotentialEquation", variable_name="Potential",
            node_model="PotentialNodeCharge", edge_model="PotentialEdgeFlux",
            time_node_model="", variable_update="log_damp")



def CreateElectronAndHoleSolutions(device, region):
  for c in ("Electrons","Holes"):
    CreateSolution(device, region, c)
    set_node_values(device=device, region=region, name=c, init_from="Intrinsic%s"%c)
    if c =="Electrons":
      expression="Potential+log(%s/n_i)*V_t"%c
    else :
      expression="Potential-log(%s/n_i)*V_t"%c
    CreateNodeModel(device, region, "QuasiFermi_%s"%c, expression)


def CreateSiliconDriftDiffusion(device, region, mu_n="mu_n_eff", mu_p="mu_p_eff",SRHMolde=True):
  CreatePE(device, region)
  CreateBernoulli(device, region)
  if SRHMolde:
    CreateSRH(device, region)
  CreateECE(device, region, mu_n)
  CreateHCE(device, region, mu_p)



def SetSchottkyContactParamaters(TherVelocities_h,Equi_Electrons,TherVelocities_e,Equi_Holes):
    set_parameter(device=device, region=region, name="TherVelocities_e",  value=TherVelocities_e)
    set_parameter(device=device, region=region, name="TherVelocities_h",  value=TherVelocities_h)
    set_parameter(device=device, region=region, name="Equi_Electrons",    value=Equi_Electrons)
    set_parameter(device=device, region=region, name="Equi_Holes",        value=Equi_Holes)

def CreateSchottkyDriftDiffusionAtContact(device, region, contact): 
    electrons_schottky_name="{0}_electrons_schottky".format(contact)
    electrons_schottky_equation="ElectronCharge*TherVelocities_e*(Electrons-Equi_Electrons)*ContactSurfaceArea/NodeVolume"
    holes_schottky_name="{0}_holes_schottky".format(contact)
    holes_schottky_equation="-ElectronCharge*TherVelocities_h*(Holes-Equi_Holes)*ContactSurfaceArea/NodeVolume"

    CreateContactNodeModel(device, contact, electrons_schottky_name, electrons_schottky_equation)
    CreateContactNodeModelDerivative(device, contact, electrons_schottky_name, electrons_schottky_equation, "Electrons")
    CreateContactNodeModel(device, contact, holes_schottky_name, holes_schottky_equation)
    CreateContactNodeModelDerivative(device, contact, holes_schottky_name, holes_schottky_equation, "Holes")

    contact_equation(device=device, contact=contact, name=ece_name, variable_name="Electrons",
                  node_model=electrons_schottky_name, edge_model="ElectronCurrent", edge_current_model="ElectronCurrent")
    contact_equation(device=device, contact=contact, name=hce_name, variable_name="Holes",
                  node_model=holes_schottky_name, edge_model="HoleCurrent", edge_current_model="HoleCurrent")

    # e_options=get_contact_equation_command(device=device, contact=contact, name=ece_name)
    # e_options['edge_model'] = "ElectronCurrent"
    # e_options['node_model'] = "{0}_electrons_schottky".format(contact)

    # h_options=get_contact_equation_command(device=device, contact=contact, name=hce_name)
    # h_options['edge_model'] = "HoleCurrent"
    # h_options['node_model'] = "{0}_holes_schottky".format(contact)

    # contact_equation(**e_options)
    # contact_equation(**h_options)

def CreateSiliconDriftDiffusionAtContact(device, region, contact, is_circuit=False): 
  '''
    Restrict electrons and holes to their equilibrium values
    Integrates current into circuit
  '''
  contact_electrons_model = "Electrons - ifelse(NetDoping > 0, {0}, n_i^2/{1})".format(celec_model, chole_model)
  contact_holes_model = "Holes - ifelse(NetDoping < 0, +{1}, +n_i^2/{0})".format(celec_model, chole_model)

  contact_electrons_name = "{0}nodeelectrons".format(contact)
  contact_holes_name = "{0}nodeholes".format(contact)

  CreateContactNodeModel(device, contact, contact_electrons_name, contact_electrons_model)
  CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_electrons_name, "Electrons"), "1")

  CreateContactNodeModel(device, contact, contact_holes_name, contact_holes_model)
  CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_holes_name, "Holes"), "1")
  
  #TODO: The simplification of the ifelse statement is time consuming
  if is_circuit:
    contact_equation(device=device, contact=contact, name=ece_name, variable_name="Electrons",  node_model=contact_electrons_name,
                         edge_current_model="ElectronCurrent", circuit_node=GetContactBiasName(contact))
    contact_equation(device=device, contact=contact, name=hce_name, variable_name="Holes",      node_model=contact_holes_name,
                         edge_current_model="HoleCurrent", circuit_node=GetContactBiasName(contact))
  else:
    contact_equation(device=device, contact=contact, name=ece_name, variable_name="Electrons",  node_model=contact_electrons_name,
                         edge_current_model="ElectronCurrent")
    contact_equation(device=device, contact=contact, name=hce_name, variable_name="Holes",      node_model=contact_holes_name,
                         edge_current_model="HoleCurrent")


def CreateContactRegionModel(device, Contact, ContactRegion):
    ### Set the contact 
    edge_model(device=device, region=ContactRegion, name="ElectricField", equation="0")
    node_model(device=device, region=ContactRegion, name="Potential", equation="%s_bias" % Contact)


#########
########Functions For Insulator:
#########

def CreateOxidePotentialOnly(device, region, update_type="default"):
  '''
    Creates Potential equation
  '''
  equation(device=device, region=region, name="PotentialEquation", variable_name="Potential",
      edge_model="PotentialEdgeFlux", variable_update=update_type)

def CreateElementOxidePotentialOnly(device, region, update_type="default"):
  '''
    Creates Element Model and Equation solution variable if not available
  '''
  element_from_edge_model(edge_model="ElectricField", device=device, region=region)
  element_from_edge_model(edge_model="ElectricField", device=device, region=region, derivative="Potential")
  PotentialElementFlux="Permittivity*dot2d(ElectricField_x, ElectricField_y, unitx, unity)"
  CreateElementModel2d(device, region, "PotentialElementFlux", PotentialElementFlux)
  CreateElementModelDerivative2d(device, region, "PotentialElementFlux", PotentialElementFlux, "Potential")

  equation(device=device, region=region, name="PotentialEquation", variable_name="Potential",
      element_model="PotentialElementFlux", variable_update=update_type)

#in the future, worry about workfunction
def CreateOxidePotentialContact(device, contact, element_contact=False, attached_to=None, offset=False):
  if attached_to==None:
    contact_bias_name  = GetContactBiasName(contact)
    set_parameter(device=device, name=contact_bias_name, value=0.0)
  else:
    contact_bias_name  = GetContactBiasName(attached_to)
  contact_model_name = GetContactNodeModelName(contact)
  if offset:
    eq = "Potential - {0} -({1}_PotentialOffset)".format(contact_bias_name, contact)
  else:
    eq = "Potential - {0}".format(contact_bias_name)
  CreateContactNodeModel(device, contact, contact_model_name, eq)
  CreateContactNodeModelDerivative(device, contact, contact_model_name, eq, "Potential")

  if element_contact:
    contact_equation(device=device, contact=contact, name="PotentialEquation", variable_name = "Potential",
                    node_model=contact_model_name, element_charge_model="PotentialElementFlux")
  else:
    contact_equation(device=device, contact=contact, name="PotentialEquation", variable_name = "Potential",
                    node_model=contact_model_name, edge_charge_model= "PotentialEdgeFlux")

def SetContactPotentilOffset(device, contact, PotentialOffset):
  print("Set %s PotentialOffset"%contact, PotentialOffset)
  set_parameter(name = "%s_PotentialOffset"%contact, value=PotentialOffset)

def GetContactPotentilOffset(device, contact):
  return get_parameter(name = "%s_PotentialOffset"%contact)

#######
######Interface Functions:
#######

def CreateSiliconOxideInterface(device, interface, type="continuous"):
  '''
    continuous potential at interface
  '''
  model_name = CreateContinuousInterfaceModel(device, interface, "Potential")
  interface_equation(device=device, interface=interface, name="PotentialEquation", variable_name="Potential", 
                      interface_model=model_name, type="continuous")

#
##TODO: similar model for silicon/silicon interface
## should use quasi-fermi potential
def CreateSiliconSiliconInterface(device, interface):
  '''
    Enforces potential, electron, and hole continuity across the interface
  '''
  CreateSiliconOxideInterface(device, interface)
  ename = CreateContinuousInterfaceModel(device, interface, "Electrons")
  interface_equation(device=device, interface=interface, name=ece_name, variable_name="Electrons", interface_model=ename, type="continuous")
  hname = CreateContinuousInterfaceModel(device, interface, "Holes")
  interface_equation(device=device, interface=interface, name=hce_name, variable_name="Holes", interface_model=hname, type="continuous")







def CheckElementValues(device, region, ElementChecklist=None, Length=5):
    if ElementChecklist != None:
        dit={j:list(get_element_model_values(device=device, region=region, name=j ))[:Length] for j in ElementChecklist}
        print(pd.DataFrame(dit)) 
        for i in ElementChecklist:
            lists=get_element_model_values(device=device, region=region, name=i)
            print("********{0:20}\tmax:{1: .6e}\t min:{2: .6e}".format(i,max(lists),min(lists)))

def CheckEdgeValues(device, region, EdgeChecklist=None, Length=5):
    if EdgeChecklist != None:
      dit={j:list(get_edge_model_values(device=device, region=region, name=j ))[:Length] for j in EdgeChecklist}
      print(pd.DataFrame(dit)) 
      for i in EdgeChecklist:
            lists=get_edge_model_values(device=device, region=region, name=i)
            print("********{0:20}\tmax:{1: .6e}\t min:{2: .6e}".format(i,max(lists),min(lists)))

            

def CreateAbsCurrentAndEfield(device, region, magnitude=False):
  ''' reate absolute value of Efield and current'''
  for i in ("ElectricField", "ElectronCurrent", "HoleCurrent", "ElementElectronCurrent", "ElementHoleCurrent", "SingleCarrierCurrent"):
    if InEdgeModelList(device, region, i) :
      if not InElementModelList(device, region, "%s_x"%(i)):
        element_from_edge_model(device=device, region=region, edge_model=i)
      if not InElementModelList(device, region, "%s_mag"%(i,)) and magnitude:
        EQ="pow({0}_x^2 + {0}_y^2 + 1e-300, 0.5)".format(i,)
        CreateElementModel2d(device, region, "%s_mag"%(i,), EQ)
        
    if InElementModelList(device, region, i) :
      if not InElementModelList(device, region, "%s_x"%(i)):
        vector_element_model(device=device, region=region, element_model=i)
      if not InElementModelList(device, region, "%s_mag"%(i,)) and magnitude:
        EQ="pow({0}_x^2 + {0}_y^2 + 1e-300, 0.5)".format(i,)
        CreateElementModel2d(device, region, "%s_mag"%(i,), EQ)

##########
######Mobility Functions:
##########

def ElectricFieldMagnitude(device, region):
  element_from_edge_model(edge_model="ElectricField", device=device, region=region)
  element_from_edge_model(edge_model="ElectricField", device=device, region=region, derivative="Potential")
  
  ElectricField_mag="pow(ElectricField_x^2 + ElectricField_y^2 + 1e-300, 0.5)"
  CreateElementModel2d(device, region, "ElectricField_mag", ElectricField_mag)
  CreateElementModelDerivative2d(device, region, "ElectricField_mag", ElectricField_mag, "Potential")
  #the , turns this into an actual tuple
  # for j in ("Potential",):
    # for i in ("@en0","@en1", "@en2"):
    #   ex="(ElectricField_x * ElectricField_x:{0}{1} + ElectricField_y * ElectricField_y:{0}{1})/ElectricField_mag".format(j, i)
    #   CreateElementModel2d(device, region, "ElectricField_mag:{0}{1}".format(j, i), ex)


def MicorSaturationMobilityParameters(device, region, V_n_Sat=1.0e11, V_p_Sat=1.0e11, LorentzEField=1.0e0):
 # velocity saturation
  set_parameter(device=device, region=region, name="V_n_Sat",       value=V_n_Sat)#um/s
  set_parameter(device=device, region=region, name="V_p_Sat",       value=V_p_Sat)#um/s
  set_parameter(device=device, region=region, name="LorentzEField", value=LorentzEField)#v/um


def CentiSaturationMobilityParameters(device, region, V_n_Sat=1.0e7, V_p_Sat=1.0e7, LorentzEField=1.0e4):
 # velocity saturation
  set_parameter(device=device, region=region, name="V_n_Sat",       value=V_n_Sat)#cm/s
  set_parameter(device=device, region=region, name="V_p_Sat",       value=V_p_Sat)#cm/s
  set_parameter(device=device, region=region, name="LorentzEField", value=LorentzEField)#v/cm

def ProjectSaturationMobility(device, region):
  for (mu, mu_eff, V_Sat, Current) in (("mu_n", "mu_n_eff", "V_n_Sat", "ElectronCurrent"), 
                                        ("mu_p", "mu_p_eff", "V_p_Sat", "HoleCurrent")):
    EFieldProName="EFieldAt{0}".format(Current)
    EFieldProEqu="abs(ElectricField_x*{0}_x+ElectricField_y*{0}_y)".format(Current)
    #*pow({0}_x^2 + {0}_y^2 + 1e-300, -0.5)
    CreateElementModel2d           (device,  region, EFieldProName, EFieldProEqu)
    CreateElementModelDerivative2d (device,  region, EFieldProName, EFieldProEqu, "Potential", "Electrons", "Holes")
    mu="({0}) / ((1 + ({0} * {1}/{2})^2)^0.5)".format(mu_eff, EFieldProName, V_Sat)
    CreateElementModel2d           (device,  region, mu, mu)
    CreateElementModelDerivative2d (device,  region, mu, mu, "Potential")


def QSSaturationMobility(device, region, ChargeType):
  if not InElementModelList(device, region, "ElectricField_mag"):
    ElectricFieldMagnitude(device, region)
  if ChargeType =="Electrons":
    (i,j,k) = ("mu_n","mu_n_eff","V_n_Sat")
  elif ChargeType =="Holes":
    (i,j,k) = ("mu_p","mu_p_eff","V_p_Sat")
  else:
    raise NameError('Please neme the charge tpye like Elertrons or Holes')
  mu="({0}) / ((1 + ({0} * ElectricField_mag/{1})^2)^0.5)".format(j,k)
  # mu="({0}) / ((1 + ({0} * ElectricField_x/{1})^2)^0.5)".format(j,k)
  CreateElementModel2d           (device,  region, i, mu)
  CreateElementModelDerivative2d (device,  region, i, mu, "Potential")

def LorentzSaturationMobility(device, region, mu_sat, LorentzEField, ChargeType):
  if ChargeType =="Electrons":
    (i,j) = ("mu_n","mu_n_eff")
  elif ChargeType =="Holes":
    (i,j) = ("mu_p","mu_p_eff")
  else:
    raise NameError('Please neme the charge tpye like Elertrons or Holes')
  mu="{mu_eff} * ({LorentzEField})^2 / ({LorentzEField}^2 + abs(ElectricField_x)^2)".format(mu_eff=j, LorentzEField=LorentzEField)
  CreateElementModel2d           (device,  region, i, mu)
  CreateElementModelDerivative2d (device,  region, i, mu, "Potential")


def ElementElectronCurrent2d(device, region):
  name="ElementElectronCurrent"
  if InEdgeModelList(device, region, name):
    delete_edge_model(device=device, region=region, name=name)
  Jn = "ElectronCharge*{0}*EdgeInverseLength*V_t*kahan3(Electrons@en1*Bern01,  Electrons@en1*vdiff,  -Electrons@en0*Bern01)".format("mu_n")
  CreateElementModel2d(device, region, name, Jn)
  vector_element_model(device=device, region=region, element_model=name)
  for i in ("Electrons", "Holes", "Potential"):
    CreateElementModelDerivative2d(device, region, name, Jn, i)
    CreateElementModelDerivative2d(device, region, "%s_x"%name, Jn, i)
    CreateElementModelDerivative2d(device, region, "%s_y"%name, Jn, i)

def ElementHoleCurrent2d(device, region):
  name="ElementHoleCurrent"
  if InEdgeModelList(device, region, name):
    delete_edge_model(device=device, region=region, name=name)
  Jp ="-ElectronCharge*{0}*EdgeInverseLength*V_t*kahan3(Holes@en1*Bern01, -Holes@en0*Bern01, -Holes@en0*vdiff)".format("mu_p")
  CreateElementModel2d(device, region, name, Jp)
  vector_element_model(device=device, region=region, element_model=name)
  for i in ("Electrons", "Holes", "Potential"):
    CreateElementModelDerivative2d(device, region, name, Jp, i)
    CreateElementModelDerivative2d(device, region, "%s_x"%name, Jp, i)
    CreateElementModelDerivative2d(device, region, "%s_y"%name, Jp, i)


def ElementElectronContinuityEquation(device, region):
  '''
    Uses element current model for equation
  '''
  if InNodeModelList(device, region,"ElectronGeneration"):
    equation(device=device, region=region, name=ece_name, variable_name="Electrons",
           time_node_model="NCharge", node_model="ElectronGeneration",
           element_model="ElementElectronCurrent", variable_update="positive")
  else:
    equation(device=device, region=region, name=ece_name, variable_name="Electrons",
           time_node_model="NCharge", element_model="ElementElectronCurrent", variable_update="positive")


def ElementHoleContinuityEquation(device, region):
  '''
    Uses element current model for equation
  '''
  if InNodeModelList(device, region,"HoleGeneration"):
    equation(device=device, region=region, name=hce_name, variable_name="Holes",
           time_node_model="PCharge", node_model="HoleGeneration",
           element_model="ElementHoleCurrent", variable_update="positive")
  else:
    equation(device=device, region=region, name=hce_name, variable_name="Holes",
           time_node_model="PCharge", element_model="ElementHoleCurrent", variable_update="positive")

def CreateElementSiliconDriftDiffusion(device, region):
  ElectricFieldMagnitude(device, region)
  # ProjectSaturationMobility(device, region)
  QSSaturationMobility(device, region)
  ElementElectronCurrent2d(device, region)
  ElementHoleCurrent2d(device, region)
  ElementElectronContinuityEquation(device, region)
  ElementHoleContinuityEquation(device, region)

#TODO: expand for circuit
def ElementContactElectronContinuityEquation(device, contact):
  '''
    Uses element current model for equation
  '''
  contact_electrons_name = "{0}nodeelectrons".format(contact)
  contact_equation(device=device, contact=contact, name=ece_name, variable_name="Electrons",
                   node_model=contact_electrons_name, element_current_model="ElementElectronCurrent")


#TODO: expand for circuit
def ElementContactHoleContinuityEquation(device, contact):
  '''
    Uses element current model for equation
  '''
  contact_holes_name = "{0}nodeholes".format(contact)
  contact_equation(device=device, contact=contact, name=hce_name, variable_name="Holes",
                    node_model=contact_holes_name,  element_current_model="ElementHoleCurrent")

def CreateElementCurrentContact(device, contact):
  ElementContactElectronContinuityEquation(device, contact)
  ElementContactHoleContinuityEquation(device, contact)


def GetSaveFileName(device,FileName):
  if InParameterList(device, "Description",):
    Description=get_parameter(device=device, name="Description")
    if not os.path.exists(Description):
      os.makedirs(Description)
    fname="%s/%s"%(Description,FileName)
    # fname="%s/%s_%s"%(Description,Description,FileName)
  else:
    fname=FileName
  return fname

def QSSaveDevice(device, file, ftype="devsim"):
  fname=GetSaveFileName(device,file)
  write_devices(file=fname, type=ftype)
  if not InDisplayHost() and not InNoneCompressHost():
    tar = tarfile.open("%s%s"%(fname,".tar.gz"),'w:gz')   
    tar.add(fname)
    tar.close()
    os.remove(fname)



def GetContactBiasName(contact):
  return "{0}_bias".format(contact)

def GetContactNodeModelName(contact):
  return "{0}nodemodel".format(contact)

def PrintCurrents(device, contact, ContactList=[]):
  '''
     print out contact currents
  '''
  # TODO add charge

  VoltageList=[]
  for j in ContactList:
    VoltageList="\t".join(i+":"+str(get_parameter(device=device, name=GetContactBiasName(i))) for i in j)
  contact_bias_name = GetContactBiasName(contact)
  voltage         = get_parameter(device=device, name=GetContactBiasName(contact))

  if InParameterList(device, "SingleCarrier"):
    SingleCarrier = get_parameter(device=device, name="SingleCarrier")
  else:
    SingleCarrier = False
  if SingleCarrier:
    SingleCarrier_current    = get_contact_current(device=device, contact=contact, equation="SingleCarrierEquation")
    print("********{0} Current\t{1}\t{2}\t{3}".format(contact, SingleCarrier_current, voltage, VoltageList))
    return [SingleCarrier_current,]
  else:
    electron_current= get_contact_current(device=device, contact=contact, equation=ece_name)
    hole_current    = get_contact_current(device=device, contact=contact, equation=hce_name)
    total_current   = electron_current + hole_current
    print("********{0} Current\t{1}\t{2}\t{3}\t{4}\t{5}".format(contact,  electron_current, hole_current, total_current, voltage, VoltageList))
    return [total_current,hole_current,electron_current]


def SetContactDopingOffset(device, contact):
  region=get_region_list(device=device, contact=contact)[0]
  Ni = get_parameter(device=device, region=region, name="n_i")
  V_t= get_parameter(device=device, region=region, name="V_t")
  net_contact_doping = get_parameter(device=device, region=region, name="%s_doping"%contact)\
                     + get_parameter(device=device, region=region, name="bulk_doping")
  contact_electrons=0.5*net_contact_doping+ 0.5*(net_contact_doping**2+4*Ni**2)**0.5
  contact_holes    =0.5*(net_contact_doping**2+4*Ni**2)**0.5-0.5*net_contact_doping
  # print()
  if net_contact_doping>0 :
    DopingOffset = V_t*math.log(contact_electrons/Ni)
  else:
    DopingOffset = -V_t*math.log(contact_holes/Ni)
  print("Set %s DopingOffset:"%contact, DopingOffset)
  set_parameter(device=device, name = "%s_DopingOffset"%contact, value=DopingOffset)

def GetContactDopingOffset(device, contact):
  return get_parameter(device=device, name = "%s_DopingOffset"%contact)



def SpecificValueInitiateRamp(device, Region, Parameter, abs_error=1e30, rel_error=1e-8, iterations=30,ratio=2):
  print("\n*******%s Initiate"%Parameter)
  OriginValue=get_parameter(device=device, region=Region, name=Parameter)
  portion=1
  while portion<=1 and portion>0.01:
    print("\n********Solve %s of %s"%(portion, Parameter))
    set_parameter(device=device, region=Region, name=Parameter,value=OriginValue*portion)
    try:
      solve(type="dc", absolute_error=abs_error, relative_error=rel_error, maximum_iterations=iterations)
    except devsim_py3.error as msg:
      print(msg)
      portion=portion/ratio
      continue
    except:
      traceback.print_exc()
      raise False
    else:
      portion=portion*ratio
    finally:
      set_parameter(device=device, region=Region, name=Parameter,value=OriginValue)
  return True

def InitialSolve(device, orders=5, abs_error=1e30, rel_error=1e-8, iterations=30):
  print("\n*******InitialSolve:")
  TempRel_error=rel_error
  while True:
    print("\n*******Solve for relative error of %s"%TempRel_error)
    try:
      solve(type="dc", absolute_error=abs_error, relative_error=TempRel_error, maximum_iterations=iterations)
    except devsim_py3.error as msg:
      if str(msg).find("Convergence failure") == 0:
        TempRel_error=TempRel_error*1e3
        continue
      else:
        raise 
    except:
      traceback.print_exc()
      raise 
    else:
      if TempRel_error<=rel_error:
        return True
      elif TempRel_error>10**orders:
        raise NameError('Increase the tolerance of rel_error range')
      TempRel_error=TempRel_error*1e-3
      continue
    

