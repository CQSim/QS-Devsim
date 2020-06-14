######## Single Carriers Functions
from devsim.python_packages.simple_dd import *
from devsim import *
sys.path.append('/home/ccc/devsim/QS_Micro_packages/')
from QSmodel_create import *

def SetMicroMeterSingleSemiParameters(device, region, 
                          T             = 300   , #K
                          rPermittivity = 11.1  ,
                          n_i           = 1.0e-2, # /um^3
                          Mobility      = 1e4   , # um/s^2
                          LifeTime      = 1e-5  , # s
                          ):
  '''
    Sets physical parameters assuming constants
  '''
  #### TODO: make T a free parameter and T dependent parameters as models
  Permittivity = rPermittivity * get_parameter(name="eps0")
  ElectronCharge = get_parameter(name = "ElectronCharge")
  Kb  = get_parameter(name="Kb")
  V_t = Kb * T/ElectronCharge
  WriteSingleSemiParameters(device, region,Permittivity, n_i, T, V_t, Mobility, LifeTime)

def WriteSingleSemiParameters(device, region,Permittivity, n_i, T, V_t, Mobility, LifeTime):
  set_parameter(device=device, name="SingleCarrier",  value=True  )
  set_parameter(device=device, region=region, name="Permittivity",   value=Permittivity  )
  set_parameter(device=device, region=region, name="T",              value=T    )
  set_parameter(device=device, region=region, name="n_i",            value=n_i  )
  # set_parameter(device=device, region=region, name="kT",             value=kT)
  set_parameter(device=device, region=region, name="V_t",            value=V_t)
  set_parameter(device=device, region=region, name="mu",             value=Mobility)
  #default SRH parameters
  set_parameter(device=device, region=region, name="taun",           value=LifeTime)



def CreateSingleCarrierSolution(device, region):
  print("Create Single Carrier Solution for %s"%region)
  '''
    Creates the physical models for a Silicon region
  '''
  CreateSolution(device, region, "SingleCharge")
  CreateNodeModel(device, region, "DopingCharge", "abs(NetDoping)")
  set_node_values(device=device, region=region, name="SingleCharge", init_from="DopingCharge")
  # edge_from_node_model(device=device, region=region, node_model="SingleCharge")


def CreateSingleCarrierPotentialEquation(device, region):
  print("Assemble Single Carrier Potential Equation for %s"%region)
  '''
    Creates the physical models for a Silicon region
  '''
  PotentialNodeCharge = "-ElectronCharge * (SingleCharge + NetDoping)"
  CreateNodeModel(device, region, "PotentialNodeCharge", PotentialNodeCharge)
  CreateNodeModelDerivative(device, region, "PotentialNodeCharge", PotentialNodeCharge, "SingleCharge")

  equation(device=device, region=region, name="PotentialEquation", variable_name="Potential",
           node_model="PotentialNodeCharge", edge_model="PotentialEdgeFlux", variable_update="log_damp")


def CreateSingleCarrierPotentialAtContact(device, region, contact, is_circuit=False, offset=False):
  print("Assemble Single Carrier Potential Equation At %s"%contact)
  '''
    Creates the potential equation at the contact
    if is_circuit is true, than use node given by GetContactBiasName
  '''
  # Means of determining contact charge
  # Same for all contacts
  if not InNodeModelList(device, region, "contactcharge_node"):
    CreateNodeModel(device, region, "contactcharge_node", "ElectronCharge*(SingleCharge + NetDoping)")
  #### TODO: This is the same as D-Field

  set_parameter(device=device, name=GetContactBiasName(contact), value=0.0)

  if offset:
    contact_model = "Potential -{0} + {1}_PotentialOffset".format(GetContactBiasName(contact), contact)
  else:
    contact_model = "Potential -{0}".format(GetContactBiasName(contact))
  
  contact_model_name = GetContactNodeModelName(contact)
  CreateContactNodeModel(device, contact, contact_model_name, contact_model)
  # Simplify it too complicated
  CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_model_name,"Potential"), "1")
  if is_circuit:
    CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_model_name, GetContactBiasName(contact)), "-1")

  if is_circuit:
    contact_equation(device=device, contact=contact, name="PotentialEquation", variable_name="Potential",
                   node_model=contact_model_name, node_charge_model="contactcharge_node", 
                   edge_charge_model="PotentialEdgeFlux", circuit_node=GetContactBiasName(contact))
  else:
    contact_equation(device=device, contact=contact, name="PotentialEquation", variable_name="Potential",
                   node_model=contact_model_name, node_charge_model="contactcharge_node", 
                   edge_charge_model="PotentialEdgeFlux")


def CreateSingleCarrierCurrent(device, region, mu):
  '''
  SingleCharge current
  '''
  #### test for requisite models here
  print("Create Single Carrier Current")
  CreateBernoulli(device, region)
  Js = "-ElectronCharge*{0}*EdgeInverseLength*V_t*kahan3(SingleCharge@n1*Bern01, -SingleCharge@n0*vdiff, -SingleCharge@n0*Bern01)".format(mu)
  CreateEdgeModel(device, region, "SingleCarrierCurrent", Js)
  for i in ("SingleCharge", "Potential"):
    CreateEdgeModelDerivatives(device, region, "SingleCarrierCurrent", Js, i)

def CreateSingleCarrierCurrentEquation(device, region):
  print("Assemble Single Carrier Current Equation")
  NCharge = "-ElectronCharge * SingleCharge"
  CreateNodeModel(device, region, "NCharge", NCharge)
  CreateNodeModelDerivative(device, region, "NCharge", NCharge, "SingleCharge")
  if InNodeModelList(device, region,"SingleChargeGeneration"):
    equation(device=device, region=region, name="SingleCarrierEquation", variable_name="SingleCharge",
      time_node_model = "NCharge", node_model="SingleChargeGeneration",
      edge_model="SingleCarrierCurrent", variable_update="positive")
  else:
    equation(device=device, region=region, name="SingleCarrierEquation", variable_name="SingleCharge",
      time_node_model = "NCharge",
      edge_model="SingleCarrierCurrent", variable_update="positive")


def CreateSingleCarrierDiffusionAtContact(device, region, contact, is_circuit=False): 
  '''
    Restrict electrons and holes to their equilibrium values
    Integrates current into circuit
  '''
  contact_Carriers_model = "SingleCharge + NetDoping"
  contact_Carriers_name = "{0}nodecarriers".format(contact)

  CreateContactNodeModel(device, contact, contact_Carriers_name, contact_Carriers_model)
  CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_Carriers_name, "SingleCharge"), "1")

  #TODO: keyword args
  if is_circuit:
    contact_equation(device=device, contact=contact, name="SingleCarrierEquation", variable_name="SingleCharge",
                         node_model=contact_Carriers_name, edge_current_model="SingleCarrierCurrent", circuit_node=GetContactBiasName(contact))
  else:
    contact_equation(device=device, contact=contact, name="SingleCarrierEquation", variable_name="SingleCharge",
                         node_model=contact_Carriers_name, edge_current_model="SingleCarrierCurrent")
