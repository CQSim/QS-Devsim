import os, socket, tarfile, traceback, math, sys
from devsim.python_packages.simple_dd import *
from devsim import *
import pandas as pd


sys.path.append('/home/ccc/devsim/QS_Micro_packages/')
from QSmodel_create import *
hostname = socket.gethostname()


def InDisplayHost():
  #### Determine whether the simulated results be displayed immediately
  HostList=("nashost","C-Twist","TurboKing","ccc-Inspiron-660s","T640U")
  return hostname in HostList

def InNoneCompressHost():
  #### Determine whether the exported files be compressed in server
  HostList=("pvd","C-Twist","nashost","ub18")#,"T640-Ubuntu","T640U"
  return hostname in HostList

#TODO: make this a class so that paramters can be changed
contactcharge_node="contactcharge_node"
PotentialEdgeFlux="PotentialEdgeFlux"
ece_name="Electron_Equation"
hce_name="Hole_____Equation"
celec_model = "(0.5*abs( NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"
chole_model = "(0.5*abs(-NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"

def SetExtendedPrecision(Threads=4,TaskSize=2048):
  set_parameter(name = "extended_solver",   value=True)
  set_parameter(name = "extended_model",    value=True)
  set_parameter(name = "extended_equation", value=True)
  set_parameter(name = "threads_available", value=Threads )
  set_parameter(name = "threads_task_size", value=TaskSize)

def SetVerbosity(verbose=True):
  if verbose==True:
    set_parameter(name="debug_level", value="verbose")
  else:
    set_parameter(name="debug_level", value="info")

def SetCentiMeterBasicParameters(T=300):
  q     = 1.6021766e-19 # Coul
  Kb    = 1.3806503e-23 # J/K
  eps0  = 8.8541878e-14 # F/cm
  V_t = Kb * T/q
  E_t = Kb * T
  set_parameter(name = "ElectronCharge",    value=q   )
  set_parameter(name = "Kb",                value=Kb  )
  set_parameter(name = "eps0",              value=eps0)
  set_parameter(name =  "E_t",              value=E_t)
  set_parameter(name =  "V_t",              value=V_t)

def SetMicroMeterBasicParameters(T=300):
  q     = 1.6021766e-19 # Coul
  Kb    = 1.3806503e-23 # J/K
  eps0  = 8.8541878e-18 # F/um

  V_t = Kb * T/q
  E_t = Kb * T
  set_parameter(name = "ElectronCharge",    value=q   )
  set_parameter(name = "Kb",                value=Kb  )
  set_parameter(name = "eps0",              value=eps0)
  set_parameter(name =  "E_t",              value=E_t)
  set_parameter(name =  "V_t",              value=V_t)

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
  E_t = Kb * T
  WriteSiliconParameters(device, region, Permittivity, n_i,n_i,n_i, T, E_t, V_t, nMobility, pMobility, nLifeTime, pLifeTime)

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
  E_t = Kb * T
  WriteSiliconParameters(device, region,Permittivity, n_i,n_i,n_i, T, E_t, V_t, nMobility, pMobility, nLifeTime, pLifeTime)


def WriteSiliconParameters(device, region, Permittivity, n_i, n1, p1, T, E_t, V_t, nMobility, pMobility, nLifeTime, pLifeTime):
  set_parameter(device=device, region=region, name="Permittivity",   value=Permittivity  )
  set_parameter(device=device, region=region, name="T",              value=T    )
  set_parameter(device=device, region=region, name="n_i",            value=n_i  )
  set_parameter(device=device, region=region, name="n1",             value=n1)
  set_parameter(device=device, region=region, name="p1",             value=p1)

  set_parameter(device=device, region=region, name="E_t",            value=E_t)
  set_parameter(device=device, region=region, name="V_t",            value=V_t)
  set_parameter(device=device, region=region, name="mu_n_eff",       value=nMobility)
  set_parameter(device=device, region=region, name="mu_p_eff",       value=pMobility)

  #default SRH parameters
  set_parameter(device=device, region=region, name="taun",           value=nLifeTime)
  set_parameter(device=device, region=region, name="taup",           value=pLifeTime)




#### Set the HOMO and LUMO's parameters for Gaussian DOS

def SetGaussianDOSCarrierParameters(device, region, TypeName, BulkDensity, Width, Center, Rate0, LocalScale, HDistance):
  ### Set the parameter for the GaussianDOSCarrier
  ### LocalScale: Localize scale of GaussianDOS
  ### HDistance: the average distance form Donorlike trap to Hole
  ### BulkDensity: the effective density of state of GaussianDOS
  ### HoleRate0: the generation and recombination rate of Traps between Hole's sites;
  ### TypeName="GaussianHomo" or "GaussianLumo"

  if TypeName=="GaussianHomo" or TypeName=="GaussianLumo":
    set_parameter(device=device, region=region, name = "{0}_BDensity".format(TypeName),   value=BulkDensity   )
    set_parameter(device=device, region=region, name = "{0}_Width".format(TypeName),      value=Width   )
    set_parameter(device=device, region=region, name = "{0}_Center".format(TypeName),     value=Center   )
    set_parameter(device=device, region=region, name = "{0}_Rate0".format(TypeName),      value=Rate0   )
    set_parameter(device=device, region=region, name = "{0}_LocalScale".format(TypeName), value=LocalScale   )  #nm
    set_parameter(device=device, region=region, name = "{0}_HDistance".format(TypeName),  value=HDistance   )   #nm
    BasicRate=Rate0 * math.exp(-2 * HDistance /LocalScale)
    set_parameter(device=device, region=region, name = "{0}BasicRate".format(TypeName),   value=BasicRate)
  else:
    NameError("Please sepecify the TypeName as GaussianHomo or GaussianLumo")


  ### Create the quasi fermi level for Gaussian Distribution Charge state
def CreateElementGaussianQuasiFermi(device, region, TypeName, CarrierName):
  # Eq= "0.5*({0}@en0+{0}@en1)".format("%sQuasiFermi"%TypeName)
  SqrtName="Sqrt%s"%CarrierName
  if TypeName in ("GaussianDLkT", "GaussianHomo"):
    Eq="V_t*igfi(min({s},{t}_BDensity)/{t}_BDensity,{t}_Width/V_t)".format(s=SqrtName, t=TypeName)
  elif TypeName in ("GaussianALkT", "GaussianLumo"):
    Eq="- V_t*igfi(min({s},{t}_BDensity)/{t}_BDensity,{t}_Width/V_t)".format(s=SqrtName, t=TypeName)
            
  ElementName="Element%sQuasiFermi"%CarrierName
  CreateElementModel2d(device, region, ElementName, Eq)
  CreateElementModelDerivative2d(device,  region, ElementName, Eq, CarrierName, "Potential")
  return ElementName

def CreateEdgeGaussianQuasiFermi(device, region, TypeName, CarrierName):
  SqrtName="Sqrt%s"%CarrierName
  if TypeName in ("GaussianDLkT", "GaussianHomo"):
    Eq='V_t*igfi(min({s},{t}_BDensity)/{t}_BDensity,{t}_Width/V_t)'.format(s=SqrtName, t=TypeName)
  elif TypeName in ("GaussianALkT", "GaussianLumo"):
    Eq="-V_t*igfi(min({s},{t}_BDensity)/{t}_BDensity,{t}_Width/V_t)".format(s=SqrtName, t=TypeName)

  EdgeName="Edge%sQuasiFermi"%CarrierName
  CreateEdgeModelAndDerivatives(device, region, EdgeName, Eq, CarrierName, "Potential")
  # CreateEdgeModelDerivatives(device,  region, EdgeName, Eq, CarrierName)
  return EdgeName

  ### Create the quasi fermi level for Gaussian Distribution Charge state
def CreateNodeGaussianQuasiFermi(device, region, TypeName, SolutionName=None):

  if SolutionName==None:
    SolutionName = "{0}Dst".format(TypeName)

  Name="%sQuasiFermi"%TypeName
  if TypeName in ("GaussianDLkT", "GaussianHomo"):
    Eq="Potential + {t}_Center + V_t*igfi(min({s},{t}_BDensity)/{t}_BDensity,{t}_Width/V_t)"\
              .format(s=SolutionName,t=TypeName)
  elif TypeName in ("GaussianALkT", "GaussianLumo"):
    Eq="Potential + {t}_Center - V_t*igfi(min({s},{t}_BDensity)/{t}_BDensity,{t}_Width/V_t)"\
              .format(s=SolutionName,t=TypeName)
  CreateNodeModelAndDerivative(device, region, Name, Eq, SolutionName)
  return Name

def CreateGaussianCarriersExchange(device, region, BasicRate, Distance):

  HigherState = "GaussianLumo"
  LowerState  = "GaussianHomo"

  HigherStateQuasiFermi = CreateNodeGaussianQuasiFermi(device, region, HigherState, "Electrons")
  LowerStateQuasiFermi  = CreateNodeGaussianQuasiFermi(device, region, LowerState, "Holes")

  EffectiveHigherHoleName ="Effective%sHoles"%HigherState
  EffectiveHigherHoleEQ   ="{h}_BDensity*gfi(({hqf}-(Potential+{h}_Center+{h}_Width^2/V_t))/V_t,{h}_Width/V_t)".format(h=HigherState,hqf=HigherStateQuasiFermi)
  CreateNodeModelAndDerivative(device, region, EffectiveHigherHoleName, EffectiveHigherHoleEQ, "Electrons", "Holes")

  EffectiveLowerElectronName ="Effective%sElectons"%LowerState
  EffectiveLowerElectronEQ ="{l}_BDensity*gfi(((Potential+{l}_Center-{l}_Width^2/V_t)-{lqf})/V_t,{l}_Width/V_t)".format(l=LowerState,lqf=LowerStateQuasiFermi)
  CreateNodeModelAndDerivative(device, region, EffectiveLowerElectronName, EffectiveLowerElectronEQ, "Electrons", "Holes")

  CdetName="GaussianCarriersExchangeCdet"
  CdetEQ="exp(({ls}_Width^2+{hs}_Width^2)/2/V_t^2+({hs}_Center-{ls}_Center)/V_t)".format(ls=LowerState,hs=HigherState)
  CreateNodeModel(device, region, CdetName, CdetEQ)

  RecombinationName = "GaussianCarrierRRate"
  RecombinationEQ   ="{r0} * exp(-{ds}/{hs}_LocalScale) * Electrons * Holes "\
                      .format(r0=BasicRate, ds=Distance, hs=HigherState)
  CreateNodeModelAndDerivative(device, region, RecombinationName, RecombinationEQ, "Electrons", "Holes")

  GenerationRateName ="GaussianCarrierGRate"
  GenerationRateEQ   ="{r0} * exp(-{ds}/{ls}_LocalScale) * {c} * {eLE} * {eHH}"\
               .format(r0=BasicRate, ds=Distance, ls=LowerState, c=CdetName, 
                      eLE=EffectiveLowerElectronName, eHH=EffectiveHigherHoleName)
  CreateNodeModelAndDerivative(device, region, GenerationRateName, GenerationRateEQ, "Electrons", "Holes")
  
  Gn = "ElectronCharge * (GaussianCarrierGRate - GaussianCarrierRRate)"

  CreateNodeModelAndDerivative(device, region, "GaussianCarrierGR", Gn, "Electrons", "Holes")
  CreateNodeModelAndDerivative(device, region, "ElectronGeneration", "GaussianCarrierGR", "Electrons", "Holes")
  CreateNodeModelAndDerivative(device, region, "HoleGeneration",     "-GaussianCarrierGR", "Electrons", "Holes")

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

def CreateSiliconNetDoping(device, region, bulk_doping):
  set_parameter(device=device, region= region, name="bulk_doping", value=bulk_doping )
  CreateNodeModel(device, region, "NetDoping", "bulk_doping")
  # CreateNodeModel(device, region, "NetDoping", "%1.15e +1e-100" % bulk_doping)
  # set_parameter(device=device, region= region, name="NetDoping", value=bulk_doping )


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

def DeleteIntrinsicNodeModelDerivatives(device, region):
    for i in ("PotentialIntrinsicCharge", 
              "IntrinsicCharge", 
              "IntrinsicElectrons",
              "IntrinsicHoles",  ):
      delete_node_model(device=device, region=region, name=i)
      delete_node_model(device=device, region=region, name="{0}:Potential".format(i))

def CreateSiliconPotentialOnlyContact(device, contact, is_circuit=False, offset=False):
  '''
    Creates the potential equation at the contact
    if is_circuit is true, than use node given by GetContactBiasName
  '''
  # Means of determining contact charge
  # Same for all contacts
  # if not InNodeModelList(device, region, "contactcharge_node"):
  #   CreateNodeModel(device, region, "contactcharge_node", "ElectronCharge*IntrinsicCharge")
  #### TODO: This is the same as D-Field
  BiasName=GetContactBiasName(contact)
  set_parameter(device=device, name=BiasName, value=0.0)

  if offset:
    contact_model = "Potential -{0} - {1}_PotentialOffset".format(BiasName, contact)
    # contact_model = "Potential -{0} + ifelse(NetDoping > 0, -V_t*log({1}/n_i), \
    # V_t*log({2}/n_i))+{3}_PotentialOffset".format(BiasName, celec_model, chole_model, contact)
  else:
    contact_model = "Potential -{0} + ifelse(NetDoping > 0, -V_t*log({1}/n_i), \
                      V_t*log({2}/n_i))".format(BiasName, celec_model, chole_model)
  
  contact_model_name = GetContactNodeModelName(contact)
  CreateContactNodeModel(device, contact, contact_model_name, contact_model)
  # Simplify it too complicated
  CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_model_name,"Potential"), "1")
  if is_circuit:
    CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_model_name,BiasName), "-1")
    contact_equation(device=device, contact=contact, name="PotentialEquation", node_model=contact_model_name,
                   edge_charge_model="PotentialEdgeFlux", circuit_node=BiasName)
  else:
    contact_equation(device=device, contact=contact, name="PotentialEquation", node_model=contact_model_name,
                   edge_charge_model="PotentialEdgeFlux",)

def CreateSRH(device, region):
  USRH="(Electrons*Holes - n_i^2)/(taup*(Electrons + n1) + taun*(Holes + p1))"
  Gn = "-ElectronCharge * USRH"
  Gp = "+ElectronCharge * USRH"

  CreateNodeModelAndDerivative(device, region, "USRH", USRH, "Electrons", "Holes")
  CreateNodeModelAndDerivative(device, region, "ElectronSRH", Gn, "Electrons", "Holes")
  CreateNodeModelAndDerivative(device, region, "HoleSRH", Gp, "Electrons", "Holes")

  CreateNodeModelAndDerivative(device, region, "ElectronGeneration", "ElectronSRH", "Electrons", "Holes")
  CreateNodeModelAndDerivative(device, region, "HoleGeneration",     "HoleSRH",     "Electrons", "Holes")


def CreateECE(device, region, mu_n, Bandtype="Silicon"):
  if Bandtype=="Silicon":
    CreateElectronCurrent(device, region, mu_n)
  NCharge = "-ElectronCharge * Electrons"
  CreateNodeModelAndDerivative(device, region, "NCharge", NCharge, "Electrons")
  if InNodeModelList(device, region,"ElectronGeneration"):
    equation(device=device, region=region, name=ece_name, variable_name="Electrons",
      time_node_model = "NCharge", node_model="ElectronGeneration",
      edge_model="ElectronCurrent", variable_update="positive")
  else:
    equation(device=device, region=region, name=ece_name, variable_name="Electrons",
      time_node_model = "NCharge",
      edge_model="ElectronCurrent", variable_update="positive")

def CreateHCE(device, region, mu_p, Bandtype="Silicon"):
  if  Bandtype=="Silicon":
    CreateHoleCurrent(device, region, mu_p)  
  PCharge = "ElectronCharge * Holes"
  CreateNodeModelAndDerivative(device, region, "PCharge", PCharge, "Holes")

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
  CreateNodeModelDerivative(device, region, "PotentialNodeCharge", pne, "Electrons", "Holes")

  equation(device=device, region=region, name="PotentialEquation", variable_name="Potential",
            node_model="PotentialNodeCharge", edge_model="PotentialEdgeFlux",
            time_node_model="", variable_update="log_damp")



def CreateElectronAndHoleSolutions(device, region, PreValue=True, Bandtype="Silicon"):
  for c in ("Electrons","Holes"):
    CreateSolution(device, region, c)
    if Bandtype=="Silicon":
      if c =="Electrons":
        expression="Potential-log(ifelse({0}<=0,1e-100,{0})/n_i)*V_t".format(c)
      if c =="Holes":
        expression="Potential+log(ifelse({0}<=0,1e-100,{0})/n_i)*V_t".format(c)
      CreateNodeModel(device, region, "QuasiFermiPotential_%s"%c, expression)
    if PreValue:
      set_node_values(device=device, region=region, name=c, init_from="Intrinsic%s"%c)

def SetInitialSolutionFromDopping(device, region, Bandtype="Silicon"):
  # Bandtype="Silicon" or "Gaussian"
  if Bandtype == "Silicon":
    celec_model = "(0.5*abs( NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"
    chole_model = "(0.5*abs(-NetDoping+(NetDoping^2 + 4 * n_i^2)^(0.5)))"
    for (name, equation, solution) in (
         ("InitialElectrons", celec_model  ,  "Electrons"  ),
         ("InitialHoles",     chole_model  ,  "Holes"  ),
         ):
      CreateNodeModel(device, region, name, equation)
      set_node_values(device=device, region=region, name=solution, init_from=name)
      delete_node_model(device=device, region=region, name=name)
  elif Bandtype == "Gaussian":
    for name in ("Electrons",  "Holes"):
      set_node_values(device=device, region=region, name=solution, init_from=name)


def CreateSiliconDriftDiffusion(device, region, mu_n="mu_n_eff", mu_p="mu_p_eff",SRHMolde=True, Bandtype="Silicon"):
  CreatePE(device, region)
  if Bandtype=="Silicon":
    CreateBernoulli(device, region)
  if SRHMolde:
    CreateSRH(device, region)
  CreateECE(device, region, mu_n, Bandtype)
  CreateHCE(device, region, mu_p, Bandtype)

def CreateSiliconDriftDiffusionAtContact(device, contact, is_circuit=False): 
  '''
    Restrict electrons and holes to their equilibrium values
    Integrates current into circuit
  '''
  contact_electrons_name = "{0}nodeelectrons".format(contact)
  contact_electrons_model = "Electrons - ifelse(NetDoping > 0, {0}, n_i^2/{1})".format(celec_model, chole_model)

  contact_holes_name = "{0}nodeholes".format(contact)
  contact_holes_model = "Holes - ifelse(NetDoping < 0, {1}, n_i^2/{0})".format(celec_model, chole_model)


  CreateContactNodeModel(device, contact, contact_electrons_name, contact_electrons_model)
  CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_electrons_name, "Electrons"), "1")

  CreateContactNodeModel(device, contact, contact_holes_name, contact_holes_model)
  CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_holes_name, "Holes"), "1")
  
  #TODO: The simplification of the ifelse statement is time consuming
  if is_circuit:
    contact_equation(device=device, contact=contact, name=ece_name, node_model=contact_electrons_name,
                         edge_current_model="ElectronCurrent", circuit_node=GetContactBiasName(contact))
    contact_equation(device=device, contact=contact, name=hce_name, node_model=contact_holes_name,
                         edge_current_model="HoleCurrent", circuit_node=GetContactBiasName(contact))
  else:
    contact_equation(device=device, contact=contact, name=ece_name, node_model=contact_electrons_name,
                         edge_current_model="ElectronCurrent")
    contact_equation(device=device, contact=contact, name=hce_name, node_model=contact_holes_name,
                         edge_current_model="HoleCurrent")




def CreateElectronAndHoleAtContact(device, contact, is_circuit=False): 
  region=get_region_list(device=device, contact=contact)[0]
  Ni = get_parameter(device=device, region=region, name="n_i")
  V_t= get_parameter(device=device, region=region, name="V_t")
  
  if InParameterList(device, "%s_doping"%contact, region=region):
    net_contact_doping = get_parameter(device=device, region=region, name="%s_doping"%contact)\
                       + get_parameter(device=device, region=region, name="bulk_doping")
  else:
    net_contact_doping = get_parameter(device=device, region=region, name="bulk_doping")
  contact_electrons= 0.5*net_contact_doping + 0.5*(net_contact_doping**2+4*Ni**2)**0.5
  contact_holes    =-0.5*net_contact_doping + 0.5*(net_contact_doping**2+4*Ni**2)**0.5
  print("Electrons and Holes at %s is %s and %s"%(contact,contact_electrons,contact_holes))

  contact_electrons_name = "{0}nodeelectrons".format(contact)
  contact_electrons_model = "Electrons - {0}".format(contact_electrons)

  contact_holes_name = "{0}nodeholes".format(contact)
  contact_holes_model = "Holes - {0}".format(contact_holes)
  CreateContactNodeModel(device, contact, contact_electrons_name, contact_electrons_model)
  CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_electrons_name, "Electrons"), "1")

  CreateContactNodeModel(device, contact, contact_holes_name, contact_holes_model)
  CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_holes_name, "Holes"), "1")
  #TODO: The simplification of the ifelse statement is time consuming
  if is_circuit:
    contact_equation(device=device, contact=contact, name=ece_name, node_model=contact_electrons_name,
                         edge_current_model="ElectronCurrent", circuit_node=GetContactBiasName(contact))
    contact_equation(device=device, contact=contact, name=hce_name, node_model=contact_holes_name,
                         edge_current_model="HoleCurrent", circuit_node=GetContactBiasName(contact))
  else:
    contact_equation(device=device, contact=contact, name=ece_name, node_model=contact_electrons_name,
                         edge_current_model="ElectronCurrent")
    contact_equation(device=device, contact=contact, name=hce_name, node_model=contact_holes_name,
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

##### in the future, worry about workfunction
def CreateOxidePotentialContact(device, contact, element_contact=False, attached_to=None, offset_contact=None, is_circuit=False):
  ##### The varable "offset_contact" is used for the situation while two side of one contact are connected with two different region
  ##### The actual value of potential offset is set by the function of "SetContactPotentilOffset"
  ##### The varable "attached_to" is used for the situation  while a long contact covers several regions. 
  if attached_to==None:
    contact_bias_name  = GetContactBiasName(contact)
    set_parameter(device=device, name=contact_bias_name, value=0.0)
  else:
    contact_bias_name  = GetContactBiasName(attached_to)
  contact_model_name = GetContactNodeModelName(contact)
  if offset_contact:
    eq = "Potential - {0} -({1}_PotentialOffset)".format(contact_bias_name, offset_contact)
  else:
    eq = "Potential - {0}".format(contact_bias_name)
  CreateContactNodeModel(device, contact, contact_model_name, eq)
  CreateContactNodeModelDerivative(device, contact, contact_model_name, eq, "Potential")

  if is_circuit:
    CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_model_name, contact_bias_name), "-1")
    if element_contact:
      # contact_edge_model(device=device, contact=contact, name="PotentialElementFlux:%s"%contact_bias_name, equation="0", display_type='vector')
      contact_equation(device=device, contact=contact, name="PotentialEquation", 
                      node_model=contact_model_name, element_charge_model="PotentialElementFlux",circuit_node=contact_bias_name)
    else:
      contact_equation(device=device, contact=contact, name="PotentialEquation", 
                      node_model=contact_model_name, edge_charge_model= "PotentialEdgeFlux", circuit_node=contact_bias_name)
  else:
    if element_contact:
      contact_equation(device=device, contact=contact, name="PotentialEquation", 
                      node_model=contact_model_name, element_charge_model="PotentialElementFlux")
    else:
      contact_equation(device=device, contact=contact, name="PotentialEquation", 
                      node_model=contact_model_name, edge_charge_model= "PotentialEdgeFlux")


#######
######Interface Functions:
#######

def CreateSiliconOxideInterface(device, interface, InterfaceType="continuous"):
  '''
    continuous potential at interface
  '''
  model_name = CreateContinuousInterfaceModel(device, interface, "Potential")
  interface_equation(device=device, interface=interface, name="PotentialEquation", 
                      interface_model=model_name, type=InterfaceType)

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



def CheckNodeValues(device, region, NodeChecklist=None, Length=5):
  if NodeChecklist != None:
    dit={j:list(get_node_model_values(device=device, region=region, name=j ))[:Length] for j in NodeChecklist}
    print(pd.DataFrame(dit)) 
    for i in NodeChecklist:
      lists=get_node_model_values(device=device, region=region, name=i)
      print("********{0:20}\tmax:{1: .6e}\t min:{2: .6e}".format(i,max(lists),min(lists)))




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

def CreateAbsEfield(device, region):
  if get_dimension(device=device)==1:
    CreateEdgeModel(device, region, "%s_mag"%("ElectricField"), "abs(ElectricField)")
    return
  if not InElementModelList(device, region, "%s_x"%("ElectricField")):
    element_from_edge_model(device=device, region=region, edge_model="ElectricField")
  if not InElementModelList(device, region, "%s_mag"%("ElectricField")):
    EQ="pow({0}_x^2 + {0}_y^2 + 1e-300, 0.5)".format("ElectricField")
    CreateElementModel2d(device, region, "%s_mag"%("ElectricField"), EQ)

def CreateAbsCurrentAndEfield(device, region, magnitude=False):
  ''' create absolute value of Efield and current'''
  if get_dimension(device=device)==1:
    return
  for i in ("ElectricField", "ElectronCurrent", "HoleCurrent", "ElementElectronCurrent", "ElementHoleCurrent", "SingleCarrierCurrent"):
    if InEdgeModelList(device, region, i) :
      if not InElementModelList(device, region, "%s_x"%(i)):
        element_from_edge_model(device=device, region=region, edge_model=i)
      if not InElementModelList(device, region, "%s_mag"%(i)) and magnitude:
        EQ="pow({0}_x^2 + {0}_y^2 + 1e-300, 0.5)".format(i)
        CreateElementModel2d(device, region, "%s_mag"%(i), EQ)
        
    if InElementModelList(device, region, i) :
      if not InElementModelList(device, region, "%s_x"%(i)):
        vector_element_model(device=device, region=region, element_model=i)
      if not InElementModelList(device, region, "%s_mag"%(i)) and magnitude:
        EQ="pow({0}_x^2 + {0}_y^2 + 1e-300, 0.5)".format(i)
        CreateElementModel2d(device, region, "%s_mag"%(i), EQ)

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


def MicorSaturationMobilityParameters(device, region, V_n_Sat=1.0e11, V_p_Sat=1.0e11, LorentzEField=1.0e-2):
 # velocity saturation
  set_parameter(device=device, region=region, name="V_n_Sat",       value=V_n_Sat)#um/s
  set_parameter(device=device, region=region, name="V_p_Sat",       value=V_p_Sat)#um/s
  set_parameter(device=device, region=region, name="LorentzEField", value=LorentzEField)#v/um


def CentiSaturationMobilityParameters(device, region, V_n_Sat=1.0e7, V_p_Sat=1.0e7, LorentzEField=1.0e4):
 # velocity saturation
  set_parameter(device=device, region=region, name="V_n_Sat",       value=V_n_Sat)#cm/s
  set_parameter(device=device, region=region, name="V_p_Sat",       value=V_p_Sat)#cm/s
  set_parameter(device=device, region=region, name="LorentzEField", value=LorentzEField)#v/cm

###### make this on a per carrier basis function
###### assume that low field model already exists, but not projected
def CreateNormalElectricFieldFromCurrentFlow(device, region, low_curr):
    dimension=get_dimension(device=device)
    if dimension != 2:
        raise ValueError("Supported in 2d only")
    if low_curr=="ElectronCurrent":
      Carrisers="Electrons"
    elif low_curr=="HoleCurrent":
      Carrisers="Holes"
    else:
      raise NameError('Please neme the current name as ElectronCurrent or HoleCurrent')

    element_from_edge_model(edge_model=low_curr, device=device, region=region)
    element_from_edge_model(edge_model=low_curr, device=device, region=region, derivative="Potential")
    element_from_edge_model(edge_model=low_curr, device=device, region=region, derivative="Electrons")
    element_from_edge_model(edge_model=low_curr, device=device, region=region, derivative="Holes")

    #### Get the current magnitude on the element
    #### do we need the derivative, since this is only scaling the direction?
    #### Worry about small denominator 1e-300 is about the limit 
    J_lf_mag="pow({0}_x^2 + {0}_y^2 + 1e-300, 0.5)".format(low_curr)
    CreateElementModel2d(device, region, "{0}_mag".format(low_curr), "{0}".format(J_lf_mag))
    for j in ("Electrons", "Holes", "Potential"):
        for i in ("@en0", "@en1", "@en2"):
            ex="({0}_x * {0}_x:{1}{2} + {0}_y * {0}_y:{1}{2})/{0}_mag".format(low_curr, j, i)
            CreateElementModel2d(device, region, "{0}_mag:{1}{2}".format(low_curr, j, i), ex)

    ### This calculates the normalized current in each direction
    for i in  ("x", "y"):
        J_norm="{0}_{1} / {0}_mag".format(low_curr, i)
        CreateElementModel2d(device, region, "{0}_norm_{1}".format(low_curr, i), J_norm)
        CreateElementModelDerivative2d(device, region, "{0}_norm_{1}".format(low_curr, i), J_norm, "Electrons", "Holes", "Potential")

    #### Get the parallel e field to current flow
    Eparallel_J="{0}_norm_x * ElectricField_x + {0}_norm_y * ElectricField_y".format(low_curr)
    CreateElementModel2d          (device, region, "Eparallel_{0}".format(low_curr), Eparallel_J)
    CreateElementModelDerivative2d(device, region, "Eparallel_{0}".format(low_curr), Eparallel_J, "Potential", "Electrons", "Holes")

    # magnitude e field
    # ElectricField_mag="pow(ElectricField_x^2 + ElectricField_y^2 + 1e-300, 0.5)"
    # CreateElementModel2d(device, region, "ElectricField_mag", ElectricField_mag)
    # #the , turns this into an actual tuple
    # for j in ("Potential",):
    #     for i in ("@en0","@en1", "@en2"):
    #         ex="(ElectricField_x * ElectricField_x:{0}{1} + ElectricField_y * ElectricField_y:{0}{1})/ElectricField_mag".format(j, i)
    #         CreateElementModel2d(device, region, "ElectricField_mag:{0}{1}".format(j, i), ex)

    # #### Get the normal e field to current flow
    # Enormal_J="pow(max(ElectricField_mag^2 - Eparallel_{0}^2,1e-300), 0.5)".format(low_curr)
    # CreateElementModel2d(device, region, "Enormal_{0}".format(low_curr), Enormal_J)
    # #CreateElementModelDerivative2d $device $region Enormal_{low_curr} {Enormal_J} Electrons Holes Potential
    # for j in ("Electrons", "Holes", "Potential"):
    #     for i in ("@en0", "@en1", "@en2"):
    #         ex="(ElectricField_mag * ElectricField_mag:{0}{1} - Eparallel_{2} * Eparallel_{2}:{0}{1})/Enormal_{2}".format(j, i, low_curr)
    #         CreateElementModel2d(device, region, "Enormal_{0}:{1}{2}".format(low_curr, j, i), ex)


def ProjectSaturationMobility(device, region, ChargeType):
  if not InElementModelList(device, region, "ElectricField_mag"):
    ElectricFieldMagnitude(device, region)
  if ChargeType =="Electrons":
    (mu, mu_eff, V_Sat, Current) = ("mu_n", "mu_n_eff", "V_n_Sat", "ElectronCurrent")
  elif ChargeType =="Holes":
    (mu, mu_eff, V_Sat, Current) = ("mu_p", "mu_p_eff", "V_p_Sat", "HoleCurrent")
  else:
    raise NameError('Please neme the charge tpye as Elertrons or Holes')
  # EFieldProName="EFieldAt{0}".format(Current)
  # EFieldProEqu="abs(ElectricField_x*{0}_x+ElectricField_y*{0}_y)/pow({0}_x^2 + {0}_y^2 + 1e-300, -0.5)".format(Current)
  # #*pow({0}_x^2 + {0}_y^2 + 1e-300, -0.5)
  # CreateElementModel2d           (device,  region, EFieldProName, EFieldProEqu)
  # CreateElementModelDerivative2d (device,  region, EFieldProName, EFieldProEqu, "Potential", "Electrons", "Holes")
  muEQ="({0}) / ((1 + ({0} * {1}/{2})^2)^0.5)".format(mu_eff, "ElectricField_x", V_Sat)
  CreateElementModel2d           (device,  region, mu, muEQ)
  CreateElementModelDerivative2d (device,  region, mu, muEQ, "Potential")

######
###### 
def QSSaturationMobility(device, region, ChargeType, EfieldModel="ElectricField_x"):
  if get_dimension() == 1:
    pass
  else:
    if not InElementModelList(device, region, "ElectricField_mag"):
      ElectricFieldMagnitude(device, region)
    if EfieldModel=="ParallelCurrent":
      CreateNormalElectricFieldFromCurrentFlow(device, region, "%sCurrent"%ChargeType)
    elif EfieldModel=="ElectricField_x" or EfieldModel=="ElectricField_mag":
      Efield=EfieldModel

    if ChargeType =="Electron":
      (name,mu_eff,V_Sat) = ("mu_n","mu_n_eff","V_n_Sat")
      if EfieldModel=="ParallelCurrent":
        Efield="Eparallel_ElectronCurrent"
    elif ChargeType =="Hole":
      (name,mu_eff,V_Sat) = ("mu_p","mu_p_eff","V_p_Sat")
      if EfieldModel=="ParallelCurrent":
        Efield="Eparallel_HoleCurrent"
    else:
      raise NameError('Please neme the charge tpye as Elertron or Hole')
    # mu="({0}) / ((1 + ({0} * ElectricField_mag/{1})^2)^0.5)".format(j,k)
    muEQ="({0}) / ((1 + ({0} * {1}/{2})^2)^0.5)".format(mu_eff, Efield, V_Sat)
    CreateElementModel2d           (device,  region, name, muEQ)
    CreateElementModelDerivative2d (device,  region, name, muEQ, "Potential", ChargeType)


def LorentzSaturationMobility(device, region, mu_sat, LorentzEField, ChargeType):
  if ChargeType =="Electrons":
    (i,j) = ("mu_n","V_n_Sat")
  elif ChargeType =="Holes":
    (i,j) = ("mu_p","V_p_Sat")
  else:
    raise NameError('Please name the charge tpye like Elertrons or Holes')
  mu="{V_Sat} * ({LorentzEField})^2 / ({LorentzEField}^2 + abs(ElectricField_x)^2)".format(V_Sat=j, LorentzEField=LorentzEField)
  CreateElementModel2d           (device,  region, i, mu)
  CreateElementModelDerivative2d (device,  region, i, mu, "Potential")


def PooleFrenkelMobilityParameters(device, region, ChargeType, ActivePotential=0, Beta=0.1, Gama=0):
  set_parameter(device=device, region=region, name="%s_muActivePotential"%ChargeType, value=ActivePotential)#v/cm
  set_parameter(device=device, region=region, name="%s_muBeta"%ChargeType, value=Beta)
  set_parameter(device=device, region=region, name="%s_muGama"%ChargeType, value=Gama)


def PooleFrenkelMobility(device, region, ChargeType):
  CreateAbsEfield(device, region)
  if ChargeType =="Electrons":
    (i,j) = ("mu_n","mu_n_eff")
  elif ChargeType =="Holes":
    (i,j) = ("mu_p","mu_p_eff")
  else:
    raise NameError('Please name the charge tpye as Elertrons or Holes')
  mu="({0})*exp(-{1}_muActivePotential/V_t)*exp(({1}_muBeta/V_t-{1}_muGama)*ElectricField_mag^0.5)".format(j, ChargeType)
  CreateElementModel2d          (device,  region, i, mu)
  CreateElementModelDerivative2d(device,  region, i, mu, "Potential")


def MacCDTMobility(device, region, CarrierName, MacCDTPre=2e-9, MacCo=1, PooleFrenkel=0.26,TrapName=None, MiniValue=1e2, TailState=1e-8):
  MuPreName="%s_MuMacCDTPre"%CarrierName
  set_parameter(device=device, region=region, name=MuPreName, value=MacCDTPre)
  set_parameter(device=device, region=region, name="%sMacCo"%CarrierName, value=MacCo)
  CreateAbsEfield(device, region)
  if CarrierName=="Holes":
    Name="mu_p"
    TypeName="GaussianHomo"
  if CarrierName=="Electrons":
    Name="mu_n"
    TypeName="GaussianLumo"

  if get_dimension(device=device)==1:
    EdgeCarrierName=CreateEdgeFromSqrtNode(device, region, CarrierName, CarrierName)
    QuasiFermi=CreateEdgeGaussianQuasiFermi(device, region, TypeName, CarrierName)
    Eq = "{pre}/(ElectronCharge*({eh}))*exp({qf}/V_t)*exp({pf}*ElectricField_mag^0.5)"\
                  .format(pre=MuPreName,eh=EdgeCarrierName,pf=PooleFrenkel,qf=QuasiFermi)
    CreateEdgeModelAndDerivatives(device, region, Name, Eq, CarrierName)
  else:
    ElementHolesName=CreateElementEdgeFromSqrtNode(device, region, CarrierName, CarrierName)
    QuasiFermi=CreateElementGaussianQuasiFermi(device, region, TypeName, CarrierName)
    if TrapName==None:
      if CarrierName=="Holes":
        if PooleFrenkel==0:
          Eq = "{pre}/(ElectronCharge*({eh}))*exp({qf}/V_t*HolesMacCo)"\
                      .format(pre=MuPreName,eh=ElementHolesName,qf=QuasiFermi)
        else:
          Eq = "{pre}/(ElectronCharge*({eh}))*exp({qf}/V_t*HolesMacCo)*exp({pf}*ElectricField_mag^0.5)"\
                      .format(pre=MuPreName,eh=ElementHolesName,pf=PooleFrenkel,qf=QuasiFermi)
      elif CarrierName=="Electrons":
        Eq = "{pre}/(ElectronCharge*({eh}))*exp(-{qf}/V_t*ElectronsMacCo)*exp({pf}*ElectricField_mag^0.5)"\
            .format(pre=MuPreName,eh=ElementHolesName,pf=PooleFrenkel,qf=QuasiFermi)

    else:
      ElementTrapName=CreateElementEdgeFromSqrtNode(device, region, TrapName)
      Eq = "{pre}*{eh}/({eh}+{tp})/(ElectronCharge*({eh}+{ts}))*exp({qf}/V_t)+{mv}"\
              .format(pre=MuPreName,eh=ElementHolesName,qf=QuasiFermi,ts=TailState,tp=ElementTrapName,mv=MiniValue)
    CreateElementModel2d(device, region, Name, Eq)
    CreateElementModelDerivative2d(device, region, Name, Eq, CarrierName)

def ShurHackMobility(device, region, Carrier, TrapName, LRatio=0.01,Exponent=1):
  ### the mobility is proportional with carrier density ratio in xoreaponding charges
  ### LRatio lowest ratio of mobility for low carrier density
  if Carrier =="Electrons":
    (i,mu0) = ("mu_n","mu_n_eff")
  elif Carrier =="Holes":
    (i,mu0) = ("mu_p","mu_p_eff")
  else:
    raise NameError('Please name the charge tpye as Elertron or Hole')

  for t in (Carrier, TrapName):
    edge_from_node_model(device=device, region=region, node_model=t)
    ElementNm="Element%s"%t
    ElementEq="({0}@n0 * {0}@n1)^0.5".format(t)
    CreateElementModel2d(device, region, ElementNm, ElementEq)
    CreateElementModelDerivative2d(device, region, ElementNm, ElementEq, t)

  mu="{u}* ({LR}+(1-{LR})*(Element{c}/(Element{c}+Element{t}))^{ep})*exp(-{c}_muActivePotential/V_t)*exp(({c}_muBeta/T-{c}_muGama)*ElectricField_mag^0.5)"\
  .format(u=mu0, c=Carrier, t=TrapName,LR=LRatio,ep=Exponent)
  CreateElementModel2d          (device,  region, i, mu)
  CreateElementModelDerivative2d(device,  region, i, mu, "Potential")



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

def ElementElectronContinuityEquation(device, region):
  '''
    Uses element current model for equation
  '''
  if not InNodeModelList(device, region,"NCharge"):
    CreateNodeModelAndDerivative(device, region, "NCharge", "-ElectronCharge * Electrons", "Electrons")
  if InNodeModelList(device, region,"ElectronGeneration"):
    equation(device=device, region=region, name=ece_name, variable_name="Electrons",
           time_node_model="NCharge", node_model="ElectronGeneration",
           element_model="ElementElectronCurrent", variable_update="positive")
  else:
    equation(device=device, region=region, name=ece_name, variable_name="Electrons",
           time_node_model="NCharge", element_model="ElementElectronCurrent", variable_update="positive")

def ElementHoleCurrent2d(device, region, OrganicCurrent=False):
  name="ElementHoleCurrent"
  if InEdgeModelList(device, region, name):
    delete_edge_model(device=device, region=region, name=name)
  if OrganicCurrent:
     # "{0}_Width"    "{0}_Center"    "{0}_BDensity"
    Jp ="ElectronCharge*ElectricField*{0}*(Holes@en1*Holes@en0)**0.5 + \
    (Holes@en0-Holes@en1)*EdgeInverseLength*{0}*(Holes@en1*Holes@en0)**0.5/({1}_BDensity*dgfidx())".format("mu_p","GaussianHomo")
  else:
    Jp ="-ElectronCharge*{0}*EdgeInverseLength*V_t*kahan3(Holes@en1*Bern01, -Holes@en0*Bern01, -Holes@en0*vdiff)".format("mu_p")
  CreateElementModel2d(device, region, name, Jp)
  vector_element_model(device=device, region=region, element_model=name)
  for i in ("Electrons", "Holes", "Potential"):
    CreateElementModelDerivative2d(device, region, name, Jp, i)
    CreateElementModelDerivative2d(device, region, "%s_x"%name, Jp, i)
    CreateElementModelDerivative2d(device, region, "%s_y"%name, Jp, i)

def ElementHoleContinuityEquation(device, region):
  '''
    Uses element current model for equation
  '''
  if not InNodeModelList(device, region,"PCharge"):
    CreateNodeModelAndDerivative(device, region, "PCharge", "ElectronCharge * Holes", "Holes")

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
  # QSSaturationMobility(device, region)
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
  contact_equation(device=device, contact=contact, name=ece_name, 
                   node_model=contact_electrons_name, element_current_model="ElementElectronCurrent")


#TODO: expand for circuit
def ElementContactHoleContinuityEquation(device, contact):
  '''
    Uses element current model for equation
  '''
  contact_holes_name = "{0}nodeholes".format(contact)
  contact_equation(device=device, contact=contact, name=hce_name, 
                    node_model=contact_holes_name,  element_current_model="ElementHoleCurrent")

def CreateElementCurrentContact(device, contact):
  ElementContactElectronContinuityEquation(device, contact)
  ElementContactHoleContinuityEquation(device, contact)


###AC Model
def CreateCapacitorCircuit(contact, attached_to=None, value=0.0, acreal=0.0, resistance=1e0):
  if attached_to:
      circuit_element(name="R_%s"%contact, n1=GetContactBiasName(attached_to), n2=1, value=resistance)
      circuit_element(name="V_%s"%contact, n1=1,         n2=0, value=0.0, acreal=acreal)
  else:
      circuit_element(name="V_%s"%contact, n1=GetContactBiasName(contact),         n2=0, value=1.0, acreal=acreal)

  # circuit_element(name="R1", n1=GetContactBiasName(contact), n2=1, value=resistance)

def CreateContactCircuit(contact, attached_to=None, value=0.0, acreal=0.0, resistance=1e0):
  if attached_to:
      circuit_element(name="R_%s"%contact, n1=GetContactBiasName(attached_to), n2=1, value=resistance)
      circuit_element(name="V_%s"%contact, n1=1,         n2=0, value=0.0, acreal=acreal)
  else:
      circuit_element(name="V_%s"%contact, n1=GetContactBiasName(contact),         n2=0, value=1.0, acreal=acreal)

  # circuit_element(name="R1", n1=GetContactBiasName(contact), n2=1, value=resistance)

def CreateTransientCircuit(contact, resistance=1e0):
  circuit_element(name="V1", n1=1,         n2=0, value=0.0)
  circuit_element(name="R1", n1=GetContactBiasName(contact), n2=1, value=resistance)


def PrintCapacitor(device, contact, frequency=[1.0,]):
  Data=[]
  if frequency==None:
    Data.append(get_contact_charge(device=device, contact=contact, equation="PotentialEquation")) 
  else:
    for Freq in frequency:
      solve(type="ac",frequency=Freq)
      Bias     =get_circuit_node_value(node=GetContactBiasName(contact), solution="ssac_real")
      ssac_real=get_circuit_node_value(node="V1.I", solution="ssac_real")
      ssac_imag=get_circuit_node_value(node="V1.I", solution="ssac_imag")
      ssac_amp=(ssac_real**2+ssac_imag**2)**0.5
      Capacitor=-ssac_amp**2/(2*math.pi*Freq*Bias*ssac_imag)
      # Capacitor=-get_circuit_node_value(node="V1.I", solution="ssac_imag")/(2*math.pi*Freq*Bias)
      print("The capacitance: {0} at frequency: {1} is: {2}".format(contact, Freq, Capacitor))
      Data.append(Capacitor)
  print(Data)
  return Data


###File Export Model
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

def ExportParameters(device, file="DeviceParameters"):
  fname=GetSaveFileName(device,file)
  with open("%s.py"%fname, "w", encoding="utf-8") as ofh:
    ofh.write('import devsim\n')
    for p in get_parameter_list():
      v=repr(get_parameter(name=p))
      print(p,v)
      ofh.write('devsim.set_parameter(name="%s", value=%s)\n' % (p, v))
    for i in get_device_list():
      for p in get_parameter_list(device=i):
        v=repr(get_parameter(device=i, name=p))
        print(p,v)
        ofh.write('devsim.set_parameter(device="%s", name="%s", value=%s)\n' % (i, p, v))
    for i in get_device_list():
      for j in get_region_list(device=i):
        for p in get_parameter_list(device=i, region=j):
          v=repr(get_parameter(device=i, region=j, name=p))
          print(p,v)
          ofh.write('devsim.set_parameter(device="%s", region="%s", name="%s", value=%s)\n' % (i, j, p, v))
    ofh.close()


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
    print("********{0} Current\t{1}\t{2}\t{3}\t{4}\t{5}".format(contact, total_current, electron_current, hole_current, voltage, VoltageList))
    return [total_current, hole_current, electron_current]


def SetContactDopingOffset(device, contact):
  region=get_region_list(device=device, contact=contact)[0]
  Ni = get_parameter(device=device, region=region, name="n_i")
  V_t= get_parameter(device=device, region=region, name="V_t")
  
  if InParameterList(device, "%s_doping"%contact, region=region):
    if InParameterList(device, "%s_doping"%contact, region=region):
      net_contact_doping = get_parameter(device=device, region=region, name="%s_doping"%contact)\
                         + get_parameter(device=device, region=region, name="bulk_doping")
    else:
      net_contact_doping = get_parameter(device=device, region=region, name="%s_doping"%contact)
  else:
    net_contact_doping = get_parameter(device=device, region=region, name="bulk_doping")
  contact_electrons= 0.5*net_contact_doping + 0.5*(net_contact_doping**2+4*Ni**2)**0.5
  contact_holes    =-0.5*net_contact_doping + 0.5*(net_contact_doping**2+4*Ni**2)**0.5
  # print()
  if net_contact_doping>0 :
    DopingOffset = V_t*math.log(contact_electrons/Ni)
  else:
    DopingOffset = -V_t*math.log(contact_holes/Ni)
  print("Set %s DopingOffset:"%contact, DopingOffset)
  set_parameter(device=device, name = "%s_DopingOffset"%contact, value=DopingOffset)


def GetContactDopingOffset(device, contact):
  return get_parameter(device=device, name = "%s_DopingOffset"%contact)

def SetContactPotentilOffset(device, contact, PotentialOffset):
  #### The value of "PotentialOffset" can be set with dopping concentration or specific number
  print("Set %s_PotentialOffset"%contact, PotentialOffset)
  set_parameter(device=device, name = "%s_PotentialOffset"%contact, value=PotentialOffset)

def GetContactPotentilOffset(device, contact):
  return get_parameter(device=device, name = "%s_PotentialOffset"%contact)

def SpecificValueInitiateRamp(device, Region, Parameter, abs_error=1e30, rel_error=1e-8, iterations=30, ratio=2):
  print("\n*******%s Initiate"%Parameter)
  OriginValue=get_parameter(device=device, region=Region, name=Parameter)
  portion=1
  while portion<=1 and portion>0.001:
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
  if portion<1:
    raise False
  return True

def SpecificValueApproach(device, Region, Parameter, SetToValue, abs_error=1e30, rel_error=1e-8, iterations=30, step_lag=0.1, min_setp=0.001, method="Line"):
  print("\n*******%s Initiate"%Parameter)
  OriginValue=get_parameter(device=device, region=Region, name=Parameter)
  PresentValue=PreviousValue=SetToValue
  step_size=step_lag
  while PresentValue!=OriginValue:
    print("\n********Solve %s of %s with the lag of %s"%(PresentValue, Parameter, step_size))
    set_parameter(device=device, region=Region, name=Parameter,value=PresentValue)
    try:
      solve(type="dc", absolute_error=abs_error, relative_error=rel_error, maximum_iterations=iterations)
    except devsim_py3.error as msg:
      print(msg)
      step_size=step_size/2
      PresentValue = PreviousValue + math.copysign(step_size, OriginValue-PresentValue)
      continue
    except:
      traceback.print_exc()
      raise False
    else:
      PreviousValue = PresentValue
      if method=="Line":
        PresentValue = PreviousValue + math.copysign(step_size, OriginValue-PresentValue)
      elif method=="Log":
        if abs(PresentValue)<abs(SetToValue):
          PresentValue = PreviousValue*(1+abs(step_size))
        else:
          PresentValue = SetToValue
    finally:
      set_parameter(device=device, region=Region, name=Parameter,value=OriginValue)
      if step_size < min_setp:
        break
  print(">>>>>>>>>>>>>>Finish the SpecificValueApproach")
  return True

def InitialSolve(device, orders=5, type="dc", abs_error=1e30, rel_error=1e-8, iterations=30, tdelta=0.0):
  print("\n*******InitialSolve:")
  TempRel_error=rel_error
  while True:
    print("\n*******Solve for relative error of %s"%TempRel_error)
    try:
      solve(type=type, absolute_error=abs_error, relative_error=TempRel_error, maximum_iterations=iterations, tdelta=tdelta)
    except devsim_py3.error as msg:
      print(msg)
      if str(msg).find("Convergence failure") == 0:
        if TempRel_error>10**orders:
          raise NameError('Increase the tolerance of rel_error range')
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
    
def NodeModelVolumeIntegraton(device, NodeDic):
  region    =NodeDic["region"]
  NodeName  =NodeDic["NodeName"]
  SampleType=NodeDic["SampleType"]
  # input([region,NodeName,SampleType])
  ### SimpleType VolSum Max Min
  volume=get_node_model_values(device=device, region=region, name="NodeVolume" )
  nodemodel=get_node_model_values(device=device, region=region, name=NodeName )
  print('****%s for %s '%(SampleType,NodeName),max(nodemodel),min(nodemodel))
  if SampleType=="VolSum":
    value=sum([x*y for x,y in zip(volume,nodemodel)])
  elif SampleType=="Max":
    value=max(nodemodel)
  elif SampleType=="Min":
    value=min(nodemodel)
  elif isinstance(SampleType,int):
    value=nodemodel[SampleType]


  # sumnodes=sum(volume)
  # input(sumnodes)
  return value


