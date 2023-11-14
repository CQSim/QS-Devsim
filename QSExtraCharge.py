from devsim.python_packages.simple_dd import *
from devsim import *
sys.path.append('/home/ccc/devsim/QS_Micro_packages/')
from QSmodel_create import *
from QSsimple_physics import *

import math

# def SetGaussianDOSTrapParameters(device, region, Type, Width, Center, BulkDensity, Interface=None,Thickness=0.05):
#   # GaussianWidth (V) width of DOS
#   # GaussianCenter (V) the value of the center of Gaussian DOS above intrinsic Fermi level 
#   BulkDensityName="%s_GaussianBulkDensity"%Type
#   set_parameter(device=device, region=region, name = "%s_GaussianWidth"%Type,   value=Width   )
#   set_parameter(device=device, region=region, name = "%s_GaussianCenter"%Type,  value=Center  )
#   set_parameter(device=device, region=region, name = BulkDensityName,  value=BulkDensity  )
#   if Interface!=None:
#     NodeName=Interface_NodeParameter(device=device, region=region, Interface=Interface, ParameterName=BulkDensityName,Thickness=Thickness)
#     set_parameter(device=device, region=region, name = "%s_GaussianBulkDensityNodeSolutionName"%Type,  value=NodeName  )
#   else:
#     NodeName=BulkDensityName

#   return NodeName



def CreateGaussianDOS(device, region, Type, Width, Center, BulkDensity, Interface=None, Thickness=0.05):
  # Degenerate Approximation Model
  # Delta (V) width of DOS
  # GaussianCenter (V) the value of the center of Gaussian DOS below intrinsic Fermi level 
  BulkDensityName=SetGaussianDOSParameters(device, region, Type, Width, Center, BulkDensity, Interface=Interface, Thickness=Thickness)
  if Type=="DonorLike":
    DerivateVariable="Holes"
    eq="{n}*0.5*(1-erf((-V_t*log({d}/n_i)-{t}_GaussianCenter)/{t}_GaussianWidth/pow(2,0.5)))".format(n=BulkDensityName, t=Type, d=DerivateVariable)
  elif Type=="AccepterLike":
    DerivateVariable="Electrons"
    eq="{n}*0.5*(1+erf((V_t*log({d}/n_i)-{t}_GaussianCenter)/{t}_GaussianWidth/pow(2,0.5)))".format(n=BulkDensityName, t=Type, d=DerivateVariable)
  else:
    raise RuntimeError("Please Set the type of _GaussianDOS as: DonorLike or AccepterLike ")
  NameOfNode="%s_GaussianDOS"%Type
  CreateNodeModel(device, region, NameOfNode , eq   )
  CreateNodeModelDerivative(device, region, NameOfNode, eq, DerivateVariable)
  # CreateNodeModel(device, region,  "{n}:{v}".format(n=NameOfNode, v=DerivateVariable) , eq1)
  print("the equation of %s_GaussianDOS is %s "%(Type,eq))
  # print(reslut)
  return NameOfNode

def FixExtraNodeCharge(device, region, Charge_Name):
  if not InNodeModelList(device, region, Charge_Name):
    raise RuntimeError("The Charge Node Model of %s is not exist"%Charge_Name)
  # ChargeDensity=get_node_model_values(device, region, Charge_Name)
  print("Fix the charge of %s"%Charge_Name)
  FixChargeName="%s_Fix"%Charge_Name
  node_solution(device=device, region=region, name=FixChargeName)
  set_node_values(device=device, region=region, name=FixChargeName, init_from=Charge_Name)
  DeleteNodeModelAndDerivates(device, region, Charge_Name)
  return FixChargeName

def InitDetrappingGaussianCharge(device, region, Charge_Name, Type):
  BulkDensityName=get_parameter(device=device, region=region, name = "%s_GaussianBulkDensityNodeSolutionName"%Type)
  if not InNodeModelList(device, region, BulkDensityName):
    raise RuntimeError("The Charge Node Model of %s is not exist"%BulkDensityName)
  set_node_values(device=device, region=region, name=BulkDensityName, init_from=Charge_Name)

def SetGaussianDOSTrapParameters(device, region, TrapName, TotalDensity, Width, Center, Rate0, LocalScale, HDistance,
                                  Interface=None,  Thickness=0.001):
  ### Set the parameter for the Donor_like traps(DLkT)
  # LocalSize: Localize scale of Donor like Trap
  # HDistance: the average distance form Donorlike trap to Hole
  # BulkDensity: the effective density of state of Donorlike trap
  # HoleRate0: the generation and recombination rate of Traps between Hole's sites;
  # TrapName="GaussianDLkT" or "GaussianALkT"


  set_parameter(device=device, region=region, name = "{0}_TDS".format(TrapName),        value=TotalDensity   )
  set_parameter(device=device, region=region, name = "{0}_Width".format(TrapName),      value=Width   )
  set_parameter(device=device, region=region, name = "{0}_Center".format(TrapName),     value=Center   )
  set_parameter(device=device, region=region, name = "{0}_Rate0".format(TrapName),      value=Rate0   )
  set_parameter(device=device, region=region, name = "{0}_LocalScale".format(TrapName), value=LocalScale   )
  set_parameter(device=device, region=region, name = "{0}_HDistance".format(TrapName),  value=HDistance   )
  ### Creat the Trap to Hole rate as basic rate
  # BasicRate="{0}_Rate0 * exp(-2 * {0}_HDistance / {0}_LocalScale)".format(TrapName)
  BasicRate=Rate0 * math.exp(-2 * HDistance /LocalScale)
  set_parameter(device=device, region=region, name="{0}BasicRate".format(TrapName),  value=BasicRate)

  ### Set the total density varing with interface distance
  if Interface==None:
    BulkDensityName="%s_TDS"%TrapName
  else:
    BulkDensityName=Interface_NodeParameter(device=device, region=region, \
                        Interface=Interface, ParameterName="%s_TDS"%TrapName,Thickness=Thickness)
  CreateNodeModel(device, region, "{0}_BDensity".format(TrapName), BulkDensityName   )


def CreateDonorLikeGaussianTrapEquation(device, region, PooleFrenkel=True,LumoRate=False):
  # Degenerate Approximation Model of Gaussian Distribution Trap state

  TrapName="GaussianDLkT"
  DerivateVariable="Holes"

  ### the name for Hole paramaters
  GaussianHomoName="GaussianHomo"
  GaussianLumoName="GaussianLumo"

  ### Solution for Gaussian Distribution Trap state charge density
  SolutionName = "{0}Dst".format(TrapName)
  CreateSolution(device, region, SolutionName)

  ## Create a node model to initate Solution with based on the intrainsic Fermi level
  ### This node name will be replaced due to different purpose
  EffectiveSolutionName="Effective%s"%SolutionName
  EffectiveSolutionEQ ="{t}_BDensity*gfi({t}_Center/V_t,{t}_Width/V_t)".format(t=TrapName)
  CreateNodeModel(device, region, EffectiveSolutionName, EffectiveSolutionEQ)

  ### Set the initial value 
  NodeNumbers=len(get_node_model_values(device=device, region=region,name="NodeVolume"))
  set_node_values(device=device, region=region, name=SolutionName, init_from=EffectiveSolutionName)
  ### 
  CreateSolution(device, region, "%s2"%SolutionName)
  # set_node_values(device=device, region=region, name="%s2"%SolutionName, values=[1e-100]*NodeNumbers)

  # HoleQuasiFermi="-log({d}/n_i)*V_t".format(d=DerivateVariable)
  # HoleQuasiFermi="HoleQuasiFermi"
  # HoleQuasiFermiEQ="{gh}_Center+V_t*igfi(ifelse({d}<{gh}_BDensity,{d},{gh}_BDensity)/{gh}_BDensity, {gh}_Width/V_t)"\
  #                     .format(d=DerivateVariable,gh=GaussianHomoName)
  # CreateNodeModel(device, region, HoleQuasiFermi, HoleQuasiFermiEQ)

  # GaussianHoleDstName="ifelse({d}<{gh}_BDensity,{d},0.999999*{gh}_BDensity)".format(d=DerivateVariable,gh=GaussianHomoName)
  GaussianHoleDstName="min({d},0.999999*{gh}_BDensity)".format(d=DerivateVariable,gh=GaussianHomoName)
  # GaussianHoleDstName="%sDst"%GaussianHomoName
  # GaussianHoleDstEQ="{gh}_BDensity*gfi(({gh}_Center-{hqf})/V_t,{gh}_Width/V_t)"\
  #                               .format(gh=GaussianHomoName, hqf=HoleQuasiFermi)
  # CreateNodeModelAndDerivative(device, region, GaussianHoleDstName, GaussianHoleDstEQ, DerivateVariable)
  if "%sQuasiFermi"%GaussianHomoName not in get_node_model_list(device=device, region=region):
    HoleQuasiFermi=CreateGaussianPositiveQuasiFermiNode(device, region, GaussianHomoName, DerivateVariable)
  else:
    HoleQuasiFermi="%sQuasiFermi"%GaussianHomoName
    
  EffectiveGaussianHomoEName="Effective%sDst"%GaussianHomoName
  EffectiveGaussianHoleDstEQ="{gh}_BDensity*gfi(({hqf}-{gh}_Center-{gh}_Width^2/V_t)/V_t,{gh}_Width/V_t)"\
                                .format(gh=GaussianHomoName, hqf=HoleQuasiFermi)
  CreateNodeModelAndDerivative(device, region, EffectiveGaussianHomoEName, EffectiveGaussianHoleDstEQ, DerivateVariable)
  
  GaussianDlkTQuasiFermi=CreateGaussianPositiveQuasiFermiNode(device, region, TrapName)

  EffectiveSolutionName="Effective%s"%SolutionName
  EffectiveSolutionEQ ="{t}_BDensity*gfi(({t}_Center-{t}_Width^2/V_t-{tqf})/V_t,{t}_Width/V_t)"\
                                .format(t=TrapName,d=DerivateVariable,s=SolutionName,tqf=GaussianDlkTQuasiFermi)
  CreateNodeModelAndDerivative(device, region, EffectiveSolutionName, EffectiveSolutionEQ, DerivateVariable, SolutionName)

  # EffectiveSolutionEQ ="{n}_BDensity*gfi(({n}_Center-{n}_Width^2/V_t-V_t*log({d}/n_i))/V_t,{n}_Width/V_t)"\
  #                               .format(n=TrapName,d=DerivateVariable,s=SolutionName)

  DLTElectons ="%sElectons"%TrapName
  DLTElectonsEQ ="{t}_BDensity*gfi(({tqf}-{t}_Center)/V_t,{t}_Width/V_t)".format(t=TrapName,tqf=GaussianDlkTQuasiFermi)
  CreateNodeModel(device, region, DLTElectons, DLTElectonsEQ)

  if PooleFrenkel==True:
    EffectiveFieldName="EffectiveField"
    EffectiveFieldEQ="ifelse(Holes>Electrons,\
                            2*V_t*NetDoping*ElectronCharge/Permittivity*log(ifelse(abs(NetDoping)*Holes>n_i^2,abs(NetDoping)*Holes/n_i^2,1))^0.5,\
                            2*V_t*NetDoping*ElectronCharge/Permittivity*log(ifelse(abs(NetDoping)>Electrons,abs(NetDoping)/Electrons,1))^0.5)"
    CreateNodeModel(device, region, EffectiveFieldName, EffectiveFieldEQ)

    PooleFrenkelName="PooleFrenkel"
    PooleFrenkelEQ="(ElectronCharge*EffectiveField/3.14/Permittivity)^0.5/V_t"
    CreateNodeModel(device, region, PooleFrenkelName, PooleFrenkelEQ)
  else:
    PooleFrenkelName="0"

  CdetName="{0}Cdet".format(TrapName)
  CdetEQ="exp(({h}_Width^2+{n}_Width^2)/2/V_t^2+({h}_Center-{n}_Center)/V_t)".format(h=GaussianHomoName,n=TrapName)
  CreateNodeModel(device, region, CdetName, CdetEQ)



  HomoGenerationRateName ="%sGRate"%TrapName
  HomoGenerationRateEQ   ="{n}BasicRate* {etd}*{gh}".format(n=TrapName,etd=DLTElectons,gh=GaussianHoleDstName)
  # HomoGenerationRateEQ   ="{n}BasicRate* ({n}_BDensity-{s})*{gh}".format(n=TrapName,gh=GaussianHoleDstName,s=SolutionName)
  # GenerationRateEQ   ="{n}BasicRate* (ifelse({n}_BDensity>{s},{n}_BDensity-{s},0)*{gh}"\
  #                       .format(n=TrapName,gh=GaussianHoleDstName,s=SolutionName)
  CreateNodeModelAndDerivative(device, region, HomoGenerationRateName, HomoGenerationRateEQ, SolutionName, DerivateVariable)
  

  HomoRecombinationName = "%sRRate"%TrapName
  # RecombinationRateEQ="{n}BasicRate * {PF} * {c} *ifelse({gh}_BDensity>{egh},{gh}_BDensity-{egh},0)*{es}"\
  HomoRecombinationEQ="{h}BasicRate * exp({PF}) * {c} *{ehn}*{es}".format(h=GaussianHomoName,PF=PooleFrenkelName,
                                    c=CdetName,ehn=EffectiveGaussianHomoEName,es=EffectiveSolutionName)

  # "ifelse({h}_BDensity>{eh},{h}_BDensity-{eh},{h}_BDensity)"
  # CreateNodeModelAndDerivative(device, region, HomoRecombinationName, RecombinationRateEQ, SolutionName, DerivateVariable)
  CreateNodeModel(device, region, HomoRecombinationName, HomoRecombinationEQ)

  if LumoRate==True:
    LumoRecombinationName = "%sLumoRRate"%TrapName
    LumoRecombinationEQ   ="{n}BasicRate * min({gL}_BDensity,Electrons)*{s}".format(n=TrapName,gL=GaussianLumoName,s=SolutionName)
    CreateNodeModel(device, region, LumoRecombinationName, LumoRecombinationEQ)
    NetGRateEQ="{HG}-{HR}-{LR}".format(HG=HomoGenerationRateName, HR=HomoRecombinationName, LR=LumoRecombinationName)
  else:
    NetGRateEQ="{HG}-{HR}".format(HG=HomoGenerationRateName, HR=HomoRecombinationName)

  NetGRateName="%sNetRate"%TrapName
  CreateNodeModel(device, region, NetGRateName, NetGRateEQ)
  # CreateNodeModelAndDerivative(device, region, NetGRateName, NetGRateEQ, SolutionName, DerivateVariable)

  # MaximumSolutionName="%sMaximum"%TrapName
  # MaximumSolutionEQ ="{t}_BDensity*gfi(({t}_Center-{hqf})/V_t,{t}_Width/V_t)".format(t=TrapName, hqf=HoleQuasiFermi)
  # CreateNodeModel(device, region, MaximumSolutionName, MaximumSolutionEQ),mt=MaximumSolutionName
  
  DensityUpdateNodeName="%sUpdate"%TrapName
  # DensityUpdateNodeEQ="ifelse({s}+{nt}*TimeLag>0,{s}+{nt}*TimeLag,1e-30)".format(nt=NetGRateName,s=SolutionName)
  DensityUpdateNodeEQ="max(min({s}+{ng}*TimeLag, {t}_BDensity*(1-1e-6)), {t}_BDensity*1e-15)"\
                          .format(t=TrapName,ng=NetGRateName,s=SolutionName)
  # DensityUpdateNodeEQ="ifelse({s}+{ng}*TimeLag >0, ifelse({s}+{ng}*TimeLag <{mt},{s}+{ng}*TimeLag,{mt}), {t}_BDensity*1e-15)"\
  CreateNodeModel(device, region, DensityUpdateNodeName, DensityUpdateNodeEQ)
  set_parameter(device=device, name="ExtraNodeChargeDic", value={"region":region, 
                                                              "SolutionName":SolutionName, 
                                                              "NodeName":DensityUpdateNodeName})
  # input(get_parameter(device=device, name="ExtraNodeChargeDic"))

  # GaussianChargeName = "%s_GaussianCharge"%TrapName
  # GaussianChargeEq   = "ElectronCharge * %s"%SolutionName
  # CreateNodeModelAndDerivative(device, region, GaussianChargeName, GaussianChargeEq, SolutionName)

  # equation(device=device, region=region, name="DonorlikeEquation", variable_name=SolutionName,
  #             time_node_model = GaussianChargeName, node_model=NetChargeGRateName,
  #             edge_model="", variable_update="positive")

  # ### Replace the Hole generation rate:
  # CreateNodeModelAndDerivative(device, region, "HoleGeneration", "HoleSRH - %s"%(NetChargeGRateName), "Electrons", "Holes", SolutionName)

def CreateGaussianTrapNode(device, region, TrapName, LumoRate=False, HomoRate=False, PooleFrenkel=True):
  # Gaussian Distribution Trap state
  # TrapName="GaussianDLkT" or "GaussianALkT"

  if TrapName=="GaussianDLkT":
    DerivateVariable = "Holes"
    CarrierStateName = "GaussianHomo"
  elif TrapName=="GaussianALkT":
    DerivateVariable = "Electrons"
    CarrierStateName = "GaussianLumo"
  ### the name for Homo and Lumo
  GaussianHomoName="GaussianHomo"
  GaussianLumoName="GaussianLumo"

  ### the name for Hole paramaters

  ### Solution for Gaussian Distribution Trap state charge density
  SolutionName = "{0}Dst".format(TrapName)
  CreateSolution(device, region, SolutionName)
  ## Create a node model to initate Solution with based on the intrainsic Fermi level
  CreateSolution(device, region, "%s2"%SolutionName)
  set_parameter(device=device, name="ExtraNodeChargeRegion", value=region) 
  if "ExtraNodeChargeList" in get_parameter_list(device=device, region=region):
    nodechargelist=get_parameter(device=device, region=region, name="ExtraNodeChargeList")
    nodechargelist.append(SolutionName)
    print(nodechargelist)

  else:
    nodechargelist=[SolutionName]
  set_parameter(device=device, region=region, name="ExtraNodeChargeList", value=nodechargelist)


  ### 

  CarrierDstName="min({d},0.999999*{gh}_BDensity)".format(d=DerivateVariable,gh=CarrierStateName)

  #### Create Quasi Fermi Level for traps and carriers
  if "%sQuasiFermi"%CarrierStateName not in get_node_model_list(device=device, region=region):
    CarrierQuasiFermi=CreateNodeGaussianQuasiFermi(device, region, CarrierStateName, DerivateVariable)
  else:
    CarrierQuasiFermi="%sQuasiFermi"%CarrierStateName

  GaussianTrapQuasiFermi=CreateNodeGaussianQuasiFermi(device, region, TrapName)

  ### Calculate the initial density through carrier quasi fermi level
  InitialSolutionName="Initial%s"%SolutionName
  if TrapName=="GaussianDLkT":
    InitialSolutionEQ ="{t}_BDensity*gfi(({cqf}-(Potential + {t}_Center))/V_t,{t}_Width/V_t)".format(t=TrapName,cqf=CarrierQuasiFermi)
  elif TrapName=="GaussianALkT":
    InitialSolutionEQ ="{t}_BDensity*gfi((Potential + {t}_Center-{cqf})/V_t,{t}_Width/V_t)".format(t=TrapName,cqf=CarrierQuasiFermi)
  CreateNodeModel(device, region, InitialSolutionName, InitialSolutionEQ)
  ### Set the initial value 
  # NodeNumbers=len(get_node_model_values(device=device, region=region,name="NodeVolume"))
  set_node_values(device=device, region=region, name=SolutionName, init_from=InitialSolutionName)
  DeleteNodeModelAndDerivates(device, region, InitialSolutionName)

  ### Create Trap and Carrier Complementary Density
  if TrapName=="GaussianDLkT":
    HigherState = TrapName
    LowerState  = CarrierStateName
    hqf         = GaussianTrapQuasiFermi
    lqf         = CarrierQuasiFermi

    HigherElectronName ="%sElectons"%TrapName
    HigherElectronEQ ="{t}_BDensity*gfi((Potential + {t}_Center-{tqf})/V_t,{t}_Width/V_t)".format(t=TrapName,tqf=GaussianTrapQuasiFermi)
    CreateNodeModel(device, region, HigherElectronName, HigherElectronEQ)

    LowerHoleName      = CarrierDstName


  elif TrapName=="GaussianALkT":
    HigherState = CarrierStateName
    LowerState  = TrapName
    hqf         = CarrierQuasiFermi
    lqf         = GaussianTrapQuasiFermi

    HigherElectronName = CarrierDstName
    LowerHoleName ="%sHoles"%TrapName
    LowerHoleEQ   ="{t}_BDensity*gfi(({tqf}- Potential - {t}_Center)/V_t,{t}_Width/V_t)".format(t=TrapName,tqf=GaussianTrapQuasiFermi)
    CreateNodeModel(device, region, LowerHoleName, LowerHoleEQ)

  EffectiveHigherHoleName ="Effective%sHoles"%HigherState
  EffectiveHigherHoleEQ   ="{h}_BDensity*gfi(({hqf}-(Potential + {h}_Center + {h}_Width^2/V_t))/V_t,{h}_Width/V_t)".format(h=HigherState,hqf=hqf)
  CreateNodeModel(device, region, EffectiveHigherHoleName, EffectiveHigherHoleEQ)

  EffectiveLowerElectronName ="Effective%sElectons"%LowerState
  EffectiveLowerElectronEQ ="{l}_BDensity*gfi(((Potential + {l}_Center - {l}_Width^2/V_t)-{lqf})/V_t,{l}_Width/V_t)".format(l=LowerState,lqf=lqf)
  CreateNodeModel(device, region, EffectiveLowerElectronName, EffectiveLowerElectronEQ)

  if PooleFrenkel==True:
    EffectiveFieldName="EffectiveField"
    EffectiveFieldEQ="ifelse(Holes>Electrons,\
                            2*V_t*NetDoping*ElectronCharge/Permittivity*log(ifelse(abs(NetDoping)*Holes>n_i^2,abs(NetDoping)*Holes/n_i^2,1))^0.5,\
                            2*V_t*NetDoping*ElectronCharge/Permittivity*log(ifelse(abs(NetDoping)>Electrons,abs(NetDoping)/Electrons,1))^0.5)"
    CreateNodeModel(device, region, EffectiveFieldName, EffectiveFieldEQ)

    PooleFrenkelName="PooleFrenkel"
    PooleFrenkelEQ="(ElectronCharge*EffectiveField/3.14/Permittivity)^0.5/V_t"
    CreateNodeModel(device, region, PooleFrenkelName, PooleFrenkelEQ)
  else:
    PooleFrenkelName="0"

  CdetName="{0}Cdet".format(TrapName)
  CdetEQ="exp(({ls}_Width^2+{hs}_Width^2)/2/V_t^2+({hs}_Center-{ls}_Center)/V_t)".format(hs=HigherState, ls=LowerState)
  CreateNodeModel(device, region, CdetName, CdetEQ)

  TrapGenerationRateName ="%sGRate"%TrapName
  TrapGenerationRateEQ   ="{n}BasicRate* {HE}*{LH}".format(n=TrapName,HE=HigherElectronName,LH=LowerHoleName)
  CreateNodeModel(device, region, TrapGenerationRateName, TrapGenerationRateEQ)
  
  TrapRecombinationName = "%sRRate"%TrapName
  TrapRecombinationEQ="{n}BasicRate * exp({PF}) * {c} * {eLE} * {eHH}"\
               .format(n=CarrierStateName,PF=PooleFrenkelName,c=CdetName, eLE=EffectiveLowerElectronName, eHH=EffectiveHigherHoleName)
  CreateNodeModel(device, region, TrapRecombinationName, TrapRecombinationEQ)


  if TrapName=="GaussianDLkT" and LumoRate==True:
    LumoRecombinationName = "%sLumoRRate"%TrapName
    LumoRecombinationEQ   ="{gL}BasicRate * min({gL}_BDensity,Electrons)*{s}".format(gL=GaussianLumoName,s=SolutionName)
    CreateNodeModel(device, region, LumoRecombinationName, LumoRecombinationEQ)
    NetGRateEQ="{GR}-{RR}-{LR}".format(GR=TrapGenerationRateName, RR=TrapRecombinationName, LR=LumoRecombinationName)
  elif TrapName== "GaussianALkT" and HomoRate==True:
    HomoRecombinationName = "%sHomoRRate"%TrapName
    HomoRecombinationEQ   ="{gH}BasicRate * min({gH}_BDensity,Holes)*{s}".format(gH=GaussianHomoName,s=SolutionName)
    CreateNodeModel(device, region, HomoRecombinationName, HomoRecombinationEQ)
    NetGRateEQ="{GR}-{RR}-{HR}".format(GR=TrapGenerationRateName, RR=TrapRecombinationName, HR=HomoRecombinationName)
  else:
    NetGRateEQ="{GR}-{RR}".format(GR=TrapGenerationRateName, RR=TrapRecombinationName)

  NetGRateName="%sNetRate"%TrapName
  CreateNodeModel(device, region, NetGRateName, NetGRateEQ)
  
  DensityUpdateNodeName="%sUpdate"%SolutionName
  # DensityUpdateNodeEQ="ifelse({s}+{nt}*TimeLag>0,{s}+{nt}*TimeLag,1e-30)".format(nt=NetGRateName,s=SolutionName)
  DensityUpdateNodeEQ="max(min({s} + {nr}*TimeLag, {t}_BDensity*(1-1e-6)), {t}_BDensity*1e-19)"\
                          .format(t=TrapName, nr=NetGRateName, s=SolutionName)
  # DensityUpdateNodeEQ="ifelse({s}+{ng}*TimeLag >0, ifelse({s}+{ng}*TimeLag <{mt},{s}+{ng}*TimeLag,{mt}), {t}_BDensity*1e-15)"\
  CreateNodeModel(device, region, DensityUpdateNodeName, DensityUpdateNodeEQ)
  #                                                             "SolutionName":SolutionName, 
  #                                                             "NodeName":DensityUpdateNodeName})

def CreateDonorLikeGaussianTrapContactEquation(device, region, contact):
  TrapName="GaussianDLkT"
  SolutionName = "{0}Dst".format(TrapName)
  contact_nodeEQ = "{0} + 0".format(SolutionName)
  contact_node_model = "{0}ContactNode".format(TrapName)
  CreateContactNodeModel(device, contact, contact_node_model, contact_nodeEQ)
  # Simplify it too complicated
  CreateContactNodeModel(device, contact, "{0}:{1}".format(contact_node_model, SolutionName), "1")
  contact_equation(device=device, contact=contact, name="DonorlikeEquation", variable_name=SolutionName,
               node_model=contact_node_model)

def CreateGaussianDOSDetrapping(device, region, ChargeModel, ChargeName, DetrappingDistance, Interface):
  FixChargeName=FixExtraCharge(device, region, ChargeName)
  DistanceName="%sDetrapDis"%FixChargeName
  set_parameter(device=device, region=region, name = DistanceName,   value=DetrappingDistance)

  DerivateVariable="Holes"
  DetrapPotentialVaryingEq="{0}*pow(abs(2*abs(bulk_doping)*E_t/Permittivity*log(abs(bulk_doping)*{1}/pow(n_i,2))),0.5)"\
                           .format(DistanceName, DerivateVariable)
  DetrapPotentialVaryingName="DetrapPotentialVarying"
  CreateNodeModel(device, region, DetrapPotentialVaryingName , DetrapPotentialVaryingEq)
  CreateNodeModelDerivative(device, region,DetrapPotentialVaryingName, DetrapPotentialVaryingEq, DerivateVariable)

  eq   ="ifelse({3}<abs(bulk_doping),{0}*0.886226925*(1+erf(({1}+{2}_GaussianCenter)/{2}_GaussianWidth)),{0})"\
        .format(FixChargeName, DetrapPotentialVaryingName, ChargeModel, DerivateVariable)
  # eqdif="ifelse({3}<abs(bulk_doping),{0}*0.886226925*derfdx(({1}+{2}_GaussianCenter)/{2}_GaussianWidth)*{1}:{3},0)".format(FixChargeName, DetrapPotentialVarying, ChargeModel, DerivateVariable)
    # eq1="{0}_GaussianBulkDensity/exp(pow(({0}_GaussianCenter+log({1}/n_i)*V_t)/{0}_GaussianWidth,2))*(V_t/{0}_GaussianWidth*{1})".format(Type, DerivateVariable)
  NameOfNode="%s_Detrapping"%ChargeName
  CreateNodeModel(device, region, NameOfNode , eq   )
  # CreateNodeModel(device, region, "%s:%s"%(NameOfNode,DerivateVariable) , eqdif)
  CreateNodeModelDerivative(device, region, NameOfNode, eq, DerivateVariable)
  return NameOfNode

def  ReplaceExtraPotentialNodeCharge(device, region, *Varables):
  ###Add extra charge to PotentialNondeCharge
  # set_parameter(device=device, name="SweepSpeed", value=SweepSpeed)
  # set_parameter(device=device, name="TimeLag",    value=BiasStep/SweepSpeed)
  
  ENCEQ="NetDoping"
  print(get_parameter(device=device, region=region, name="ExtraNodeChargeList"))
  for i in get_parameter(device=device, region=region, name="ExtraNodeChargeList"):
    if i=="GaussianALkTDst":
      ENCEQ="%s - %s"%(ENCEQ, i)
    elif i=="GaussianDLkTDst":
      ENCEQ="%s + %s"%(ENCEQ, i)
  print("Extra Node Charge in %s in %s "%(region,ENCEQ))
  CreateNodeModel(device, region, "ExtraChargeNode", ENCEQ)

  print("**********Replace the Potential Node Charge")
  DeleteNodeModelAndDerivates(device, region, "PotentialNodeCharge")
  pne = "-ElectronCharge*kahan3(Holes, -Electrons, ExtraChargeNode)"
  CreateNodeModelAndDerivative(device, region, "PotentialNodeCharge", pne, "Potential", "Electrons", "Holes")

  for var in Varables:
    # CreateNodeModelDerivative(device, region, "ExtraCharge", ExtraCharge, var)
    CreateNodeModelDerivative(device, region, "PotentialNodeCharge", pne, var)

def InterfaceCharge_NodeModel(device, region, Interface, Charge_Type, Thickness):
  #get the nodes on Interface
  ChargeName="%s_At_%s"%(Charge_Type,Interface,)
  length="exp(-y/{t})".format(t=Thickness)
  # length="exp(-y/{t})*if(y<={t}*10,1)".format(t=Thickness)
  eq2="{l}*{c}".format(l=length, c=Charge_Type)
  NodeModel=CreateNodeModel(device, region, ChargeName , eq2)
  print(NodeModel)
  for v in ("Potential", "Electrons", "Holes"):
    result=CreateNodeModel(device, region, "{m}:{v}".format(m=ChargeName, v=v),
                         "{l}*{c}:{v}".format(l=length, c=Charge_Type, v=v))
    print(result)
  return ChargeName

def Interface_NodeParameter(device, region, Interface, ParameterName, Thickness):
  #get the nodes on Interface
  TotalDensity=get_parameter(device=device, region=region, name = ParameterName)
  BulkDensity=TotalDensity/Thickness
  NodeName  ="%s_At_%sNode"%(ParameterName,Interface)
  # length="exp(-y/{t})*if(y<={t}*10,1)".format(t=Thickness)
  Eq="{d:1.15e}*exp(-y/{t})".format(d=BulkDensity, t=Thickness)
  NodeModel=CreateNodeModel(device, region, NodeName , Eq)

  NodeSolutionName="%s_At_%s"%(ParameterName,Interface)
  node_solution(device=device, region=region, name=NodeSolutionName)
  set_node_values(device=device, region=region, name=NodeSolutionName, init_from=NodeName)
  DeleteNodeModelAndDerivates(device, region, NodeName)
  return NodeSolutionName

def DeleteNodeModelAndDerivates(device, region, name):
  if InNodeModelList(device, region, name): 
    delete_node_model(device=device, region=region, name=name)
    for i in get_node_model_list(device=device, region=region):
      if i.startswith("%s:"%name):
        delete_node_model(device=device, region=region, name=i)
  else:
    raise ValueError("Could not find the NodeModel named of: %s"%name)
  
def ExpotionalInterfaceCharge_NodeModel(device, Interface, region, Charge_Type, Distance):
  pass
