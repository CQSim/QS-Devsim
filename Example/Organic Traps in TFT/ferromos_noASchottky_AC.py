import devsim
import sys, numpy, os,  tarfile

sys.path.append('../../')
from QSsimple_physics import *
from QSPlotSweep import *
from QSFerro import *
from QSramp import *
from QSExtraCharge import *
from QSSchottkyContact import *

#set extended 128bit for device
SetExtendedPrecision()
# SetVerbosity(verbose=True)
SetMicroMeterBasicParameters()
import create_mos

device          = "mos_topContact"
oxide_region    = "oxide"
silicon_region  = "bulk"
air_region      = "air"
silicon_regions = (silicon_region,)#None
regions         = (silicon_region, oxide_region,  air_region)
interfaces      = ("bulk_air","bulk_oxide")
FerroRegion     = None #oxide_region
Gaussian_region = silicon_region
DLT_region      = "bulk"# "bulk"
ALT_region      = None# "bulk"
bulk_doping     = 1e-2
SweepSpeed      = 6.479 #V/s
#resistance      = 2e5
Description     = "Spd%0.2f"%(SweepSpeed)

for i in regions:
  CreatePotentialAndFlux(device, i)

mu=4.81e5
# LRatio=1e-3
# Exponent=0.4
# Description="%sMu%0.0e"%(Description,mu)
for i in silicon_regions:
  CreateSiliconNetDoping(device, i, bulk_doping)
  SetMicroMeterSiliconParameters(device, i, 300, rPermittivity=3.5  , pMobility=mu, nMobility=mu)
  CreateSiliconPotentialOnly(device, i)

Insulators={"air":1.0, "oxide":4}
for i in Insulators:
  print(i,Insulators[i])
  SetOxideParameters(device, i, rPermittivity=Insulators[i])
  CreateOxidePotentialOnly(device, i, "log_damp")

# 
### Set up contacts
for c in ("ground","gate"):
  CreateOxidePotentialContact(device, c)

contacts=("source", "drain" )
for i in contacts:
  SetContactDopingOffset(device, i)
  CreateSiliconPotentialOnlyContact(device, i)

attached_contacts=(("attached_source","source"),  ("attached_drain", "drain"))
for (i,j) in attached_contacts:
  print("Dopping offset for %s is: "%i, GetContactDopingOffset(device, j))
  SetContactPotentilOffset(device, i, GetContactDopingOffset(device, j))
  CreateOxidePotentialContact(device, i, attached_to=j, offset_contact=i)

for i in interfaces:
  CreateSiliconOxideInterface(device, i)
#, type="hybrid"

InitialSolve(device,rel_error=1e-15)
for i in silicon_regions:
  CreateElectronAndHoleSolutions(device, i)
  CreateSiliconDriftDiffusion(device, i)
  DeleteIntrinsicNodeModelDerivatives(device, i)
  
# QSSaveDevice(device, file="SInitialSolve.tec",ftype="tecplot")

for c in ("source", "drain"):
  tmp = get_region_list(device=device, contact=c)
  r = tmp[0]
  CreateSchottkyDriftDiffusionAtContact(device, c)

for (ac,c) in attached_contacts:
  CreateSiliconPotentialOnlyContact(device, c, offset=True)
  CreateOxidePotentialContact(device, ac, attached_to=c, offset_contact=c)

SchottkyContactOffset=-0.67
Description="%sPoff%s"%(Description,SchottkyContactOffset)

lag=math.copysign(0.2, SchottkyContactOffset-GetContactDopingOffset(device, "source"))
for setp in numpy.arange(GetContactDopingOffset(device, "source"),  SchottkyContactOffset, lag):
  print("Change the Schottky contact potential offset to :",setp)
  for (i,j) in attached_contacts:
    SetContactPotentilOffset(device, i, float(setp))
    SetContactPotentilOffset(device, j, float(setp))
    tmp = get_region_list(device=device, contact=j)
    r = tmp[0]
    SetSchottkyContactParamaters(device, r, j)
  InitialSolve(device,rel_error=1e-15)
  PrintCurrents(device,j)
  # QSSaveDevice(device, file="InitialSolve%s.tec"%setp,ftype="tecplot")
InitialSolve(device,rel_error=1e-15)
# input("aa")


if True:
  HomoTotalDensity= 1.2e9 # um^-3
  HomoWidth       = 0.11  # eV
  HomoCenter      = 0.8   # eV
  HomoRate0       = 2.31347e9  # s^-1
  HomoLocalScale  = 0.29  # nm
  HomoDistance    = 5.526  # nm
  SetGaussianDOSCarrierParameters(device, Gaussian_region, TypeName="GaussianHomo", 
                                      BulkDensity=HomoTotalDensity, Width=HomoWidth, Center=HomoCenter, 
                                      Rate0=HomoRate0, LocalScale=HomoLocalScale, HDistance=HomoDistance)
  Description="%sHn%0.1ew%0.2fc%0.2fr%0.0e"%(Description, HomoTotalDensity, HomoWidth, HomoCenter, HomoRate0)

  LumoTotalDensity= HomoTotalDensity
  LumoWidth       = 0.11  # eV
  LumoCenter      = -0.8  # eV
  LumoRate0       = 1e8   # s^-1
  LumoLocalScale  = 0.29  # nm
  LumoDistance    = 3.163  # nm
  SetGaussianDOSCarrierParameters(device, Gaussian_region, "GaussianLumo", LumoTotalDensity, LumoWidth, LumoCenter, LumoRate0, LumoLocalScale, LumoDistance)
  # QuasiFermi=CreateNodeGaussianQuasiFermi(device, Gaussian_region, "GaussianHomo", "Holes")
if False:
  CreateGaussianCarriersExchange(device, Gaussian_region, 1e1)
  Description="GasRG%s"%(Description)

GaussianMobility=True
for r in ("bulk",):
  DeleteEdgeModelAndDerivatives(device, r, "HoleCurrent")
  # DeleteEdgeModelAndDerivatives(device, r, "ElectronCurrent")
  if GaussianMobility:
    MacCDTPre=3e1  #2.5e-5 S/cm =2.5e-5 S
    PooleFrenkel=0
    MacCo=1
    # TailState=1e-6,Ts%0.0e TailState
    Description="%sMacCDT%0.0eco%spf%s"%(Description,MacCDTPre,MacCo,PooleFrenkel)
    MacCDTMobility(device, r, "Holes", MacCDTPre=MacCDTPre, MacCo=MacCo, PooleFrenkel=PooleFrenkel)
  else:
    ActivePotential=0.1
    Beta=0.01
    Gama=0.1
    Description="%sPFu%2.0eA%sB%sG%s"%(Description,mu,ActivePotential,Beta,Gama)
    PooleFrenkelMobilityParameters(device, r, "Holes", ActivePotential=ActivePotential, Beta=Beta, Gama=Gama)
    PooleFrenkelMobility(device, r, "Holes")

  # ElementElectronCurrent2d(device, r)
  # ElementElectronContinuityEquation(device, r)
  ElementHoleCurrent2d(device, r)
  ElementHoleContinuityEquation(device, r)

for contact in ( "drain", "source"):
  CreateSchottkyElementHolesContactEquation(device, contact) 
  # CreateSchottkyElementElectronsContactEquation(device, contact)

InitialSolve(device, rel_error=1e-15)
PrintCurrents(device, "drain")
# input('bb')
# SetVerbosity(verbose=True)
# Gaussian_region=[]
# VolumeIntegrateList=[]
VolumeIntegrateList=[{"region":"bulk", "NodeName":"Holes",           "SampleType":134},
                     {"region":"bulk", "NodeName":"Electrons",       "SampleType":134},]


NodeChecklist=[]

if DLT_region:
#   # input('bb')
  # DLTChargeModel="DonorLike"
  DLT="GaussianDLkT"
  DLTTotalDensity = 1.6e5
  DLTWidth        = 0.03
  DLTCenter       = 0.545
  DLTRate0        = 2.31347e9
  HDistance       = 5.526 
  DLTLocalScale   = 0.29
  DLTThickness    = 5e-3
  LumoRate        = False
  PooleFrenkel    = True
  TapInterface    ="bulk_oxide"
  Charge_Name = "{0}Dst".format(DLT)
  NodeChecklist=NodeChecklist+[Charge_Name,"GaussianDLkTNetRate","GaussianDLkTGRate","GaussianDLkTRRate",'GaussianDLkTElectons',
  'Holes','EffectiveGaussianHomoElectons','EffectiveGaussianDLkTHoles']##,,TDSName

  Description="%sDltN%0.0ew%0.3fc%0.3fr%0.0el%sPf%iLr%iThc%0.0e"%(Description, DLTTotalDensity, DLTWidth, DLTCenter, DLTRate0, DLTLocalScale, PooleFrenkel, LumoRate, DLTThickness)
  SetGaussianDOSTrapParameters(device, DLT_region, TrapName=DLT, TotalDensity=DLTTotalDensity, Width=DLTWidth, Center=DLTCenter, 
                                Rate0=DLTRate0, LocalScale=DLTLocalScale, HDistance=HDistance, 
                                Interface=TapInterface, Thickness=DLTThickness)
  CreateGaussianTrapNode(device, DLT_region, TrapName=DLT, LumoRate=LumoRate,  PooleFrenkel=PooleFrenkel)
  ReplaceExtraPotentialNodeCharge(device, DLT_region)
  # MiniValue=1e4
  # Description="%sMinimu%0.0e"%(Description,MiniValue)
  # MacCDTMobility(device, DLT_region, "Holes", MacCDTPre=MacCDTPre, TailState=TailState, TrapName="%sDst"%DLT, MiniValue=MiniValue)
  VolumeIntegrateList=VolumeIntegrateList+[
                       {"region":"bulk", "NodeName":Charge_Name,       "SampleType":134},
                       {"region":"bulk", "NodeName":"%sGRate"%DLT,     "SampleType":134},
                       {"region":"bulk", "NodeName":"%sRRate"%DLT,     "SampleType":134},
                       {"region":"bulk", "NodeName":"%sNetRate"%DLT,   "SampleType":134},
                       {"region":"bulk", "NodeName":"EffectiveGaussianHomoElectons",     "SampleType":134},
                       {"region":"bulk", "NodeName":"EffectiveGaussianDLkTHoles",   "SampleType":134},
                       ]

if ALT_region:
  ALT="GaussianALkT"
  ALTTotalDensity = 3e4
  ALTWidth        = 0.13
  ALTCenter       = -0.4
  ALTThickness    = 3e-3
  HomoRate        = True
  DLTLocalScale   = 0.00029
  HDistance       = 0.005
  PooleFrenkel    = True
  TapInterface    ="bulk_oxide"

  SetGaussianDOSTrapParameters(device, ALT_region, TrapName=ALT, Width=ALTWidth, Center=ALTCenter, TotalDensity=ALTTotalDensity, 
                                LocalScale=DLTLocalScale, HDistance=HDistance, Rate0=Rate0, 
                                Interface=TapInterface, Thickness=ALTThickness)
  CreateGaussianTrapNode(device, ALT_region, TrapName=ALT, HomoRate=HomoRate, PooleFrenkel=PooleFrenkel)
  ReplaceExtraPotentialNodeCharge(device, ALT_region)
  Description="%sAltN%0.0ew%sc%sr%0.0ePf%iHr%iThc%0.0e"%(Description, ALTTotalDensity, ALTWidth, ALTCenter, Rate0, PooleFrenkel, HomoRate, ALTThickness)

  InitialSolve(device, rel_error=1e-15)

  TDSName ="{0}_BDensity".format(ALT)
  Charge_Name = "{0}Dst".format(ALT)
  NodeChecklist=NodeChecklist+[Charge_Name,"GaussianALkTNetRate","GaussianALkTGRate","GaussianALkTRRate",'GaussianALkTHoles',
  'Holes','EffectiveGaussianALkTElectons','EffectiveGaussianLumoHoles']##,,TDSName
  # input(NodeChecklist)
  CheckNodeValues(device, r, NodeChecklist=NodeChecklist, Length=15)
  VolumeIntegrateList=VolumeIntegrateList+[
                       {"region":"bulk", "NodeName":Charge_Name,       "SampleType":134},
                       {"region":"bulk", "NodeName":"%sGRate"%ALT,     "SampleType":134},
                       {"region":"bulk", "NodeName":"%sRRate"%ALT,     "SampleType":134},
                       {"region":"bulk", "NodeName":"%sNetRate"%ALT,   "SampleType":134},
                       ]
set_parameter(device=device, name="SweepSpeed", value=SweepSpeed)
if DLT_region!=None or ALT_region!=None:
  for i in range(-12,2):
    Step_size=SweepSpeed*10**(i/2)
    UpdateExtraNodeCharge(device, Step_size, NodeChecklist)
    InitialSolve(device, rel_error=1e-15)

if FerroRegion:
  SaturationPolarization =7.749e-14
  RemanentPolarization   =7.130e-14
  Description="%sPs%0.0e"%(Description,SaturationPolarization)
  print("*******Creat Ferro and Solve")
  MicroMeterFerroParameters(device, FerroRegion, SaturationPolarization=SaturationPolarization, RemanentPolarization=RemanentPolarization) 
  CreateFerroRegion(device, FerroRegion,  update_type="log_damp")
  CreateOxidePotentialContact(device,"gate", element_contact=True)
  SpecificValueInitiateRamp(device, FerroRegion, "SaturationPolarization", abs_error=1e30, rel_error=1e-15, iterations=30)

for r in regions:
  CreateAbsCurrentAndEfield(device, r, magnitude=True)

ElementChecklist=("SweepDirectionY","CoerciveSignY","OldPCoefficientY","NewPCoefficientY", "NumeratorTanhY","DenominatorTanhY","PolarizationY","PreElectricField_y","ElectricField_y")


if FerroRegion!=None:
  Coersive_bias  = 15.00
  Other_bias     = -5 
  # Star_bias      = max(Other_bias,0)+3*Coersive_bias
  # End_bias       = min(Other_bias,0)-3*Coersive_bias 
  Star_bias      = min(Other_bias,0)-3*Coersive_bias
  End_bias       = max(Other_bias,0)+3*Coersive_bias 
  # End_bias       = -abs(Star_bias)/Star_bias * Coersive_bias + Other_bias/2
  Description="Fer%s_D%sG%s~%s"%(Description, Other_bias, Star_bias, End_bias)
  SweepSettings=( 
                  ("drain", Other_bias, "pre",  0.5, False, False,True),
                  ("gate",  Star_bias,  "star", 0.5, False, True, True ),
                  ("gate",  End_bias,   "forw", 0.2, False, True, True ),
                  ("gate",  Star_bias,  "back", 0.2, True,  True, True ),
                  # ("drain", -Other_bias, "output", 0.1, True, True),
                )
else:
  Other_bias = -0.12 #V -0.024 0.048
  For_bias = -100
  Back_bias = -For_bias/2
  Description="%s_D%sG%s~%s"%(Description, Other_bias, For_bias, Back_bias)
  step = 0.2
  SweepSettings=( 
                ("drain", Other_bias,  "pre",  step, False, False, True, True,"dc"),
                ("gate",  For_bias, "forw", step, False, True, True, True, "dc" ),
                ("gate",  Back_bias, "back", step, False, True, True, True, "dc"),
              )

print(SweepSettings)
print("Description=",Description)
set_parameter(device=device, name="Description", value=Description)
QSSaveDevice(device, file="InitialSolve.tec",ftype="tecplot")
ExportParameters(device)

SweepModel=["Current", "NodeCharge"]   ###, "Capacitor"
ChargeContacts=["attached_source","attached_drain","gate"]
CurrentContacts=["drain"]

           ###{"region":"bulk","NodeName":"Electrons"},
for (SweepContact, Bias, Lable, step_limit, SaveZero, SaveFinal, PlotResults, ChargeUpdate, solvetype) in SweepSettings:
  if solvetype=="transient_bdf1":
    # input("AA")
    for contact, attached_to, is_circuit in (("gate", "Vgate", True), ):
      CreateContactCircuit(contact, attached_to=attached_to, resistance=resistance)
      CreateOxidePotentialContact(device, contact, attached_to=attached_to, is_circuit=is_circuit)
      last_bias=devsim.get_parameter(device=device, name=GetContactBiasName(contact))
      circuit_alter(name="V_%s"%contact, value=last_bias)

    
    for contact, attached_contact, DyBias, is_circuit in (("source", "attached_source", "Vsource", True), ("drain", "attached_drain", "Vdrain", True)):
      contact_bias=devsim.get_parameter(device=device, name=GetContactBiasName(contact))
      CreateContactCircuit(contact, attached_to=DyBias, resistance=1)
      CreateContactCircuit(attached_contact, attached_to=DyBias, resistance=1)
      CreateOxidePotentialContact(device, contact, attached_to=DyBias, offset_contact=contact, is_circuit=is_circuit)
      CreateOxidePotentialContact(device, attached_contact, attached_to=DyBias, offset_contact=contact, is_circuit=is_circuit)
      circuit_alter(name="V_%s"%contact, value=contact_bias)

  if DLT_region or ALT_region:
    if Lable=="forw1":
      SweepSpeed=20
    if Lable=="forw2":
      SweepSpeed=10
    if Lable=="forw3":
      SweepSpeed=5
    if Lable=="forw4":
      SweepSpeed=2
    if Lable=="forw5":
      SweepSpeed=1
    set_parameter(device=device, name="SweepSpeed", value=SweepSpeed)
  print(SweepContact, Bias, Lable, step_limit, SaveZero, SaveFinal)
  # input()
  Msg=VoltagePlotSweep(device, solvetype = solvetype, SweepContact=SweepContact, End_bias=Bias, 
                        FerroRegion=FerroRegion, SweepModel=SweepModel,
                        CurrentContacts=CurrentContacts,ChargeContacts=ChargeContacts,
                        VolumeIntegrateList=VolumeIntegrateList,NodeChecklist=NodeChecklist,
                        iterations=50, step_limit=step_limit, min_step = 0.001, rel_error=1e-10, Lable="%s_%2.2e"%(Lable,SweepSpeed),
                        CapacitorContacts=["gate"],  
                        ElementChecklist=ElementChecklist,
                        SaveAll=False, PlotResults = True, ExtraNodeChargeUpdate=ChargeUpdate)
# , PlotResults = False,SaveFinal=SaveFinal, SaveZero=SaveZero,SaveCurrentVariation=1.01,frequency=frequency,

print("Finished")
