from devsim import *
import sys, devsim
import os,  tarfile

sys.path.append('/home/ccc/devsim/QS_Micro_packages/')
sys.path.append('../')
from QSsimple_physics import *
from QSPlotSweep import *
from QSFerro import *
from QSramp import *


#set extended 128bit for device
SetExtendedPrecision()
SetMicroMeterBasicParameters()

import create_mos

device          = "mos_BotContact"
silicon_regions = ("bulk",)
oxide_region    = "oxide"
air_region      = "air"
oxide_regions   = ("oxide",)
regions         = ("bulk", "oxide", air_region)
interfaces      = ("bulk_oxide", "bulk_air")
FerroRegion     = oxide_region

for i in regions:
  CreatePotentialAndFlux(device, i)
mu=1e5
for i in silicon_regions:
  SetMicroMeterSiliconParameters(device, i, 300, pMobility=mu, nMobility=mu)
  CreateSiliconPotentialOnly(device, i)

for (c,r,k) in (("ground",air_region,1.0),
                ("gate",oxide_region,3.9)):
  SetOxideParameters(device, r, rPermittivity=k)
  CreateOxidePotentialOnly(device, r, "log_damp")
  CreateOxidePotentialContact(device, c)

contacts=("source", "drain" )
for c in contacts:
  tmp = get_region_list(device=device, contact=c)
  r = tmp[0]
  print("%s %s" % (r, c))
  SetContactDopingOffset(device, c)
  CreateSiliconPotentialOnlyContact(device, c)

attached_contacts=(("attached_source","source"),  ("attached_drain", "drain"))
for (i,j) in attached_contacts:
  SetContactPotentilOffset(device, j, GetContactDopingOffset(device, j))
  CreateOxidePotentialContact(device, i, attached_to=j, offset_contact=j)

for i in interfaces:
  CreateSiliconOxideInterface(device, i)
#, type="hybrid"


InitialSolve(device,rel_error=1e-9)
# input("aa")

for i in silicon_regions:
  CreateElectronAndHoleSolutions(device, i)
  CreateSiliconDriftDiffusion(device, i)

for c in ("source", "drain"):
  tmp = get_region_list(device=device, contact=c)
  r = tmp[0]
  CreateSiliconDriftDiffusionAtContact(device, c)

V_Sat=1.0e11

for r in ("bulk",):
  DeleteEdgeModelAndDerivatives(device, r, "HoleCurrent")
  ElectricFieldMagnitude(device, r)
  MicorSaturationMobilityParameters(device, r,V_n_Sat=V_Sat, V_p_Sat=V_Sat)
  QSSaturationMobility(device, r, "Hole")
  # LorentzSaturationMobility(device, r, "mu_n_Sat", "mu_n_eff",  "LorentzEField")
  ElementHoleCurrent2d(device, r)
  ElementHoleContinuityEquation(device, r)

for contact in ( "drain", "source"):
  ElementContactHoleContinuityEquation(device, contact) 

InitialSolve(device, rel_error=1e-9)
PrintCurrents(device, "drain")
# input('bb')


if FerroRegion!=None:
  print("*******Creat Ferro and Solve")
  MicroMeterFerroParameters(device, FerroRegion) 
  CreateFerroRegion(device, FerroRegion,  update_type="log_damp")
  CreateOxidePotentialContact(device,"gate", element_contact=True)
  for (i,j) in attached_contacts:
    SetContactPotentilOffset(device, i, GetContactDopingOffset(device, j))
    CreateOxidePotentialContact(device, i, attached_to=j, offset_contact=j, element_contact=True)
  SpecificValueInitiateRamp(device, FerroRegion, "SaturationPolarization", abs_error=1e30, rel_error=1e-9, iterations=30)

# input('cc')

for r in regions:
  CreateAbsCurrentAndEfield(device, r, magnitude=True)
# input('aa')

ElementChecklist=()#("SweepDirectionY","CoerciveSignY","OldPCoefficientY","NewPCoefficientY", "NumeratorTanhY","DenominatorTanhY","PolarizationY","PreElectricField_y","ElectricField_y")

Coersive_bias  = 10.00
Other_bias     = -5 
#Star_bias      = max(Other_bias,0)+3*Coersive_bias
Star_bias      = min(Other_bias,0)-3*Coersive_bias

End_bias       = 16 #ther_bias,0)-3*Coersive_bias
# -abs(Star_bias)/Star_bias * Coersive_bias +Other_bias/2
#    
Description="1e-16V_Sat%0.2e_mu%0.2e_D%sG%s~%s"%(V_Sat, mu, Other_bias, Star_bias, End_bias )
# Description=""
print("Description=",Description)
set_parameter(device=device, name="Description", value=Description)
QSSaveDevice(device, file="InitialSolve.tec",ftype="tecplot")
set_parameter(device=device, name="SweepSpeed", value=10)

ChargeContacts=["attached_source","attached_drain","gate"]
SweepSettings=(   ("drain", Other_bias, "pre",  0.5, False, False),
                  ("gate",  Star_bias,  "star", 0.5, False, True ),
                  ("gate",  End_bias,   "forw", 0.1, True,  True ),
                  ("gate",  Star_bias,  "back", 0.1, True,  True ),
              )
print(SweepSettings)

for (SweepContact, Bias, Lable, step_limit, SaveZero, SaveFinal) in SweepSettings:
  Msg=VoltagePlotSweep(device, SweepContact=SweepContact, CurrentContacts = ["drain"],End_bias=Bias, FerroRegion=FerroRegion, ChargeContacts=ChargeContacts,\
                        iterations=30, step_limit=step_limit, min_step = 0.005, rel_error=1e-18, Lable=Lable,\
                        ElementChecklist=ElementChecklist, SaveAll=False, PlotResults = True)
# , PlotResults = False,SaveCurrentVariation=1.01, SaveZero=SaveZero, SaveFinal=SaveFinal,SaveCurrentVariation=1.01
print("Finished")



