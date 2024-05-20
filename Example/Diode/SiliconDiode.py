import devsim,sys

sys.path.append('../../')
from QSsimple_physics import *
from QSPlotSweep import *
import QSdiode_common
Description=""
#set extended 128bit for device
SetExtendedPrecision()
# SetVerbosity(verbose=True)
SetCentiMeterBasicParameters()

device  = "GaussianDiode"
region  = "bulk"
contacts = ("bot", "top")

Length=2e-5
refinescale=1e-8
DDoping=ADoping=2e18
Description="%sD%0.0e"%(Description,DDoping)
QSdiode_common.CreateMesh(device, region, Length,refinescale)
QSdiode_common.SetNetDoping(device=device, region=region, DonorDoping=DDoping, AccepterDoping=ADoping, Length=Length)
CreatePotentialAndFlux(device, region)

mun=mup=1e3
taun=taup=1e-5
Description="%sMu%0.0e"%(Description,mun)

SetCentiMeterSiliconParameters(device, region, 300, pMobility=mup, nMobility=mun, nLifeTime=taun, pLifeTime=taup)
CreateSiliconPotentialOnly(device, region)

for c in contacts:
  CreateSiliconPotentialOnlyContact(device, c)

InitialSolve(device,rel_error=1e-10)

CreateElectronAndHoleSolutions(device, region)
CreateSiliconDriftDiffusion(device, region)
DeleteIntrinsicNodeModelDerivatives(device, region)
  
for c in contacts:
  CreateSiliconDriftDiffusionAtContact(device, c)
InitialSolve(device,rel_error=1e-6)

## Set models for monitor during simulation
DeviceMonitorList=[ {"region":region, "ModelType":"NodeModel", "ylogScale":False, "ModelNames":["Electrons","Holes"] },]#,"USRH"
DeviceMonitorList.append({"region":region, "ModelType":"NodeModel", "ylogScale":False, "ModelNames":["ElectronGeneration","HoleGeneration"] })
DeviceMonitorList.append({"region":region, "ModelType":"NodeModel", "ylogScale":False, "ModelNames":["Potential","QuasiFermiPotential_Electrons","QuasiFermiPotential_Holes"] })

QSSaveDevice(device, file="InitialSolve.dev",ftype="devsim")

InitialSolve(device, orders=10, rel_error=1e-6)

for c in contacts:
  PrintCurrents(device, c)

VolumeIntegrateList=[]
NodeChecklist=[]

CreateAbsCurrentAndEfield(device, region, magnitude=True)

if True:
  Start_bias = -1
  Stop_bias  = 3
  Description="%sT%s~%s"%(Description, Start_bias,Stop_bias)
  step = 0.1
  SweepSettings=( 
                ("top", Start_bias,   "pre",  step, False, False,True),
                # ("top",  0,  "star", step, False, True, True ),
                # ("top",  1,  "star", step, False, True, True ),
                # ("top",  2,  "star", step, False, True, True ),
                ("top",  Stop_bias,  "star", step, False, True, True ),
             )
  
print(SweepSettings)
print("Description=",Description)
set_parameter(device=device, name="Description", value=Description)
QSSaveDevice(device, file="InitialSolve.tec",ftype="tecplot")
ExportParameters(device)
# input("aa")
SweepModel=["Current", "NodeCharge"]   ###, "Capacitor"
ChargeContacts=[]
CurrentContacts=["top"]

frequency=ElementChecklist=[]

           ###{"region":"bulk","NodeName":"Electrons"},
for (SweepContact, Bias, Lable, step_limit, SaveZero, SaveFinal,PlotResults) in SweepSettings:
  print(SweepContact, Bias, Lable, step_limit, SaveZero, SaveFinal)
  Msg=VoltagePlotSweep(device, SweepContact=SweepContact, End_bias=Bias,  SweepModel=SweepModel,
                        CurrentContacts=CurrentContacts,ChargeContacts=ChargeContacts,
                        VolumeIntegrateList=VolumeIntegrateList,DeviceMonitorList=DeviceMonitorList,NodeChecklist=NodeChecklist,
                        iterations=30, step_limit=step_limit, rel_error=1e-10, Lable=Lable,
                        ElementChecklist=ElementChecklist, SaveAll=True, PlotResults = True)
print("Finished")
