from devsim import *
import sys, devsim
import os,  tarfile

sys.path.append('../../')
from QSsimple_physics import *
from QSPlotSweep import *
from QSFerro import *
from QSramp import *


#set extended 128bit for device
SetExtendedPrecision()
SetMicroMeterBasicParameters()

Description="Lw1000-t300Hw1000-t200Slope100nm"
device="test"
region="bulk"

# input("netx")

create_gmsh_mesh(file='cap-2step.msh', mesh=device)
add_gmsh_region    (mesh=device, gmsh_name="oxide",   region=region, material="Oxide"  )
add_gmsh_contact   (mesh=device, gmsh_name="top_contact",  region=region, name="top", material="metal")
add_gmsh_contact   (mesh=device, gmsh_name="bot_contact", region=region, name="bot", material="metal")

finalize_mesh(mesh=device)
create_device(mesh=device, device=device)

MicroMeterFerroParameters(device, region)
CreateFerroRegion(device, region)

# input("aa")
### Set the contact 
#for (contact,ContactRegion) in (("top","m1"), ("bot", "m2")):
#    CreateContactRegionModel(device, Contact=contact, ContactRegion=ContactRegion)
    
### Contact models and equations
for contact in ("top", "bot"):
    CreateFerroContactEquation(device, contact)

set_parameter(device=device, name="bot_bias", value=0.0)
set_parameter(device=device, name="top_bias", value=0)
solve(type="dc", absolute_error=1e-6, relative_error=1e-4, maximum_iterations=30)
set_parameter(device=device, name="SweepSpeed",        value=10)

set_parameter(device=device, name="Description", value=Description)
write_devices(file="capacitor.dev", type="devsim")
write_devices(file="capacitor.tec", type="tecplot")

# set_parameter(device=device, region=region, name="InitialCoefficient", value=1)
# set_parameter(device=device, region=region, name="Iterate", value=False)
# set_element_values(device=device, region=region, name="OldPCoefficientX", init_from="StartCoefficient")
# set_element_values(device=device, region=region, name="OldPCoefficientY", init_from="StartCoefficient")

# set_parameter(device=device, region=region, name="StepByStep", value=True)
    #,"CoerciveSignY","PreElectricField_y", "ElectricField_y",
    #,"PolarizationY","Polarization", "OldPCoefficientY","NewPCoefficientY","NewPCoefficientY","OldCoerciveSignY","CoerciveSignY","PreElectricField_y","ElectricField_y"
Checklist=()#("CoerciveSignY","OldPCoefficientY","NewPCoefficientY","NewPCoefficientY", "NumeratorTanhY","DenominatorTanhY","PolarizationY","PreElectricField_y","ElectricField_y")
#"SweepDirectionY",
V=60
for SV in (V,-V, V):
	VoltagePlotSweep(device, 
                    Lable           = "",
                    SweepModel      = ["Capacitor"],
                    solvetype       = "dc",
                    SweepContact    = 'top',
                    CurrentContacts = [],
                    ChargeContacts  = [],
                    CapacitorContacts    = ['top'],
                    End_bias        = SV,
                    step_limit      = 0.4,
                    min_step        = 0.01, 
                    rel_error       = 1e-8, 
                    abs_error       = 1e30,
                    iterations      = 30,
                    frequency       = None,
                    FerroRegion     = region,
                    ElementChecklist= Checklist,
                    SaveAll         = True,
                    SaveData        = True,
                    SaveFinal       = True,
                    PlotResults     = True,)

#input("finished")
