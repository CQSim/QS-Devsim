# Copyright 2013 Devsim LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from devsim import *
sys.path.append('/home/ccc/devsim/QS_Micro_packages/')
sys.path.append('../')
from QSsimple_physics import *
device="mos_BotContact"
create_gmsh_mesh(file='mos_BotContact-intLay2.msh', mesh=device)

device_width        =   1.0 #um
gate_width          =   0.8 #um

bulk_thickness      =   4e-2 #um
air_thickness       =   4e-1 #um
oxide_thickness     =   2e-1 #um
diffusion_thickness =   5e-3 #um

x_diffusion_decay   =   5e-3 #um
y_diffusion_decay   =   5e-3 #um

refine_spacing      =   1e-2 #um

x_bulk_left =    0.0
x_bulk_right =   x_bulk_left + device_width
x_center =       0.5 * (x_bulk_left + x_bulk_right)
x_gate_left =    x_center - 0.5 * (gate_width)+x_diffusion_decay
x_gate_right =   x_center + 0.5 * (gate_width)-x_diffusion_decay
x_device_left =  x_bulk_left - air_thickness
x_device_right = x_bulk_right + air_thickness

y_bulk_top =       0.0
y_oxide_top =    y_bulk_top - oxide_thickness
y_oxide_mid    =    0.5 * (y_oxide_top + y_bulk_top)
y_bulk_bottom=    y_bulk_top + bulk_thickness
y_bulk_mid    =    0.5 * (y_bulk_top + y_bulk_bottom)
y_device_bottom=y_bulk_bottom + air_thickness
y_air_top    =    y_bulk_bottom + air_thickness
y_diffusion    =    y_bulk_top + diffusion_thickness


add_gmsh_region    (mesh=device, gmsh_name="bulk",    region="bulk",  material="Silicon")
add_gmsh_region    (mesh=device, gmsh_name="oxide",   region="oxide", material="Oxide"  )
add_gmsh_region    (mesh=device, gmsh_name="air",    region="air",  material="Air")
add_gmsh_contact   (mesh=device, gmsh_name="drain_contact",  region="bulk", name="drain", material="metal")
add_gmsh_contact   (mesh=device, gmsh_name="source_contact", region="bulk", name="source", material="metal")
add_gmsh_contact   (mesh=device, gmsh_name="body_contact",   region="air", name="ground", material="metal")

add_gmsh_contact   (mesh=device, gmsh_name="gate_contact",   region="oxide", name="gate", material="metal")
add_gmsh_contact   (mesh=device, gmsh_name="attached_source_contact", region="oxide", name="attached_source", material="metal")
add_gmsh_contact   (mesh=device, gmsh_name="attached_drain_contact",  region="oxide", name="attached_drain",  material="metal")

add_gmsh_interface (mesh=device, gmsh_name="bulk_air_interface", region0="bulk", region1="air", name="bulk_air")
add_gmsh_interface (mesh=device, gmsh_name="bulk_oxide_interface", region0="bulk", region1="oxide", name="bulk_oxide")
finalize_mesh(mesh=device)
create_device(mesh=device, device=device)

#p doping
bulk_doping     =   1e3    #1um^3
body_doping     =   1e7    #1um^3
#n doping
drain_doping    =   -1e7    #1um^3
source_doping   =   -1e7    #1um^3    
gate_doping     =   -1e8    #1um^3

SetDopingParameters(device,  "bulk",
                                    bulk_doping  =bulk_doping, 
                                    body_doping  =body_doping, 
                                    drain_doping =drain_doping,
                                    source_doping=source_doping,)

#### all variable substitutions are immediate, since they are locked into the mesh
mydict = {}
mydict["drain_doping"] = drain_doping
mydict["body_doping"] = body_doping
mydict["gate_doping"] = gate_doping
mydict["source_doping"] = source_doping
mydict["bulk_doping"] = bulk_doping
mydict["x_gate_left"] = x_gate_left
mydict["x_gate_right"] = x_gate_right
mydict["x_diffusion_decay"] = x_diffusion_decay
mydict["y_diffusion"] = y_diffusion
mydict["y_bulk_bottom"] = y_bulk_bottom
mydict["y_diffusion_decay"] = y_diffusion_decay

# node_model(name="NetDoping",    device=device, region="gate", equation="%(gate_doping)1.15e + 1" % (mydict))

node_model(name="DrainDoping",  device=device, region="bulk", equation="0.25*%(drain_doping)1.15e*erfc((x-%(x_gate_left)1.15e)/%(x_diffusion_decay)1.15e)*erfc(-(%(y_diffusion)1.15e-y)/%(y_diffusion_decay)1.15e)" % mydict)
node_model(name="SourceDoping", device=device, region="bulk", equation="0.25*%(source_doping)1.15e*erfc(-(x-%(x_gate_right)1.15e)/%(x_diffusion_decay)1.15e)*erfc(-(%(y_diffusion)1.15e-y)/%(y_diffusion_decay)1.15e)" % mydict)
# node_model(name="BodyDoping",   device=device, region="bulk", equation="0.5*%(body_doping)1.15e*erfc(-(y-%(y_bulk_bottom)1.15e)/%(y_diffusion_decay)1.15e)" % mydict)
node_model(name="NetDoping",    device=device, region="bulk", equation="DrainDoping + SourceDoping  +%(bulk_doping)1.15e + 1" % mydict)


# write_devices(file="doping.tec",ftype="tecplot")

