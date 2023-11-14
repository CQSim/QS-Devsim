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
#sys.path.append('../')
from QSsimple_physics import *
create_gmsh_mesh(file='mos_topContact.msh', mesh='mos_topContact')
device="mos_topContact"

# device_width        =   1.0 #um
# gate_width          =   0.8 #um

# bulk_thickness      =   4e-2 #um
# air_thickness       =   4e-1 #um
# oxide_thickness     =   2e-1 #um
# diffusion_thickness =   5e-3 #um

# x_diffusion_decay   =   5e-3 #um
# y_diffusion_decay   =   5e-3 #um

# refine_spacing      =   1e-2 #um

# x_bulk_left =    0.0
# x_bulk_right =   x_bulk_left + device_width
# x_center =       0.5 * (x_bulk_left + x_bulk_right)
# x_gate_left =    x_center - 0.5 * (gate_width)+x_diffusion_decay
# x_gate_right =   x_center + 0.5 * (gate_width)-x_diffusion_decay
# x_device_left =  x_bulk_left - air_thickness
# x_device_right = x_bulk_right + air_thickness

# y_bulk_top =       0.0
# y_oxide_top =    y_bulk_top - oxide_thickness
# y_oxide_mid    =    0.5 * (y_oxide_top + y_bulk_top)
# y_bulk_bottom=    y_bulk_top + bulk_thickness
# y_bulk_mid    =    0.5 * (y_bulk_top + y_bulk_bottom)
# y_device_bottom=y_bulk_bottom + air_thickness
# y_air_top    =    y_bulk_bottom + air_thickness
# y_diffusion    =    y_bulk_bottom - diffusion_thickness


add_gmsh_region    (mesh="mos_topContact", gmsh_name="bulk",    region="bulk",  material="Silicon")
add_gmsh_region    (mesh="mos_topContact", gmsh_name="oxide",   region="oxide", material="Oxide"  )
add_gmsh_region    (mesh="mos_topContact", gmsh_name="air",    region="air",  material="Air")
add_gmsh_contact   (mesh="mos_topContact", gmsh_name="drain_contact",  region="bulk", name="drain", material="metal")
add_gmsh_contact   (mesh="mos_topContact", gmsh_name="source_contact", region="bulk", name="source", material="metal")
add_gmsh_contact   (mesh="mos_topContact", gmsh_name="body_contact",   region="air", name="ground", material="metal")

add_gmsh_contact   (mesh="mos_topContact", gmsh_name="gate_contact",   region="oxide", name="gate", material="metal")
add_gmsh_contact   (mesh="mos_topContact", gmsh_name="attached_source_contact", region="air", name="attached_source", material="metal")
add_gmsh_contact   (mesh="mos_topContact", gmsh_name="attached_drain_contact",  region="air", name="attached_drain",  material="metal")

add_gmsh_interface (mesh="mos_topContact", gmsh_name="bulk_air_interface", region0="bulk", region1="air", name="bulk_air")
add_gmsh_interface (mesh="mos_topContact", gmsh_name="bulk_oxide_interface", region0="bulk", region1="oxide", name="bulk_oxide")
finalize_mesh(mesh="mos_topContact")
create_device(mesh="mos_topContact", device=device)


