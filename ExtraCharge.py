from devsim.python_packages.simple_dd import *
from devsim import *
sys.path.append('/home/ccc/devsim/QS_Micro_packages/')
from QSmodel_create import *


def SetGaussianDOSParameters(device, region, Type, Width, Center, TotalDensity):
  # GaussianWidth (V) width of DOS
  # GaussianCenter (V) the value of the center of Gaussian DOS above intrinsic Fermi level 

  set_parameter(device=device, region=region, name = "%s_GaussianWidth"%Type,   value=Width   )
  set_parameter(device=device, region=region, name = "%s_GaussianCenter"%Type,  value=Center  )
  set_parameter(device=device, region=region, name = "%s_GaussianTotalDensity"%Type,  value=TotalDensity  )

def CreateGaussianDOS(device, region, Type, Width, Center, TotalDensity):
  # Degenerate Approximation Model
  # Delta (V) width of DOS
  # GaussianCenter (V) the value of the center of Gaussian DOS below intrinsic Fermi level 
  SetGaussianDOSParameters(device, region, Type, Width, Center, TotalDensity)
  if Type=="DonorLike":
    DravateVariable="Holes"
    eq="{0}_GaussianTotalDensity*0.5*(1+erf((log({1}/n_i)*V_t+{0}_GaussianCenter)/{0}_GaussianWidth/1.414))".format(Type, DravateVariable)
  elif Type=="AccepterLike":
    DravateVariable="Electrons"
    eq="{0}_GaussianTotalDensity*0.5*(1+erf((log({1}/n_i)*V_t-{0}_GaussianCenter)/{0}_GaussianWidth/1.414))".format(Type, DravateVariable)
  else:
    raise RuntimeError("Please Set the type of _GaussianDOS as: DonorLike or AccepterLike ")
  print("the equation of %s_GaussianDOS is %s "%(Type,eq))
  NameOfNode="%s_GaussianDOS"%Type
  CreateNodeModel(device, region, NameOfNode , eq   )
  CreateNodeModelDerivative(device, region, NameOfNode, eq, DravateVariable)
  eq2=" log({1}/n_i)*V_t-{0}_GaussianCenter".format(Type, "Electrons")
  CreateNodeModel(device, region, "QuaisLag" , eq2   )
  return NameOfNode

def CreateGaussianDOS_PotentialEquation(device, region, interface=None,Thickness=0.05):
  Charge_Type="GaussianDOS"
  CreateExtraCharge_PotentialEquation(device, region, Charge_Type, interface=interface,Thickness=Thickness)


def CreateExtraCharge_PotentialEquation(device, region, Charge_Type, interface=None,Thickness=0.05):
  DeleteNodeModelAndDerivates(device, region, "PotentialNodeCharge")

  DonorLikeCharge="DonorLike_%s"%(Charge_Type)
  AccpterlikeCharge="AccepterLike_%s"%(Charge_Type)
  if interface!=None:
    if InNodeModelList(device, region, DonorLikeCharge):
      DonorLikeCharge=InterfaceCharge_NodeModel(device=device, interface=interface, region=region, Charge_Type=DonorLikeCharge,Thickness=Thickness)
    if InNodeModelList(device, region, AccpterlikeCharge):
      AccpterlikeCharge=InterfaceCharge_NodeModel(device=device, interface=interface, region=region, Charge_Type=AccpterlikeCharge,Thickness=Thickness)

  if InNodeModelList(device, region, DonorLikeCharge) and InNodeModelList(device, region, AccpterlikeCharge):
    print("Assemble {0} and {1} {2} Charge in total charge of Potential Equation".format(DonorLikeCharge,AccpterlikeCharge,Charge_Type))
    ExtraCharge="{0}-{1}".format(DonorLikeCharge,AccpterlikeCharge)
  elif InNodeModelList(device, region, DonorLikeCharge):
    print("Only Assemble Donor like %s Charge in Potential Equation"%(Charge_Type,))
    ExtraCharge="{0}".format(DonorLikeCharge)
  elif InNodeModelList(device, region, AccpterlikeCharge):
    print("Only Assemble Accepter like %s Charge in Potential Equation"%(Charge_Type,))
    ExtraCharge="-{0}".format(AccpterlikeCharge)
  else:
    raise RuntimeError("Could not find any type of extar charge as %s"%Charge_Type,)
  CreateNodeModel(device, region, "ExtraCharge", ExtraCharge)
  CreateNodeModelDerivative(device, region, "ExtraCharge", ExtraCharge, "Potential", "Electrons", "Holes")
  pne = "-ElectronCharge*kahan4(Holes, -Electrons, NetDoping, ExtraCharge)"
  CreateNodeModel(device, region, "PotentialNodeCharge", pne)
  CreateNodeModelDerivative(device, region, "PotentialNodeCharge", pne, "Potential", "Electrons", "Holes")

  equation(device=device, region=region, name="PotentialEquation", variable_name="Potential",
            node_model="PotentialNodeCharge", edge_model="PotentialEdgeFlux",
            time_node_model="", variable_update="log_damp")

def InterfaceCharge_NodeModel(device, interface, region, Charge_Type, Thickness):
  #get the nodes on interface
  # for i, j in enumerate(get_region_list(device=device, interface=interface)):
  #   if j == region:
  #     eq1 = "node_index@r%d" %(i,)
  #     break
  # interface_model(device=device, interface=interface, name="interface_node_indexes", equation=eq1)
  # interface_node_indexes = get_interface_model_values(device=device, interface=interface, name="interface_node_indexes")
  # node_index=get_node_model_values (device=device, region=region, name="node_index")
  # Atinterface=[int(i in interface_node_indexes) for i in node_index]
  # NameAtInterface="At{0}_{1}".format(interface,region)
  # node_solution(device=device, region=region, name=NameAtInterface)
  # set_node_values(device=device, region=region, name=NameAtInterface, values=Atinterface)

  ChargeName="%s_At_%s"%(Charge_Type,interface,)
  
  eq2="{0}*exp(-y/{1})".format(Charge_Type,Thickness)
  # eq2="{0}*{1}*SurfaceArea/NodeVolume".format(NameAtInterface,Charge_Type)
  CreateNodeModel(device, region, ChargeName , eq2)
  CreateNodeModelDerivative(device, region, ChargeName, eq2, "Potential", "Electrons", "Holes")
  return ChargeName

def ExpotionalInterfaceCharge_NodeModel(device, interface, region, Charge_Type, Distance):
  pass