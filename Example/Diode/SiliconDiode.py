#! /home/ccc/miniconda3/bin/python
#词句指定本文件执行的解释程序，也就是python所在的路径，指定后可以在终端中直接执行本文件。本局必须为首局

# 导入必备模块，“devsim”即使所使用的Devsim求解器
import devsim,sys
#将QS-Devsim模块所在路径导入，然后才能导入具体模块
sys.path.append('/home/ccc/devsim/QS_Micro_packages/')
#导入本文件所在路径（当前路径），然后才能导入当前路径下的模块
sys.path.append('./')

# 导入QS-Devsim模块中的QSsimple_physics（基本物理模型）库中的所有函数
from QSsimple_physics import *
# 导入QS-Devsim模块中的QSPlotSweep（扫描电压模块，包含I_V和C-V曲线扫描）库中的所有函数
from QSPlotSweep import *
# 导入QS-Devsim模块中的QSdiode_common（二极管相关模型）库中的所有函数
import QSdiode_common
# Description用于器件特征描述，Devsim中具体器件为一个对象，Description字符窜将用于器件的Description属性，然后创建文件夹用于保存器件的状态文件和扫描I_V和C-V曲线结果文件。
Description=""

# 设置使用128位浮点型进行计算，对部分器件而言，关台电流依然重要，常规64位无法精确模拟，因此需要使用更高精度的128位，但目前计算机是64位的，128位是linux系统通过软件实现，计算效率较低
#set extended 128bit for device
SetExtendedPrecision()

# 设置使用厘米为长度基本单位的物理参数，包含：元电荷电量、玻尔兹曼常数、真空介电常数
SetCentiMeterBasicParameters()

# device为器件名称，用于创建Devsim中的device对象
device  = "GaussianDiode"
# region 区域名称，一个器件由多个区域构成，每个区域是同一种的材料，具有相同的物理特性，
#如：氧化物区域属于介电材料，仅需要计算电势和电场，而半导体区域则需要再加上电子和空穴的浓度
region  = "bulk"
#电极名称，二极管属于两端器件，仅有两个电极：顶电极“top”和底电极“bot”，FET则有三个“源-source、漏-drain、栅-gate”或四个电极“外加一个背底材料，体电极-body”，
contacts = ("bot", "top")

##下面将使用QSdiode_common模块中的函数创建器件几何结构
#器件长度，由于二极管可以简化为一维进行模拟，因此仅需要指定器件的长度，由于前面使用SetCentiMeterBasicParameters()函数，所以两个维度为长度基本单位为：厘米
Length=2e-5
# 二极管器件中网格划分的最小点距，在PN结处各种物理属性变化较大，此处的点距需要设置比较精细，才能计算出符合实际的结果
refinescale=1e-8
#设置P区和N区的参杂浓度
DDoping=ADoping=2e18
#器件的Description属性，如上
Description="%sD%0.0e"%(Description,DDoping)

#调用QSdiode_common模块中CreateMesh函数创建器件结构
QSdiode_common.CreateMesh(device, region, Length,refinescale)

#调用QSdiode_common模块中SetNetDoping函数设置器件参杂浓度，此时会创建两个节点变量（NodeModel），用于存储给体和受体的参杂浓度信息
QSdiode_common.SetNetDoping(device=device, region=region, DonorDoping=DDoping, AccepterDoping=ADoping, Length=Length)

#创建电势节点变量（NodeModel），同时根据此变量创建电场类矢量（此处电场用edgemodel，一维情况与矢量相同，二维和三维情况下，矢量（element_model）需要笛卡尔坐标系，edgemodel不需要）
CreatePotentialAndFlux(device, region)

#设置电子和空穴的迁移率
mun=mup=1e3
##设置电子和空穴的寿命
taun=taup=1e-5
#器件的Description属性，如上
Description="%sMu%0.0e"%(Description,mun)

# 使用厘米为长度基本单位设置半导体物理参数，包含：温度、玻电子和空穴的迁移率、子和空穴的寿命
SetCentiMeterSiliconParameters(device, region, 300, pMobility=mup, nMobility=mun, nLifeTime=taun, pLifeTime=taup)
# 以0为器件电势，设置空穴和电子初始浓度，用作迭代计算的初值
CreateSiliconPotentialOnly(device, region)

# 设置初始态电极边界条件（仅有电势一个变量的边界条件），此处使用欧姆接触
for c in contacts:
  CreateSiliconPotentialOnlyContact(device, c)

#根据电极的欧姆接触和参杂浓度计算平衡态电势分布（此时器件中仅有电势的高斯定理方程）
InitialSolve(device,rel_error=1e-10)

#创建电子和空穴节点变量（NodeModel）
CreateElectronAndHoleSolutions(device, region)
#创建电子和空穴的edgemodel（类矢量）和对应的连续性方程
CreateSiliconDriftDiffusion(device, region)
#删除计算初始平衡态用到的变量
DeleteIntrinsicNodeModelDerivatives(device, region)
  
#设置电子和空穴的欧姆接触边界条件
for c in contacts:
  CreateSiliconDriftDiffusionAtContact(device, c)

#求解初始条件
InitialSolve(device,rel_error=1e-6)


#设置需要在图像中监视的变量，以字典形式设置，需指定区域、类型、坐标方式、模型列表
## Set models for monitor during simulation
DeviceMonitorList=[ {"region":region, "ModelType":"NodeModel", "ylogScale":False, "ModelNames":["Electrons","Holes"] },]#,"USRH"
DeviceMonitorList.append({"region":region, "ModelType":"NodeModel", "ylogScale":False, "ModelNames":["ElectronGeneration","HoleGeneration"] })
DeviceMonitorList.append({"region":region, "ModelType":"NodeModel", "ylogScale":False, "ModelNames":["Potential","QuasiFermiPotential_Electrons","QuasiFermiPotential_Holes"] })

#保存器件数据到当前路径，此时还未添加Description属性，不会创建新文件夹
QSSaveDevice(device, file="InitialSolve.dev",ftype="devsim")

#求解器件
InitialSolve(device, orders=10, rel_error=1e-6)

#输出器件当前contacts列表中电极的电流
for c in contacts:
  PrintCurrents(device, c)

#二维和三维情况下，由电场类矢量（edgemodel）创建电场矢量（element_model）
CreateAbsCurrentAndEfield(device, region, magnitude=True)

if True:
  Start_bias = -1 #第一次扫描电压的终止值
  Stop_bias  = 3 #第二次扫描电压的终止值
  Description="%sT%s~%s"%(Description, Start_bias,Stop_bias)
  step = 0.1 ##描电压的步长
  SweepSettings=( #扫描电压的每一步参数,包含：电极、终止电压、标签（用于导出文件标记）、步长、零电压是否必须经过、最后状态是否保存
                ("top", Start_bias,   "pre",  step, False, False,True),#第一次扫描电压的参数
                ("top",  Stop_bias,  "star", step, False, True, True ),#第二次扫描电压的参数
             )
# 输出扫描参数和Description
print(SweepSettings)
print("Description=",Description)
#设置Description属性
set_parameter(device=device, name="Description", value=Description)
#保存器件数据到当前路径，此时已经添加Description属性，会创建新文件夹
QSSaveDevice(device, file="InitialSolve.tec",ftype="tecplot")
#输出device的当前所有参数
ExportParameters(device)
#电压扫描模式，即需要输出的物理信息，包含电流、节点数据，可以包含Capacitor
SweepModel=["Current", "NodeCharge"]   ###, "Capacitor"
#扫描时需要导出电荷量的电极
ChargeContacts=[]
#扫描时电极
CurrentContacts=["top"]

#执行扫描过程
for (SweepContact, Bias, Lable, step_limit, SaveZero, SaveFinal,PlotResults) in SweepSettings:
  print(SweepContact, Bias, Lable, step_limit, SaveZero, SaveFinal)
  Msg=VoltagePlotSweep(device, SweepContact=SweepContact, End_bias=Bias,  SweepModel=SweepModel,
                        CurrentContacts=CurrentContacts,ChargeContacts=ChargeContacts,
                        VolumeIntegrateList=VolumeIntegrateList,DeviceMonitorList=DeviceMonitorList,NodeChecklist=NodeChecklist,
                        iterations=30, step_limit=step_limit, rel_error=1e-10, Lable=Lable,
                        ElementChecklist=ElementChecklist, SaveAll=True, PlotResults = True)
print("Finished")
