import math,sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from devsim import *
sys.path.append('/home/ccc/devsim/QS_Micro_packages/')
from QSsimple_physics import *

def CentiMeterFerroParameters(device, region, 
                        SaturationPolarization  =   12e-6, #C/cm^2
                        RemanentPolarization    =   10e-6, #C/cm^2
                        CoerciveField           =   5.0e5, #V/cm
                        eps_Ferro               =   12,
                        Iterate                 =   True,
                        FunctionModel           =   "tanh"):
    if RemanentPolarization>SaturationPolarization or RemanentPolarization<0:
        raise NameError('RemanentPolarization could not larger than SaturationPolarization or negetive')
    WriteFerroParameters(device, region,  SaturationPolarization,RemanentPolarization,CoerciveField, eps_Ferro ,Iterate, FunctionModel)

def MicroMeterFerroParameters(device, region, 
                        SaturationPolarization  =   12e-14, #C/um^2
                        RemanentPolarization    =   10e-14, #C/um^2
                        CoerciveField           =   5e1, #V/um
                        eps_Ferro               =   12,
                        Iterate                 =   True, 
                        FunctionModel           =   "tanh"):
    if RemanentPolarization>SaturationPolarization or RemanentPolarization<0:
        raise NameError('RemanentPolarization could not larger than SaturationPolarization or negetive')
    WriteFerroParameters(device, region,  SaturationPolarization,RemanentPolarization,CoerciveField, eps_Ferro ,Iterate, FunctionModel)

def WriteFerroParameters(device, region,  SaturationPolarization,RemanentPolarization,CoerciveField, eps_Ferro ,Iterate, FunctionModel):
    Omega=math.log((SaturationPolarization+RemanentPolarization)/(SaturationPolarization-RemanentPolarization))/CoerciveField/2.0
    if Iterate :
        if FunctionModel=="erf":
            InitialCoefficient=1.0/(math.erf(Omega*CoerciveField)+1.0)
        elif FunctionModel=="tanh":
            InitialCoefficient=1.0/(math.tanh(Omega*CoerciveField)+1.0)
    else:
        InitialCoefficient=0.9
    Permittivity = eps_Ferro * get_parameter(name="eps0")
    set_parameter(device=device, region=region, name="Permittivity",            value=Permittivity  )
    set_parameter(device=device, region=region, name="SaturationPolarization",  value=SaturationPolarization)
    set_parameter(device=device, region=region, name="RemanentPolarization",    value=RemanentPolarization)
    set_parameter(device=device, region=region, name="CoerciveField",           value=CoerciveField)
    set_parameter(device=device, region=region, name="Omega",                   value=Omega)
    set_parameter(device=device, region=region, name="InitialCoefficient",      value=InitialCoefficient)
    set_parameter(device=device, region=region, name="Iterate",                 value=Iterate)
    set_parameter(device=device, region=region, name="StepByStep",              value=False)
    set_parameter(device=device, region=region, name="RampPolarization",        value=False)
    print(get_dimension(device=device))
    if get_dimension(device=device)==2:
        DimList=("X","Y")
    elif get_dimension(device=device)==3:
        DimList=("X","Y","Z")
    set_parameter(device=device, region=region, name="DimList", value=DimList)



def CreateFerroRegion(device, region, update_type="default", FunctionModel="tanh"):
    ###Import the math function of tanh and it differential function
    # symdiff(expr="declare(tanh(x))")
    # symdiff(expr="define(tanh(x), 1-pow(tanh(x),2))")
    # register_function(name="tanh", nargs=1, procedure=math.tanh)

    ###
    ### Create the Potential solution variable

    if not InNodeModelList(device, region, "Potential"):
        node_solution(device=device, region=region, name="Potential")
        ### Creates the Potential@n0 and Potential@n1 edge model
        edge_from_node_model(device=device, region=region, node_model="Potential")

    if not InEdgeModelList(device, region, "ElectricField"):
        ### Electric field on each edge, as well as its derivatives with respect to the potential at each node
        ### Then Previous ElectricField
        edge_model(device=device, region=region, name="ElectricField",   equation="(Potential@n0 - Potential@n1)*EdgeInverseLength")
        edge_model(device=device, region=region, name="ElectricField:Potential@n0",  equation="EdgeInverseLength")
        edge_model(device=device, region=region, name="ElectricField:Potential@n1",  equation="-EdgeInverseLength")

    if not InElementModelList(device, region, "ElectricField_x"):
        element_from_edge_model(device=device, region=region,edge_model="ElectricField")
        element_from_edge_model(device=device, region=region,edge_model="ElectricField",derivative="Potential")

    DimList=get_parameter(device=device, region=region, name="DimList")

    edge_solution(device=device, region=region, name="PreElectricField")
    element_from_edge_model(device=device, region=region,edge_model="PreElectricField")

    #Create element_edge_solutions to store coefficient so on:
    element_model(device=device, region=region, name="StartCoefficient",equation="InitialCoefficient")
    ElementEdgeNumbers=len(get_element_model_values(device=device, region=region,name="ElementEdgeCouple"))

    for dim in DimList:
        element_solution(device=device, region=region, name="OldPCoefficient{0}".format(dim))
        element_solution(device=device, region=region, name="CoerciveSign{0}".format(dim))
        element_solution(device=device, region=region, name="OldCoerciveSign{0}".format(dim))

        ###get the data numbers in edge_mode and element_model
        #Set the initial value of cofficients and coercive signs
        set_element_values(device=device, region=region, name="OldPCoefficient{0}".format(dim), init_from="StartCoefficient")
        set_element_values(device=device, region=region, name="CoerciveSign{0}".format(dim),    values=[1]*ElementEdgeNumbers)
        set_element_values(device=device, region=region, name="OldCoerciveSign{0}".format(dim), init_from="CoerciveSign{0}".format(dim))
        #Prepare edge_models for iteration of edge_solutions
        VaryFieldDirection="ifelse(ElectricField_{0}>=PreElectricField_{0},1,-1)".format(dim.lower())
        element_model(device=device, region=region, name="VaryFieldDirection{0}".format(dim), equation=VaryFieldDirection)
        SweepDirection="ifelse(ElectricField_{1}==PreElectricField_{1},OldCoerciveSign{0},VaryFieldDirection{0})".format(dim, dim.lower())
        element_model(device=device, region=region, name="SweepDirection{0}".format(dim), equation=SweepDirection)

        # #Modify the coefficient while it large than "1"
        # edge_model(device=device, region=region, name="JudgePCoefficient",  equation="ifelse(PolorizationCoefficient>1,1,PolorizationCoefficient)")
        NumeratorTanh="(OldPCoefficient{0}*({2}(Omega*(PreElectricField_{1}-CoerciveSign{0}*CoerciveField))-CoerciveSign{0})+CoerciveSign{0}-SweepDirection{0})".format(dim, dim.lower(), FunctionModel)
        element_model(device=device, region=region, name="NumeratorTanh{0}".format(dim), equation=NumeratorTanh)  
        DenominatorTanh="({2}(Omega*(PreElectricField_{1}-SweepDirection{0}*CoerciveField))-SweepDirection{0})".format(dim, dim.lower(), FunctionModel)
        element_model(device=device, region=region, name="DenominatorTanh{0}".format(dim), equation=DenominatorTanh)  
        NewPCoefficient="ifelse(abs(DenominatorTanh{0})<2e-16, 1.0, NumeratorTanh{0}/ifelse(abs(DenominatorTanh{0})<2e-16, 1.0, DenominatorTanh{0}))".format(dim, dim.lower())
        element_model(device=device, region=region, name="NewPCoefficient{0}".format(dim), equation=NewPCoefficient)

        ### get the unitx and unity edge models adapted to the element edge.
        # element_model(device=device, region=region, name="ElementUnit{0}", equation="unitx")
        ##Ferroelectric Polarization and its derivatives with respect to nodes
        Polarization="SaturationPolarization*(OldPCoefficient{0}*{2}(Omega*(ElectricField_{1}-CoerciveField*CoerciveSign{0}))+CoerciveSign{0}*(1-OldPCoefficient{0}))".format(dim, dim.lower(), FunctionModel)
        CreateElementModel2d(device, region, "Polarization{0}".format(dim), Polarization)
        CreateElementModelDerivative2d(device, region, "Polarization{0}".format(dim), Polarization, "Potential")
    if get_dimension(device=device)==2 :
        PotentialElementFlux="dot2d(PolarizationX+Permittivity*ElectricField_x, PolarizationY+Permittivity*ElectricField_y, unitx, unity)"
    elif get_dimension(device=device)==3 :
        PotentialElementFlux="dot3d(PolarizationX+Permittivity*ElectricField_x, PolarizationY+Permittivity*ElectricField_y, PolarizationZ+Permittivity*ElectricField_z, unitx, unity, unitz)"
    CreateElementModel2d(device, region, "PotentialElementFlux", PotentialElementFlux)
    CreateElementModelDerivative2d(device, region, "PotentialElementFlux", PotentialElementFlux, "Potential")
    #
    ### Create the bulk equation
    equation(device=device, region=region, name="PotentialEquation", variable_name="Potential", 
             element_model="PotentialElementFlux", variable_update=update_type)

# def CreateSemiFerroInterface(device, interface):
#   '''
#     continuous potential at interface
#   '''
#   model_name = CreateContinuousInterfaceModel(device, interface, "Potential")
#   interface_equation(device=device, interface=interface, name="PotentialEquation", variable_name="Potential", interface_model=model_name, type="continuous")


def CreateFerroContactEquation(device, contact, attached_to=None):
    ### Contact models and equations--edge_charge_model="PotentialElementFlux" ,
    CreateOxidePotentialContact(device, contact, element_contact=True)



def FerroRegionIterate(device, region, ElementChecklist=None):
    Iterate=get_parameter(device=device, region=region, name="Iterate")
    DimList=get_parameter(device=device, region=region, name="DimList")
    CheckElementValues(device, region, ElementChecklist)
    for dim in DimList:
        set_element_values(device=device, region=region, name="OldCoerciveSign{0}".format(dim), init_from="CoerciveSignX")
        if Iterate:
            set_element_values(device=device, region=region, name="OldPCoefficient{0}".format(dim), init_from="NewPCoefficient{0}".format(dim))
        set_element_values(device=device, region=region, name="CoerciveSign{0}".format(dim),    init_from="SweepDirection{0}".format(dim))
    set_edge_values(device=device, region=region, name="PreElectricField", init_from="ElectricField")
    # CheckElementValues(device, region, ElementChecklist)



def FerroCapacitorPlotSweep(device, region,
                SweepContact,
                ChargeContact,
                CurrentContact=None,
                End_bias=1.0,
                step_limit=0.1,
                min_step=0.001, 
                rel_error=1e-8, 
                abs_error=1e30,
                iterations=30,
                SaveAs=None,
                Checklist=None):
    
    ### Create the chart for display
    x = y1 = y2 = []
    plt.ion()  # interactive mode on
    fig=plt.figure(num=SweepContact, figsize=(10, 5))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    Plot1Curve, = ax1.plot(x, y1)  # plot the data and specify the 2d line
    Plot2Curve, = ax1.plot(x, y2)
    ax1.set_xlabel('%s Voltage'% SweepContact)
    ax1.set_ylabel('%s Charge'% ChargeContact)

    # if SweepContact=="gate" : ax1.set_yscale('log')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    text = plt.text(0.5, 0.1, 'Voltage', transform=ax1.transAxes,
                    fontsize=10, verticalalignment='top', bbox=props)
        
    # input("pause")

    Expand=2.0
    step_size=step_limit
    last_bias=get_parameter(device=device, name=GetContactBiasName(SweepContact))
    y1_new = get_contact_charge(device=device, contact=ChargeContact, equation="PotentialEquation")
    StoreData=[last_bias,y1_new]
    y1 = np.append(y1, y1_new)
    x = np.append(x, last_bias)

    if (last_bias < End_bias):
        step_sign=1
    else:
        step_sign=-1
    # last_bias=start_bias
    StepByStep=get_parameter(device=device, region=region, name="StepByStep")
    while(abs(last_bias - End_bias) > min_step):
        if step_size>step_limit :step_size=step_limit
        # input("aa")
        next_bias=last_bias + step_sign * step_size
        if next_bias < End_bias:
            next_step_sign=1
        else:
            next_step_sign=-1
        if next_step_sign != step_sign:
            next_bias=End_bias
        print("Setting The next_bias:%s, End_bias:%s,step_size:%s" % (next_bias,End_bias,step_size,))
        try:
            # set_element_values(device=device, region=region, name="OldCoerciveSignX", values="CoerciveSignX")
            # set_element_values(device=device, region=region, name="OldCoerciveSignY", values="CoerciveSignY")
            set_parameter(device=device, name=GetContactBiasName(SweepContact), value=next_bias)
            solve(type="dc", absolute_error=abs_error, relative_error=rel_error, maximum_iterations=iterations)
            print("Success. Set at bias:%s, last_bias:%s, step_size:%s, step_sign:%s"%(next_bias,last_bias,step_size,step_sign))
            FerroRegionIterate(device, region, ElementChecklist=Checklist)
            y1_new = get_contact_charge(device=device, contact=ChargeContact, equation="PotentialEquation")
            if StepByStep:
                input("bb")
            print(next_bias, y1_new, get_contact_charge(device=device, contact="bot", equation="PotentialEquation"))
            StoreData.append([next_bias,y1_new])
            y1 = np.append(y1, y1_new)
            x = np.append(x, next_bias)
            text.set_text("%s:%.2fV" % (SweepContact,next_bias))
            Plot1Curve.set_xdata(x)
            Plot1Curve.set_ydata(y1)
            ax1.relim()  # renew the data limits
            ax1.autoscale_view(True, True, True)  # rescale plot view
            plt.pause(0.001)
            step_size=step_size*Expand
        except error as msg:
            if str(msg).find("Convergence failure") != 0:
                print(msg,"%s:%.2fV " % (SweepContact,next_bias))
                return "_Fatal"
            set_parameter(device=device, name=GetContactBiasName(SweepContact), value=last_bias)
            step_size *= 1/Expand
            print("Failure at %sV!!!, Setting new step size :%s, "%(next_bias,step_size))
            if step_size < min_step:
                print("Min step size too small at %s:"%min_step)
                CheckElementValues(device, region, Checklist)
                return "_Reduce Interval"
            continue
        last_bias=next_bias
    text.set_text("")
    #export results
    if SaveAs != None:
        SaveAs= "%s_%s %.2fV.csv"%(SaveAs,SweepContact,End_bias)
        StoreData=np.array(StoreData)
        dataframe=pd.DataFrame({"Voltage":x,"contact_charge":y})
        dataframe.to_csv(SaveAs,index=True,sep=",",line_terminator="\r\n")
    return ""



