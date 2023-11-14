import math,sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from devsim import *

from QSsimple_physics import *

def CentiMeterFerroParameters(device, region, 
                        SaturationPolarization  =   7.749e-6, #C/cm^2
                        RemanentPolarization    =   7.130e-6, #C/cm^2
                        CoerciveField           =   5.0e5, #V/cm
                        eps_Ferro               =   12,
                        FunctionModel           =   "tanh"):
    if RemanentPolarization>SaturationPolarization or RemanentPolarization<0:
        raise NameError('RemanentPolarization could not larger than SaturationPolarization or negetive')
    WriteFerroParameters(device, region,  SaturationPolarization,RemanentPolarization,CoerciveField, eps_Ferro, FunctionModel)

def MicroMeterFerroParameters(device, region, 
                        SaturationPolarization  =   7.749e-14, #C/um^2
                        RemanentPolarization    =   7.130e-14, #C/um^2
                        CoerciveField           =   5e1, #V/um
                        eps_Ferro               =   12,
                        FunctionModel           =   "tanh"):
    if RemanentPolarization>SaturationPolarization or RemanentPolarization<0:
        raise NameError('RemanentPolarization could not larger than SaturationPolarization or negetive')
    WriteFerroParameters(device, region,  SaturationPolarization,RemanentPolarization,CoerciveField, eps_Ferro, FunctionModel)

def WriteFerroParameters(device, region,  SaturationPolarization,RemanentPolarization,CoerciveField, eps_Ferro, FunctionModel):
    Omega=math.log((SaturationPolarization+RemanentPolarization)/(SaturationPolarization-RemanentPolarization))/CoerciveField/2.0
    if FunctionModel=="erf":
        InitialCoefficient=1.0/(math.erf(Omega*CoerciveField)+1.0)
    elif FunctionModel=="tanh":
        InitialCoefficient=1.0/(math.tanh(Omega*CoerciveField)+1.0)
    else:
        InitialCoefficient=0.9
    Permittivity = eps_Ferro * get_parameter(name="eps0")
    set_parameter(device=device, region=region, name="FerroModel",              value="Empirical")
    set_parameter(device=device, region=region, name="Permittivity",            value=Permittivity  )
    set_parameter(device=device, region=region, name="SaturationPolarization",  value=SaturationPolarization)
    set_parameter(device=device, region=region, name="RemanentPolarization",    value=RemanentPolarization)
    set_parameter(device=device, region=region, name="CoerciveField",           value=CoerciveField)
    set_parameter(device=device, region=region, name="Omega",                   value=Omega)
    set_parameter(device=device, region=region, name="InitialCoefficient",      value=InitialCoefficient)
    set_parameter(device=device, region=region, name="StepByStep",              value=False)
    set_parameter(device=device, region=region, name="RampPolarization",        value=False)
    CreateDimensionLable(device)



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

    edge_solution(device=device, region=region, name="PreElectricField")


    #Create element_edge_solutions to store coefficient so on:
    if get_dimension(device=device)>1 :
        DimList=get_parameter(device=device, region=region, name="DimList")
        element_model(device=device, region=region, name="StartCoefficient",equation="InitialCoefficient")
        ElementEdgeNumbers=len(get_element_model_values(device=device, region=region,name="ElementEdgeCouple"))
        element_from_edge_model(device=device, region=region,edge_model="PreElectricField")
        if not InElementModelList(device, region, "ElectricField_x"):
            element_from_edge_model(device=device, region=region,edge_model="ElectricField")
            element_from_edge_model(device=device, region=region,edge_model="ElectricField",derivative="Potential")

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
            CreateElementModel2d(device, region, "Polarization_{0}".format(dim), Polarization)
            CreateElementModelDerivative2d(device, region, "Polarization_{0}".format(dim), Polarization, "Potential")
        if get_dimension(device=device)==2 :
            PotentialElementFlux="dot2d(Polarization_x+Permittivity*ElectricField_x, Polarization_y+Permittivity*ElectricField_y, unitx, unity)"
        elif get_dimension(device=device)==3 :
            PotentialElementFlux="dot3d(Polarization_x+Permittivity*ElectricField_x, Polarization_y+Permittivity*ElectricField_y, PolarizationZ+Permittivity*ElectricField_z, unitx, unity, unitz)"
        CreateElementModel2d(device, region, "PotentialElementFlux", PotentialElementFlux)
        CreateElementModelDerivative2d(device, region, "PotentialElementFlux", PotentialElementFlux, "Potential")

        ### Create the bulk equation
        equation(device=device, region=region, name="PotentialEquation", variable_name="Potential", 
                        element_model="PotentialElementFlux", variable_update=update_type)

    elif get_dimension(device=device)==1 :
        edge_model(device=device, region=region, name="StartCoefficient",equation="InitialCoefficient")
        EdgeNumbers=len(get_edge_model_values(device=device, region=region,name="EdgeCouple"))

        edge_solution(device=device, region=region, name="OldPCoefficient")
        edge_solution(device=device, region=region, name="CoerciveSign")
        edge_solution(device=device, region=region, name="OldCoerciveSign")

        ###get the data numbers in edge_mode and edge_model
        #Set the initial value of cofficients and coercive signs
        set_edge_values(device=device, region=region, name="OldPCoefficient", init_from="StartCoefficient")
        set_edge_values(device=device, region=region, name="CoerciveSign",    values=[1]*EdgeNumbers)
        set_edge_values(device=device, region=region, name="OldCoerciveSign", init_from="CoerciveSign")
        #Prepare edge_models for iteration of edge_solutions
        VaryFieldDirection="ifelse(ElectricField>=PreElectricField,1,-1)"
        edge_model(device=device, region=region, name="VaryFieldDirection", equation=VaryFieldDirection)
        SweepDirection="ifelse(ElectricField==PreElectricField,OldCoerciveSign,VaryFieldDirection)"
        edge_model(device=device, region=region, name="SweepDirection", equation=SweepDirection)

        # #Modify the coefficient while it large than "1"
        # edge_model(device=device, region=region, name="JudgePCoefficient",  equation="ifelse(PolorizationCoefficient>1,1,PolorizationCoefficient)")
        NumeratorTanh="(OldPCoefficient*({0}(Omega*(PreElectricField-CoerciveSign*CoerciveField))-CoerciveSign)+CoerciveSign-SweepDirection)".format(FunctionModel)
        edge_model(device=device, region=region, name="NumeratorTanh", equation=NumeratorTanh)  
        DenominatorTanh="({0}(Omega*(PreElectricField-SweepDirection*CoerciveField))-SweepDirection)".format(FunctionModel)
        edge_model(device=device, region=region, name="DenominatorTanh", equation=DenominatorTanh)  
        NewPCoefficient="ifelse(abs(DenominatorTanh)<2e-16, 1.0, NumeratorTanh/ifelse(abs(DenominatorTanh)<2e-16, 1.0, DenominatorTanh))"
        edge_model(device=device, region=region, name="NewPCoefficient", equation=NewPCoefficient)

        ### get the unitx and unity edge models adapted to the edge edge.
        # edge_model(device=device, region=region, name="edgeUnit", equation="unitx")
        ##Ferroelectric Polarization and its derivatives with respect to nodes
        Polarization="SaturationPolarization*(OldPCoefficient*{0}(Omega*(ElectricField-CoerciveField*CoerciveSign))+CoerciveSign*(1-OldPCoefficient))".format(FunctionModel)
        CreateEdgeModel(device, region, "Polarization", Polarization)
        CreateEdgeModelDerivatives(device, region, "Polarization", Polarization, "Potential")
        PotentialEdgeFlux="Polarization+Permittivity*ElectricField"   
        CreateEdgeModel(device, region, "PotentialEdgeFlux", PotentialEdgeFlux)
        CreateEdgeModelDerivatives(device, region, "PotentialEdgeFlux", PotentialEdgeFlux, "Potential")
 
        ### Create the bulk equation
        equation(device=device, region=region, name="PotentialEquation", variable_name="Potential", 
                 edge_model="PotentialEdgeFlux", variable_update=update_type)

# def CreateSemiFerroInterface(device, interface):
#   '''
#     continuous potential at interface
#   '''
#   model_name = CreateContinuousInterfaceModel(device, interface, "Potential")
#   interface_equation(device=device, interface=interface, name="PotentialEquation", variable_name="Potential", interface_model=model_name, type="continuous")


def CreateFerroContactEquation(device, contact, attached_to=None,  is_circuit=False):
    ### Contact models and equations--edge_charge_model="PotentialElementFlux" ,
    if get_dimension(device=device)==1:
        CreateOxidePotentialContact(device, contact, is_circuit=is_circuit)
    elif get_dimension(device=device)>1:
        CreateOxidePotentialContact(device, contact, element_contact=True,  is_circuit=is_circuit)

def FerroRegionIterate(device, region, ElementChecklist=None):

    if not InParameterList(device, "FerroModel", region=region):
        raise NameError('The region of %s is not Ferro'%region)

    if get_parameter(device=device, region=region, name="FerroModel") == "GLFerro":
        CheckElementValues(device, region, ElementChecklist)
        if get_dimension(device=device)>1:
            for dim in get_parameter(device=device, region=region, name="DimList"):
                set_element_values(device=device, region=region, name="PrePolarization_{0}".format(dim), init_from="Polarization_{0}".format(dim))
        elif  get_dimension(device=device)==1:
            set_edge_values(device=device, region=region, name="PrePolarization", init_from="Polarization")

    elif get_parameter(device=device, region=region, name="FerroModel") == "Empirical":
        CheckElementValues(device, region, ElementChecklist)
        if get_dimension(device=device)>1:
            for dim in get_parameter(device=device, region=region, name="DimList"):
                set_element_values(device=device, region=region, name="OldCoerciveSign{0}".format(dim), init_from="CoerciveSign{0}".format(dim))
                set_element_values(device=device, region=region, name="OldPCoefficient{0}".format(dim), init_from="NewPCoefficient{0}".format(dim))
                set_element_values(device=device, region=region, name="CoerciveSign{0}".format(dim),    init_from="SweepDirection{0}".format(dim))
        elif  get_dimension(device=device)==1:
            set_edge_values(device=device, region=region, name="OldCoerciveSign", init_from="CoerciveSign")
            set_edge_values(device=device, region=region, name="OldPCoefficient", init_from="NewPCoefficient")
            set_edge_values(device=device, region=region, name="CoerciveSign",    init_from="SweepDirection")

        set_edge_values(device=device, region=region, name="PreElectricField", init_from="ElectricField")
    # CheckElementValues(device, region, ElementChecklist)


def GinzburgLandauFerroParameters(device, region, 
                        GLFerro_rho = 0, #
                        GLFerro_alpha2 = 0*1e4 , # cm/F to um/F
                        GLFerro_alpha4 = 0*1e20, # cm5/FC2 to um5/FC2
                        ):
    ####Physicis Model form :
    ####Physical Cause and Impact of Negative Capacitance Effect in Ferroelectric P(VDF-TrFE) Gate Stack and Its Application to Landau Transistor
    set_parameter(device=device, region=region, name="GLFerro_rho",     value=GLFerro_rho  )
    set_parameter(device=device, region=region, name="GLFerro_alpha2",  value=GLFerro_alpha2)
    set_parameter(device=device, region=region, name="GLFerro_alpha4",  value=GLFerro_alpha4)
    set_parameter(device=device, region=region, name="FerroModel",      value="GLFerro")
    CreateDimensionLable(device)

def CreateGinzburgLandauFerroRegion(device, region, update_type="default"):
    ####Physicis Model form :
    ####Physical Cause and Impact of Negative Capacitance Effect in Ferroelectric P(VDF-TrFE) Gate Stack and Its Application to Landau Transistor

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

    edge_solution(device=device, region=region, name="PreElectricField")

    if get_dimension(device=device) > 1 :
        if not InElementModelList(device, region, "ElectricField_x"):
            element_from_edge_model(device=device, region=region,edge_model="ElectricField")
            element_from_edge_model(device=device, region=region,edge_model="ElectricField",derivative="Potential")

        DimList=get_parameter(device=device, name="DimList")
        element_solution(device=device, region=region, name="PrePolarization_x")
        element_solution(device=device, region=region, name="PrePolarization_y")
        ElementEdgeNumbers=len(get_element_model_values(device=device, region=region,name="ElementEdgeCouple"))  #get the number of element models
        set_element_values(device=device, region=region, name="PrePolarization_x", values=[0]*ElementEdgeNumbers)
        set_element_values(device=device, region=region, name="PrePolarization_y", values=[0]*ElementEdgeNumbers) #-3.0e-14
        for dim in DimList:
            # input(dim)
            Conterdim= "y" if dim == "x" else "x"
            Polarization="PrePolarization_{0} + tdelta * GLFerro_rho * (ElectricField_{0} - 2*PrePolarization_{0}*GLFerro_alpha2 - 4*PrePolarization_{0}*GLFerro_alpha4*(PrePolarization_{0}^2 + PrePolarization_{1}^2))".format(dim, Conterdim)
            CreateElementModel2d(device, region, "Polarization_{0}".format(dim), Polarization)
            CreateElementModelDerivative2d(device, region, "Polarization_{0}".format(dim), Polarization, "Potential")
        if get_dimension(device=device)==2 :
            PotentialElementFlux="dot2d(Polarization_x+Permittivity*ElectricField_x, Polarization_y+Permittivity*ElectricField_y, unitx, unity)"
        elif get_dimension(device=device)==3 :
            PotentialElementFlux="dot3d(Polarization_x+Permittivity*ElectricField_x, Polarization_y+Permittivity*ElectricField_y, PolarizationZ+Permittivity*ElectricField_z, unitx, unity, unitz)"
        CreateElementModel2d(device, region, "PotentialElementFlux", PotentialElementFlux)
        CreateElementModelDerivative2d(device, region, "PotentialElementFlux", PotentialElementFlux, "Potential")
        ### Create the bulk equation
        equation(device=device, region=region, name="PotentialEquation", variable_name="Potential", 
                        element_model="PotentialElementFlux", variable_update=update_type)

    elif get_dimension(device=device)==1 :
        EdgeNumbers=len(get_edge_model_values(device=device, region=region,name="EdgeCouple"))
        edge_solution(device=device, region=region, name="PrePolarization")
        set_edge_values(device=device, region=region, name="PrePolarization", values=[0]*EdgeNumbers)
        Polarization="PrePolarization + tdelta * GLFerro_rho * (ElectricField - 4*PrePolarization^3*GLFerro_alpha4 - 2*PrePolarization*GLFerro_alpha2)"
        CreateEdgeModel(device, region, "Polarization", Polarization)
        CreateEdgeModelDerivatives(device, region, "Polarization", Polarization, "Potential")
        PotentialEdgeFlux="Polarization+Permittivity*ElectricField"   
        CreateEdgeModel(device, region, "PotentialEdgeFlux", PotentialEdgeFlux)
        CreateEdgeModelDerivatives(device, region, "PotentialEdgeFlux", PotentialEdgeFlux, "Potential")
        ### Create the bulk equation
        equation(device=device, region=region, name="PotentialEquation", variable_name="Potential", 
                 edge_model="PotentialEdgeFlux", variable_update=update_type)

def PolarizationRamp(device, FerroRegion, abs_error=1e30, rel_error=1e-8, iterations=30):
    print("\n*******Polarization Ramp")
    if InParameterList(device, "FerroModel", region):
        DimList=get_parameter(device=device, region=region, name="DimList")
        if get_parameter(device=device, region=region, name="FerroModel")=="GLFerro":
            for dim in DimList:
                PrePolarization=get_parameter(device=device, region=FerroRegion, name="PrePolarization_{0}".format(dim))
            # set_parameter(device=device, region=FerroRegion, name="",value=*i/8)
            try:
                solve(type="dc", absolute_error=abs_error, relative_error=rel_error, maximum_iterations=iterations)
            except:
                traceback.print_exc()
                return False
            else:
                pass
        if get_parameter(device=device, region=region, name="FerroModel")=="Empirical":
            SaturationPolarization=get_parameter(device=device, region=FerroRegion, name="SaturationPolarization")
            for i in range(4,9):
                print("\n********Solve %s/8 Polarization"%i)
                set_parameter(device=device, region=FerroRegion, name="SaturationPolarization",value=SaturationPolarization*i/8)
                try:
                    solve(type="dc", absolute_error=abs_error, relative_error=rel_error, maximum_iterations=iterations)
                except:
                    traceback.print_exc()
                    return False
                else:
                    pass
                finally:
                    set_parameter(device=device, region=FerroRegion, name="SaturationPolarization",value=SaturationPolarization)
    return True

def CapacitorPlotSweep(device, region,
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
                Checklist=None,
                frequency=1e0):
    
    ### Create the chart for display
    x = y1 = y2 = []
    plt.ion()  # interactive mode on
    fig=plt.figure(num=SweepContact, figsize=(10, 5))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    Plot1Curve, = ax1.plot(x, y1)  # plot the data and specify the 2d line
    Plot2Curve, = ax2.plot(x, y2)
    ax1.set_xlabel('%s Voltage'% SweepContact)
    ax1.set_ylabel('%s Charge'% ChargeContact)
    ax2.set_xlabel('%s Voltage'% SweepContact)
    ax2.set_ylabel('%s Capancitance'% ChargeContact)
    # if SweepContact=="gate" : ax1.set_yscale('log')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    text1 = plt.text(0.5, 0.1, 'Voltage', transform=ax1.transAxes,
                    fontsize=10, verticalalignment='top', bbox=props)
    text2 = plt.text(0.5, 0.1, 'Voltage', transform=ax2.transAxes,
                    fontsize=10, verticalalignment='top', bbox=props)
        
    # input("pause")

    Expand=2.0
    step_size=step_limit
    last_bias=get_parameter(device=device, name=GetContactBiasName(SweepContact))
    y1_new = get_contact_charge(device=device, contact=ChargeContact, equation="PotentialEquation")
    y1 = np.append(y1, y1_new)
    x = np.append(x, last_bias)

    solve(type="ac",frequency=frequency)
    y2_new=PrintCapacitor(device, SweepContact, frequency=frequency)
    y2 = np.append(y2, y2_new)
    StoreData=[last_bias,y1_new,y2_new]

    if (last_bias < End_bias):
        step_sign=1
    else:
        step_sign=-1
    # last_bias=start_bias
    if InElementModelList(device, region, "StepByStep"):
        StepByStep=get_parameter(device=device, region=region, name="StepByStep")
    else:
        StepByStep=False
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
        set_parameter(device=device, name=GetContactBiasName(SweepContact), value=next_bias)
        try:
            solve(type="dc", absolute_error=abs_error, relative_error=rel_error, maximum_iterations=iterations)
            print("Success. Set at bias:%s, last_bias:%s, step_size:%s, step_sign:%s"%(next_bias,last_bias,step_size,step_sign))
            print("Parameter_list:", get_parameter_list(device=device, region=region))
            if InParameterList(device, "FerroIterate", region):
                FerroRegionIterate(device, region, ElementChecklist=Checklist)
            y1_new = get_contact_charge(device=device, contact=ChargeContact, equation="PotentialEquation")
            if StepByStep:
                input("bb")
            print(next_bias, y1_new, get_contact_charge(device=device, contact="bot", equation="PotentialEquation"))
            
            solve(type="ac",frequency=frequency)
            y2_new=PrintCapacitor(device, SweepContact, frequency=frequency)
            # input("aa")

            x = np.append(x, next_bias)
            y1 = np.append(y1, y1_new)
            text1.set_text("%s:%.2fV" % (SweepContact,next_bias))
            Plot1Curve.set_xdata(x)
            Plot1Curve.set_ydata(y1)
            ax1.relim()  # renew the data limits
            ax1.autoscale_view(True, True, True)  # rescale plot view

            y2 = np.append(y2, y2_new)
            Plot2Curve.set_xdata(x)
            Plot2Curve.set_ydata(y2)
            ax2.relim()  # renew the data limits
            ax2.autoscale_view(True, True, True)  # rescale plot view
            plt.pause(0.00001)

            StoreData.append([next_bias,y1_new,y2_new])
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
    text1.set_text("")
    #export results
    if SaveAs != None:
        SaveAs= "%s_%s %.2fV.csv"%(SaveAs,SweepContact,End_bias)
        StoreData=np.array(StoreData)
        dataframe=pd.DataFrame({"Voltage":x,"contact_charge":y})
        dataframe.to_csv(SaveAs,index=True,sep=",",line_terminator="\r\n")
    return ""
