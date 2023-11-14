import time,math,sys,traceback,os.path, matplotlib
import devsim
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from QSsimple_physics import *
from QSmodel_create import *
from QSFerro import FerroRegionIterate



def BackGroundPlot():
    matplotlib.use('Agg')
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_colwidth',100)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)

def ForGroundPlot():
    matplotlib.use('Qt5Agg')
    # matplotlib.get_backend()

# if not InDisplayHost():
#     BackGroundPlot()

def CreatePlotFigures(name, SubFigureNums=1, figsize=(10, 5), labels=[["Voltage","Current"], ["Voltage","Current"]]):
    ### Create the chart for display and return the Figure and subfigures object as list
    SubFigures=[]
    plt.ion()  # interactive mode on
    Figure=plt.figure(num=name, figsize=figsize)
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    lable=None #plt.text(0.1, 0.1, 'Voltage', fontsize=10, verticalalignment='bottom', bbox=props)
    for i in range(SubFigureNums):
        ax = Figure.add_subplot(1,SubFigureNums,i+1)
        ax.set_xlabel(labels[i][0])
        ax.set_ylabel(labels[i][1])
        SubFigures.append(ax)
    return Figure,SubFigures


def CreateMainFigure(name, figsize=(10, 8)):
    Figure=plt.figure(num=name, figsize=figsize)
    # Figure.show()
    return Figure

def CreateSubPlot(Figure, Row, Column, num, axeslabel=["Voltage","Current"],xlog=False,ylog=False, title= None):
    # Figure.get_axes()[0].get_gridspec().ncols
    # Figure.get_axes()[0].get_gridspec().nrows
    # input((Row, Column, num, num))
    for i in Figure.get_axes():
        # print(i.get_subplotspec().get_geometry() ,(Row, Column, num, num))
        # input()
        #if i.get_gridspec().nrows==Row and i.get_gridspec().ncols==Column and 
        if i.get_subplotspec().get_geometry() ==(Row, Column, num-1, num-1):# 
            return i
    SubPlot = Figure.add_subplot(Row, Column, num, title= None)
    SubPlot.set_xlabel(axeslabel[0])
    SubPlot.set_ylabel(axeslabel[1])
    if xlog==True:
        SubPlot.set_xscale('log')
    if ylog==True:
        SubPlot.set_yscale('log')
    return SubPlot

def RefreshSubPlot(SubFigure,Curves,lists,Title=None,ylogScale=False):
    t0=time.time()
    SubFigure.set_title(Title)
    for i in range(int(len(lists)/2)):
        # input(lists[2*i])
        xdata=lists[2*i]
        # Curves[i].set_xdata(lists[2*i])
        if ylogScale:
            ydata=[abs(po) for po in lists[2*i+1]]
            # Curves[i].set_ydata([abs(po) for po in lists[2*i+1]])
        else:
            ydata=lists[2*i+1]
            # Curves[i].set_ydata(lists[2*i+1])
        Curves[i].set_data(xdata,ydata)
    SubFigure.relim()  # renew the data limits
    SubFigure.autoscale_view(True, True, True)  # rescale plot view
    if InDisplayHost():
        plt.pause(0.0001)
    print("Time consumption of matplot refresh: ",time.time()-t0)

def CreatePlotCurves(Figure, SubFigure, lists=([0.1,0.2],[0.1,0.2]), ylogScale=False,labels=[]):
    ##create curves and reture the curve Object
    Curves=[]
    CurveNum=int(len(lists)/2)
    LineNum =max(CurveNum,len(labels))
    for i in range(LineNum):
        if i<CurveNum :
            curve=(lists[i], lists[i+1])
        else:
            curve=(lists[0], lists[1])
        if ylogScale :
            curve=(curve[0], [abs(po)+1e-100 for po in curve[1]])
            SubFigure.set_yscale('log')
        PlotCurve, = SubFigure.plot(np.array(curve[0]),np.array(curve[1]),linestyle=(i%4*2,(10,i%4*2)))
        # print(PlotCurve)
        Curves.append(PlotCurve)
    if labels:
        SubFigure.legend(handles=Curves,labels=labels, loc='best')
    SubFigure.relim()  # renew the data limits
    SubFigure.autoscale_view(True, True, True)  # rescale plot view
    return Curves

def SavePlotFigures(Figure, FileName, device=None, FileTypes=["pdf"]):
    if device!=None:
        fname=GetSaveFileName(device,FileName)
    else:
        fname="Test"
    for t in FileTypes:
        Figure.savefig("%s.%s"%(fname, t))
    # plt.close(Figure)

def RemovePlotFigures(Figure, FileName, device=None, FileTypes=["svg","pdf"],):
    if device!=None:
        fname=GetSaveFileName(device,FileName)
    else:
        fname="Test"
    for t in FileTypes:
        if os.path.exists("%s.%s"%(fname, t)):
            os.remove("%s.%s"%(fname, t))


def UpdateExtraNodeCharge(device, Step_size=0.1, NodeChecklist=None, Length=10):

    if not InParameterList(device, "ExtraNodeChargeRegion"):
        raise NameError("There isn't ExtraNodeCharge Region for iterate")
    region=get_parameter(device=device, name="ExtraNodeChargeRegion")

    if not InParameterList(device, "ExtraNodeChargeList", region):
        raise NameError("There isn't ExtraNodeChargeList in Region of %s for iterate"%region)
    print("update the ExtraNodeChargeDic")

    if InParameterList(device, "SweepSpeed"):
        TimeLag=Step_size/get_parameter(device=device, name="SweepSpeed")
        set_parameter(device=device, name="TimeLag", value=TimeLag)
        # input(TimeLag)
    else:
        TimeLag=0.0
    # input([region,SolutionName,NodeName])
    for SolutionName in get_parameter(device=device, region=region, name="ExtraNodeChargeList"):
        NodeModelName = "%sUpdate"%SolutionName
        SolutionName2 = "%s2"%SolutionName
        # input(NodeChecklist)
        CheckNodeValues(device, region, NodeChecklist=NodeChecklist, Length=Length)
        set_node_values(device=device, region=region, name=SolutionName2, init_from=NodeModelName)
        set_node_values(device=device, region=region, name=SolutionName,  init_from=SolutionName2)
        CheckNodeValues(device, region, NodeChecklist=NodeChecklist, Length=Length)

    # InitialSolve(device, rel_error=1e-15,tdelta=TimeLag)

def VoltagePlotSweep(device, 
                    SweepModel      = ["Current"],
                    solvetype       = "dc",
                    SweepContact    = None,
                    YlogScale       = True,
                    CurrentContacts = [],
                    ChargeContacts  = [],
                    End_bias        = 1.0,
                    step_limit      = 0.1,
                    min_step        = 0.01, 
                    rel_error       = 1e-8, 
                    abs_error       = 1e30,
                    iterations      = 30,
                    FerroRegion     = None,
                    ElementChecklist= None,
                    NodeChecklist   = None,
                    Lable           = '',
                    LowBodyBias     = False,
                    CapacitorContacts    = [],
                    frequency       = [1.0,],
                    VolumeIntegrateList  = [],
                    DeviceMonitorList    = [],
                    SaveAll         = False,
                    SaveData        = True,
                    SaveFinal       = True,
                    SaveZero        = False,
                    SaveCurrentVariation=0,
                    SaveRange       = [0,0],
                    PlotResults     = True,
                    ExtraNodeChargeUpdate= True):
    #### For the single carrier model
    if InParameterList(device, "SingleCarrier"):
        SingleCarrier=get_parameter(device=device, name="SingleCarrier")
    else:
        SingleCarrier=False
    CurrentWidth=0
    if "Current" in SweepModel:
        if SingleCarrier:
            CurrentWidth=1
        else:
            CurrentWidth=3

    #### The Sweep Model decide which kind of result to be expoted: "Current" "Capacitor"
    CurrentTitles=[]
    CapacitorTitles=[]
    NodeChargeTitles=[]
    ChargeContactsTitles=[]
    if "Current" in SweepModel:
        if CurrentContacts==[]:
            raise NameError("Please Specify the Current Contact")
        else:
            for C in CurrentContacts:
                for i in ("Current","HoleCurrent","ElectronCurrent"):
                    CurrentTitles.append("%s_%s"%(C,i))
    if "Capacitor" in SweepModel:
        if CapacitorContacts==[] :
            CapacitorContacts = [SweepContact]
        for C in CapacitorContacts:
            if frequency==None:
                CapacitorTitles.append(C)
            else:
                for Freq in frequency:
                    CapacitorTitles.append("%sCap%0.0e"%(C,Freq))
    for c in ChargeContacts:
        ChargeContactsTitles.append(c)
    if "NodeCharge" in SweepModel:
        for Vl in VolumeIntegrateList:
            NodeChargeTitles.append("%s_%s_%s"%(Vl["region"],Vl["NodeName"],Vl["SampleType"]))

    titlelist=["Voltage"]+CurrentTitles+CapacitorTitles+ChargeContactsTitles+NodeChargeTitles

    DataWidth=len(titlelist)

    print("*******Tiles for store data: ",titlelist)
    StoreData= pd.DataFrame(columns = titlelist)
    # input(StoreData.shape)

    if step_limit<=0:
        raise NameError("Please Specify a positive step_limit")

    if SweepContact :  #and solvetype == "dc"
        if SweepContact == "drain":
            CounterContact = "gate"
        elif SweepContact =="gate":
            CounterContact = "drain"
            # YlogScale=True
        elif SweepContact  =="top":
            CounterContact  ="bot"
            # YlogScale=True
        elif SweepContact  =="bot":
            CounterContact  ="top"
        Start_bias      = devsim.get_parameter(device=device, name=GetContactBiasName(SweepContact))
        Counter_bias    = devsim.get_parameter(device=device, name=GetContactBiasName(CounterContact))
    # if solvetype == "transient_bdf1":
    #     Start_bias      = devsim.circuit_alter(name="V1", value=next_bias)
    Msg=""

    if InParameterList(device, "Description"):
        Description=devsim.get_parameter(device=device, name="Description")
    else:
        Description=""
    FileStr="%s_%.2fV_%s_%.2fV~%.2fV%s"%(CounterContact, Counter_bias, SweepContact, Start_bias, End_bias, Lable)


    ## Prepare the matplotlib for graphic display or export
    if not InDisplayHost()  or not PlotResults:
        BackGroundPlot()
    PlotLables=[]
    SubPlotColumns=max(len(CurrentContacts),len(CapacitorContacts))
    SubPlotRows=len(SweepModel)
    Figure=CreateMainFigure(SweepContact, figsize=(8*SubPlotColumns,6*SubPlotRows))
    Figure.subplots_adjust(left=0.2, right=0.9, top=0.90, bottom=0.15)
    
    if "Current" in SweepModel:
        Current_SubPlots=[]
        Current_Lines=[]
        row=SweepModel.index("Current")
        for i,C in enumerate(CurrentContacts):
            Current_SubPlots.append(CreateSubPlot(Figure, SubPlotRows, SubPlotColumns, row*SubPlotColumns+i+1,
                                                    axeslabel=("%s Voltage(V)"%SweepContact,"%s Current(A)"%C)))
            Current_Lines.append(CreatePlotCurves(Figure,Current_SubPlots[i],ylogScale=YlogScale,
                                                    labels=[CurrentTitles[i*CurrentWidth]]))

    if "Capacitor" in SweepModel:
        Capacitor_SubPlots=[]
        Capacitor_Lines=[]
        if frequency==None: 
            FrequencyLen=1
        else:
            FrequencyLen=len(frequency)
        row=SweepModel.index("Capacitor")
        for i,C in enumerate(CapacitorContacts):
            Capacitor_SubPlots.append(CreateSubPlot(Figure,SubPlotRows,SubPlotColumns,row*SubPlotColumns+i+1,
                                                    axeslabel=("%s Voltage(V)"%SweepContact,"%s Charge(C)"%C)))
            Capacitor_Lines.append(CreatePlotCurves(Figure,Capacitor_SubPlots[i],
                labels= CapacitorTitles[i*FrequencyLen:i*FrequencyLen+FrequencyLen]))

    if "NodeCharge" in SweepModel:
        NodeCharge_SubPlots=[]
        NodeCharge_Lines=[]
        row=SweepModel.index("NodeCharge")
        NodeCharge_SubPlots.append(CreateSubPlot(Figure,SubPlotRows,SubPlotColumns,row*SubPlotColumns+1,
                axeslabel=("%s Voltage(V)"%SweepContact,"Value")))
        NodeCharge_Lines.append(CreatePlotCurves(Figure,NodeCharge_SubPlots[0],ylogScale=True,
                                                    labels=NodeChargeTitles))
    if DeviceMonitorList!=[]:
        PositionX=get_node_model_values(device=device, region=DeviceMonitorList[0]["region"], name="x" )
        MonitorData=pd.DataFrame(map(int,PositionX), columns = ["Position"],)
        MonitorFigure=CreateMainFigure(device, figsize=(8*0.8,5*len(DeviceMonitorList)*0.8))
        MonitorFigure.subplots_adjust(left=0.2, right=0.9, top=0.95, bottom=0.1)
        Monitor_SubPlots=[]
        Monitor_Lines=[]
        for i,m in enumerate(DeviceMonitorList):
            Monitor_SubPlots.append(
                CreateSubPlot(MonitorFigure,len(DeviceMonitorList),1,i+1, axeslabel=("PositionX","Value")))
            Monitor_Lines.append(CreatePlotCurves(MonitorFigure, Monitor_SubPlots[i], ylogScale=m["ylogScale"], labels= m["ModelNames"]))

    last_bias=Start_bias
    Expand=2.0
    SuccessTimes=0
    step_size=step_limit
    Plot=True
    SaveTail=True
    FileName="%sV" % (FileStr)
    if solvetype=="transient_bdf1":
        InitialSolve(device, type="transient_dc", rel_error=rel_error)
    while( SaveTail):
        if Plot :
            # text="%s:%.2f$V %s:%.2f$V" % (SweepContact,last_bias,CounterContact,Counter_bias)
            if SweepContact in ("gate","drain") and "Current" in SweepModel:
                PrintCurrents(device,"source")
            Step_Data=[last_bias]
            if "Current" in SweepModel:
                for C in CurrentContacts:
                    Step_Data=Step_Data+PrintCurrents(device,C)
            if "Capacitor" in SweepModel:
                for C in CapacitorContacts:
                    Step_Data=Step_Data+PrintCapacitor(device, C, frequency=frequency)
            if ChargeContacts!= []:
                for chargeC in ChargeContacts:
                    Step_Data.append(devsim.get_contact_charge(device=device, contact=chargeC, equation="PotentialEquation") )
            if "NodeCharge" in SweepModel:
                for Vl in VolumeIntegrateList:
                    Step_Data.append(NodeModelVolumeIntegraton(device, Vl))
            print("Data:",last_bias, Step_Data)
            StoreData.loc[len(StoreData)] = Step_Data

            # Title="%s %s:%.2f$V %s:%.2f$V" % (Description, SweepContact,last_bias,CounterContact,Counter_bias)
            # if PlotResults and SaveTail:
            if PlotResults:
                plt.figure(SweepContact)
            if "Current" in SweepModel:
                for i,C in enumerate(CurrentContacts):
                    RefreshSubPlot(Current_SubPlots[i], Current_Lines[i], 
                        (list(StoreData["Voltage"]),list(StoreData["%s_Current"%C])),  ylogScale=YlogScale)
            if "Capacitor" in SweepModel:
                for i,C in enumerate(CapacitorContacts):
                    CapacitorLineData=[]
                    if frequency==None:
                        CapacitorLineData.append(list(StoreData["Voltage"]))
                        CapacitorLineData.append(list(StoreData[C]))
                    else:
                        for Freq in frequency:
                            CapacitorLineData.append(list(StoreData["Voltage"]))
                            CapacitorLineData.append(list(StoreData["%sCap%0.0e"%(C,Freq)]))
                    RefreshSubPlot(Capacitor_SubPlots[i], Capacitor_Lines[i], CapacitorLineData)
            if "NodeCharge" in SweepModel:
                NodeChargeLineData=[]
                for i,C in enumerate(VolumeIntegrateList):
                    NodeChargeLineData.append(list(StoreData["Voltage"]))
                    NodeChargeLineData.append(list(StoreData[NodeChargeTitles[i]]))
                RefreshSubPlot(NodeCharge_SubPlots[0], NodeCharge_Lines[0], NodeChargeLineData,  ylogScale=True)
            
            if DeviceMonitorList!=[]:
                plt.figure(device)
                for i,l in enumerate(DeviceMonitorList):
                    MonitorLineData=[]
                    PositionX=get_node_model_values(device=device, region=l["region"], name="x" )
                    for m in l["ModelNames"]:
                        MonitorLineData.append(PositionX)
                        MonitorValue=get_node_model_values(device=device, region=l["region"], name=m)
                        MonitorLineData.append(MonitorValue)
                        MonitorData.insert(MonitorData.shape[1],"%s%s%0.6sV"%(Lable,m,last_bias), MonitorValue)
                    RefreshSubPlot(Monitor_SubPlots[i], Monitor_Lines[i], MonitorLineData,  ylogScale=l["ylogScale"])
                if SaveAll:
                    SavePlotFigures(MonitorFigure, FileName+"monitor", device)

            if not PlotResults or not InDisplayHost() :
                RemovePlotFigures(Figure, FileName, device)
            FileName="%s~%.2fV" % (FileStr,last_bias)
            if not PlotResults or not InDisplayHost() :
                SavePlotFigures(Figure, FileName, device)

            DataLength = len(StoreData)
            # input(StoreData)
            # CurrentlineData=StoreData["%s_Current"%CurrentContacts[0]]
            CurrentlineData=StoreData["Voltage"]
            if DataLength >=2 and ( SaveCurrentVariation>1 and \
                abs(math.log(abs(CurrentlineData[DataLength-1]/CurrentlineData[DataLength-2]),SaveCurrentVariation))>=1 \
                    or last_bias >= min(SaveRange) and last_bias <= max(SaveRange) \
                    or SaveAll):
                QSSaveDevice(device, file=FileName+".tec", ftype="tecplot")
            if last_bias==0 and Start_bias!=0 and SaveZero or last_bias==End_bias:
                # QSSaveDevice(device, file=FileName+".tec", ftype="tecplot")
                QSSaveDevice(device, file=FileName+".dev", ftype="devsim")
                ExportParameters(device, file=FileName)

        if abs(last_bias - End_bias) < min_step:
            # input("last_bias:%s,End_bias:%s"%(last_bias,End_bias)) or Msg in ("_DIVG","_FATAL","_Reduce Interval")
            break

        if step_size>step_limit :
            step_size=step_limit
        next_bias=last_bias + math.copysign(step_size, End_bias-Start_bias)
        if next_bias*last_bias<0 and SaveZero:
            next_bias=0
            ###set the bias to record the zero state
        if Start_bias<End_bias and next_bias>End_bias or Start_bias>End_bias and next_bias<End_bias:
            next_bias=End_bias
            print("********setting to the End_bias %e" % (End_bias))
        print(("\n********%s Last:%eV, \tNext:%eV, \tEnd:%eV. \tStep:%s \tCounterContact:%s %sV") % (SweepContact, last_bias, next_bias, End_bias, step_size, CounterContact, Counter_bias))
        
        if InContactList(device, "body") and SweepContact == "drain" and LowBodyBias:
            lowbias=devsim.get_parameter(device=device, name=GetContactBiasName("drain"))
            devsim.set_parameter(device=device, name=GetContactBiasName("body"), value=min(lowbias,0))
        # CheckElementValues(device, "bulk", ElementChecklist, Length=5)
        #####
        ### Set time inteval for Ferroelectric polarization evole...
        if FerroRegion != None:
            print("*******FerroRegionIterate")
            tdelta=step_size/get_parameter(device=device, name="SweepSpeed")
            set_parameter(device=device, name="tdelta", value=tdelta)
            FerroRegionIterate(device, FerroRegion, ElementChecklist=ElementChecklist)

        if InParameterList(device, "ExtraNodeChargeRegion") and ExtraNodeChargeUpdate:
            UpdateExtraNodeCharge(device,step_size,NodeChecklist)
            
        try:
            if solvetype=="transient_bdf1":
                if SweepContact=="gate":
                    circuit_alter(name="V_%s"%SweepContact, value=next_bias)
                else:
                    devsim.set_parameter(device=device, name=GetContactBiasName(SweepContact), value=next_bias)
                tdelta=step_size/get_parameter(device=device, name="SweepSpeed")
            else:
                # circuit_alter(name="V1", value=next_bias)
                devsim.set_parameter(device=device, name=GetContactBiasName(SweepContact), value=next_bias)
                tdelta=0
            devsim.solve(type=solvetype, absolute_error=abs_error, relative_error=rel_error, charge_error =1e-1,  maximum_iterations=iterations, tdelta=tdelta)
        except devsim.error as msg:
            SuccessTimes=0
            print(msg,"at %s:%.2fV \t%s:%.2fV \tstep_size:%sV" % (SweepContact,next_bias,CounterContact,Counter_bias,step_size))
            Plot=False
            # if FerroRegion != None:
            #     if devsim.get_parameter(device=device, region=FerroRegion, name="RampPolarization") and step_size<0.02:
            #         if PolarizationRamp(device, FerroRegion, abs_error=abs_error, rel_error=rel_error, iterations=iterations):
            #             next_bias=last_bias
            #             continue
            #     CheckElementValues(device, FerroRegion, ElementChecklist=ElementChecklist)
            # input(devsim.get_parameter(device=device, name=GetContactBiasName(SweepContact)))
            if str(msg).find("Convergence failure") != 0:
                Msg= "_FATAL"
                SaveTail=False
                break
            if step_size < min_step:
                print("Min step size too small at %s:"%min_step)
                Msg= "_Reduce_Interval"
                break
            step_size *= 1/Expand
            devsim.set_parameter(device=device, name=GetContactBiasName(SweepContact), value=last_bias)
            print("********Failure!!! Setting new step size :", step_size)
        except:
            devsim.set_parameter(device=device, name=GetContactBiasName(SweepContact), value=last_bias)
            SuccessTimes=0
            print("Unexpected error:")#[0]
            traceback.print_exc()
            Msg= "_UNERR"
            break
        else:
            Plot=True
            SuccessTimes+=1
            if SuccessTimes>=4:
                step_size=step_size*Expand
                SuccessTimes=0
            last_bias=next_bias


    Final_bias=get_parameter(device=device, name=GetContactBiasName(SweepContact))
        
    # if Msg!='':
    FileName="%s~%.2fV%s"%(FileStr, Final_bias, Msg)
    QSSaveDevice(device, file= FileName+".dev" , ftype="devsim")
    QSSaveDevice(device, file= FileName+".tec" , ftype="tecplot")

    if SaveData:
        # dataframe=pd.DataFrame({titlelist[i]:StoreData[:,i] for i in range(0,DataWidth)})
        StoreData.to_csv("%s.csv"%GetSaveFileName(device,FileName), index=True, sep=",", lineterminator="\r\n")
        if DeviceMonitorList!=[]:
            MonitorData.to_csv("%sMonitor.csv"%GetSaveFileName(device,FileName), index=True, sep=",", lineterminator="\r\n")
    SavePlotFigures(Figure, FileName, device)
    ExportParameters(device, file=FileName)
    if DeviceMonitorList!=[]:
        SavePlotFigures(MonitorFigure, FileName+"monitor", device)

    return Msg


def CapacitorSweep(device, region,
                SweepContact,
                ChargeContact,
                CurrentContacts=None,
                End_bias=1.0,
                step_limit=0.1,
                min_step=0.001, 
                rel_error=1e-8, 
                abs_error=1e30,
                iterations=30,
                SaveAs=None,
                Checklist=None,
                frequency=1e10):
    
    ### Create the chart for display
    x = y1 = y2 = []

    # input("pause")

    Expand=2.0
    step_size=step_limit
    last_bias=get_parameter(device=device, name=GetContactBiasName(SweepContact))
    y1_new = get_contact_charge(device=device, contact=ChargeContact, equation="PotentialEquation")
    y1 = np.append(y1, y1_new)
    x = np.append(x, last_bias)

    solve(type="ac",frequency=frequency)
    y2_new=get_circuit_node_value(node="V1.I", solution="ssac_imag")/ (-2*math.pi)
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
            y2_new=get_circuit_node_value(node="V1.I", solution="ssac_imag")/ (-2*math.pi)
            print("The capacitance of {0} at {1}V is: {2}".format(SweepContact, next_bias, y2_new),"\n")
            # input("aa")

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

def PlotTransient(device, 
                VoltateContacts = [],
                CurrentContacts = [],
                ChargeContacts  = [],
                CapacitorContacts= [],
                TimeModel       = "Log",
                gamma           = 1,
                minT_Step       = 1e-15,
                maxTime         = 1e-5,
                rel_error       = 1e-8, 
                abs_error       = 1e30,
                iterations      = 30,
                FerroRegion     = None,
                ElementChecklist= None,
                Lable           = '',
                SaveAll         = False,
                SaveData        = True,
                SaveFinal       = True,
                SaveZero        = False,
                SaveCurrentVariation=0,
                SaveRange       = [0,0],
                PlotResults     = True,
                NodeChecklist =[]
                ):
    if InParameterList(device, "Description"):
        Description=devsim.get_parameter(device=device, name="Description")
    else:
        Description=""
    FileName=""
    # for c,v in VoltateContacts:
    #     set_parameter(device=device, name=GetContactBiasName(c), value=v)
    #     FileName="%s%s%s"%(FileName,c,v)
    #tdelta is the time step
    tdelta=minT_Step
    t=0
    time=[]
    list_g=[]
    list_d=[]
    
    ylogScale=True
    # devsim.solve(type="transient_dc", absolute_error=1.0, relative_error=1e-14, maximum_iterations=3)
    FalseTimes=0
    SubPlotColumns=1
    SubPlotRows=2
    Figure=CreateMainFigure("Transient", figsize=(8*SubPlotColumns*0.8,5*SubPlotRows*0.8))
    Figure.subplots_adjust(left=0.2, right=0.9, top=0.96, bottom=0.1)
    SubPlots=CreateSubPlot(Figure,SubPlotRows, 1, 1, axeslabel=("Time","Current"),xlog=False,ylog=False)
    SubPlots2=CreateSubPlot(Figure,SubPlotRows, 1, 2, axeslabel=("Time","Current"),xlog=False,ylog=True)
    # GLines=CreatePlotCurves(Figure,SubPlots,ylogScale=ylogScale, labels=("I_g","I_d"))
    DLines=CreatePlotCurves(Figure, SubPlots, ylogScale=False, labels=["I_d"])
    DLines2=CreatePlotCurves(Figure, SubPlots2, ylogScale=True, labels=["I_d"])

    # devsim.solve(type="transient_dc", absolute_error=abs_error, relative_error=rel_error, 
    #                     maximum_iterations=iterations, tdelta=tdelta, charge_error=1e10)
    
    circuit_alter(name="V1", param="value", value=-60)
    #doing 10000 time steps
    while True:
        # doing backward euler time integration
        try:
            devsim.solve(type="dc", absolute_error=abs_error, relative_error=rel_error, 
                                    maximum_iterations=iterations, tdelta=tdelta, charge_error=1e10)
        # "dcop" is the correct solution for time integration, as well as dc
        except devsim.error as msg:
            FalseTimes+=1
            print("********Failure!!! Setting new time lag :", tdelta)
            if str(msg).find("Convergence failure") != 0:
                Msg= "_FATAL"
                SaveTail=False
                break
            if FalseTimes >4:
                print("time lag  size too small at %s:"%tdelta)
                Msg= "_Reduce_timelag"
                break
            tdelta=tdelta/5
        else:
            t += tdelta
            # accumulate time for visualization

        # QSSaveDevice(device, file= "%sT%0.0e_%0.0e.tec"%(FileName,t,tdelta) , ftype="tecplot")

        current_drain=PrintCurrents(device,"drain")
        # current_gate=get_circuit_node_value(node="V1.I", solution="dcop")
        print(t,tdelta,current_drain)
 
        time.append(t)
        # list_g.append(current_gate)
        list_d.append(current_drain[0])

        if FerroRegion != None:
            FerroRegionIterate(device, FerroRegion, ElementChecklist=ElementChecklist)
        if InParameterList(device, "ExtraNodeChargeRegion"):
            UpdateExtraNodeCharge(device,tdelta,NodeChecklist)

        # if tdelta/minT_Step==0:
        # plotting the result with matplotlib
        # RefreshSubPlot(SubPlots, DLines, (time,list_g, time,list_d),  ylogScale=ylogScale) #
        RefreshSubPlot(SubPlots, DLines, (time,list_d),  ylogScale=False)
        RefreshSubPlot(SubPlots2, DLines2, (time,list_d),  ylogScale=True)
        if tdelta<10:
            tdelta=tdelta*1.2
        else:
            tdelta=10
        if t>maxTime:
            break
            
    SavePlotFigures(Figure, "Transient", device)
