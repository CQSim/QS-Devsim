import time,math,sys,traceback,os.path
import devsim, matplotlib
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
    
if not InDisplayHost():
    BackGroundPlot()

def CreatePlotFigures(name, SubNums=1, figsize=(10, 5), labels=[["Voltage","Current"], ["Voltage","Current"]]):
    ### Create the chart for display and return the Figure and subfigures object as list
    SubFigures=[]
    plt.ion()  # interactive mode on
    Figure=plt.figure(num=name, figsize=figsize)
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    lable=None #plt.text(0.1, 0.1, 'Voltage', fontsize=10, verticalalignment='bottom', bbox=props)
    for i in range(SubNums):
        ax = Figure.add_subplot(1,SubNums,i+1)
        ax.set_xlabel(labels[i][0])
        ax.set_ylabel(labels[i][1])
        SubFigures.append(ax)
    return Figure,SubFigures
    
def SavePlotFigures(Figure, FileName, device, FileType="svg",):
    fname=GetSaveFileName(device,FileName)
    Figure.savefig("%s.%s"%(fname, FileType))
    Figure.savefig("%s.%s"%(fname, "pdf"))
    # plt.close(Figure)

def CreatePlotCurves(Figure, PlotSubFigure, lists=([1,2],[1,2]), ylogScale=False):
    ##create curves and reture the curve Object
    Curves=[]
    for i in range(int(len(lists)/2)):
        if ylogScale :
            PlotCurve, = PlotSubFigure.plot(lists[i], [abs(po) for po in lists[i+1]])
            PlotSubFigure.set_yscale('log')
        else:
            PlotCurve, = PlotSubFigure.plot(lists[i], lists[i+1])
        Curves.append(PlotCurve)
    PlotSubFigure.relim()  # renew the data limits
    PlotSubFigure.autoscale_view(True, True, True)  # rescale plot view
    return Curves

def RefreshSubPlot(PlotSubFigure,Curves,lists,Title=None,ylogScale=False):
    t0=time.clock()
    PlotSubFigure.set_title(Title)
    for i in range(len(Curves)):
        Curves[i].set_xdata(lists[2*i  ])
        if ylogScale:
            Curves[i].set_ydata([abs(po) for po in lists[2*i+1]])
        else:
            Curves[i].set_ydata(lists[2*i+1])
    PlotSubFigure.relim()  # renew the data limits
    PlotSubFigure.autoscale_view(True, True, True)  # rescale plot view
    if InDisplayHost():
        plt.pause(0.01)
    print("Time consumption of matplot refresh: ",time.clock()-t0)



def PolarizationRamp(device, FerroRegion, abs_error=1e30, rel_error=1e-8, iterations=30):
    print("\n*******Polarization Ramp")
    SaturationPolarization=devsim.get_parameter(device=device, region=FerroRegion, name="SaturationPolarization")
    for i in range(4,9):
        print("\n********Solve %s/8 Polarization"%i)
        devsim.set_parameter(device=device, region=FerroRegion, name="SaturationPolarization",value=SaturationPolarization*i/8)
        try:
            devsim.solve(type="dc", absolute_error=abs_error, relative_error=rel_error, maximum_iterations=iterations)
        except:
            traceback.print_exc()
            return False
        else:
            pass
        finally:
            devsim.set_parameter(device=device, region=FerroRegion, name="SaturationPolarization",value=SaturationPolarization)
    return True


def PlotSweep(device, 
                SweepContact="gate",
                CurrentContact="drain",
                ChargeContacts=[],
                End_bias    = 1.0,
                step_limit  = 0.1,
                min_step    = 0.01, 
                rel_error   = 1e-8, 
                abs_error   = 1e30,
                iterations  = 30,
                FerroRegion = None,
                ElementChecklist= None,
                Lable       = '',
                LowBodyBias = False,
                SaveAll     = False,
                SaveCurrent = True,
                SaveFinal   = True,
                SaveZero    = False,
                SaveCurrentVariation=0,
                SaveRange=[0,0],
                PlotResults = True):
    ylogScale=False
    if SweepContact == "drain":
        CounterContact = "gate"
    elif SweepContact =="gate":
        CounterContact = "drain"
        ylogScale=True
    elif SweepContact  =="top":
        CounterContact  ="bot"

    Counter_bias    = devsim.get_parameter(device=device, name=GetContactBiasName(CounterContact))
    Start_bias      = devsim.get_parameter(device=device, name=GetContactBiasName(SweepContact))

    Msg=""

    if InParameterList(device, "Description"):
        Description=devsim.get_parameter(device=device, name="Description")
    else:
        Description=""

    if InParameterList(device, "SingleCarrier"):
        SingleCarrier=get_parameter(device=device, name="SingleCarrier")
    else:
        SingleCarrier=False
    if SingleCarrier:
        DataLength=2+len(ChargeContacts)
    else:
        DataLength=4+len(ChargeContacts)

    StoreData= np.empty(shape=[0, DataLength])

    FileStr="%s_%.2fV_%s_%.2fV~%.2fV%s"%(CounterContact, Counter_bias, SweepContact, Start_bias, End_bias, Lable)

    if not PlotResults:
        BackGroundPlot()

    Figure,SubFigures=CreatePlotFigures(SweepContact, labels=(("%s Voltage"%SweepContact,"%s Current"%CurrentContact,),))
    if InDisplayHost() and PlotResults:
        lines0=CreatePlotCurves(Figure,SubFigures[0],ylogScale=ylogScale)

    last_bias=Start_bias
    Expand=2.0
    SuccessTimes=0
    step_size=step_limit
    Plot=True
    SaveTail=True
    FileName="%sV" % (FileStr)

    while( SaveTail):
        if Plot :
            text="%s:%.2f$V %s:%.2f$V" % (SweepContact,last_bias,CounterContact,Counter_bias)
            Current_Data=PrintCurrents(device,CurrentContact)
            if CurrentContact in ("source","drain"):
                PrintCurrents(device,"source")
            if ChargeContacts!= []:
                for chargeC in ChargeContacts:
                    Current_Data.append(devsim.get_contact_charge(device=device, contact=chargeC, equation="PotentialEquation") )
            print("Data:",last_bias, Current_Data)
            StoreData=np.concatenate((StoreData,[[last_bias]+Current_Data]))
            Title="%s %s:%.2f$V %s:%.2f$V" % (Description, SweepContact,last_bias,CounterContact,Counter_bias)
            if InDisplayHost() and PlotResults:
                RefreshSubPlot(SubFigures[0], lines0, (StoreData[:,0],StoreData[:,1]), Title=Title, ylogScale=ylogScale)
            DateWeith = StoreData.shape[0]
            if DateWeith >=2 and ( SaveCurrentVariation>1 and abs(math.log(abs(StoreData[DateWeith-1][1]/StoreData[DateWeith-2][1]),SaveCurrentVariation))>=1 \
                    or last_bias >= min(SaveRange) and last_bias <= max(SaveRange) or SaveAll):
                QSSaveDevice(device, file=FileName+".tec", ftype="tecplot")

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
        devsim.set_parameter(device=device, name=GetContactBiasName(SweepContact), value=next_bias)
        if InContactList(device, "body") and SweepContact == "drain" and LowBodyBias:
            lowbias=devsim.get_parameter(device=device, name=GetContactBiasName("drain"))
            devsim.set_parameter(device=device, name=GetContactBiasName("body"), value=min(lowbias,0))
        # CheckElementValues(device, "bulk", ElementChecklist, Length=5)
        try:
            devsim.solve(type="dc", absolute_error=abs_error, relative_error=rel_error, maximum_iterations=iterations)
        except devsim.error as msg:
            SuccessTimes=0
            print(msg,"at %s:%.2fV \t%s:%.2fV \tstep_size:%sV" % (SweepContact,next_bias,CounterContact,Counter_bias,step_size))
            Plot=False
            if FerroRegion != None:
                if devsim.get_parameter(device=device, region=FerroRegion, name="RampPolarization") and step_size<0.02:
                    if PolarizationRamp(device, FerroRegion, abs_error=abs_error, rel_error=rel_error, iterations=iterations):
                        next_bias=last_bias
                        continue
                CheckElementValues(device, FerroRegion, ElementChecklist=ElementChecklist)
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
            if FerroRegion != None:
                FerroRegionIterate(device, FerroRegion, ElementChecklist=ElementChecklist)
            FileName="%s~%.2fV" % (FileStr,next_bias)
            if next_bias==0 and Start_bias!=0 and SaveZero:
                QSSaveDevice(device, file=FileName+".tec", ftype="tecplot")
                # QSSaveDevice(device, file=FileName+".devsim", ftype="devsim")
            last_bias=next_bias


    Final_bias=get_parameter(device=device, name=GetContactBiasName(SweepContact))
        
    if Msg!='':
        FileName="%s%s~%.2fV"%(FileStr, Msg, Final_bias)

    if SaveFinal :
    #and not os.path.isfile('./'+FileName+".dev")
        QSSaveDevice(device, file= FileName+".dev" , ftype="devsim")
#    if SaveFinal and not os.path.isfile('./'+FileName+".tec"):
        QSSaveDevice(device, file= FileName+".tec" , ftype="tecplot")

    titlelist=tuple(["Voltage","totalCurrent","HoleCurrent","ElectronCurrent",]+ChargeContacts)
    print(titlelist)

    if SaveCurrent:
        dataframe=pd.DataFrame({titlelist[i]:StoreData[:,i] for i in range(0,DataLength)})
        dataframe.to_csv("%s.csv"%GetSaveFileName(device,FileName), index=True, sep=",", line_terminator="\r\n")

    if InDisplayHost() and PlotResults:
        RefreshSubPlot(SubFigures[0], lines0, (StoreData[:,0],StoreData[:,1]), Title=Title, ylogScale=ylogScale)
    else:
        CreatePlotCurves(Figure,SubFigures[0], (StoreData[:,0],StoreData[:,1]), ylogScale)
    SavePlotFigures(Figure, FileName, device)

    return Msg
