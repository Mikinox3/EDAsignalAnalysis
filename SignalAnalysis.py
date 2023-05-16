
from cmath import sqrt
from pickle import TRUE
import numpy as np
import neurokit2 as nk
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import simps

#prepare the raw EDA to calculate parameters later--> filter + calculate SCR, SCL component...
def EdaAnalysis(Edalst,SamplingFreq,doPlots,type):

    #Clean raw EDA and calculate SCL ,SCR, SCR onset, SCR peaks etc...  
    eda_FP,info=nk.eda_process(Edalst,sampling_rate=SamplingFreq) 
    

    #filtering the peaks again 
    currentPeaks = info["SCR_Peaks"]
    currentAmps = eda_FP["SCR_Amplitude"].values[currentPeaks]
    # keep only peaks occuring 1s after stim onset 
    keep1 = currentPeaks >=2000
    # keep only peaks over 0.04 microSiemens & under 1 microSiemens
    keep2 = currentAmps >= 0.04 
    keep3 =currentAmps<=1
  
    
    keeping = keep1*keep2*keep3
  
    for key in info:
        if not isinstance(info[key],int):
            info[key]=info[key][keeping] 
   
    peak_signal = nk.signal.signal_formatpeaks(info, desired_length=len(eda_FP["EDA_Phasic"].values),
                                     peak_indices=info["SCR_Peaks"], other_indices=info["SCR_Recovery"])

    signalsTMP = pd.DataFrame({"EDA_Raw": eda_FP["EDA_Raw"].values, "EDA_Clean": eda_FP["EDA_Clean"].values, 
    "EDA_Tonic": eda_FP["EDA_Tonic"].values, "EDA_Phasic": eda_FP["EDA_Phasic"].values})
    eda_FPNew = pd.concat([signalsTMP, peak_signal], axis=1)
    #eda_FPNew.to_csv('C:/Users\mathi/Documents/Visual code/M1 Projet Rembau/testnew{0}.csv'.format(j)) --> to save img data in csv --> add j in function parameters
  
    #plot the cleaned phasic component of the signal with the detected peaks
    features=[info["SCR_Onsets"],info["SCR_Peaks"],info["SCR_Recovery"]]
    if doPlots:
        plot = nk.events_plot(features,eda_FP["EDA_Phasic"],color=["red","blue","orange"])
        plt.title("EDA phasic")
        plt.ylabel('Conductance (in microS)')
        plt.xlabel('Time (in seconds)')
        plt.show()

    if type=='Img':
    # save data for img with scr peak 
        if eda_FPNew["SCR_Peaks"].values.sum()!=0: # to minimize the qtty of data
            return(eda_FPNew["SCR_Peaks"],eda_FPNew["SCR_Amplitude"],eda_FP["SCR_RiseTime"],eda_FP["SCR_RecoveryTime"])
        else:
            return(0,0,0,0)
    else :
        return(eda_FPNew)
        

# save in one df the peaks parameters of one epoch before analysis
def PrepDataImg(EdaIP,SamplingFreq,doPlots):
    type='Img'
    IPSCRpeaks=list()
    IPSCRamp=list()
    IPSCRrise=list()
    IPSCRrecov=list()
    IP=list()
    
    for j in range(8):
        IP=EdaIP.iat[j,1] #IMGPOS one by one
        resultIP=EdaAnalysis(IP,SamplingFreq,doPlots,type)

        IPSCRpeaks=np.append(IPSCRpeaks,resultIP[0])
        IPSCRamp=np.append(IPSCRamp,resultIP[1])
        IPSCRrise=np.append(IPSCRrise,resultIP[2])
        IPSCRrecov=np.append(IPSCRrecov,resultIP[3])
        #resultIP=sa.EdaAnalysis(IP,SamplingFreq,doPlots,j) # --> create csv files with data if function unlocked in signal analysis.py (only for img)

    IPvalues=pd.DataFrame(list(zip(IPSCRpeaks, IPSCRamp,IPSCRrise,IPSCRrecov)),columns=['SCR_Peaks','SCR_Amplitude','SCR_RiseTime','SCR_RecoveryTime'])
    return(IPvalues)
    #IPvalues.to_csv('C:/Users\mathi/Documents/Visual code/M1 Projet Rembau/IP_values.csv')

#save in one df the peaks + general parameters for one epoch before analysis
def PrepDataFilm(EdaFilmPos,SamplingFreq,doPlots):
    type='Film'
    FP=list()

    for x in EdaFilmPos :
        FP=np.append(FP,x)
    resultP=EdaAnalysis(FP,SamplingFreq,doPlots,type) 
    #resultP.to_csv(r'C:\Users\mathi\Documents\Visual code\M1 Projet Rembau\edaFP.csv') #save in a csv file the data
    return(resultP)
   
#Save in one df the general parameters of the signal for one epoch and analyse the parameters
def ImgAnalysis (EdaIP,SamplingFreq) :
    IP =list()
    IPphasic=list()
    IPClean=list()
    PhasicEda = {}

    for j in range(8):
            IP=EdaIP.iat[j,1] #IMGPOS one by one
            resultIP=EdaAnalysis(IP,SamplingFreq, doPlots=False,type="I")
            IPphasic=np.append(IPphasic,resultIP['EDA_Phasic'])
            IPClean=np.append(IPClean,resultIP['EDA_Clean'])
            PhasicEda["SCR{0}".format(j)]=()
            PhasicEda["SCR{0}".format(j)]=np.append(PhasicEda["SCR{0}".format(j)],resultIP['EDA_Phasic'])

    data_items=PhasicEda.items()
    data_list=list(data_items)
    PhasicEda=pd.DataFrame(data_list)
            
    IPphasic=pd.DataFrame(IPphasic)
    vals={}
    #calculate phasic signal parameters
    vals["Phasic_Mean"] = (np.nansum(IPphasic[0].values)) / len(IPphasic)
    vals["Phasic_Std"] = (np.nanstd(IPphasic[0].values)) 
    vals["Phasic_Var"] = (np.nanvar(IPphasic[0].values)) 

    #calculate spectral power value
    Symp=nk.eda_sympathetic(IPClean,sampling_rate=2000,show=TRUE)
    vals["SympAct"]=Symp["EDA_SympN"]

    #calculate AUC 
    AUC=()
    AUC2=()
    for m in range(8):
        y = abs(np.array(PhasicEda.iat[m,1]))
        y2= np.square(np.array(PhasicEda.iat[m,1]))

        # Compute the area using the composite Simpson's rule.
        area = simps(y, dx=1)
        AUC=np.append(AUC,area)

        area2 = simps(y2, dx=1)
        AUC2=abs(np.append(AUC2,area2))


    AUC = sum(AUC)
    APSC= sum(AUC2)/len(AUC2)
    RMSC = sqrt(APSC)
    vals["AUC"]=AUC.real
    vals["RMSC"]=RMSC.real

    eda_vals = pd.DataFrame.from_dict(vals, orient="index").T.add_prefix("SCR_")
    print(eda_vals)
    return(eda_vals)
    
# analyse peaks features for img + film and analyse general signal parameters for film 
def moreAnalysis (data,what):
    vals = {}
    if isinstance(data, pd.DataFrame):

        if what=="Film" or what=="Img": 

            peaks_cols = [col for col in data.columns if "SCR_Peaks" in col]
            if len(peaks_cols) == 1:
                #calculate the number of peaks
                vals["Peaks_N"] = (
                data[peaks_cols[0]].values.sum()

                )

            amp_cols = [col for col in data.columns if "SCR_Amplitude" in col]
            ### to get an array with the Amplitudes without the 0s ####
            Ampltd =[]
            w = data[amp_cols[0]].values
            i=0
            for x in w :
                if w[i] !=0 :
                    Ampltd=np.append(Ampltd,x)
                i=i+1

            if len(amp_cols) == 1:
                vals["Peaks_Amplitude_Sum"] = (
                np.nansum(data[amp_cols[0]].values)
                )
                vals["Peaks_Amplitude_Mean"] = (
                np.nansum(data[amp_cols[0]].values) / data[peaks_cols[0]].values.sum()
                )
                vals["Peaks_Amplitude_Std"] = (
                np.nanstd(Ampltd) 
                )
                vals["Peaks_Amplitude_Var"]=(
                np.nanvar(Ampltd) 
                )  
                vals["Magnitude"] = (
                np.nansum(data[amp_cols[0]].values) / 8
                )

            riseTime_cols = [col for col in data.columns if "SCR_RiseTime" in col]
            if len(riseTime_cols) == 1:
                vals["Peaks_RiseTime_Sum"] = (data[riseTime_cols[0]].values.sum())

                vals["Peaks_RiseTime_Mean"] = (np.nansum(data[riseTime_cols[0]].values) / data[peaks_cols[0]].values.sum())
            
            if what=="Film" :    
                Phasic_cols = [col for col in data.columns if "EDA_Phasic" in col]
                if len(Phasic_cols) == 1:

                    vals["Phasic_Mean"] = (np.nansum(data[Phasic_cols[0]].values)) / len(data[Phasic_cols[0]].values)
                    vals["Phasic_Std"] = (np.nanstd(data[Phasic_cols[0]].values)) 
                    vals["Phasic_Var"] = (np.nanvar(data[Phasic_cols[0]].values))

                # calculate spectral power
                Symp=nk.eda_sympathetic(data["EDA_Clean"],sampling_rate=2000,show=TRUE)
                vals["SympAct"]=Symp["EDA_SympN"]
                
                #calculate AUC and signal energy parameters
                AUC=()
                AUC2=()
                y = abs(np.array(data["EDA_Phasic"]))
                y2= np.square(np.array(data["EDA_Phasic"]))

                # Compute the area using the composite Simpson's rule.
                area = simps(y, dx=1)
                AUC=abs(np.append(AUC,area))

                area2 = simps(y2, dx=1)
                AUC2=abs(np.append(AUC2,area2))

                AUC= sum(AUC)
                APSC= sum(AUC2)/len(AUC2)
                RMSC = sqrt(APSC)
                vals["AUC"]=AUC.real
                vals["RMSC"]=RMSC.real

            eda_vals = pd.DataFrame.from_dict(vals, orient="index").T.add_prefix("SCR_")
   
                
        if what =="SCL" :       
            Tonic_cols = [col for col in data.columns if "EDA_Tonic" in col]
            if len(Tonic_cols) == 1:
                vals["Tonic_Mean"] = (np.nansum(data[Tonic_cols[0]].values)) / len(data[Tonic_cols[0]].values)

            eda_vals = pd.DataFrame.from_dict(vals, orient="index").T.add_prefix("SCL_")

        print(eda_vals)
        return(eda_vals)


#calculate mezn SCL for rest 1 or rest 2 period + plot the signal
def RestSCL(datax,bounds,SamplingFreq,doPlots) :

    what='SCL'
    EdaRest1={"EDA":[],"Time":[]}
    x=0
    i= (bounds.iat[1,0]-40000) #40000 lines = 20 seconds --> time window for SCL mean calculation
    while i <= int(bounds.iat[1,0]) :
        EdaRest1["EDA"]=np.append( EdaRest1["EDA"],datax.iat[i,1])
        EdaRest1["Time"]=np.append( EdaRest1["Time"],datax.iat[i,0])
        i=i+1

    EdaRest1=pd.DataFrame.from_dict(EdaRest1,orient='index')
    EdaRest1= EdaRest1.T
    R1=list()
    for x in EdaRest1['EDA']:
        R1=np.append(R1,x)
    eda_R1,info=nk.eda_process(R1,sampling_rate=SamplingFreq)

    if doPlots :
        plt.plot(EdaRest1['Time'], eda_R1['EDA_Tonic'], color="blue", label='Raw EDA')
        plt.title("Tonic EDA : Rest 1 sequence, last 20 sec")
        plt.ylabel('Conductance (in microS)')
        plt.xlabel('Time (in seconds)')
        plt.show()
    
    REDAvals=moreAnalysis(eda_R1,what)
    return(pd.DataFrame(REDAvals))




'''
EDA_power=nk.signal_power(eda_FP["EDA_Clean"],frequency_band=[(0.045,0.15)],sampling_rate=SamplingFreq,show=True)
print(EDA_power) #power spectrum values --> doing the same thing as eda_sympathetic but more general function
+ check wavelets with time_frquency function
'''


 