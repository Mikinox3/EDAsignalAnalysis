import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import Cut_the_file as cf
import SignalAnalysis as sa
import csv

########################################
## stage 0 : Data loading & general file details
######################################### 

# path infos
path = 'C:/xxx/xxx/'
fileToLookAt = 'xxxx_biopac.tsv'
doAllInFolder = True
doPlots = False
fileExtension = '.tsv'

# files details
SamplingFreq = 2000 # in Hz
TimeColumn = 0
TriggerMultipliers = [1, 2, 4, 8, 16, 32, 64] #c notre bin to dec
#detectors = Detectors(SamplingFreq)

## Stage 0.2 : list of the files to process

filesToProcess = []
if doAllInFolder:  # all the files 
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(fileExtension):
                filesToProcess.append(file)
else:  # one file only
    filesToProcess.append(path+fileToLookAt)


bigData={}
i=0

## Start : go through selected files 
for aFile in filesToProcess:

    print("Doing File : " + aFile)

    ## Stage 0.3 parameters selection
    _, currentFileExtention = os.path.splitext(aFile)
    if currentFileExtention == ".txt":
        headerSize = 20
        timeInSeconds = False
        #EDAColumn = ?? no idea 1??
        TriggerColumns = [2, 3, 4, 5, 6, 7, 8]
    elif currentFileExtention == ".tsv":
        headerSize = 0 
        timeInSeconds = True
        EDAColumn = 1
        TriggerColumns = [5, 6, 7, 8, 9, 10, 11]
    else:
        headerSize = 0
        timeInSeconds = True
        EDAColumn = 1
        TriggerColumns = [5, 6, 7, 8, 9, 10, 11]

    ########################################
    ## Stage 1 : read the file
    ########################################

    df = pd.read_csv(aFile, sep='\t', header=headerSize)

    ########################################
    ## Stage 2 : triggers extraction & other features 
    ########################################

    print("\t Extracting Triggers")
    triggers = np.sum(df.iloc[:,TriggerColumns] * TriggerMultipliers/5, axis=1)
    triggers = triggers.astype('int')
    EDA = df.iloc[:,EDAColumn]
    times = df.iloc[:, TimeColumn]

    # Time converter : min to sec
    if not timeInSeconds:
        times = times * 60

    if doPlots:
        fig, axs = plt.subplots(2)
        fig.suptitle('Trigger Locations')
        axs[1].plot(times, EDA, color="blue", label='Raw EDA')
        axs[1].legend(loc='upper right')
        axs[1].set_ylabel('Amplitude (in microS)')
        axs[1].set_xlabel('Time (in seconds)')

        axs[0].plot(times,triggers,label='Triggers')
        axs[0].legend(loc='upper right')
        axs[0].set_ylabel('Trigger Idx')
        axs[0].set_xlabel('Time (in seconds)')

        # set useblit True on gtkagg for enhanced performance
        
        plt.show()

    ######################################## 
    ## Stage 3 : Cut the file (rest 1 / img / rest 2 / film) 
    ########################################

    frame = { 'Time': times, 'EDA': EDA,'Triggers' :triggers }
    datax = pd.DataFrame(frame)
    
    #calculate SCL mean for the choosen sequence
    choix = 1 #Rest1
    #choix = 2 #Img seq
    #choix = 3 #Rest2
    #choix = 4 #Film seq

    analysis=cf.Cut_file(datax,choix,doPlots) #cut the signal in 4 parts (and save EDA by valence in img and film part), plot signals with triggers if doPlots

    Eda_by_epochs = (analysis[2]) #dictionnary with all EDA by epochs (img+film)
    

    ########################################
    # Stage 4 : Data Analysis 
    ########################################

    # Stage 4.0 : Rest Analysis

    if choix==1 or choix==3:
        # if choix 1 or 3 selected, calculate mean SCL for rest period
        bounds=(analysis[0],analysis[1])
        bounds=pd.DataFrame(bounds) #idx for the EDA part choosen with choix
        REDAvals=sa.RestSCL(datax,bounds,SamplingFreq,doPlots) #dataframe with meanSCL for the choosen rest period
    else :
        REDAvals=pd.DataFrame()

    # Stage 4.1 :Img Analysis

    ###IMG POS###
    EdaIP=Eda_by_epochs['ImgPos']
    IPvalues=sa.PrepDataImg(EdaIP,SamplingFreq,doPlots)
    if sum(IPvalues["SCR_Peaks"].values) !=0 : #if peaks detected
        IPEDAvals=sa.moreAnalysis(IPvalues,"Img")
    else :
        print("no peaks found here")
        IPEDAvals=pd.DataFrame()

    IPEDAvals2=sa.ImgAnalysis(EdaIP,SamplingFreq)
   
    ###IMG NEU###
    EdaINeu=Eda_by_epochs['ImgNeu']
    INEUvalues=sa.PrepDataImg(EdaINeu,SamplingFreq,doPlots)
    if sum(INEUvalues["SCR_Peaks"].values) !=0 :
        INeuEDAvals=sa.moreAnalysis(INEUvalues,"Img")
    else :
        print("no peaks found here")
        INeuEDAvals=pd.DataFrame()

    INeuEDAvals2=sa.ImgAnalysis(EdaINeu,SamplingFreq)

    ###IMG NEG###
    EdaINeg=Eda_by_epochs['ImgNeg']
    INEGvalues=sa.PrepDataImg(EdaINeg,SamplingFreq,doPlots)
    if sum(INEGvalues["SCR_Peaks"].values) !=0 :
        INegEDAvals=sa.moreAnalysis(INEGvalues,"Img")
    else :
        print("no peaks found here")
        INegEDAvals=pd.DataFrame()

    INegEDAvals2=sa.ImgAnalysis(EdaINeg,SamplingFreq)

    # Stage 4.2 : Film Analysis

    # Prepdata
    EdaFilmPos=Eda_by_epochs['FilmPos']
    EdaFilmNeu=Eda_by_epochs['FilmNeu']
    EdaFilmNeg=Eda_by_epochs['FilmNeg']

    if len(Eda_by_epochs['FilmPos'])>3 : #if film
        print("Film Pos analysis...")
        FP_values=sa.PrepDataFilm(EdaFilmPos,SamplingFreq,doPlots)  
        FPEDAvals=sa.moreAnalysis(FP_values,"Film")

        print("Film Neu analysis...")
        FNeu_values=sa.PrepDataFilm(EdaFilmNeu,SamplingFreq,doPlots)  
        FNeuEDAvals=sa.moreAnalysis(FNeu_values,"Film")

        print("Film Neg analysis...Almost done champ !")
        FNeg_values=sa.PrepDataFilm(EdaFilmNeg,SamplingFreq,doPlots)  
        FNegEDAvals=sa.moreAnalysis(FNeg_values,"Film")
    else:
        FPEDAvals=pd.DataFrame()
        FNegEDAvals=pd.DataFrame()
        FNeuEDAvals=pd.DataFrame()

    ########################################
    # Stage 5 : Data Export 
    ########################################

    if len(IPEDAvals) == 0 : # to avoid errors if no peaks detected
        IPEDAvals = pd.DataFrame({'SCR_Peaks_N':[0],'SCR_Peaks_Amplitude_Sum':[np.NaN],
        'SCR_Peaks_Amplitude_Mean':[np.NaN],"SCR_Peaks_Amplitude_Std":[np.NaN],
        "SCR_Peaks_Amplitude_Var":[np.NaN],"SCR_Magnitude": [np.NaN],
        "SCR_Peaks_RiseTime_Sum":[np.NaN],"SCR_Peaks_RiseTime_Mean":[np.NaN]})
        
    
    IPos=pd.concat([IPEDAvals,IPEDAvals2],axis=1,join='inner')
    print(IPos)


    if len(INeuEDAvals) == 0 : # to avoid errors if no peaks detected
        INeuEDAvals = pd.DataFrame({'SCR_Peaks_N':[0],'SCR_Peaks_Amplitude_Sum':[np.NaN],
        'SCR_Peaks_Amplitude_Mean':[np.NaN],"SCR_Peaks_Amplitude_Std":[np.NaN],
        "SCR_Peaks_Amplitude_Var":[np.NaN],"SCR_Magnitude": [np.NaN],
        "SCR_Peaks_RiseTime_Sum":[np.NaN],"SCR_Peaks_RiseTime_Mean":[np.NaN]})
    INeu=pd.concat([INeuEDAvals,INeuEDAvals2],axis=1,join='inner')
    print(INeu)

    if len(INegEDAvals) == 0 : # to avoid errors if no peaks detected
        INegEDAvals = pd.DataFrame({'SCR_Peaks_N':[0],'SCR_Peaks_Amplitude_Sum':[np.NaN],
        'SCR_Peaks_Amplitude_Mean':[np.NaN],"SCR_Peaks_Amplitude_Std":[np.NaN],
        "SCR_Peaks_Amplitude_Var":[np.NaN],"SCR_Magnitude": [np.NaN],
        "SCR_Peaks_RiseTime_Sum":[np.NaN],"SCR_Peaks_RiseTime_Mean":[np.NaN]})

    INeg=pd.concat([INegEDAvals,INegEDAvals2],axis=1,join='inner')
    print(INeg)

    name=aFile.split()
    n=pd.DataFrame({'File ':[name[len(name)-1]]})
    Total=pd.concat([n,REDAvals,IPos,INeu,INeg,FPEDAvals,FNeuEDAvals,FNegEDAvals]) #df with all parameters

    #save data in a dict
    bigData[i]=Total
    i=i+1
    
bigData_list=list(bigData.values())
print (bigData_list)

columns= ["File","SCL_Tonic_Mean","SCR_Peaks_N","SCR_Peaks_Amplitude_Sum","SCR_Peaks_Amplitude_Mean",
"SCR_Peaks_Amplitude_Std","SCR_Peaks_Amplitude_Var","SCR_Magnitude","SCR_Peaks_RiseTime_Sum",
"SCR_Peaks_RiseTime_Mean","SCR_Phasic_Mean","SCR_Phasic_Std","SCR_Phasic_Var","SCR_SympAct","SCR_AUC","SCR_RMSC"]

with open ('datafinal.csv','w') as csvfile :
    writer = csv.writer(csvfile,dialect=csv.excel)
    writer.writerow(columns)

    for x in bigData_list :
        writer.writerows(x.values.tolist())




