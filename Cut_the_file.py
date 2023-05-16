
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

###BIOPAC Trigger Codes
# -- CODE IMAGE
# IMGTYPE:-IMGScrambled--IMG--IMGScrambled--COTation
##NEG--71--72--73--74 ## MID--81--82--83--84 ## POS--91--92--93--94
# -- CODE FILM
# FILMTYPE:-AV--FILM--AP--COT
##NEG--101--102--103--104 ## MID--111--112--113--114 ## POS--121--122--123--124
# -- CODE REPOS
##REPOS1--65 ## REPOS2--66/65 ## CALIBRATION--67 ## VALIDATION--68


def getIndexes(dfObj, value):
    ''' Get index position of value in dataframe i.e. dfObj.'''
    listOfPos = list()
    # Get bool dataframe with True at positions where the given value exists
    result = dfObj.isin([value])
    # Get list of columns that contains the value
    seriesObj = result.any()
    columnNames = list(seriesObj[seriesObj == True].index)
    # Iterate over list of columns and fetch the rows indexes where value exists
    for col in columnNames:
        rows = list(result[col][result[col] == True].index)
        for row in rows:
            listOfPos.append((row))
    # Return a list of tuples indicating the positions of value in the dataframe
    return listOfPos

def Cut_file(datax,choix,doPlots) :

    #-----------B. Peak REPOS 1 first & last value detection---------------
    Rest_time = [] #times for peaks 65 or 66
    i= 0
    for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
            if x ==66 or x==65:
                Rest_time = np.append(Rest_time,datax.iat[i,0])
            i=i+1
            

    Rest1_FstIdx=getIndexes(datax,Rest_time[0])
    Rest1_Lastidx = getIndexes(datax,Rest_time[199]) 
    # each trigger 200 lines, 1 line= 0.0005s, 120000 lines = 1 min 


    #-----------C. First Img peak detection (either = 71 /81 or 91)------------

    Img_time = [] #times when img's first trigger
    i= 0
    for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
            if x ==71 or x==81 or x==91:
                Img_time = np.append(Img_time,datax.iat[i,0])
            i=i+1

    Img_FstIdx=getIndexes(datax,Img_time[0]) #idx of the first img's trigger

    #-----------C. 1. Control : presence of unexpected peaks between Repos1 & img1-------

    i = Rest1_Lastidx[0]+1
    PeakError=[]

    while i < Img_FstIdx[0] :
        if datax.iat[i,2]!=0 :
            PeakError = np.append(PeakError,datax.iat[i,2])
        i=i+1

    if len(PeakError)==0:
        print("u good")
    else :
        print("unexpected triggers found : ",PeakError)

    #-----------D. Last img peak detection (either = 74 /84 or 94)---------------------
    LastImg_time = []
    i= 0
    for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
            if x ==74 or x==84 or x==94:
                LastImg_time = np.append(LastImg_time,datax.iat[i,0])
            i=i+1

    Img_LastIdx=getIndexes(datax,LastImg_time[len(LastImg_time)-1]) #idx of the last img's trigger

    #-----------D. 1. Control : presence of unexpected peak values in img sequence----
    i = Img_FstIdx[0]+1
    PeakError=[]

    while i < Img_LastIdx[0] :
        if datax.iat[i,2]!=0 and datax.iat[i,2]!=71 and datax.iat[i,2]!=72 and  datax.iat[i,2]!=73 and  datax.iat[i,2]!=74 and  datax.iat[i,2]!=81 and  datax.iat[i,2]!=82 and  datax.iat[i,2]!=83 and  datax.iat[i,2]!=84 and  datax.iat[i,2]!=91 and datax.iat[i,2]!=92 and datax.iat[i,2]!=93 and datax.iat[i,2]!=94 :
            PeakError = np.append(PeakError,datax.iat[i,2])
        i=i+1

    if len(PeakError)==0:
        print("u good")
    else :
        print("error")
        print(PeakError,"found in the Img serie")

    #-----------E. Peak REPOS 2 first & last value detection--------------------

    if len(Rest_time)>200 :  #if the record stops after img trial we skip film part (rest2 and film abs)
        Rest2_FstIdx=getIndexes(datax,Rest_time[200])
        Rest2_Lastidx=getIndexes(datax,Rest_time[len(Rest_time)-1])

        #-----------F. First Film peak detection (either = 101 /111 or 121)---------------
        Film_time = []
        i= 0
        for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
                if x ==101 or x==111 or x==121:
                    Film_time = np.append(Film_time,datax.iat[i,0])
                i=i+1

        if len(Film_time)!=0 :
            Film_FstIdx=getIndexes(datax,Film_time[0]) #idx of the first film trigger
        
            #-----------F. 1. Control : presence of unexpected peaks between Repos2 & film1--------------------
        
            i=Rest2_Lastidx[0]+1
            PeakError=[]
            
            while i < Film_FstIdx[0] :
                if datax.iat[i,2]!=0 :
                    PeakError = np.append(PeakError,datax.iat[i,2])
                i=i+1

            if len(PeakError)==0:
                print("u good")
            else :
                print("error")


            #-----------G. Last film peak detection (either = 104 /114 or 124)------------------
            LastFilm_time = []
            i= 0
            for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
                    if x ==104 or x==114 or x==124:
                        LastFilm_time = np.append(LastFilm_time,datax.iat[i,0])
                    i=i+1

            Film_LastIdx=getIndexes(datax,LastFilm_time[len(LastFilm_time)-1]) #idx of the last film trigger

            #-----------G. 1. Control : presence of unexpected peak values in film sequence-------
            i = Film_FstIdx[0]+1
            PeakError=[]

            while i < Film_LastIdx[0] :
                if datax.iat[i,2]!=0 and datax.iat[i,2]!=101 and datax.iat[i,2]!=102 and  datax.iat[i,2]!=103 and  datax.iat[i,2]!=104 and  datax.iat[i,2]!=111 and  datax.iat[i,2]!=112 and  datax.iat[i,2]!=113 and  datax.iat[i,2]!=114 and  datax.iat[i,2]!=121 and datax.iat[i,2]!=122 and datax.iat[i,2]!=123 and datax.iat[i,2]!=124 :
                    PeakError = np.append(PeakError,datax.iat[i,2])
                i=i+1

            if len(PeakError)==0:
                print("u good")
            else :
                print("error")
                print(PeakError,"found in the Img serie")
        else :
            print("no film")   
    else :
        print("no film")


    #---------- A1. Legend SAVE EDA BY IMG TYPE---------------------------
    #Pos = dataX's times corresponding to peaks 92 
    #True92 = df with the 8 idx only
    #Endpos = dataX's times corresponding to peaks 93 + 20000 lines
    #True93 = df with the 8 idx only 
    #Image = df containing 2 cols for the 2 x 8 peaks (92 & 93+20000)
    #EDAImgPos = df with 8 cols --> 1 col = Eda for one positive img (from 92 to 93 peak)
    #The same structure is used for negative img (72/73) and neutral img (82/83)


    #######################################################################################################
    #---------- A2. EDA FOR POSITIVE IMG (92 peaks in a dataframe) ----------------
    Pos = []
    i= 0
    for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
            if x ==92:
                Pos = np.append(Pos,datax.iat[i,0])
            i=i+1

    datalist_pos = list
    i=0
    datalist_pos = Pos[0]
    while i < len(Pos)-1 :
        if Pos[i+1]-Pos[i]>1:
            datalist_pos=np.append(datalist_pos,Pos[i+1])
        i=i+1

    dataframe_pos=pd.DataFrame(datalist_pos) #df with times for 92 peaks

    i=0
    true92=()
    for i in range(8):
        true92=np.append(true92,getIndexes(datax,dataframe_pos.iat[i,0]))   
    
    true92=pd.DataFrame(true92) # df with idx when 92 peaks occur

   

    #----------- save the 93 peaks in a df-----------------------
    EndPos = []
    i= 16000 # 93 line + 16000 lines to add 8 s of EDA recording (total : 6 sec of stim + 8 sec for latency)
    #i=0
    for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
            if x ==93:
                EndPos = np.append(EndPos,datax.iat[i,0])
            i=i+1


    datalist_Endpos = list
    i=0
    datalist_Endpos = EndPos[0]
    while i < len(EndPos)-1 :
        if EndPos[i+1]-EndPos[i]>1:
            datalist_Endpos=np.append(datalist_Endpos,EndPos[i+1])
        i=i+1


    dataframe_Endpos=pd.DataFrame(datalist_Endpos) #df with times for 93 peaks + 10sec

    i=0
    true93=()
    for i in range(8):
        true93=np.append(true93,getIndexes(datax,dataframe_Endpos.iat[i,0]))   
    
    true93=pd.DataFrame(true93) # df with 93 peaks +10sec idx

    #----------------Combine dataframes----------------------
    Image = pd.concat([true92,true93],axis=1,join='inner')
    print("92||93 peaks idx : ")
    print(Image)

    #----------------EDA between 92 & 93 peaks (Pos Img time)-------------------
    EdaImgPos={}    
    for x in range(8):
        EdaImgPos["EDA{0}".format(x)]=()
        i= int(Image.iat[x,0])
        while i <= int(Image.iat[x,1]) :
            EdaImgPos["EDA{0}".format(x)]=np.append(EdaImgPos["EDA{0}".format(x)],datax.iat[i,1])
            i=i+1

    data_items=EdaImgPos.items()
    data_list=list(data_items)
    EdaImgPos=pd.DataFrame(data_list)

    print("EDA pour l'époque Img Pos :")
    print(EdaImgPos)

    #################################################################################################

    #---------- A3. EDA FOR NEUTRAL IMG (82 peaks in a dataframe) ----------------
    Neu = []
    i= 0
    for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
            if x ==82:
                Neu = np.append(Neu,datax.iat[i,0])
            i=i+1

    datalist_Neu = list
    i=0
    datalist_Neu = Neu[0]
    while i < len(Neu)-1 :
        if Neu[i+1]-Neu[i]>1:
            datalist_Neu=np.append(datalist_Neu,Neu[i+1])
        i=i+1

    dataframe_Neu=pd.DataFrame(datalist_Neu) #df with times for 82 peaks

    i=0
    true82=()
    for i in range(8):
        true82=np.append(true82,getIndexes(datax,dataframe_Neu.iat[i,0]))   
    
    true82=pd.DataFrame(true82) # df with idx when 82 peaks occur

    #----------- save the 83 peaks in a df-----------------------
    EndNeu = []
    i= 16000
    #i=0
    for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
            if x ==83:
                EndNeu = np.append(EndNeu,datax.iat[i,0])
            i=i+1

    datalist_EndNeu = list
    i=0
    datalist_EndNeu = EndNeu[0]
    while i < len(EndNeu)-1 :
        if EndNeu[i+1]-EndNeu[i]>1:
            datalist_EndNeu=np.append(datalist_EndNeu,EndNeu[i+1])
        i=i+1

    dataframe_EndNeu=pd.DataFrame(datalist_EndNeu) 
    i=0
    true83=()
    for i in range(8):
        true83=np.append(true83,getIndexes(datax,dataframe_EndNeu.iat[i,0]))   
    
    true83=pd.DataFrame(true83) #df with times for 83 peaks+10sec

    #----------------Combine dataframes----------------------
    Image2 = pd.concat([true82,true83],axis=1,join='inner')
    print("82||83 peaks idx : ")
    print(Image2)

    #----------------EDA between 82 & 83 peaks (Neutral Img )-------------------
    EdaImgNeu={}    
    for x in range(8):
        EdaImgNeu["EDA{0}".format(x)]=()
        i= int(Image2.iat[x,0])
        while i <= int(Image2.iat[x,1]) :
            EdaImgNeu["EDA{0}".format(x)]=np.append(EdaImgNeu["EDA{0}".format(x)],datax.iat[i,1])
            i=i+1

    data_items=EdaImgNeu.items()
    data_list=list(data_items)
    EdaImgNeu=pd.DataFrame(data_list)

    print("EDA pour l'époque Img Neu:")
    print(EdaImgNeu)

    #################################################################################################
    #---------- A4. EDA FOR NEGATIVE IMG (72 peaks in a dataframe) ----------------
    Neg = []
    i= 0
    for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
            if x ==72:
                Neg = np.append(Neg,datax.iat[i,0])
            i=i+1

    datalist_Neg = list
    i=0
    datalist_Neg = Neg[0]
    while i < len(Neg)-1 :
        if Neg[i+1]-Neg[i]>1:
            datalist_Neg=np.append(datalist_Neg,Neg[i+1])
        i=i+1

    dataframe_Neg=pd.DataFrame(datalist_Neg) #df with times for 72 peaks

    i=0
    true72=()
    for i in range(8):
        true72=np.append(true72,getIndexes(datax,dataframe_Neg.iat[i,0]))   
    
    true72=pd.DataFrame(true72) #df with idx of the 72 peaks

    #----------- save the 73 peaks in a df-----------------------
    EndNeg = []
    i= 16000
    #i=0
    for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
            if x ==73:
                EndNeg = np.append(EndNeg,datax.iat[i,0])
            i=i+1


    datalist_EndNeg = list
    i=0
    datalist_EndNeg = EndNeg[0]
    while i < len(EndNeg)-1 :
        if EndNeg[i+1]-EndNeg[i]>1:
            datalist_EndNeg=np.append(datalist_EndNeg,EndNeg[i+1])
        i=i+1


    dataframe_EndNeg=pd.DataFrame(datalist_EndNeg) #df with times for 73 peaks+10sec

    i=0
    true73=()
    for i in range(8):
        true73=np.append(true73,getIndexes(datax,dataframe_EndNeg.iat[i,0]))   
    
    true73=pd.DataFrame(true73) #df with idx of the 73 peaks


    #----------------Combine dataframes----------------------
    Image3 = pd.concat([true72,true73],axis=1,join='inner')
    print("72||73 peaks idx:")
    print(Image3)

    #----------------EDA between 72 & 73 peaks (Neg Img )-------------------
    EdaImgNeg={}    
    for x in range(8):
        EdaImgNeg["EDA{0}".format(x)]=()
        i= int(Image3.iat[x,0])
        while i <= int(Image3.iat[x,1]) :
            EdaImgNeg["EDA{0}".format(x)]=np.append(EdaImgNeg["EDA{0}".format(x)],datax.iat[i,1])
            i=i+1

    data_items=EdaImgNeg.items()
    data_list=list(data_items)
    EdaImgNeg=pd.DataFrame(data_list)

    print("EDA de l'époque Img NEg:")
    print(EdaImgNeg)

    #################################################################################################
    #---------- A2. Legend SAVE EDA BY FILM TYPE---------------------------
    #PosF = dataX's times corresponding to peaks 122 
    #True122 = df with the 6 idx only
    #EndposF = dataX's times corresponding to peaks 123 
    #True123 = df with the 6 idx only 
    #Film = df containing 2 cols for the 2 x 6 peaks (122/123)
    #EDAFilmPos = df with 2 cols = Eda for positive film sequence  + time
    #The same structure is used for negative film (112/113) and neutral img (102/103)


    #---------- AA2. EDA FOR POSITIVE FILM (122 peaks in a dataframe) ----------------
    PosF = []
    i= 0
    for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
            if x ==122:
                PosF = np.append(PosF,datax.iat[i,0])
            i=i+1

    if len(PosF)!=0 :
        datalist_posF = list
        i=0
        datalist_posF = PosF[0]

        while i < len(PosF)-1 :
            if PosF[i+1]-PosF[i]>1:
                datalist_posF=np.append(datalist_posF,PosF[i+1])
            i=i+1

        dataframe_posF=pd.DataFrame(datalist_posF) #df with times for 122 peaks
        Peaks122=dataframe_posF[0:3].values.tolist()

        i=0
        true122=()
        for i in range(6):
            true122=np.append(true122,getIndexes(datax,dataframe_posF.iat[i,0]))   
    
        true122=pd.DataFrame(true122) # df with 122 peaks idx 
        
        #----------- save the 123 peak in a df-----------------------
        EndPosF = []
        i= 0
        for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
                if x ==123:
                    EndPosF = np.append(EndPosF,datax.iat[i,0])
                i=i+1

        true123=()
        true123=getIndexes(datax,EndPosF[0])
        true123=pd.DataFrame(true123) # df with 123 peaks idx

        Film = { 'Pstart': true122.iat[0,0], 'Pend': true123[0] }
        Film = pd.DataFrame(Film)
        print("122||123 peaks idx:")
        print(Film)

        #----------------EDA between 122 & 123 peaks (Pos Film time)--------------------

        EdaFilmPos={"EDA":[],"Time":[]}
        x=0
        i= int(Film.iat[x,0])
        while i <= int(Film.iat[x,1]) :
            EdaFilmPos["EDA"]=np.append(EdaFilmPos["EDA"],datax.iat[i,1])
            EdaFilmPos["Time"]=np.append(EdaFilmPos["Time"],datax.iat[i,0])
            i=i+1

        EdaFilmPos=pd.DataFrame.from_dict(EdaFilmPos,orient='index')
        EdaFilmPos=EdaFilmPos.T
        print("EDA de l'époque Film Pos:")
        print(EdaFilmPos)
        
        #plot the Film Pos sequence with the triggers
        if doPlots :
            plt.plot(EdaFilmPos['Time'], EdaFilmPos['EDA'], color="blue", label='Raw EDA')
            plt.axvline(x=Peaks122[0])
            plt.axvline(x=Peaks122[1])
            plt.axvline(x=Peaks122[2])
            plt.title("Raw positive film sequence with triggers")
            plt.ylabel('Conductance(in microS)')
            plt.xlabel('Time (in seconds)')
            plt.show()

        ########################################################################################""
        #---------- AA3. EDA FOR NEUTRAL FILM (112 peaks in a dataframe) ----------------
        NeuF = []
        i= 0
        for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
                if x ==112:
                    NeuF = np.append(NeuF,datax.iat[i,0])
                i=i+1

        datalist_NeuF = list
        i=0
        datalist_NeuF = NeuF[0]
        while i < len(NeuF)-1 :
            if NeuF[i+1]-NeuF[i]>1:
                datalist_NeuF=np.append(datalist_NeuF,NeuF[i+1])
            i=i+1

        dataframe_NeuF=pd.DataFrame(datalist_NeuF) #df with times for 112 peaks
        Peaks112=dataframe_NeuF[0:3].values.tolist()

        i=0
        true112=()
        for i in range(6):
            true112=np.append(true112,getIndexes(datax,dataframe_NeuF.iat[i,0]))   
    
        true112=pd.DataFrame(true112) # df with 112 peaks idx 

        #----------- save the 113 peaks in a df-----------------------
        EndNeuF = []
        i= 0
        for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
                if x ==113:
                    EndNeuF = np.append(EndNeuF,datax.iat[i,0])
                i=i+1

        true113=()
        true113=getIndexes(datax,EndNeuF[0])
        true113=pd.DataFrame(true113) # df with 113 peaks idx
        
        Film2 = { 'Pstart': true112.iat[0,0], 'Pend': true113[0] }
        Film2 = pd.DataFrame(Film2)
        print("112||113 peaks idx:")
        print(Film2)

        #----------------EDA between 112 & 113 peaks (Neu Film time)--------------------

      
        EdaFilmNeu={"EDA":[],"Time":[]}
        x=0
        i= int(Film2.iat[x,0])
        while i <= int(Film2.iat[x,1]) :
            EdaFilmNeu["EDA"]=np.append(EdaFilmNeu["EDA"],datax.iat[i,1])
            EdaFilmNeu["Time"]=np.append(EdaFilmNeu["Time"],datax.iat[i,0])
            i=i+1

        EdaFilmNeu=pd.DataFrame.from_dict(EdaFilmNeu,orient='index')
        EdaFilmNeu=EdaFilmNeu.T
        print("EDA de l'époque Film Neu:")
        print(EdaFilmNeu)

        #plot the Film Neu sequence with the triggers
        if doPlots :
            plt.plot(EdaFilmNeu['Time'], EdaFilmNeu['EDA'], color="blue", label='Raw EDA')
            plt.axvline(x=Peaks112[0])
            plt.axvline(x=Peaks112[1])
            plt.axvline(x=Peaks112[2])
            plt.title("Raw Neutral film sequence with triggers")
            plt.ylabel('Conductance(in microS)')
            plt.xlabel('Time (in seconds)')
            plt.show()

        ###############################################################
        #---------- AA4. EDA FOR NEGATIVE FILM (102 peaks in a dataframe) ----------------
        NegF = []
        i= 0
        for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
                if x ==102:
                    NegF = np.append(NegF,datax.iat[i,0])
                i=i+1

        datalist_NegF = list
        i=0
        datalist_NegF = NegF[0]
        while i < len(NegF)-1 :
            if NegF[i+1]-NegF[i]>1:
                datalist_NegF=np.append(datalist_NegF,NegF[i+1])
            i=i+1
        dataframe_NegF=pd.DataFrame(datalist_NegF) #df with times for 102 peaks
        Peaks102=dataframe_NegF[0:3].values.tolist()

        i=0
        true102=()
        for i in range(6):
            true102=np.append(true102,getIndexes(datax,dataframe_NegF.iat[i,0]))   
    
        true102=pd.DataFrame(true102) # df with 102 peaks idx 

        #----------- save the 103 peaks in a df-----------------------
        EndNegF = []
        i= 0
        for x in datax["Triggers"] : # x prend successivement chacune des valeurs de datax trig
                if x ==103:
                    EndNegF = np.append(EndNegF,datax.iat[i,0])
                i=i+1

        true103=()
        true103=getIndexes(datax,EndNegF[0])
        true103=pd.DataFrame(true103) # df with 103 peaks idx

        Film3 = { 'Pstart': true102.iat[0,0], 'Pend': true103[0] }
        Film3 = pd.DataFrame(Film3)
        print("102||103 peaks idx:")
        print(Film3)

        #----------------EDA between 102 & 103 peaks (Neg Film time)--------------------

        EdaFilmNeg={"EDA":[],"Time":[]}
        x=0
        i= int(Film3.iat[x,0])
        while i <= int(Film3.iat[x,1]) :
            EdaFilmNeg["EDA"]=np.append(EdaFilmNeg["EDA"],datax.iat[i,1])
            EdaFilmNeg["Time"]=np.append(EdaFilmNeg["Time"],datax.iat[i,0])
            i=i+1

        EdaFilmNeg=pd.DataFrame.from_dict(EdaFilmNeg,orient='index')
        EdaFilmNeg=EdaFilmNeg.T
        print("EDA de l'époque Film Neg:")
        print(EdaFilmNeg)
        
        #plot the Film Neg sequence with the triggers
        if doPlots :
            plt.plot(EdaFilmNeg['Time'], EdaFilmNeg['EDA'], color="blue", label='Raw EDA')
            plt.axvline(x=Peaks102[0])
            plt.axvline(x=Peaks102[1])
            plt.axvline(x=Peaks102[2])
            plt.title("Raw Negative film sequence with triggers")
            plt.ylabel('Conductance (in microS)')
            plt.xlabel('Time (in seconds)')
            plt.show()
        
    else :
        print("no film")
        Film=0
        Film2=0
        Film3=0
        EdaFilmPos={"EDA":[],"Time":[]}
        EdaFilmNeu={"EDA":[],"Time":[]}
        EdaFilmNeg={"EDA":[],"Time":[]}
        EdaFilmPos["EDA"]=[0,1]
        EdaFilmNeu["EDA"]=[0,1]
        EdaFilmNeg["EDA"]=[0,1]


    #creating a dict containing EDA for img by epoch + film
    EDA_by_Epochs={'ImgPos':EdaImgPos,'ImgNeu':EdaImgNeu,'ImgNeg':EdaImgNeg,'FilmPos':EdaFilmPos["EDA"],'FilmNeu':EdaFilmNeu["EDA"],'FilmNeg':EdaFilmNeg["EDA"]}
    
    #---------Plot Raw Img seq with the triggers---------
    if doPlots:
        plt.plot(datax['Time'][Img_FstIdx[0]:Img_LastIdx[0]], datax['EDA'][Img_FstIdx[0]:Img_LastIdx[0]], color="Red", label='Raw EDA')
        for m in datalist_pos:
            plt.axvline(x=m,color="green")
        for m in datalist_Neu :
            plt.axvline(x=m,color="orange")
        for m in datalist_Neg:
            plt.axvline(x=m)
        plt.title("Raw Img sequence with triggers")
        plt.ylabel('Conductance (in microS)')
        plt.xlabel('Time (in seconds)')
        plt.legend(loc='upper right')
        plt.show()

    ######################################
    #-------------- BOUNDS FOR THE SCL MEAN------------------------------
    if choix == 1 :
        bnd1 = Rest_time[0] #in seconds
        bnd2 = Img_time[0]
        print("choix 1 : repos")

    if choix == 2 : 
        bnd1 = Img_time[0]
        bnd2 = LastImg_time[len(LastImg_time)-1]
        print("choix 2 : image")

    if choix == 3 :
        if len(PosF)!=0 :
            bnd1 = Rest_time[200]
            bnd2 = bnd1 = Film_time[0]
        else:
            bnd1=0
            bnd2=1   
        print("choix 3 : repos 2")

    if choix == 4 : 
        if len(PosF)!=0 :
            bnd1 = Film_time[0]
            bnd2 = LastFilm_time[len(LastFilm_time)-1]
        else:
            bnd1=0
            bnd2=1
        print ("choix 4 : film")

    bnd1=getIndexes(datax,bnd1)
    bnd2=getIndexes(datax,bnd2)
    
    return(bnd1,bnd2,EDA_by_Epochs)

