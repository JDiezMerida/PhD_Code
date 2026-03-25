# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 15:54:43 2017

Converts the hdf5 data given by Labber to a txt file with all the varible names, units and values. 
The txt file is saved in the same folder as the original hdf5 file.

This version gives the data for all Entries in a single column for each of the variables measured.
In the case of megasweeps (2 variable sweeps), the first column is the slow sweeping variable (outer loop) 
and the second the fast changing one (inner loop). 
The first two rows are the variables and the units. The third row is a long comment with all the details
In case of sweeps with more dimensions, the same logic as megasweeps apply. Slowest variable
will always appear first.

@authors: Aamir Ali & Jaime Diez Mérida @LDQM group at ICFO. 
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import Labber
import numpy as np
import os
import string
import pandas as pd

def hdfTest(s): return '.hdf5' in s and not '.hdf5.lock' in s

def txtThere(s): 
    return s.replace('.hdf5', '.txt') #and os.stat(s.replace('.hdf5', '.txt')).st_size!=0
''' 


'''
def dataExportTxt(filepath, filename):

    f = Labber.LogFile(filepath+'\\'+filename)
  
    logChannels = f.getLogChannels()
    stepChannels = f.getStepChannels()
    settingsValues=f.getChannelValuesAsDict()
    
    
    nEntries = f.getNumberOfEntries()
    nLogChannels = len(logChannels)
    nStepChannels = len(stepChannels)

    
    ###This saves the names and units of all the variables to put as header
    
    names=[]
    units=[]
    if nStepChannels<=2:
        for i in range(len(stepChannels)):
            # Save first the big loop and then the inner loop
            names.append(stepChannels[-1-i]['name'])
            units.append(stepChannels[-1-i]['unit'])
    else:
        names.append(stepChannels[1]['name'])
        names.append(stepChannels[0]['name'])
        units.append(stepChannels[1]['unit'])
        units.append(stepChannels[0]['unit'])
    for i in range(0,len(logChannels)):
        names.append(logChannels[i]['name'])
    for i in range(0,len(logChannels)):
        units.append(logChannels[i]['unit'])
        
  
        
    header_txt=''
    
    ### First line is the actual neame of the columns 
    
    for i in range(0,len(names)):
        header_txt=header_txt+names[i]+'\t'
    header_txt=header_txt+'\n'
    
    
    ### Second line are the units 
    
    for i in range(0,len(names)):
        header_txt=header_txt+units[i]+'\t' 
        
    header_txt=header_txt+'\n'
    
      
    ### Third line is the comment. All the info 
    if nStepChannels==1:
        header_txt=header_txt+'%linesweep\t'  
        scan_type=0
    if nStepChannels==2:
        header_txt=header_txt+'%megasweep\t'  
        scan_type=1
    if nStepChannels==3:
        header_txt=header_txt+'%megasweep\t'  
        if len(stepChannels[2]['values'])==1:
            header_txt=header_txt+stepChannels[2]['name']+str(stepChannels[2]['values'][0])+'\t'
            scan_type=1
        else:
            scan_type=2
   # else:
    #    scan_type=0

    if nStepChannels==4:
        header_txt=header_txt+'%megasweep\t'  
        
        if len(stepChannels[2]['values'])==1:
            scan_type=1
            header_txt=header_txt+stepChannels[2]['name']+str(stepChannels[2]['values'][0])+'\t'
            header_txt=header_txt+stepChannels[3]['name']+str(stepChannels[3]['values'][0])+'\t'
        elif len(stepChannels[3]['values'])==1:
            scan_type=2
            header_txt=header_txt+stepChannels[3]['name']+str(stepChannels[3]['values'][0])+'\t'
        else:
            scan_type=3
    if nStepChannels==5:
        header_txt=header_txt+'%megasweep\t'  
        if len(stepChannels[2]['values'])==1:
            scan_type=1
            for i in range(2,5):
                header_txt=header_txt+stepChannels[i]['name']+str(stepChannels[i]['values'][0])+'\t'
              
        elif len(stepChannels[3]['values'])==1:
            scan_type=2
            for i in range(3,5):
                header_txt=header_txt+stepChannels[i]['name']+str(stepChannels[i]['values'][0])+'\t'
        elif len(stepChannels[4]['values'])==1:
            scan_type=3
            header_txt=header_txt+stepChannels[4]['name']+str(stepChannels[4]['values'][0])+'\t'
        else:
            scan_type=4
    if nStepChannels==6:
        header_txt=header_txt+'%megasweep\t'  
        if len(stepChannels[2]['values'])==1:
            scan_type=1
            for i in range(2,6):
                header_txt=header_txt+stepChannels[i]['name']+str(stepChannels[i]['values'][0])+'\t'
              
        elif len(stepChannels[3]['values'])==1:
            scan_type=2
            for i in range(3,6):
                header_txt=header_txt+stepChannels[i]['name']+str(stepChannels[i]['values'][0])+'\t'
        elif len(stepChannels[4]['values'])==1:
            scan_type=3
            for i in range(4,6):
                header_txt=header_txt+stepChannels[i]['name']+str(stepChannels[i]['values'][0])+'\t'
        else:
            scan_type=4
    
    
    for i in range(len([*settingsValues])):
        header_txt=header_txt+[*settingsValues][i]+':'+str(settingsValues[[*settingsValues][i]])+'\t'
        '''
        if i in [5,10,15,20,25,30]:
            header_txt=header_txt+'\t'
        '''
    header_txt=header_txt+'\n'   
    
    
     
############################# Header ##################################
    
    Header = '"'
    
    Header = Header + 'User: ' + f.getUser()
    #Header = Header + 'Tags: ' + f.getTags()
    Header = Header + '\nProject: ' + f.getProject()
    Header = Header + '\nComment: ' + f.getComment()
    
    Header = Header + '\nNumber of entries: '+ str(nEntries)
    
    Header = Header + '\nStep channels:\n'
    for channel in stepChannels:
        Header = Header + channel['name'] + ', unit: ' + channel['unit'] + ' and with values: ' + str(channel['values']) + '\n'
    
    Header = Header + '\nLog channels:\n'
    for channel in logChannels:
        Header = Header + channel['name'] + '\n'
        
    Header = Header + '"\n\n'
    
############################# Sub-header ##################################
   
    #Header = Header + 'Long name: ' 
    for i in range(nEntries):
        Header = Header + stepChannels[0]['name'] + '\t'
        for j in range(nLogChannels):
            Header = Header + str(logChannels[j]['name']) + '\t'
    
    Header = Header + '\n'
    
    #Header = Header + 'Unit: '
    for i in range(nEntries):
        Header = Header + stepChannels[0]['unit'] + '\t' 
        for j in range(nLogChannels):
            Header = Header + logChannels[j]['unit'] + '\t'
    
    
    #for range  
    #for i in range(nEntry):
    #    for j in 
    #    Header = Header + stepChannels[i]['unit'] + '\t'
    #    for j in range(nLogChannels):
    #        Header = Header + '\t' + logChannels[i]['unit']
    
    
    
############################# Data ##################################
    
    
    if  logChannels[0]['vector'] == False:
        
        ### Record the data. Differentiate between linesweep or megasweep
        ### If only one sweep variable, it is a linesweep
        
        
        if scan_type==0:
            Data = []
            for i in range(nEntries):
                Data.append(stepChannels[0]['values'])  
                for j in range(nLogChannels):
                    Data.append(f.getData(name = j, entry = i))
            Data = np.asarray(Data)
            Data = np.transpose(Data)
            #df=pd.DataFrame(Data,columns=names)
            #df.to_csv(filepath+'\\'+filename[:-4]+'txt', header=names, index=None, sep='\t')
            #np.savetxt(filepath+'\\'+filename[:-5]+'.txt',Data, delimiter='\t')
            
            ### Appended files somehow do not save the information correctly. This will handle the error so 
            ### the rest of the data can still be extracted
            
            try: 
                np.savetxt(filepath+'\\'+filename[:-5]+'.txt',Data, delimiter='\t', header = header_txt, comments='')        
                
            except TypeError as e:        
                print(e)
                print(filename+' is somehow corrupted. Maybe it was appended. ')
        
        elif scan_type==1:
            
            ### Make a Data empty array with the right size to put every
            ### variable in a column of the size of the number of entries
            
            # number of StepChannels is 2
            Data = np.zeros((2+nLogChannels,len(f.getData(0,0))*nEntries))
            
            ### Use counter to record different entries one after each other
            ### Second sweep has the size of the data for each entry. The data from
            ### slow sweep is recorded as single value, now I want to make an array
            ### of the same size of the fast sweeping variable with the data from the
            ### slow sweeping vaiable. Secondsweep is used for that. 
            counter=0
            secondsweep=np.zeros(len(f.getData(0,0))*nEntries)
            
            for i in range(nEntries):
                ### First column is the slow sweep variable
                ### Second column the fast sweep variable
                ### The rest is the data as recorded by Labber 
                try:
                    #Data[nStepChannels-1,counter:counter+len(f.getData(0,0))]=stepChannels[nStepChannels-2]['values']
                    Data[1,counter:counter+len(f.getData(0,0))]=stepChannels[0]['values']
                    secondsweep[counter:counter+len(f.getData(0,0))]=stepChannels[1]['values'][i]
                    Data[0,counter:counter+len(f.getData(0,0))]=secondsweep[counter:counter+len(f.getData(0,0))]
                    
                    for j in range(nLogChannels):
                      Data[j+2,counter:counter+len(f.getData(0,0))]= f.getData(name = j, entry = i)
                except IndexError as e:
                    print(e)
                    print(filename+' is somehow corrupted. Maybe it was appended. ')
                    continue
               
                counter=counter+len(f.getData(0,0))
            Data = np.asarray(Data)
            Data = np.transpose(Data)
            #df=pd.DataFrame(Data,columns=names)
            #df.to_csv(filepath+'\\'+filename[:-4]+'txt', header=names, index=None, sep='\t')
            #np.savetxt(filepath+'\\'+filename[:-5]+'.txt',Data, delimiter='\t')
            np.savetxt(filepath+'\\'+filename[:-5]+'.txt',Data, delimiter='\t', header = header_txt, comments='')

        elif scan_type==2:
            ## IN case of a triple sweep the software will create different files for the different thrid variable scans creating megasweep scans
            thirdvariable_name=stepChannels[2].get('name')
            thirdvariable_name.replace(' ','')
            
            firstvariable_values=stepChannels[0].get('values')
            secondvariable_values=stepChannels[1].get('values')
            thirdvariable_values=stepChannels[2].get('values')

            size_1=len(firstvariable_values)
            size_2=len(secondvariable_values)
            size_3=len(thirdvariable_values)
            
            ## Sweep of the first variable. 
            total_counter=0
            counter_triple=1
            
            
                
            for i in range(size_3):
                
                if size_2*(counter_triple)>nEntries:
                    
                    if nEntries<(total_counter):
                        size_2=nEntries-total_counter
                        print(size_2)
                      
                    
                if total_counter>=nEntries:
                        print('scan broke too early ')
                        break
                ## Now here we can just resolve it as a megasweep
                counter=0
                Data = np.zeros((2+nLogChannels,size_1*size_2))
                secondsweep=np.zeros((size_1*size_2))
                #print(i)
                for j in range(size_2):
                    
                    if total_counter>=nEntries:
                        print('scan broke too early ')
                        break
                    ### First column is the slow sweep variable
                    ### Second column the fast sweep variable
                    ### The rest is the data as recorded by Labber 
                    
                    Data[1,counter:counter+size_1]=firstvariable_values
                    
                    secondsweep[counter:+counter+size_1]=secondvariable_values[j]
                    Data[0,counter:counter+size_1]=secondsweep[counter:counter+size_1]
                    
                    
                    for k in range(nLogChannels):
                        try:
                            Data[k+2,counter:counter+size_1]= f.getData(name = k, entry = j+size_2*(counter_triple-1))
                        except ValueError:
                            print('Is the measurement broke?')
                            break
                    counter=counter+size_1
                    total_counter=total_counter+1
                    
                ### Now convert the data to right format to extract to a text file
                ### np.savetxt seems like the most efficient way to do this, better than pandas
                Data = np.asarray(Data)
                Data = np.transpose(Data)
                #df=pd.DataFrame(Data,columns=names)
                #df.to_csv(filepath+'\\'+filename[:-4]+'txt', header=names, index=None, sep='\t')
                #np.savetxt(filepath+'\\'+filename[:-5]+'.txt',Data, delimiter='\t')
                currentthirdvariablevalue=str(thirdvariable_values[i])
                savecurrentvalue=currentthirdvariablevalue.replace('.','p')
                np.savetxt(filepath+'\\'+filename[:-5]+savecurrentvalue+'.txt',Data, delimiter='\t', header = header_txt, comments='')
                counter_triple=counter_triple+1
                
        elif scan_type==3:
            ## IN case of a cuadruple sweep the software will create different files for the different thrid/fourth variable scans 
            ##creating megasweep scans
            
            thirdvariable_name=stepChannels[2].get('name')
            thirdvariable_name=thirdvariable_name.replace(' ','')
            fourthvariable_name=stepChannels[3].get('name')
            fourthvariable_name=fourthvariable_name.replace(' ','')
            
            firstvariable_values=stepChannels[0].get('values')
            secondvariable_values=stepChannels[1].get('values')
            thirdvariable_values=stepChannels[2].get('values')
            fourthvariable_values=stepChannels[3].get('values')

            size_1=len(firstvariable_values)
            size_2=len(secondvariable_values)
            size_3=len(thirdvariable_values)
            size_4=len(fourthvariable_values)
            
            ## Sweep of the first variable. 
            total_counter=0
            counter_triple=1
            counter_cuadrup=1
            
            for l in range(size_4):
                if size_3*size_2*(counter_cuadrup)>nEntries:
                    size_3=nEntries-total_counter
                    print(size_3)
                    
                if total_counter>=nEntries:
                    print(filename+' broke too early')
                    break
                for i in range(size_3):
                
                    ## Now here we can just resolve it as a megasweep
                    counter=0
                    Data = np.zeros((nStepChannels-1+nLogChannels,size_1*size_2))
                    secondsweep=np.zeros((size_1*size_2))
                    #print(i)
                    for j in range(size_2):
                        
                        if total_counter>=nEntries:
                            print(filename+' broke too early ')
                            break
                        ### First column is the slow sweep variable
                        ### Second column the fast sweep variable
                        ### The rest is the data as recorded by Labber 
                        
                        Data[1,counter:counter+size_1]=firstvariable_values
                        
                        secondsweep[counter:+counter+size_1]=secondvariable_values[j]
                        Data[0,counter:counter+size_1]=secondsweep[counter:counter+size_1]
                        
                        
                        for k in range(nLogChannels):
                            try:
                                Data[k+2,counter:counter+size_1]= f.getData(name = k, entry = j+size_2*(counter_triple-1))
                            except ValueError:
                                print('Is the measurement broken?')
                                break
                        counter=counter+size_1
                        total_counter=total_counter+1
                    
                ### Now convert the data to right format to extract to a text file
                ### np.savetxt seems like the most efficient way to do this, better than pandas
                    Data = np.asarray(Data)
                    Data = np.transpose(Data)
                    #df=pd.DataFrame(Data,columns=names)
                    #df.to_csv(filepath+'\\'+filename[:-4]+'txt', header=names, index=None, sep='\t')
                    #np.savetxt(filepath+'\\'+filename[:-5]+'.txt',Data, delimiter='\t')
                    currentthirdvariablevalue=str(thirdvariable_values[i])
                    savecurrentthirdvalue=currentthirdvariablevalue.replace('.','p')
                    currentfourthvariablevalue=str(fourthvariable_values[l])
                    savecurrentfourthvalue=currentfourthvariablevalue.replace('.','p')
                    
                    np.savetxt(filepath+'\\'+filename[:-5]+'_'+thirdvariable_name+'-'+savecurrentthirdvalue+'_'+fourthvariable_name+'-'+savecurrentfourthvalue+'.txt',Data, delimiter='\t', header = header_txt, comments='')
                    counter_triple=counter_triple+1
                counter_cuadrup=counter_cuadrup+1
                
                
        '''
    
        
        ### Now convert the data to right format to extract to a text file
        ### np.savetxt seems like the most efficient way to do this, better than pandas
        Data = np.asarray(Data)
        Data = np.transpose(Data)
        #df=pd.DataFrame(Data,columns=names)
        #df.to_csv(filepath+'\\'+filename[:-4]+'txt', header=names, index=None, sep='\t')
        #np.savetxt(filepath+'\\'+filename[:-5]+'.txt',Data, delimiter='\t')
        np.savetxt(filepath+'\\'+filename[:-5]+'.txt',Data, delimiter='\t', header = header_txt, comments='')
        '''
    else:
        
        Data = []
        for i in range(nEntries):
            Z=f.getTraceXY(entry=i)
            Data.append(Z[0])
            Data.append(Z[1])
            
        Data = np.asarray(Data)
        Data = np.transpose(Data)
        #df=pd.DataFrame(Data,columns=names)
        #df.to_csv(filepath+'\\'+filename[:-4]+'txt', header=names, index=None, sep='\t')
        #np.savetxt(filepath+'\\'+filename[:-5]+'.txt',Data, delimiter='\t')
        np.savetxt(filepath+'\\'+filename[:-5]+'.txt',Data, delimiter='\t', header = header_txt, comments='')
    return scan_type




       
            
filepath = os.getcwd()
#filepath='//z-sv-l2smb01.ad.physik.uni-muenchen.de/ls-efetov/User/Roop/LMU/Data/Kiutra/TBG36/Data_0312/'
filepath='C:/jdiez_local/LDQM/a. Data/FabPaperData/P-TB55/Original/Fraunhofer/Good_one/'
allFiles=os.listdir(filepath)
exceptFiles=[]

list = filter(hdfTest,allFiles)

counter=0
for s in list:
    sTxt=txtThere(s)
    print(sTxt)
    exceptFiles.append(sTxt)
    if exceptFiles[counter] not in allFiles or os.stat(filepath+'/'+ sTxt).st_size==0:
        dataExportTxt(filepath, s)
    counter=counter+1  
    