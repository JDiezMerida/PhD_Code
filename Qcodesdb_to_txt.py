#%% 
'''
Basic code to export the data from the Qcodes database into txt files.
This is useful to be able to plot the data with the regular plotting code,
without having to load the database every time.
The code is written in a way that it checks if the txt file already exists, and if it does, it does not export it again. This is useful to avoid overwriting the files and to save time. 
The code also checks if the scan is broken, and if it is, it exports only the good data, excluding the broken part. This is useful to avoid plotting broken data and to save time.
The code also writes a log file with the scan id, the name of the file, and the status of the scan (good, broken, etc). This is useful to keep track of the scans and to know which ones are good and which ones are not.
Ju1 was the name of the sample used as an example here, but it can be easily modified to any other dataset. 
'''

import numpy as np 
import matplotlib.pyplot as plt 
import os
from pyparsing import line

#from qcodes import load_by_id, Parameter, ScaledParameter
#from qcodes.dataset.plotting import plot_by_id
#from qcodes.dataset.data_export import get_data_by_id, load_by_id, reshape_2D_data
from qcodes.dataset.data_set import DataSet
#from qcodes.interactive_widget import experiments_widget



class Object():
    pass

def fileExistTest(s): 
    return '.txt' in s 

### the double gate map is ID 175 
def QcodesExportTxt(dbpath,scan_id,savepath,savename):
    print(scan_id)
    writelog=str(scan_id)+'\t'
    allFiles=os.listdir(savepath)
    #print(allFiles)
    filename='ID'+str(scan_id)+'_'+savename+'.txt'
    if filename not in allFiles:# or os.stat(filename).st_size==0:
        try:

            tb32_dataset=DataSet(dbpath,scan_id)
            '''
            The data is nested as dictionaries inside of a dictionary 
            so it is a double dictionary sort to say 
            first get the data
            '''
            double_gate=tb32_dataset.get_parameter_data()
            '''
            now write the names of the variables, but this will not 
            record the measurement sets 
            '''
            meas_variables=[] 
            set_variables=[]
            variables=[]
            units=[]
            for i in double_gate.keys():
                meas_variables.append(str(i))
            '''
            now access one of the dictionaries to get the set variables
            and make a new parameter reordering these variables 
            For the set variables in a megasweep the outer variable
            comes first and the inner variable comes second, so I reverse
            them below 
            '''
            for i in double_gate[meas_variables[0]].keys():
                if str(i) not in meas_variables:
                    set_variables.append(str(i))
                    variables.append(str(i))
            
            #set_variables.reverse()
            '''In the megasweeps, the variables are inverted in Qcodes
                Outer variables goes second and inner goes first, so I invert them'''
            if double_gate[meas_variables[0]][set_variables[0]][0]==double_gate[meas_variables[0]][set_variables[0]][1]:
                pass

            else:
                variables.reverse()
                set_variables.reverse()

            '''Put the type of scan. Not very useful, just for consistency'''
            if len(set_variables)==1:

                type_of_scan= '% linesweep'
            elif len(set_variables)==2:

                type_of_scan= '% megasweep '
            '''
            Create a new variables to store the names but with 
            the swept variables as the 1st and 2nd column
            '''
            for i in meas_variables:
                variables.append(i)

            '''Save the units to put in the txt file'''
            for i in range(len(variables)):
                units.append(tb32_dataset.get_parameters()[i].unit)

            '''Arrange the data, here it is already all done'''
            data_dg=[]  # regular data, just a line
            real_data=[]    # data in the form of a matrix, used to know if the scan broke
            good_data_dg=[] # regular data, excluding the broken last scan
            count=0
            for i in set_variables:
                extract_data=double_gate[meas_variables[0]][i]
                data_dg.append(extract_data)
                real_data.append([])
                good_data_dg.append([])
                count=count+1
            for i in meas_variables:
                extract_data=double_gate[i][i]
                data_dg.append(extract_data) 
                real_data.append([])
                good_data_dg.append([])
            
            ### MISSING SOMETHING TO FIX BROKEN DATA### 

            '''
            Now I want to have the data such that it looks as in a map 
            The first variable should give you the shape of the data
            to know where the sweeps break 
            data is going up 
            '''
            ### THIS WAS ACTUALLY NOT NEEDED!! THIS WILL BE DONE 
            ### AFTER LOADING THE TXT FOR PLOTTING, ANYWAY I CAN KEEP
            ### IT FOR NOW, OTHERWISE JUST COMMENT IT OUT ###
            
            first_variable=data_dg[0]
            first_cut=1 # just to avoid any problem 
            for i in range(len(first_variable)):
                if first_variable[i]!=first_variable[i-1] and i!=0:
                    first_cut=i
                    break
            sweepcount=  int(len(first_variable)/first_cut) 
            if len(first_variable)/first_cut>sweepcount:
                sweepcount= sweepcount   +1    # Only add one if it is broken
            
            ### REAL DATA LIKE THIS ALREADY GETS RID OF THE BROKEN SCAN###
            for j in range(len(variables)):
                for i in range(sweepcount):
                    real_data[j].append(data_dg[j][first_cut*(i):(i+1)*first_cut])
                        #break 
            
            ### DETECT IF THE SCAN IS BROKEN  ###

            if len(set_variables)==2:
                total_data_len=len(data_dg[0])
                first_data_len=len(real_data[0][0])
                last_scan_len=len(real_data[0][-1])
                print(total_data_len)
                print(last_scan_len)
                if last_scan_len!=first_data_len:
                    print('Scan with id '+str(scan_id)+' broke too early')
                    writelog=writelog+'Scan broke too early \t'
                    for i in range(len(variables)):
                        good_data_dg[i]=data_dg[i][0:total_data_len-last_scan_len]
                else:
                    writelog=writelog+'Scan is a map and good \t'
                    for i in range(len(variables)):
                        good_data_dg[i]=data_dg[i]

            else: 
                writelog=writelog+'Scan is a line and good \t'
                good_data_dg=data_dg
            
            ### THIS JUST MAKES THE APPENDED DATA INTO A NUMPY ARRAY TO SAVE EASIER###
            good_data_dg=np.asanyarray(good_data_dg)
            good_data_dg=np.transpose(good_data_dg)

            ### THE HEADER SHOULD MATCH THE ONE OF THE ICFO CODE ###
            ### THIS ENSURES THAT I CAN JUST PLOT IT WITH THE REGULAR CODE ###
            header_txt=''
            for i in range(0,len(variables)):
                    header_txt=header_txt+variables[i]+'\t' 
                    
            header_txt=header_txt+'\n'
            for i in range(0,len(variables)):
                    header_txt=header_txt+units[i]+'\t' 
            header_txt=header_txt+'\n'
            header_txt=header_txt+type_of_scan+'\n'      
            '''
            Now make it into a textfile. 
            Do it only if the file does not already exist
            '''    
            #filename=filenames[scan]
            
            np.savetxt(savepath+filename,
                    good_data_dg, delimiter='\t', header = header_txt, comments='')
        except IndexError as e:
            writelog=writelog+'Scan is broken'
            print(filename+' broken' + ' error: '+str(e))

        with open(savepath+'convert_log.txt','a') as f:
                f.write(writelog+'\n')
    else:
        print(filename+' already exists.')
    
    
        
    return 

# %%
## Round 1 of loading, begining to 20/06
'''
Choose the folder and the database name that you want to open
'''
general_path="C:/jdiez_local/LDQM/a. Data/IST_Austria/"
path_to_db=general_path+"Ju1.db"

'''
Write the items which you want to load, you could translate all
but probably many of the scans are not actually important, so 
there is no need to load them 
'''


#items_to_load1=np.arange(1,485)
savefolder='Ju1/'
savepath=general_path+savefolder

# Take the filenames from the already written ones 
with open(general_path+'Ju1_scan_names.txt') as f:
    lines = f.readlines()
filesList=[i.strip() for i in lines] 
filesList=[i.split('\t') for i in filesList] 
filenames1=[]
for i in range(len(filesList)):
    filenames1.append(filesList[i][1])
items_to_load1=[]
for i in range(len(filesList)):
    items_to_load1.append(int(filesList[i][0]))
#print(items_to_load1)
#allFiles=os.listdir(filepath+sample)
#list = filter(fileExistTest,allFiles)

counter=0
for scan in range(len(items_to_load1)):
    QcodesExportTxt(path_to_db,items_to_load1[scan],
                    savepath,filenames1[scan]) 
##%%


# %%
counter=0
scan=52
QcodesExportTxt(path_to_db,items_to_load1[scan],
                    savepath,filenames1[scan]) 

# %%
