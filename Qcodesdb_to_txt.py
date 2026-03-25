#%% 
import numpy as np 
import matplotlib.pyplot as plt 
import os
from pyparsing import line

from qcodes import load_by_id, Parameter, ScaledParameter
from qcodes.dataset.plotting import plot_by_id
from qcodes.dataset.data_export import get_data_by_id, load_by_id, reshape_2D_data
from qcodes.dataset.data_set import DataSet
from qcodes.interactive_widget import experiments_widget



class Object():
    pass

def fileExistTest(s): 
    return '.txt' in s 

a=1
### the double gate map is ID 175 
def QcodesExportTxt(dbpath,scan_id,savepath,savename):
    allFiles=os.listdir(savepath)
    #print(allFiles)
    filename=savename+'_ID'+str(scan_id)+'.txt'
    if filename not in allFiles:# or os.stat(filename).st_size==0:
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
        if double_gate[meas_variables[0]][set_variables[0]][0]==double_gate[meas_variables[0]][set_variables[0]][1]:
            pass
            #print(scan_id)
        else:
            variables.reverse()
            set_variables.reverse()
            #print('reversed '+str(scan_id))
        if len(set_variables)==1:
            #linesweep
            type_of_scan= '% linesweep'
        elif len(set_variables)==2:
            #megasweep
            type_of_scan= '% megasweep '
        '''
        Create a new variables to store the names but with 
        the swept variables as the 1st and 2nd column
        '''
        for i in meas_variables:
            variables.append(i)
        #print(variables)
        for i in range(len(variables)):
            units.append(tb32_dataset.get_parameters()[i].unit)

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
        for i in range(len(first_variable)):
            if first_variable[i]!=first_variable[i-1] and i!=0:
                first_cut=i
                print(first_cut)
                break
        sweepcount=  int(len(first_variable)/first_cut)        
        for j in range(len(variables)):
            #print(j)
            for i in range(sweepcount):
                real_data[j].append(data_dg[j][first_cut*(i):(i+1)*first_cut])
                    #break 
        ### DETECT IF THE SCAN IS BROKEN ###
        real_data_len=len(data_dg)
        last_scan_len=len(real_data[0][sweepcount-1])
        if len(real_data[0][sweepcount-1]!=len(real_data[0][0])):
            print('Scan with id '+str(scan_id)+' broke too early')
            for i in range(len(variables)):
                good_data_dg[i]=data_dg[i][0:real_data_len-last_scan_len]
        else:
            good_data_dg[i]=data_dg[i]


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
    else:
        print(filename+' already exists.')
    return 
# %%
## Round 1 of loading, begining to 20/06
'''
Choose the folder and the database name that you want to open
'''
general_path="C:/Users/jaime/PhD/a. Data/IST_Austria/"
path_to_db=general_path+"tbg32.db"

'''
Write the items which you want to load, you could translate all
but probably many of the scans are not actually important, so 
there is no need to load them 
'''
#filepath = os.getcwd()

#exceptFiles=[]
### Initial measurements I:19-16, Vxx-10-11
items_to_load1=[62, 128,
                175, 186,209,    ##First measurements 27/05,29/05 
                251, 269,       ## Change current leads I:49-43 (fraunhofer)
                272,            ## Change current 28-11, JJ line 
                302,            # back gate vs DC
                312,313,314,316, # Double gate  
                329,343,344,345, #Fraunhofers
                353,354,              #Double gate complete
                368,370 ,          # Fraunhfoers agains at -580 mV
                401,402,            # Vbg and JJ line after Dacs
                434,441, 454,455,  ## JJ line after changing DACs         
                459,460,466,467,
                #468,477,478,
                    ]
savefolder='TB32/'
savepath=general_path+savefolder
filenames1=['220525_I-19-16_Vxx-10-11_mTVbg',   # Cooldown
            '220526_I-19-16_Vxx-10-11_mIVVbg', # Iv vs gates
            '220527_I-19-16_Vxx-10-11_mVbgVtg',# double gate 
            '220529_I-19-16_Vxx-10-11_mVbgIdc',# find SC?
            '220529_I-19-16_Vxx-10-11_Vbg(-590-Vtgdiv3p3)_mVbgVtg', # Jj line
            '220530_I-49-43_Vxx-30-10_Vbg(-580)_mBxIdc', # Fraunhofer 1
            '220531_I-49-43_Vxx-30-10_Vbg(-580)_mBxIdc_focus', #Fraunhofer focus
            '220531_I-28-11_Vxx-30-10_Vbg(-590-Vtgdiv3p3)_mVbg-Idc',
            '220601_I-50-43_Vxx-30-10_mVbg-Idc',
            '220602_I-50-43_Vxx-30-10_mVbgVtg_n1hfocus_1',# double gate 
            '220602_I-50-43_Vxx-30-10_mVbgVtg_n1hfocus_2',# double gate 
            '220602_I-50-43_Vxx-30-10_mVbgVtg_n1hfocus_3',# double gate 
            '220602_I-50-43_Vxx-30-10_mVbgVtg_n1hfocus_4',# double gate 
            '220602_I-50-43_Vxx-30-10_Vbgn580mV_mBxIdc',
            '220603_I-50-43_Vxx-30-10_Vbgn580mV_mBxIdc',
            '220604_I-50-43_Vxx-30-10_Vbgn570mV_mBxIdc',
            '220604_I-50-43_Vxx-30-10_Vbg(-590-Vtgdiv3p3)_Vtg-520mV_mBxIdc',
            '220609_I-50-43_Vxx-30-10_mVbgVtg_Vtg3ton0p5',# double gate 
            '220609_I-50-43_Vxx-30-10_mVbgVtg_Vtgn1to3',
            '220611_I-50-43_Vxx-30-10_Vbgn580mV_mBxIdc',
            '220611_I-50-43_Vxx-30-10_Vbgn580mV_mBxIdc_focus',
            '220614_I-50-43_Vxx-30-10_Vbg(-590-Vtgdiv3p3)_mVtg-Idc',
            '220614_I-50-43_Vxx-30-10_mVbgIdc',
            '220616_I-50-43_Vxx-30-10_Vbg(-600-Vtgdiv3p3)_mVtg-Idc',
            '220616_I-50-43_Vxx-30-10_Vbg(-580-Vtgdiv3p3)_mVtg-Idc',
            '220617_I-50-43_Vxx-30-10_Vbg(-590-Vtgdiv3p3)_mVtg-Idc_1hfocus_1',
            '220617_I-50-43_Vxx-30-10_Vbg(-590-Vtgdiv3p3)_mVtg-Idc_1hfocus_2',
            '220617_I-50-43_Vxx-30-10_Vbg(-590)_mBxIdc',
            '220617_I-50-43_Vxx-30-10_Vbg(-590)_mBxIdc_focus',
            '220618_I-50-43_Vxx-30-10_Vbg(-600-Vtgdiv3p3)_Vtg(-310mV)_mBxIdc',
            '220619_I-50-43_Vxx-30-10_Vbg(-600-Vtgdiv3p3)_Vtg(-300mV)_mBxIdc',
            '220620_I-50-43_Vxx-30-10_Vbg(-600-Vtgdiv3p3)_Vtg(-1p7V)_mBxIdc',
            '220620_I-50-43_Vxx-30-10_Vbg(-590-Vtgdiv3p3)_Vtg(60mV)_mBxIdc',
            '220620_I-50-43_Vxx-30-10_Vbg(-590-Vtgdiv3p3)_Vtg(100mV)_mBxIdc',
            
            ]
#allFiles=os.listdir(filepath+sample)
#list = filter(fileExistTest,allFiles)

counter=0
for scan in range(len(items_to_load1)):
    QcodesExportTxt(path_to_db,items_to_load1[scan],
                    savepath,filenames1[scan])
    counter=counter+1  
#%%

# %%
##Round two of loading 21/06
'''
Choose the folder and the database name that you want to open
'''
general_path="C:/Users/jaime/PhD/a. Data/IST_Austria/"
path_to_db=general_path+"tbg32.db"

'''
Write the items which you want to load, you could translate all
but probably many of the scans are not actually important, so 
there is no need to load them 
'''
#filepath = os.getcwd()

#exceptFiles=[]
### Initial measurements I:19-16, Vxx-10-11
items_to_load2=[472,477,478,
                486, 488, 
                489,490,491,492,493
                    ]
savefolder='TB32/'
savepath=general_path+savefolder
filenames2=['220620_I-50-43_Vxx-30-10_Vbg(-600mV-Vtgdiv3p3)_Vtg(-1p7V)_mBxIdc',
            '220620_I-50-43_Vxx-30-10_Vbg(-590mV-Vtgdiv3p3)_Vtg(60mV)_mBxIdc',
            '220620_I-50-43_Vxx-30-10_Vbg(-590mV-Vtgdiv3p3)_Vtg(100mV)_mBxIdc',

            '220620_I-50-43_Vxx-30-10_Vbg(-590mV-Vtgdiv3p3)_Idc(0V)_mBxVtg',
            '220620_I-50-43_Vxx-30-10_Vbg(-590mV-Vtgdiv3p3)_mB-hyst_Idc',

            '220621_I-50-43_Vxx-30-10_Vbg(-590mV-Vtgdiv3p3)_Vtg(20mV)_mBxIdc',
            '220621_I-50-43_Vxx-30-10_Vbg(-590mV-Vtgdiv3p3)_Vtg(40mV)_mBxIdc',
            '220621_I-50-43_Vxx-30-10_Vbg(-590mV-Vtgdiv3p3)_Vtg(60mV)_mBxIdc',
            '220621_I-50-43_Vxx-30-10_Vbg(-590mV-Vtgdiv3p3)_Vtg(160mV)_mBxIdc',
            '220621_I-50-43_Vxx-30-10_Vbg(-590mV-Vtgdiv3p3)_Vtg(120mV)_mBxIdc',
            ]
#allFiles=os.listdir(filepath+sample)
#list = filter(fileExistTest,allFiles)

counter=0
for scan in range(len(items_to_load2)):
    QcodesExportTxt(path_to_db,items_to_load2[scan],
                    savepath,filenames2[scan])
    counter=counter+1  

#%%
### ROund three of translating, this should be the script to follow from now on.
### Now I just need the ID number and the name of the file is written in the txt file

general_path="C:/Users/jaime/PhD/a. Data/IST_Austria/"
path_to_db=general_path+"tbg32.db"

'''
Write the items which you want to load, you could translate all
but probably many of the scans are not actually important, so 
there is no need to load them 
'''
#filepath = os.getcwd()

#exceptFiles=[]
### Initial measurements I:19-16, Vxx-10-11
items_to_load3=[496, 497,
            499,500, 
            501,502, 503]
savefolder='TB32/'
savepath=general_path+savefolder

filenames3=[]
with open(general_path+'TB32_scan_names.txt','r') as f:
    lines=f.readlines()
for i in range(len(lines)):
    split_line=lines[i].split()
    if float(split_line[0]) in items_to_load3:
        filenames3.append(split_line[1])
print(filenames3)

counter=0
for scan in range(len(items_to_load3)):
    QcodesExportTxt(path_to_db,items_to_load3[scan],
                    savepath,filenames3[scan])
    counter=counter+1  


#%%
filenames3=[]
with open(general_path+'TB32_scan_names.txt','r') as f:
    lines=f.readlines()
for i in range(len(lines)):
    split_line=lines[i].split()
    if float(split_line[0]) in items_to_load3:
        filenames3.append(split_line[1])
print(filenames3)
#%% 
#%%
'''
filenames=filenames1+filenames2
# open a txt and write the names I have just created
## this is actually just an exercise to apply it in the actual measurement code 
try:
    with open(general_path+'TB32_IST.txt') as f:
        lines = f.readlines()
        f.close()
except: 
    FileNotFoundError
with open(general_path+'TB32_IST.txt', 'a') as f:
    for i in range(len(filenames)):
        if filenames[i]+'\n' not in lines:
            f.write(filenames[i])
            f.write('\n')
    f.close()
# %%
with open(general_path+'TB32_IST.txt') as f:
    lines = f.readlines()
# %%
lines[0][:-1]
# %%
'''