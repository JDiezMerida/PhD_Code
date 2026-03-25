# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 16:40:13 2020
File to plot data imported from Labber with the exportDataToText_v8 file 
This includes codes to plot linescans and maps, as well as some interactive codes to quickly compare data.
@author: Jaime Diez Mérida@LDQM-LS-Efetov group at ICFO-LMU. 
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import matplotlib
import time 

#import matplotlib.cm as cm
#import math as math
#from scipy import signal
#from scipy.fftpack import fft, ifft
#from scipy.signal import savgol_filter
#from collections import Counter
#from ipywidgets import interactive
#import importlib
import exportDataToText_v6
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

orig_cmap = matplotlib.cm.RdBu_r

#%% PARSING AND LOADING DATA FOR ANALYSIS

'''
Dummy class. Necessary to be able to add attributes to the class object later
on.
'''
class Object():
    pass
'''
The ones I always use are plot_lines_fromfiles and plot_map_fromfiles
The others are either old versions or not really used 
'''

def plot_lines_fromfiles(file,x_variable=0,plotvar=False):
    '''
    Read and load the data from linescan files
    Input
    ---
    File: filename including the directory
    x_variable: in principle it will be the first column, but it can be 
    changed to other colum
    plotvar: True will plot all the variables, False will not plot 
    
    Return
    ---
    objdata: an Object with all the data
    variables: variable names
    filename: filename without the directory and the .txt ending
    '''
    
    ## First parse the data and convert it into an object, a np.array 
    ## and an array with names
    objdata,rawdata, variables=parse_data(file)
    
    ## Extract the filename from the directory name
    filename=file.split('/')[-1]
    filename=filename[:-4]
    
    ## Plot all variables, in case plotvar=True
    if plotvar:
        for k in range(len(variables)-1):
            plot_linescans_raw(objdata,variables[x_variable],variables[k],filename)

    return  objdata,variables,filename

'''
Code to plot raw linescans from the data. No conversion of units. 
In order to use it in a Jupyter logbook typical code to copy and paste
in the cell would be:
 
directory = ' ' #general directory
folder=' ' #specific folder of the data 
directory= directory+folder

files=['TBGTB5_4 probe_11_7_9_8_3K_bias10nA_sen200uV_tconst300ms_gatetest1.txt',
       'TBGTB5_4 probe_11_7_9_8_3K_bias10nA_sen200uV_tconst300ms_gatetest2.txt',
       ]

obj = {} #stores all the data
raw ={} #store the xaxis data for plotting
variables={} #store the data to plot in the right units
name={}
# scan type = 0

for i in range(len(files)):
    file=directory+files[i]
    obj['data'+str(i)], variables['data'+str(i)],name['data'+str(i)]=LDQM_dataplot.plot_lines_fromfiles(file)
 
'''

def parse_data(fname,skiprows=3):
    '''
    Cleans the data of the txt file so that it can be converted into an object
    It stores the data as np.array and the headers in an array. Then uses this
    to create the data object. 
    '''
    with open(fname) as myfile:
        
        head = [next(myfile) for x in range(skiprows)]

        head_names = head[0].replace(' ','')
        head_names = head_names.replace('-','')
        #scan_type=head[2]
        #head_unit = head_names.replace('\n','')
        data = np.loadtxt(fname,delimiter = '\t',skiprows=skiprows)
        names = head_names.split("\t");
    #print(head_names)
    # Parse the data into a class with attributes    
    fdata = Object()  
    for i in range(len(names)-1):

        setattr(fdata, names[i], data[:,i])

    return fdata,data,names


def plot_linescans_raw(fdata,x_axis='',y_axis='',file=''):
    '''
    Plots all the variables of a file when loading the data 
    In case the data is somehow wrong it will give an error, this should not
    happen ofter as there are error corrections in the parsing steps as well.
    
    '''
    #file =dataset
    try: 
        title=file 
        #split_title=dataset.split('/')[-1].split('.')
        #for i in range(len(split_title)-1):
        #    title=title+' '+split_title[i]
        plt.figure(figsize=(6,4))
        plt.plot(getattr(fdata, x_axis),getattr(fdata, y_axis),'b')
        plt.title(title)
        plt.xlabel(x_axis)
        plt.ylabel(y_axis)
        plt.grid() 
        plt.ticklabel_format(axis='both', style='sci', scilimits=(-3,3))
        plt.tight_layout()
        plt.show()
    except AttributeError:
        print ('Data corrupted')
    return 


def plot_map_fromfiles(file,plotvar=False):
    '''
    Read and load the data from 2D map files
    Input
    ---
    File: filename including the directory
    plotvar: True will plot all the variables, False will not plot 
    
    Return
    ---
    objdata_mapshape: an Object with all the data
    map_variables: variable names
    map_rawdata: np.array with all the data
    '''
    
    ## Parse data same way as for the linescans
    map_objdata,map_rawdata, map_variables=parse_data(file)
    
    filename=file.split('/')[-1]
    
    ## Organize and plot the data of the maps
    for k in range(len(map_variables)-1):
        #if scan_type=='megasweep\n':
        mapdata=map_rawdata
        mapnames=map_variables
    
        objdata_mapshape,objdata_mapshape_real=plot_maps_v3(mapdata,mapnames,map_variables[1],map_variables[0],map_variables[k],filename,plotvar)

    #return objdata_mapshape,map_variables,objdata_mapshape_real
    return objdata_mapshape,map_variables,map_rawdata

'''
Code to plot raw maps

Again to use it in Jupyter, typical code would be:
    
B_directory = 'C:/jdiez_local/LDQM/a. Data/'
B_folder='/TB32/LandauFan/'
B_directory= B_directory+B_folder

## 17-18 and 12-13 are the JJ
## 8-17 is the other side of 11-12 which is the other good SC
## 18-19 is the other side of 13-16, which is the one with the fractional insulators 

B_files=['J-TB32_Landaufan_35mK_1nA_I-11-16_Rxx1-17-18_Rxx4-12-13_Rxx3-12-17_mBVbg.txt', # contact in the JJ
         'J-TB32_LandauFan_5nA_I-9-20_Rxx1-8-17_Rxx4-18-19_Rxx3-17-18_mBVBG.txt', # all the other contacts
         'J-TB32_LandauFan_5nA_I-9-20_Rxx1-8-11_Rxx4-18-13_Rxx3-19-16_Rxys_mBVBG.txt' # Rxy of all the other contacts
         
]

B_obj = {} #stores all the data
B_raw ={} #store the xaxis data for plotting
B_variables={} #store the data to plot in the right units

# All IV curves
# scan type = 0
for i in range(len(B_files)):
    file=B_directory+B_files[i]
    B_obj['data'+str(i)], B_variables['data'+str(i)],B_raw['data'+str(i)]=LDQM_dataplot.plot_map_fromfiles(file)

'''

def plot_maps_v3(mapdata,mapnames,x_axis='',y_axis='',z_axis='',file='',plotvar=True):
    '''
    Convert the 2D map data into the correct array form and plots it if plotvar=True
    It uses plt.pcolormesh, which gives the right aspect ratio when plotting
    Input
    ---
    mapdata: raw data extracted from parse_data
    mapnames: variables used in the scans 
    axis_names: these will be written in the plots 
    file: used to give the title of the plot
    plotvar: decides whether there is plotting or not
    
    Return
    final_map_data: Object with all the data organized
    final_map_data_real: Object with all the data as it comes from the raw data
    --- 
        
    '''
    ## First convert to right structure. The original data is a single long array 
    ## This code breaks it into the actual sweeps, so that we can build a map 
    ## The code uses the "slow variable" to create breaking points where the sweeps
    ## stop and start the next one
    
    breakpoint_vec=[]
    for i in range(len(mapdata[:,0])):
        if mapdata[:,1][i]==mapdata[:,1][0]:
            breakpoint_vec.append(i)
    sweep_number=len(breakpoint_vec) # how many sweeps were in the data? 
    variables=len(mapdata[0])
    break_point=int(len(mapdata[:,0])/sweep_number)
    
    ## Create the Objects to store the data 
    final_map_data=Object()
    final_map_data_real=Object()
    
    ## Store the data into np.arrays of the right size and finally add the 
    ## data to the object organized with the variable names 
    
    for i in range(variables):
        image_shape=np.zeros((len(breakpoint_vec),break_point))
        image_shape_real=np.zeros((len(breakpoint_vec),break_point))
        for j in range(sweep_number):
            
            ## This one has everything as in the data
            if j==sweep_number-1:
                image_shape_real[j]=mapdata[breakpoint_vec[j]:,i]
            else:
                image_shape_real[j]=mapdata[breakpoint_vec[j]:breakpoint_vec[j+1],i]
               
            ## This inverts the order of the arrays to make plotting easier,
            ## This is the one that will be stored in the plot_map_fromfiles
            if j==sweep_number-1:
                image_shape[sweep_number-1-j]=mapdata[breakpoint_vec[j]:,i]
            else:
                image_shape[sweep_number-1-j]=mapdata[breakpoint_vec[j]:breakpoint_vec[j+1],i]
            
        setattr(final_map_data,mapnames[i],image_shape)
        setattr(final_map_data_real,mapnames[i],image_shape_real)
        
    ## Whether to plot the data or not 
    if plotvar==True:
        plt.figure(figsize=(10,7))
    
        x=getattr(final_map_data,x_axis)
        y=getattr(final_map_data,y_axis)
        z=getattr(final_map_data,z_axis)
    
        plt.pcolormesh(x,y,z,cmap='RdBu_r')
        title=file
        plt.xlabel(x_axis)        
        plt.ylabel(y_axis)
        cb=plt.colorbar()
        cb.set_label(z_axis,fontsize=12)
        cb.ax.tick_params(labelsize=12)
        plt.title(title)
        plt.show()
        
    return final_map_data,final_map_data_real

#%% INTERACTIVE LINESCANS AND MAPS FOR THE JUPYTER 

'''
These are different codes to create interactive plots. They can be nice to quickly
compare data, but in the end they are a bit complicated to use,
 because making them interactive means they are rather cumbersome 
I rarely use them anymore. 

'''

def int_singlelinescans(x_data,y_data, dev_name, x_name,y_name,
                        scan_type, label_x='',label_y='',title_label='', log_bool=False,
                        save=False, conduc=False, eh=False,file='',hd=0,font=12):
    
    #scan_name=['Back gate sweep ','Top gate sweep ', 'Field sweep ','IV ',]
    #x_axis= [r'$Back\ Gate\ [V]$',r'$Top\ Gate\ [V]$',r'$Field\ [T]$',r'$Bias\ [mV]$']
    
    plt.figure(figsize=(10,7))
    cond=np.zeros_like(y_data)
    cond_eh=np.zeros_like(y_data)
    for i in range(len(y_data)):
        cond[i] = 1/(y_data[i])  #units in mS
        cond_eh[i] =cond[i]*1E-3/(2*3.877E-5)
    if conduc==False and eh==False:
        plt.plot(x_data,y_data,'b')

    elif conduc==True and eh==False:
        plt.plot(x_data,cond,'b')

    elif eh==True:
        plt.plot(x_data,cond_eh,'b')

    if title_label=='':
        plt.title(dev_name+' '+file,fontsize=font )
    else:
        plt.title(title_label,fontsize=font )
    if label_x!='':
        plt.xlabel(label_x,fontsize=font)
        plt.ylabel(label_y,fontsize=font)
    else:
        plt.xlabel(x_name,fontsize=font)
        plt.ylabel(y_name,fontsize=font)
    if log_bool:
        plt.yscale('log')
    else:
        plt.ticklabel_format(axis='both', style='sci', scilimits=(-3,3))
    #plt.ylabel(r'$dV/dI\ [k\Omega]$')
    plt.grid()
    plt.tight_layout()
    if save:
        if hd==0:
            plt.savefig(dev_name+title_label +'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300)
           
        elif hd==1:
            plt.savefig(dev_name+title_label +'_'+time.strftime("%H%M%d%m%y")+'.png', dpi=600)
        elif hd==2:
            plt.savefig(dev_name+title_label +'_'+time.strftime("%H%M%d%m%y")+'.pdf')
    plt.show()
    return

'''
For the maps copy this code

%matplotlib inline

#x_dict=xIV
dev_names = ['','TBGTB5 ']
def interactive_maps(dataset=1,plot_data=2,x_data=0,y_data=1,lockin=1, scan_type=0,
                          labels=False,save=False, conduc=False, eh=False,extent_bool=False):
    
    #Choose the device you want to plot
    scan_names=TB5_variables_map['data'+str(dataset)+'sr'+str(lockin)]
    data_obj=TB5_objdata_mapshape['data'+str(dataset)+'sr'+str(lockin)]

    #x_axis=getattr(data_obj,scan_names[x_data])
    #y_axis=getattr(data_obj,scan_names[y_data])
    #z_axis=getattr(data_obj,scan_names[plot_data])
    x_name=scan_names[x_data]
    y_name=scan_names[y_data]
    z_name=scan_names[plot_data]
    #dev_norm=norm_dict['data'+str(dataset)+'sr'+str(lockin)]
    int_maps(data_obj,z_name,x_name, y_name,dev_names[lockin],
            lockin,scan_type,labels ,save, conduc, eh,extent_bool)
    
    return 
interactive(interactive_maps,dataset=(0,30,1), plot_data=(0,15,1), x_data=(0,15,1),y_data=(0,15,1),lockin=(0,2,1), scan_type=(0,6,1))

''' 

def int_maps(data, z_name,x_name,y_name, I_name, dev_name,min_value=0,max_value=1,
                        x_start=0,x_finish=1000,y_start=0, y_finish=1000, 
                        label_x='',label_y='', title_label='',save_label='',
                        norm_current=False,norm_quantum=False,
                        conduc=False, eh=False,
                        log_bool=False, abs_bool=False,
                        extent_bool=False,
                        minmax_bool=False,cut_bool=False,save=False, file='',font=15,color='RdBu_r',I_mult=1,
                        colorshift=0.5, orig_cmap=orig_cmap,hd=0):
    #print(data.VxxR)
    
    mycolor=shiftedColorMap(colorshift,orig_cmap)
    
    # Prepare the data to plot 
    
    x_data=getattr(data,x_name)
    y_data=getattr(data,y_name)
    z_data=getattr(data,z_name)

    # Choose whether to plot the whole data or not 
       
    if cut_bool:
        
        
        if y_finish-y_start>=len(z_data):
            y_finish=len(z_data)-1
        if x_finish>len(z_data[0]):
            x_finish=len(z_data[0])-1
        
        
        
        if extent_bool:
            extent=[getattr(data,x_name)[0,x_start],
            getattr(data,x_name)[0,x_finish],
            getattr(data,y_name)[y_finish,0],
            getattr(data,y_name)[y_start,0]]
        else:
            extent=[x_start,x_finish,  #x_axis
               y_start,y_finish]   
            
        plot_data=np.zeros((y_finish-y_start,x_finish-x_start))
        current_data=np.zeros_like(plot_data)
        plot_x_data=np.zeros_like(plot_data)
        plot_y_data=np.zeros_like(plot_data)
        count=y_start
        
        plot_x_data=x_data[y_start:y_finish,x_start:x_finish]
        plot_y_data=y_data[y_start:y_finish,x_start:x_finish]
        plot_data=z_data[y_start:y_finish,x_start:x_finish]
        current_data=getattr(data,I_name)[y_start:y_finish,x_start:x_finish]*I_mult
        '''
        for i in range(y_finish-y_start):
            
            plot_data[i]=z_data[count,x_start:x_finish]
            current_data[i]=getattr(data,I_name)[count,x_start:x_finish]*I_mult
            plot_x_data[i]=x_data[count,x_start:x_finish]
            plot_y_data=y_data[count,x_start:x_finish]
            
            count=count+1
        '''
    else:
        plot_data=z_data
        current_data=getattr(data,I_name)*I_mult
        plot_x_data=x_data
        plot_y_data=y_data
        
        
        if extent_bool:
            extent=[getattr(data,x_name)[0,0],
            getattr(data,x_name)[0,-1],
            getattr(data,y_name)[-1,0],
            getattr(data,y_name)[0,0]]
        else:
            extent=[0,len(z_data[0]),  #x_axis
               0,len(z_data)]   
        
    # Which type of data, raw, resistance, conductance

            
    cond=np.zeros_like(plot_data)
    cond_eh=np.zeros_like(plot_data)
    norm_Rxx=np.zeros_like(plot_data)
    norm_Rxx_he=np.zeros_like(plot_data)
    #log_data=np.zeros_like(plot_data)
    
        
    for i in range(len(plot_data)):
        for j in range(len(plot_data[i])):
            norm_Rxx[i,j]=plot_data[i,j]/current_data[i,j]/1000 # kOhms
            norm_Rxx_he[i,j]=norm_Rxx[i,j]/25.8128 #h/e2
            cond[i,j] = 1/(norm_Rxx[i,j])  #1/kOhm
            cond_eh[i,j] = 1/norm_Rxx_he[i,j] # e^2/h
    
        
    # Now make the figure 
    
    plt.figure(figsize=(12,7))
    
    
    
    if conduc==True:
        plot_data=cond
        z_name = 'Conductance (mS)' 
        
    elif eh==True:
        plot_data=cond_eh
        z_name = 'Normalized Conductance (Units of $e^2/h$)' 
        
    elif norm_current==True:
        plot_data=norm_Rxx
        z_name = 'Resistance (k$\Omega$)' 
           
    elif norm_quantum==True:
        plot_data=norm_Rxx_he
        z_name = 'Normalized Resistance (Units of $h/e^2$)' 
    
    if log_bool:
        log_data=np.log(plot_data)
        plot_data=log_data
        
    if abs_bool:
        abs_data=np.abs(plot_data)
        plot_data=abs_data
        
    if minmax_bool:
        plt.pcolormesh(plot_x_data,plot_y_data,plot_data,cmap=mycolor,vmin=min_value,vmax=max_value)
        #plt.imshow(plot_data,color,aspect='auto',extent=extent,vmin=min_value,vmax=max_value)
    else:
        plt.pcolormesh(plot_x_data,plot_y_data,plot_data,cmap=mycolor)
        #plt.imshow(plot_data,color,aspect='auto',extent=extent)
    
    cb=plt.colorbar()
    cb.set_label(z_name,fontsize=font)
    cb.ax.tick_params(labelsize=font)
    if title_label=='':
        plt.title(file)
    else:
        plt.title(title_label)
    #split_title=dataset.split('/')[-1].split('.')
    #for i in range(len(split_title)-1):
    #    title=title+' '+split_title[i]
    if label_x=='':
        plt.xlabel(x_name,fontsize=font)
        plt.ylabel(y_name,fontsize=font)    
    else:
        plt.xlabel(label_x,fontsize=font)
        plt.ylabel(label_y,fontsize=font)   
    plt.xticks(fontsize=font,rotation=0)
    plt.yticks(fontsize=font,rotation=0)
    #plt.xlabel(r'$B\ (T)$',fontsize=12)
    
    if save:
        if save_label=='':
            plt.savefig(file +'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300)
        else:
            if hd==0:
                plt.savefig(dev_name+save_label +'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300)
               
            elif hd==1:
                plt.savefig(dev_name+save_label +'_'+time.strftime("%H%M%d%m%y")+'.png', dpi=600)
            elif hd==2:
                plt.savefig(dev_name+save_label +'_'+time.strftime("%H%M%d%m%y")+'.pdf')
                
    plt.show()
    return 

'''
This makes the interactive linescans from the maps. It allows you to access
different lines of the maps, both in the vertical and horizontal direction
Again with scan_type it should ensure that the correct map and correct units 
are plotted 
The code to run such a structure is shown below:
'''
#Example code
'''
%matplotlib inline


dev_name = 'TBGTB4 '

def interactive_map_linescans(dataset=0,plot_data=2,x_data=1,y_data=0, 
                        scan_type=0,
                        x_st=0 ,x_ext=20, x_int=1, x_offset=0, 
                        y_st=0,y_ext=10,y_int=1,y_offset=0,
                        plot_x=True,plot_y=False, 
                        label_x='',label_y='', title_label='',
                        save=False, conduc=False, eh=False,legend=False,interlocked=False):

    #Choose the device you want to plot
    scan_names=TB4_map_variables['data'+str(dataset)+'sr1']
    data_obj=TB4_objdata_mapshape_real['data'+str(dataset)+'sr1']

    #x_axis=getattr(data_obj,scan_names[x_data])
    #y_axis=getattr(data_obj,scan_names[y_data])
    #z_axis=getattr(data_obj,scan_names[plot_data])
    x_name=scan_names[x_data]
    y_name=scan_names[y_data]
    z_name=scan_names[plot_data]
    
    LDQM_dataplot.make_int_linescans(data_obj,z_name,x_name, y_name,dev_name,
            scan_type, 
            x_st ,x_ext, x_int, x_offset, 
            y_st,y_ext,y_int,y_offset,                         
            plot_x,plot_y, 
            label_x,label_y, title_label,
            save, conduc, eh,legend,interlocked)
    
    
interactive(interactive_map_linescans,dataset=(0,30,1), plot_data=(0,15,1), 
            x_data=(0,15,1),y_data=(0,15,1), scan_type=(0,6,1),
            x_st=(0,1000,1),x_ext=(0,1000,1),x_int=(0,1000,1),x_offset=(0,1,0.00001),
            y_st=(0,1000,1),y_ext=(0,1000,1),y_int=(0,1000,1),y_offset=(0,1,0.00001),
            )
'''


# Create a new colormap


def make_int_linescans(file,plot_axis,data,z_name,x_name,y_name,I_name,dev_name, 
                       plot_start=0,plot_extend=10,plot_interval=1,plot_offset=0,
                       label_int=1,
                       label_z='',label_x='', title_label='',
                       legend_label='',
                       x_min=0,x_max=0,y_min=0,y_max=0, line_marker='',zero_line=False,
                       hor_point=0,x_figsize=10,y_figsize=7,
                       save=False, conduc=False, eh=False,norm_current=False,
                       norm_quantum=False,log_bool=False,abs_bool=False,legend=True,
                       background=False,x_conversion=1,
                       expand=10, I_mult=1, font=15, fontleg=15,
                       hd=0, colormap='nipy_spectral',color_st=0,color_end=0.9,
                       appended=False,interlocked=False, symmetrize=False,
                       backg_line=0):
    
    
    ### Make the arrays of data for horizontal and vertical from the image
    
    
    if appended==False:
        y_axis=getattr(data,y_name)[:,0]
    else:
        y_axis=getattr(data,y_name)[:][0]
        
    x_axis=getattr(data,x_name)[0]*x_conversion
    x_size=len(x_axis)
    y_size=len(y_axis)
    z_axis=getattr(data,z_name)
    current_data=np.zeros_like(z_axis)
    current_data=getattr(data,I_name)*I_mult
        
    
   ### Choose which lines you wish to plot. There should be a max limit of lines to avoid errors that
    ### you want to plot over the limit of the array
    ### start is where it starts and ext(extension) adds how many more lines
    ### If the end is further than the real end then we put a limit on that. This has to do in the 
    ### plotting stage
    
    if plot_axis==0: #plot_x
        plot_size=y_size
        legend_name=y_name
        plot_xaxis=x_axis
        plot_interval_axis=y_axis
        plot_xlabel=label_x
        x_name=x_name
        
    elif plot_axis==1: #plot_y
        plot_size=x_size
        legend_name=x_name
        plot_xaxis=y_axis
        plot_interval_axis=x_axis
        plot_xlabel=label_x
        x_name=y_name

    plot_end, Rxx,Rxx_he,cond, cond_eh, plot_lines, plot_colors, plot_legend_label =int_plot_dataprep(z_axis,current_data,
                                                                                x_size, 
                                                                               y_size, plot_size, 
                                                                               plot_start, plot_extend,plot_interval,
                                                                               legend_label, legend_name,colormap
                                                                               ,color_st,color_end,appended)
    
      
    ### Make the horizontal plot

    if conduc==False and eh==False and norm_current==False and norm_quantum==False:
        plot_yaxis=z_axis
        
    elif conduc==False and eh==False and norm_current==True and norm_quantum==False:
        plot_yaxis=Rxx
        
    elif conduc==False and eh==False and norm_quantum==True:
        plot_yaxis=Rxx_he
        
    elif conduc==True and eh==False and norm_quantum==True:
        plot_yaxis=cond
        
    elif eh==True:
        plot_yaxis=cond_eh
        
    
        
    if abs_bool:
        abs_data=np.abs(plot_yaxis)
        plot_yaxis=abs_data
    '''   
    if log_bool:
        log_data=np.log(plot_yaxis)
        plot_yaxis=log_data
     ''' 
    
        
    if plot_axis==0:
        int_horplot(plot_xaxis, plot_yaxis, plot_start, plot_end, plot_lines, 
                    plot_interval_axis, label_int, plot_offset,
                    plot_legend_label, plot_colors ,expand, interlocked,
                    symmetrize, background, backg_line,line_marker,
                    x_figsize,y_figsize)
    
    elif plot_axis==1:
        int_vertplot(plot_xaxis, plot_yaxis, plot_start, plot_end, plot_lines, 
                    plot_interval_axis, label_int, plot_offset,
                    plot_legend_label ,plot_colors, expand, interlocked,line_marker,
                    x_figsize,y_figsize)
    if log_bool:
        plt.yscale('log')
    ### Add the rest of the details  
    if zero_line:
        plt.plot(plot_xaxis,plot_xaxis*hor_point,'grey',linestyle='dashed',linewidth=0.5)
        
    int_plot_details(file,dev_name, title_label, plot_offset, plot_xlabel, label_z, 
                     x_min, x_max, y_min, y_max,
                     x_name, z_name, save, conduc, eh, norm_current,
                     norm_quantum,legend,
                     expand, font, fontleg, hd)
    return







def int_plot_dataprep(z_axis,current_data, x_size, y_size, plot_size, start, extension, interval,
                      legend_label, legend_name,colormap,color_st,color_end,appended):

    end=start+extension
    if end>=plot_size:
        end=plot_size
        print('Index'+ str(end)+' out of range')

    plot_lines=np.arange(start,end,interval)

    ###Create a conductance array to plot conductance . Both in mS and in e^2/h units

    cond=np.zeros_like(z_axis)
    cond_eh=np.zeros_like(z_axis)
    Rxx=np.zeros_like(z_axis)
    Rxx_he=np.zeros_like(z_axis)
    #log_data=np.zeros_like(plot_data)
    if appended==False:
        for i in range(y_size):
            for j in range(x_size):
                Rxx[i,j]=z_axis[i,j]/current_data[i,j]/1000 # kOhms
                Rxx_he[i,j]=Rxx[i,j]/25.812 #h/e2
                cond[i,j] = 1/(Rxx[i,j])  #1/kOhm
                cond_eh[i,j] = 1/Rxx_he[i,j] # e^2/h
    else:
        for i in range(y_size):
            for j in range(x_size):
                Rxx[i][j]=z_axis[i][j]/current_data[i][j]/1000 # kOhms
                Rxx_he[i][j]=Rxx[i][j]/25.812 #h/e2
                cond[i][j] = 1/(Rxx[i][j])  #1/kOhm
                cond_eh[i][j] = 1/Rxx_he[i][j] # e^2/h

    ###For the colorprogression of the linescans. The idea is to have a progressing color code
    ###with always the same limits, to make it clear in which direciton the plotting advances
    
    plot_num=len(plot_lines)
    

    nipyspbig = plt.get_cmap(colormap, 512)
    newcmp = ListedColormap(nipyspbig(np.linspace(color_st, color_end, 256)))
    
    cmap = plt.get_cmap(newcmp)
    
    plot_colors = [cmap(i) for i in np.linspace(0, 1, plot_num)]

    
    #Prepare all the labels to make good plots

    if legend_label=='':
        legend_label=legend_name
        
    return end,Rxx,Rxx_he,cond,cond_eh,plot_lines,plot_colors,legend_label


def int_plot_details(file,dev_name,title_label,offset,label_x,label_y,x_min=0,x_max=0,y_min=0,y_max=0,
                       x_name='',z_name='', save=False, conduc=False, eh=False,
                       norm_current=False,norm_quantum=False,legend=False,expand=10,font=15,fontleg=15,hd=0):
    
    '''Details are inputted here, maybe they can also be moved to a different function''' 
        
    if x_max!=0:  
        plt.xlim(x_min,x_max)
    if y_max!=0:
        plt.ylim(y_min,y_max)
        
    if title_label=='':
        title_label=file
    if offset==0.0:
        title=dev_name +title_label
    elif offset!=0.0:
        title=dev_name+title_label+ ' offseted by ' +str(offset)

            
    plt.title(title,fontsize=font)
        
    if label_y=='':
        plt.xlabel(x_name,fontsize=font)
        #plt.locator_params(axis='x', nbins=15)
        plt.ylabel(z_name,fontsize=font)
       
        if norm_current== True and norm_quantum==False and conduc==False and eh==False:
            plt.ylabel(('R '+z_name+r'$ [k\Omega]$'),fontsize=font) 
            
        elif conduc==True and eh==False:
            plt.ylabel('G '+z_name+' mS',fontsize=font)
        elif eh==True:
            plt.ylabel(('G '+z_name+r'$ [e^2/h]$'),fontsize=font)   
        
        elif norm_quantum==True and conduc==False and eh==False:
            plt.ylabel(('R '+z_name+r'$ [h/e^2]$'),fontsize=font)
    else:
        plt.xlabel(label_x,fontsize=font)
        #plt.locator_params(axis='x', nbins=15)
        plt.ylabel(label_y,fontsize=font)
        
    plt.tick_params(axis='both',which='major',labelsize=font)
        
    if legend:
        plt.legend(loc='upper center', bbox_to_anchor=(1, 1), shadow=True, ncol=1,fontsize=fontleg)
    #plt.grid()
    if save:
        if hd==0:
            plt.savefig(title+time.strftime("%H%M%S%d%m%y")+'.png', dpi=300)
        elif hd==1:
            plt.savefig(title+time.strftime("%H%M%S%d%m%y")+'.png', dpi=600)
        elif hd==2:
            plt.savefig(title+time.strftime("%H%M%S%d%m%y")+'.pdf')
    

    plt.show()    

''' 
Plotting for interactive linescans. Code usedin the above lines for doing all
the plots
'''


def int_horplot(x, y, plot_start, plot_end, plot_lines, lines_label ,label_int, 
                offset, labels, colors,
                expand, interlocked, symmetrize, background, backg_line,line_marker,
                x_figsize, y_figsize):
    
    
    plt.figure(figsize=(x_figsize,y_figsize+expand*offset))
    lines_to_label=np.arange(plot_start,plot_end,label_int)
    count =0
    
    for i in plot_lines:
        if interlocked:

            if symmetrize:
                if i in lines_to_label:  
                    plt.plot(x[plot_start:plot_end],(y[i,plot_start:plot_end]+np.flip(y[i,plot_start:plot_end]))/2+offset*count, label=labels+'{0:.2f}'.format(lines_label[i]),color=colors[count],marker=line_marker)
                else:
                    plt.plot(x[plot_start:plot_end],(y[i,plot_start:plot_end]+np.flip(y[i,plot_start:plot_end]))/2+offset*count,color=colors[count],marker=line_marker)
                count=count+1
            elif background:
                if i in lines_to_label:  
                    plt.plot(x[plot_start:plot_end],y[i,plot_start:plot_end]-y[backg_line,plot_start:plot_end]+offset*count, label=labels+'{0:.2f}'.format(lines_label[i]),color=colors[count],marker=line_marker)
                else:
                    plt.plot(x[plot_start:plot_end],y[i,plot_start:plot_end]-y[backg_line,plot_start:plot_end]+offset*count,color=colors[count],marker=line_marker)
                count=count+1
            else:
                if i in lines_to_label:  
                    plt.plot(x[plot_start:plot_end],y[i,plot_start:plot_end]+offset*count, label=labels+'{0:.2f}'.format(lines_label[i]),color=colors[count],marker=line_marker)
                else:
                    plt.plot(x[plot_start:plot_end],y[i,plot_start:plot_end]+offset*count,color=colors[count],marker=line_marker)
                count=count+1
                
        else:

            if symmetrize:
                if i in lines_to_label:  
                    plt.plot(x,(y[i,:]+np.flip(y[i,:]))/2+offset*count, label=labels+'{0:.2f}'.format(lines_label[i]),color=colors[count],marker=line_marker)
                else:
                    plt.plot(x,(y[i,:]+np.flip(y[i,:]))/2+offset*count,color=colors[count],marker=line_marker)
                count=count+1
            elif background:
                
                if i in lines_to_label:  
                    plt.plot(x[:],y[i,:]-y[backg_line,:]+offset*count, label=labels+'{0:.2f}'.format(lines_label[i]),color=colors[count],marker=line_marker)
                else:
                    plt.plot(x[:],y[i,:]-y[backg_line,:]+offset*count, color=colors[count],marker=line_marker)
                count=count+1
            else:
                if i in lines_to_label:                    
                    plt.plot(x,y[i,:]+offset*count, label=labels+'{0:.2f}'.format(lines_label[i]),color=colors[count],marker=line_marker)
                else:
                    plt.plot(x,y[i,:]+offset*count,color=colors[count])
                count=count+1
    
    return


def int_vertplot(x, y, init, end, lines_to_plot, lines_label, 
                 label_int, offset, labels, colors, expand, interlocked,line_marker,
                 x_figsize,y_figsize):
    
    plt.figure(figsize=(x_figsize,y_figsize+expand*offset))
    
    lines_to_label=np.arange(init,end,label_int)
    
    count =0
    
    for i in lines_to_plot:
        if interlocked:
            if i in lines_to_label: 
                plt.plot(x[init:end],y[init:end,i]+offset*count, label=labels+'{0:.2f}'.format(lines_label[i]),color=colors[count],marker=line_marker)
            else:
                plt.plot(x[init:end],y[init:end,i]+offset*count,color=colors[count],marker=line_marker)
            count=count+1
        else:
            if i in lines_to_label: 
                plt.plot(x,y[:,i]+offset*count, label=labels+'{0:.2f}'.format(lines_label[i]),color=colors[count],marker=line_marker)
            else:
                plt.plot(x,y[:,i]+offset*count, color=colors[count],marker=line_marker)
                
            count=count+1
    
    return


def compare_int_linescans(filenames,data_dict,z_name,x_name,leg_name,I_name,dev_name, 
                       plot_start=0,plot_extend=10,plot_interval=1,offset=0,
                       label_int=1,
                       label_x='',label_y='', title_label='',
                       legend_label='',
                       x_min=0,x_max=0,y_min=0,y_max=0,
                       name_st=0,name_end=5,
                       save=False, conduc=False, eh=False,norm_current=False,
                       norm_quantum=False,log_bool=False,abs_bool=False,legend=True,
                       background=False,
                       expand=10, colors='jet',I_mult=1, font=15, fontleg=15,
                       hd=0, interlocked=False, symmetrize=False,
                       backg_line=0):
    
    count=0
    if plot_extend>len(data_dict):
        plot_extend=len(data_dict)
        
    plt.figure(figsize=(10,7+expand*offset))
    

    for i in np.arange(plot_start,plot_start+plot_extend,plot_interval):
        data=getattr(data_dict['data'+str(i)],z_name)
        
        x_ax=getattr(data_dict['data'+str(i)],x_name)
        #lines_label=getattr(data_dict['data'+str(i)],leg_name)
        current=getattr(data_dict['data'+str(i)],I_name)
        
        #if legend_label=='':
        legend_label=filenames[i][name_st:name_end]
        
        norm_Rxx=np.zeros_like(data)
        norm_Rxx_he=np.zeros_like(data)
        cond=np.zeros_like(data)
        cond_eh=np.zeros_like(data)
        
        
        for j in range(len(data)):
            norm_Rxx[j]=data[j]/current[j]/1000 # R in kOhm
            norm_Rxx_he[j]=norm_Rxx[j]/25.8128 #h/e2
            cond[j] = 1/(norm_Rxx[j])  #1/kOhm
            cond_eh[j] = 1/norm_Rxx_he[j]
            
        if norm_current==False and norm_quantum==False and conduc==False and eh==False:
            plt.plot(x_ax,data+offset*i,label=legend_label,color=colors[count])
                
        elif norm_current== True and norm_quantum==False and conduc==False and eh==False:
            plt.plot(x_ax,norm_Rxx+offset*i,label=legend_label ,color=colors[count])
                
        elif norm_quantum==True and conduc==False and eh==False:
            plt.plot(x_ax,norm_Rxx_he+offset*i,label=legend_label ,color=colors[count])
                
        elif conduc==True and eh==False:
            plt.plot(x_ax,cond+offset*i,label=legend_label ,color=colors[count])
        elif eh==True:
            plt.plot(x_ax,cond_eh+offset*i,label=legend_label ,color=colors[count])
                
        count=count+1*plot_interval
            
    if offset== 0.0:
        plt.title(title_label,fontsize=font)
    else:
        plt.title(title_label+' offseted by '+str(offset),fontsize=font)
            
    plt.tick_params(axis='both',which='major',labelsize=font)
    
    if label_x=='':
        plt.xlabel(x_name,fontsize=font)
    else:
        plt.xlabel(label_x,fontsize=font) 
    if label_y=='':
        #plt.locator_params(axis='x', nbins=15)
        plt.ylabel(z_name,fontsize=font)
        
        if norm_current== True and norm_quantum==False and conduc==False and eh==False:
            plt.ylabel(('R '+z_name+r'$ [k\Omega]$'),fontsize=font) 
            
        elif conduc==True and eh==False:
            plt.ylabel('G '+z_name+' mS',fontsize=font)
        elif eh==True:
            plt.ylabel(('G '+z_name+r'$ [e^2/h]$'),fontsize=font)   
        
        elif norm_quantum==True and conduc==False and eh==False:
            plt.ylabel(('R '+z_name+r'$ [h/e^2]$'),fontsize=font) 
    else:
        plt.xlabel(label_x,fontsize=font)
        #plt.locator_params(axis='x', nbins=15)
        plt.ylabel(label_y,fontsize=font)
        
    
    if x_max!=0:  
        plt.xlim(x_min,x_max)
    if y_max!=0:
        plt.ylim(y_min,y_max)
    
    plt.grid()
        
    if legend:
        plt.legend(loc='right',fontsize=fontleg)
            
    if save:
        plt.savefig(title_label++str(offset)+time.strftime("%H%M%S%d%m%y")+'.png', dpi=300)
    plt.show() 
    
    return

def shiftedColorMap(midpoint=0.5,cmap=orig_cmap, start=0,  stop=1.0,
                    division=258,name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, division)
    
    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, int(division/2), endpoint=False), 
        np.linspace(midpoint, 1.0, int(division/2), endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    #plt.register_cmap(cmap=newcmap)

    return newcmap

#%% OLD CODE
'''
Different iterations used to create the final code of the plot_lines and maps
'''

def export_datafigure(fig_data,column_raw_len,save_name,header_txt):
    savedata=np.zeros((len(fig_data),column_raw_len))
    for j in range(len(fig_data)):
        count=0
        for i in range(len(fig_data[0])):
            savedata[j][count*(len(fig_data[0][0])):count*(len(fig_data[0][0]))+len(fig_data[0][0])]=fig_data[j][i]
            count=count+1
    savedata=np.transpose(savedata)

    np.savetxt(save_name,savedata,delimiter='\t',header=header_txt)
    return 

def files_from_folder(directory, keyword):
    file_list = glob.glob(os.path.join(os.getcwd(),directory, keyword))
    for i in range(len(file_list)):
        file_list[i]= file_list[i].replace('\\','/')
        file_list[i]=file_list[i].split('/')
        file_list[i]=file_list[i][-1]

    return file_list

def know_scan_type(fname,skiprows=3):
    
    with open(fname) as myfile:
        
        head = [next(myfile) for x in range(skiprows)]
        if head[2][0:10]=='%linesweep':
            scan_type=0
        elif head[2][0:10]=='%megasweep':
            scan_type=1
        else:
            scan_type=2
            
    return scan_type

def scan_from_map_to_line(fname):
    
    filein=fname
    fileout=filein
    
    f = open(filein,'r')
    filedata = f.read()
    f.close()
    
    newdata = filedata.replace("megasweep","linesweep")
    
    f = open(fileout,'w')
    f.write(newdata)
    f.close()
    

#def clean_force_sweep_back(file,variable_to_sweep):
    
    
def fix_force_swept(file,forced_sweep_column=1,skiprows=3):
    
    with open (file) as myfile:
        head = [next(myfile) for x in range(3)]
    
        head_names = head[0].replace(' ','')
        head_names = head_names.replace('-','')
        #scan_type=head[2]
        #head_unit = head_names.replace('\n','')
        data = np.loadtxt(file,delimiter = '\t',skiprows=skiprows)
        head_txt=head[0]+head[1]+head[2]

    ##Find the things to delete

    counter=0
    
    for i in range(2,len(data[:,forced_sweep_column])):
        if data[i,forced_sweep_column]-data[i-forced_sweep_column,forced_sweep_column]<data[1,forced_sweep_column]-data[0,forced_sweep_column]:
            stop=i
            
            break
    for i in range(stop,len(data[:,forced_sweep_column])):
        if data[i,forced_sweep_column]-data[i+1,forced_sweep_column]<data[stop,forced_sweep_column]-data[stop-1,forced_sweep_column]:
            start=i
            break    
    break_point=0
    for i in range(2,len(data[:,0])):
        if data[i,0]-data[i-1,0]!=data[1,0]-data[0,0]:
            break_point=break_point+1
    break_point=break_point
    
    ### Make the final data
    data_real=[]
    real_length=len(data[:,0])-(break_point+1)*(start-stop)
    each_length=int(real_length/(break_point+1))
    
    data_real_small=np.zeros((real_length))
    
    for i in range(len(data[0])):
        data_real.append(data_real_small)
        
    for i in range(len(data[0])):
        counter=0
        counter2=0
        data_real_small=np.zeros((real_length))
        for j in range(break_point+1):
            data_real_small[counter2*(each_length):each_length+counter2*(each_length)]=data[counter:counter+each_length,i]
            counter=(j+1)*(start+1)
            counter2=counter2+1
        
        data_real[i]=data_real_small
        
    ## Save the data again in a text file 
    Data = np.transpose(data_real)
    
    np.savetxt(file[:-4]+'_fixed.txt',Data, delimiter='\t', header = head_txt, comments='')

def parse_and_plot(directory):
    
    files=files_from_folder(directory, '*hdf5*')
    files_txt=files_from_folder(directory, '*txt*')


    line_objdata = {} #stores all the data
    line_rawdata ={} #store the xaxis data for plotting
    line_variables={} #store the data to plot in the right units
    map_objdata = {} #stores all the data
    map_rawdata ={} #store the xaxis data for plotting
    map_variables={} #store the data to plot in the right units
    objdata_mapshape={}
    objdata_mapshape_real={}

    map_counter=0
    line_counter=0
    
    for i in range(0,len(files)):
    
        filetxt=files[i][:-5]+'.txt'
        
        if filetxt not in files_txt or os.stat(directory+filetxt).st_size==0:
            scan_type=exportDataToText_v6.dataExportTxt(directory,files[i])
        else: 
            scan_type=know_scan_type(directory+filetxt)
    
        #print(files[i])
        if scan_type==0:
            print(line_counter)
            line_objdata['data'+str(line_counter)],line_rawdata['data'+str(line_counter)], line_variables['data'+str(line_counter)]=parse_data(directory+filetxt)
    
            for k in range(len(line_variables['data'+str(line_counter)])-1):
                plot_linescans_raw(line_objdata['data'+str(line_counter)],line_variables['data'+str(line_counter)][0],line_variables['data'+str(line_counter)][k],files[i][:-5])
            line_counter=line_counter+1
            
        elif scan_type==1:
            map_objdata['data'+str(map_counter)],map_rawdata['data'+str(map_counter)], map_variables['data'+str(map_counter)]=parse_data(directory+filetxt)
    
            for k in range(len(map_variables['data'+str(map_counter)])-1):
            #if scan_type=='megasweep\n':
                mapdata=map_rawdata['data'+str(map_counter)]
                mapnames=map_variables['data'+str(map_counter)]
    
                objdata_mapshape['data'+str(map_counter)],objdata_mapshape_real['data'+str(map_counter)]=plot_maps_v3(mapdata,mapnames,map_variables['data'+str(map_counter)][1],map_variables['data'+str(map_counter)][0],map_variables['data'+str(map_counter)][k],files[i][:-5])
            map_counter=map_counter+1
        elif scan_type==2:
            print('three inner loops, need to rethink the data shape')

    return line_objdata,line_rawdata,line_variables,map_objdata,map_rawdata,map_variables,objdata_mapshape,objdata_mapshape_real

def parse_and_plot_lines(directory):
    
    files=files_from_folder(directory, '*hdf5*')
    files_txt=files_from_folder(directory, '*txt*')


    line_objdata = {} #stores all the data
    line_rawdata ={} #store the xaxis data for plotting
    line_variables={} #store the data to plot in the right units
    filenames={}

    line_counter=0
    for i in range(0,len(files)):
        try:
            filename=files[i][:-5]
            filetxt=files[i][:-5]+'.txt'
            filenames['data'+str(line_counter)]=filename
            if filetxt not in files_txt or os.stat(directory+filetxt).st_size==0:
                scan_type=exportDataToText_v6.dataExportTxt(directory,files[i])
            else: 
                scan_type=know_scan_type(directory+filetxt)
        
            #print(files[i])
            if scan_type==0:
                print(line_counter)
                line_objdata['data'+str(line_counter)],line_rawdata['data'+str(line_counter)], line_variables['data'+str(line_counter)]=parse_data(directory+filetxt)
        
                for k in range(len(line_variables['data'+str(line_counter)])-1):
                    plot_linescans_raw(line_objdata['data'+str(line_counter)],line_variables['data'+str(line_counter)][0],line_variables['data'+str(line_counter)][k],filename)
                line_counter=line_counter+1
                
                
            elif scan_type==1:
                continue
            elif scan_type==2:
                print('three inner loops, need to rethink the data shape')
                
        except FileNotFoundError or IndexError:
            print(files[i]+'cannot be extracted. Check if it might be empty or data is corrupted in some way')

    return  line_objdata,line_rawdata,line_variables,filenames

def parse_and_plot_maps(directory):
    
    files=files_from_folder(directory, '*hdf5*')
    files_txt=files_from_folder(directory, '*txt*')



    map_objdata = {} #stores all the data
    map_rawdata ={} #store the xaxis data for plotting
    map_variables={} #store the data to plot in the right units
    objdata_mapshape={}
    objdata_mapshape_real={}

    map_counter=0

    
    for i in range(0,len(files)):
        try:
            filetxt=files[i][:-5]+'.txt'
            
            if filetxt not in files_txt or os.stat(directory+filetxt).st_size==0:
                scan_type=exportDataToText_v6.dataExportTxt(directory,files[i])
            else: 
                scan_type=know_scan_type(directory+filetxt)
        
            #print(files[i])
            if scan_type==0:
                continue
                
            elif scan_type==1:
                map_objdata['data'+str(map_counter)],map_rawdata['data'+str(map_counter)], map_variables['data'+str(map_counter)]=parse_data(directory+filetxt)

                if map_rawdata['data'+str(map_counter)][0,0]==map_rawdata['data'+str(map_counter)][-1,0]:
                    print(files[i]+' is actually a linescan. It will be rewritten as linesweep so it does not appear anymore in the maps')
                    scan_from_map_to_line(directory+filetxt)
                    for k in range(len(map_variables['data'+str(map_counter)])-1):
                        plot_linescans_raw(map_objdata['data'+str(map_counter)],map_variables['data'+str(map_counter)][1],map_variables['data'+str(map_counter)][k],files[i][:-5])
                else:  
                    for k in range(len(map_variables['data'+str(map_counter)])-1):
                #if scan_type=='megasweep\n':
                        mapdata=map_rawdata['data'+str(map_counter)]
                        mapnames=map_variables['data'+str(map_counter)]
                        objdata_mapshape['data'+str(map_counter)],objdata_mapshape_real['data'+str(map_counter)]=plot_maps_v3(mapdata,mapnames,map_variables['data'+str(map_counter)][1],map_variables['data'+str(map_counter)][0],map_variables['data'+str(map_counter)][k],files[i][:-5])

                map_counter=map_counter+1
            elif scan_type==2:
                print('three inner loops, need to rethink the data shape')
                
        except FileNotFoundError or IndexError as e:
            print(str(e)+' '+files[i]+' cannot be extracted. Check if it might be empty or data is corrupted in some way')
    return objdata_mapshape,map_variables,objdata_mapshape_real

def parse_and_plot_map_fromfiles(directory,file):

    filetxt=file[:-5]+'.txt'
        
    if filetxt not in file or os.stat(directory+filetxt).st_size==0:
        scan_type=exportDataToText_v6.dataExportTxt(directory,file)
    else: 
        scan_type=know_scan_type(directory+filetxt)
        
    if scan_type==0:
        print(file+'is a linescan')
    elif scan_type==1:
        map_objdata,map_rawdata, map_variables=parse_data(file)
        
        for k in range(len(map_variables)-1):
            #if scan_type=='megasweep\n':
            mapdata=map_rawdata
            mapnames=map_variables
        
            objdata_mapshape,objdata_mapshape_real=plot_maps_v3(mapdata,mapnames,map_variables[1],map_variables[0],map_variables[k],file)
    else: 
        print('Data type not understood, check your data format')

    return objdata_mapshape,map_variables,objdata_mapshape_real

    
