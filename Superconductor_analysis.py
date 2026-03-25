# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 12:30:31 2021

Functions for extracting superconducting critical current from minimums and maximums in Fraunhofer pattern plots in 
MATBG devices. It also includes an analysis to extract the critical current from the derivative of Ic.
In the end there are some basic simulations of Fraunhofer patterns with different magnetic signals.
@author: Jaime Diez @LDQM-LS-Efetov group at ICFO-LMU
"""

import numpy as np
import matplotlib.pyplot as plt
import importlib
import glob
import os
import matplotlib.cm as cm
import math as math
from scipy import signal
from scipy.fftpack import fft, ifft
from scipy.signal import savgol_filter
from collections import Counter
from ipywidgets import interactive
import time 
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib
from homemade_functions import runningMeanFast 


savefolder='C:/Users/jdiez/LDQM/a. Data/Ju1/JJ_BF_2021/'



def get_Ic_minimum(x,y,z,limit,start,end,middle,z_max=5,x_lim=50,y_lim=100,save=False,font=15,window=7,
                   use_min=False,limit2=0,plotting=True,correct_middle=False,normalize=False,markers=5,
                   mycolor='RdBu_r'):
    
    '''
    Get the Ic from a value given as a limit from the minimum, you would define
    Ic as the point in current at which the resistance is no longer 0, giving
    this 0 a margin of error. For Ju1 this is a value of 100 Ohms. The extraction
    needs to take into account also the fact that at higher fields, the SC
    is not longer 0 resistance, so there will be a factor which accounts
    for the increase of zero value resistance. 
    x, y and z: variables which make the Fraunhofer pattern 
    start, end: limits the position of the extraction of the Ic, 
                useful in case there are several transitions in the data 
    middle: zero current position of the dVdI
            The positive values are taken between middle and end while
            the negative values are taken between start and middle. 
    z_maz, x_lim, y_lim: boundaries of the plots 
    normalize: will plot divided by the maximum Ic 
    markers: size of the markers 
    '''
    
    #Extract it from a minimum value
    Ic_pos=np.zeros_like(x[:,0])
    Ic_neg=np.zeros_like(Ic_pos)
    middle_large=np.zeros_like(Ic_pos)
    
    count=0
    
    # Get the given Ic and save it in the predefined values 
    if use_min==False:
        for i in range(len(x)):
            min_value=np.amin(z[i,start:end])
            for j in range(middle,end):
                if z[i,j]>=min_value+limit:
                    Ic_pos[count]=x[i,j]
                    break
            for j in range(middle,start,-1):
                if z[i,j]>=min_value+limit:
                    Ic_neg[count]=x[i,j]
                    break
            count=count+1
    else:
        if correct_middle:
            for i in range(len(x)):
                min_value=np.amin(z[i,start:end])
                for j in range(start,end):
                    if z[i,j]==min_value:
                        middle=j
                        middle_large[count]=x[i,middle]
                        break
                for j in range(middle,end):
                    if z[i,j]>=min_value+limit:
                        if x[i,middle]>0:
                            Ic_pos[count]=x[i,j]-x[i,middle]
                        else:
                            Ic_pos[count]=x[i,j]+x[i,middle]
                        break
                for j in range(middle,start,-1):
                    if z[i,j]>=min_value+limit2:
                        if x[i,middle]>0:
                            Ic_neg[count]=x[i,j]+x[i,middle]
                        else:
                            Ic_neg[count]=x[i,j]-x[i,middle]
                        break
                count=count+1
        else:
            for i in range(len(x)):
                middle_range=[]
                min_value=np.amin(z[i,start:end])
                for j in range(start,end):
                    try:
                        if z[i,j]==min_value:
                            middle_range.append(j)
                    except IndexError:
                        continue
                #print(middle_range)
                if len(middle_range)>2:
                    middle=int((middle_range[0]+middle_range[-1])/2)
                    #middle=j
                    middle_large[count]=x[i,middle]
                else:
                    middle=middle_range[0]
                    middle_large[count]=x[i,middle]
                    #break
                for j in range(middle,end):
                    if z[i,j]>=min_value+limit:
                        Ic_pos[count]=x[i,j]
                        break
                for j in range(middle,start,-1):
                    if z[i,j]>=min_value+limit2:
                        Ic_neg[count]=x[i,j]
                        break
                count=count+1
    if normalize:
        for i in range(len(x)):
            Ic_pos=Ic_pos/np.amax(Ic_pos)
            Ic_neg=Ic_neg/(np.amax(-Ic_neg))
    #print(middle_large)
    D_Ic=Ic_pos-abs(Ic_neg) 

    Ic_pos_sm = savgol_filter(Ic_pos, window, 3) # window size 51, polynomial order 3
    Ic_neg_sm = savgol_filter(Ic_neg, window, 3) # window size 51, polynomial order 3

    if plotting:
        # Now plot the different figures 

        ### Orginal map with the fitting 

        fig = plt.figure(figsize=(12,6))
        ax = plt.subplot(111)
        im = plt.pcolormesh(y,x,z,cmap=mycolor,vmax=z_max)

        #ax.plot(y[:,0],Ic_pos_sm,'k--',label='$I_c^+(B)$',alpha=0.7)
        #ax.plot(y[:,0],Ic_neg_sm,'k--',label='$|I_c^-(B)|$',alpha=0.7)

        ax.plot(y[:,0],Ic_pos,'k--',label='$I_c^+(B)$',alpha=0.7)
        ax.plot(y[:,0],Ic_neg,'b--',label='$|I_c^-(B)|$',alpha=0.7)
        #ax.plot(y[:,0],middle_large,'k--',label='$|I_c^-(B)|$',alpha=0.7)
        
        ax.set_ylabel('$I_{DC}$ (nA)',fontsize=font)
        ax.set_xlabel('$B$ (mT)',fontsize=font)

        #ax.set_xticks(ticks)
        #ax.set_yticks(ticks)

        ax.tick_params(axis='both', which='major', labelsize=font)

        ### Add the colorbar 
        cbaxes = fig.add_axes([0.75, 0.90, 0.15,0.015]) #left,bottom,width,heigth
        cb = plt.colorbar(im,shrink=0.3,cax=cbaxes,orientation='horizontal') 

        #cb=plt.colorbar(pad=0)
        cb.set_label('$dV_{xx}/dI$ (k$\Omega$)',fontsize=font,position=(-0.5,1), labelpad=-15,rotation=0)
        cbaxes.xaxis.tick_top()
        cb.ax.tick_params(labelsize=font)

        ax.set_ylim(-y_lim,y_lim)
        ax.set_xlim(-x_lim,x_lim)
        #ax.set_title('I:7-18 (1 nA), $dV_{xx}/dI$: 5-1, $V_{BG}:-1.12$ V, $V_{TGs}:-147$ mV',fontsize=font,position=(0.4,0.92),color='k')
        #ax.annotate('B ', xy=(30,-80),
        #            xytext=(50, -80),
        #            arrowprops=dict(facecolor='black',width=0.5,headwidth=9, shrink=0.05),
        #            horizontalalignment='right', verticalalignment='top',fontsize=18
        #            )
        plt.show()


        plt.figure(figsize=(8,5))
        plt.plot(y[:,0],Ic_pos,'bo--',label='$I_c^+(B)$',markersize=markers)
        plt.plot(y[:,0],-Ic_neg,'ko--',label='$|I_c^-(B)|$',markersize=markers)
        plt.legend(loc='upper center', bbox_to_anchor=(0.2, 1), shadow=False, ncol=1,fontsize=font,frameon=False)
        plt.xlim(-x_lim,x_lim)
        plt.xlabel('$B$ (mT)',fontsize=font)
        plt.ylabel('$I_c$ (nA)',fontsize=font)
        plt.tick_params(axis='both',which='major',labelsize=font,direction='in')
        plt.show()

        plt.figure(figsize=(8,5))
        plt.plot(y[:,0],Ic_pos_sm,'bo--',label='$I_c^+(B)$',markersize=markers)
        plt.plot(y[:,0],-Ic_neg_sm,'ko--',label='$|I_c^-(B)|$',markersize=markers)
        plt.title('Compare the smooth data',fontsize=font)
        plt.legend(loc='upper center', bbox_to_anchor=(0.2, 1), shadow=False, ncol=1,fontsize=font,frameon=False)
        plt.xlim(-x_lim,x_lim)
        plt.xlabel('$B$ (mT)',fontsize=font)
        plt.ylabel('$I_c$ (nA)',fontsize=font)
        plt.tick_params(axis='both',which='major',labelsize=font,direction='in')
        plt.show()

        plt.figure(figsize=(8,5))
        plt.plot(y[:,0],Ic_pos,'bo--',label='$I_c^+(B)$',markersize=markers)
        plt.plot(-y[:,0],-Ic_neg,'ko--',label='$|I_c^-(-B)|$',markersize=markers)
        plt.title('Invert one of them',fontsize=font)
        plt.legend(loc='upper center', bbox_to_anchor=(0.2, 1), shadow=False, ncol=1,fontsize=font,frameon=False)
        plt.xlim(-x_lim,x_lim)
        plt.xlabel('$B$ (mT)',fontsize=font)
        plt.ylabel('$I_c$ (nA)',fontsize=font)
        plt.tick_params(axis='both',which='major',labelsize=font,direction='in')
        plt.show()

        plt.figure(figsize=(8,5))
        plt.plot(y[:,0],Ic_pos_sm,'bo--',label='$I_c^+(B)$',markersize=markers)
        plt.plot(-y[:,0],-Ic_neg_sm,'ko--',label='$|I_c^-(-B)|$',markersize=markers)
        plt.title('Invert one of them',fontsize=font)
        plt.legend(loc='upper center', bbox_to_anchor=(0.2, 1), shadow=False, ncol=1,fontsize=font,frameon=False)
        plt.xlim(-x_lim,x_lim)
        plt.xlabel('$B$ (mT)',fontsize=font)
        plt.ylabel('$I_c$ (nA)',fontsize=font)
        plt.tick_params(axis='both',which='major',labelsize=font,direction='in')
        plt.show()

    return Ic_pos,Ic_neg,D_Ic

def get_Ic_maximums(x,y,z,middle,limits,max_mult=5,mycolor='RdBu_r',
                   zmax=5):
    '''
    Extract the position of the maximums of the dVdI. 
    x,y and z make a Fraunhofer pattern
    middle: value of the 0 current
    limits: limit the max and minimum considered currents. Useful if there
            are multiple transitions 
    max_mult: random number used to see if the shape of the extracted data
              could martch the actual position of the real Ic 
    '''
    Icnegmax=[]
    Icposmax=[]
    Icmaxypos=[]
    Icmaxyneg=[]

    start=limits
    end=-1*limits
    for i in range(len(x[:,0])):
        Icnegmax.append(np.amax(z[i,start:middle]))
        Icposmax.append(np.amax(z[i,middle:end]))
    for i in range(len(x)):
        for j in range(len(x[i])):
            if z[i][j]==Icposmax[i]:
                Icmaxypos.append(y[i,j])
                continue
            elif z[i][j]==Icnegmax[i]:
                Icmaxyneg.append(y[i,j])
                continue
                
    plt.figure(figsize=(10,5))
    plt.pcolormesh(x,y,z,cmap=mycolor,vmax=zmax)
    plt.plot(x[:,0],np.asarray(Icposmax)*max_mult,'b-',)
    plt.plot(x[:,0],np.asarray(Icnegmax)*-max_mult,'r-')
    plt.colorbar()
    plt.show()
    return Icposmax, Icnegmax, Icmaxypos,Icmaxyneg 


def get_Ic_minimum_good(x,y,z,limit,start,end,middle,z_max=5,x_lim=50,y_lim=100,save=False,font=15,window=7,
                   use_min=False,limit2=0,limit3=0,plotting=True,correct_middle=False,normalize=False,
                   smooth=False,mycolor='RdBu_r',markers=5,interp=True,save_format='png'):
    
    #Extract it from a minimum value
    Ic_pos=np.zeros_like(x[:,0])
    Ic_neg=np.zeros_like(Ic_pos)
    middle_large=np.zeros_like(Ic_pos)
    
    if interp==True:
        new_x,new_y,new_z=map_interpolation(x,y,z,plot_map=False)
        x=new_x
        y=new_y
        z=new_z
        
    count=0
    
    # Get the given Ic and save it in the predefined values 
    if use_min==False:
        for i in range(len(x)):
            min_value=np.amin(z[i,start:end])
            for j in range(middle,end):
                if z[i,j]>=min_value+limit:
                    Ic_pos[count]=x[i,j]
                    break
            for j in range(middle,start,-1):
                if z[i,j]>=min_value+limit:
                    Ic_neg[count]=x[i,j]
                    break
            count=count+1
    else:
        if correct_middle:
            for i in range(len(x)):
                min_value=np.amin(z[i,start:end])
                for j in range(start,end):
                    if z[i,j]==min_value:
                        middle=j
                        middle_large[count]=x[i,middle]
                        break
                for j in range(middle,end):
                    '''
                    if min_value>limit3:
                        Ic_pos[count]=x[i,middle]
                    '''
                    if z[i,j]>=min_value+limit:
                        if x[i,middle]>0:
                            Ic_pos[count]=x[i,j]-x[i,middle]
                        else:
                            Ic_pos[count]=x[i,j]+x[i,middle]
                        break
                for j in range(middle,start,-1):
                    '''
                    if min_value>limit3:
                        Ic_neg[count]=x[i,middle]
                    '''
                    if z[i,j]>=min_value+limit2:
                        if x[i,middle]>0:
                            Ic_neg[count]=x[i,j]+x[i,middle]
                        else:
                            Ic_neg[count]=x[i,j]-x[i,middle]
                        break
                count=count+1
        else:
            for i in range(len(x)):
                middle_range=[]
                min_value=np.amin(z[i,start:end])
                for j in range(start,end):
                    try:
                        if z[i,j]==min_value:
                            middle_range.append(j)
                    except IndexError:
                        continue
                #print(middle_range)
                if len(middle_range)>2:
                    middle=int((middle_range[0]+middle_range[-1])/2)
                    #middle=j
                    middle_large[count]=x[i,middle]
                else:
                    middle=middle_range[0]
                    middle_large[count]=x[i,middle]
                    #break
                for j in range(middle,end):
                    '''
                    if min_value>limit3:
                        Ic_pos[count]=x[i,middle]
                    '''
                    if z[i,j]>=min_value+limit:
                        Ic_pos[count]=x[i,j]
                        break
                for j in range(middle,start,-1):
                    '''
                    if min_value>limit3:
                        Ic_neg[count]=x[i,middle]
                    '''
                    if z[i,j]>=min_value+limit2:
                        Ic_neg[count]=x[i,j]
                        break
                count=count+1
    if normalize:
        for i in range(len(x)):
            Ic_pos=Ic_pos/np.amax(Ic_pos)
            Ic_neg=Ic_neg/(np.amax(-Ic_neg))
    #print(middle_large)
    D_Ic=Ic_pos-abs(Ic_neg) 
    

    if plotting:
        # Now plot the different figures 

        ### Orginal map with the fitting 

        fig = plt.figure(figsize=(12,6))
        ax = plt.subplot(111)
        im = plt.pcolormesh(y,x,z,cmap=mycolor,vmax=z_max)

        #ax.plot(y[:,0],Ic_pos_sm,'k--',label='$I_c^+(B)$',alpha=0.7)
        #ax.plot(y[:,0],Ic_neg_sm,'k--',label='$|I_c^-(B)|$',alpha=0.7)

        ax.plot(y[:,0],Ic_pos,'g--',label='$I_c^+(B)$',alpha=0.7)
        ax.plot(y[:,0],Ic_neg,'k--',label='$|I_c^-(B)|$',alpha=0.7)
        #ax.plot(y[:,0],middle_large,'k--',label='$|I_c^-(B)|$',alpha=0.7)
        
        ax.set_ylabel('$I_{DC}$ (nA)',fontsize=font)
        ax.set_xlabel('$B$ (mT)',fontsize=font)

        #ax.set_xticks(ticks)
        #ax.set_yticks(ticks)

        ax.tick_params(axis='both', which='major', labelsize=font)

        ### Add the colorbar 
        cbaxes = fig.add_axes([0.75, 0.90, 0.15,0.015]) #left,bottom,width,heigth
        cb = plt.colorbar(im,shrink=0.3,cax=cbaxes,orientation='horizontal') 

        #cb=plt.colorbar(pad=0)
        cb.set_label('$dV_{xx}/dI$ (k$\Omega$)',fontsize=font,position=(-0.5,1), labelpad=-15,rotation=0)
        cbaxes.xaxis.tick_top()
        cb.ax.tick_params(labelsize=font)

        ax.set_ylim(-y_lim,y_lim)
        ax.set_xlim(-x_lim,x_lim)
        #ax.set_title('I:7-18 (1 nA), $dV_{xx}/dI$: 5-1, $V_{BG}:-1.12$ V, $V_{TGs}:-147$ mV',fontsize=font,position=(0.4,0.92),color='k')
        #ax.annotate('B ', xy=(30,-80),
        #            xytext=(50, -80),
        #            arrowprops=dict(facecolor='black',width=0.5,headwidth=9, shrink=0.05),
        #            horizontalalignment='right', verticalalignment='top',fontsize=18
        #            )
        plt.show()


        plt.figure(figsize=(8,5))
        plt.plot(y[:,0],Ic_pos,'bo--',label='$I_c^+(B)$',markersize=markers)
        plt.plot(y[:,0],-Ic_neg,'ko--',label='$|I_c^-(B)|$',markersize=markers)
        plt.legend(loc='upper center', bbox_to_anchor=(0.2, 1), shadow=False, ncol=1,fontsize=font,frameon=False)
        plt.xlim(-x_lim,x_lim)
        plt.xlabel('$B$ (mT)',fontsize=font)
        plt.ylabel('$I_c$ (nA)',fontsize=font)
        plt.tick_params(axis='both',which='major',labelsize=font,direction='in')
        plt.show()
        
        plt.figure(figsize=(8,5))
        plt.plot(y[:,0],Ic_pos,'bo--',label='$I_c^+(B)$',markersize=markers)
        plt.plot(-y[:,0],-Ic_neg,'ko--',label='$|I_c^-(-B)|$',markersize=markers)
        plt.title('Invert one of them',fontsize=font)
        plt.legend(loc='upper center', bbox_to_anchor=(0.2, 1), shadow=False, ncol=1,fontsize=font,frameon=False)
        plt.xlim(-x_lim,x_lim)
        plt.xlabel('$B$ (mT)',fontsize=font)
        plt.ylabel('$I_c$ (nA)',fontsize=font)
        plt.tick_params(axis='both',which='major',labelsize=font,direction='in')
        plt.show()

    if smooth:  
        #D_Icm=abs(Ic_neg) -Ic_pos
    
        Ic_pos_sm = savgol_filter(Ic_pos, window, 3) # window size 51, polynomial order 3
        Ic_neg_sm = savgol_filter(Ic_neg, window, 3) # window size 51, polynomial order 3
        
        plt.figure(figsize=(8,5))
        plt.plot(y[:,0],Ic_pos_sm,'bo--',label='$I_c^+(B)$',markersize=markers)
        plt.plot(y[:,0],-Ic_neg_sm,'ko--',label='$|I_c^-(B)|$',markersize=markers)
        plt.title('Compare the smooth data',fontsize=font)
        plt.legend(loc='upper center', bbox_to_anchor=(0.2, 1), shadow=False, ncol=1,fontsize=font,frameon=False)
        plt.xlim(-x_lim,x_lim)
        plt.xlabel('$B$ (mT)',fontsize=font)
        plt.ylabel('$I_c$ (nA)',fontsize=font)
        plt.tick_params(axis='both',which='major',labelsize=font,direction='in')
        plt.show()
        
        plt.figure(figsize=(8,5))
        plt.plot(y[:,0],Ic_pos_sm,'go--',label='$I_c^+(B)$',markersize=markers)
        plt.plot(-y[:,0],-Ic_neg_sm,'ro--',label='$|I_c^-(-B)|$',markersize=markers)
        plt.title('Invert one of them',fontsize=font)
        plt.legend(loc='upper center', bbox_to_anchor=(0.2, 1), shadow=False, ncol=1,fontsize=font,frameon=False)
        plt.xlim(-x_lim,x_lim)
        plt.xlabel('$B$ (mT)',fontsize=font)
        plt.ylabel('$I_c$ (nA)',fontsize=font)
        plt.tick_params(axis='both',which='major',labelsize=font,direction='in')
        plt.show()
        
        
    return Ic_pos,Ic_neg,D_Ic

def map_interpolation(x,y,z,interpx_size=1001, interpy_size=601,
                      interp_y=False,mycolor='RdBu_r',plot_map=True):
    ''' Interpolate the data to get many more points 
        The interpolation in Idc direction is useful because when
        extracting the exact point we need more points, otherwise
        the code will jump to the next close point, which in general
        will give errors later, as in 'bad data'
        The interpolation has to go in ascending x order 
        Interp y is usually False because it does not add anything
    '''
    new_x=np.linspace(x[0,0],x[0,-1],interpx_size)
    new_data=np.zeros((len(z),len(new_x))) 
    new_y_large=np.zeros((len(z),len(new_x)))
    new_x_large=np.zeros((len(z),len(new_x)))
    
    if x[0,0]<x[0,-1]: #Ascending order
        for i in range(len(z)):
            new_data[i]=np.interp(new_x,x[i],z[i])
            new_y_large[i]=y[i,0]

    else:
        for i in range(len(z)):
            new_data[i]=np.interp(np.flip(new_x),np.flip(x[i]),np.flip(z[i]))
            new_y_large[i]=y[i,0]

        new_data=np.flip(new_data) 
        
    for i in range(len(new_data)):
        new_x_large[i]=new_x  
    ## Now interpolate the field 
    if interp_y: 
        
        new_y=np.linspace(y[-1,0],y[0,0],interpy_size)
        
        new_dataB=np.zeros((len(new_y),len(new_x))) 
        new_y_largeB=np.zeros((len(new_y),len(new_x)))
        new_x_largeB=np.zeros((len(new_y),len(new_x)))
        for i in range(0,len(new_data[0])):
            new_dataB[:,i]=np.interp(new_y,new_y_large[:,i],new_data[:,i])
        
        for i in range(0,len(new_data[0])):
            #new_dataB[:,i]=np.flip(new_dataB[:,i])
            new_y_largeB[:,i]=new_y
        
        for i in range(len(new_dataB)):
            new_x_largeB[i]=new_x
    else:
        new_x_largeB=new_x_large
        new_y_largeB=new_y_large
        new_dataB=new_data

    if plot_map:
        plt.figure(figsize=(5,5))
        plt.pcolormesh(new_x_largeB,new_y_largeB,new_dataB,cmap=mycolor)
        plt.show()    
    
    return new_x_largeB,new_y_largeB,new_dataB

def parula_cmap():
    cm_data = [[0.2081, 0.1663, 0.5292], [0.2116238095, 0.1897809524, 0.5776761905], 
     [0.212252381, 0.2137714286, 0.6269714286], [0.2081, 0.2386, 0.6770857143], 
     [0.1959047619, 0.2644571429, 0.7279], [0.1707285714, 0.2919380952, 
      0.779247619], [0.1252714286, 0.3242428571, 0.8302714286], 
     [0.0591333333, 0.3598333333, 0.8683333333], [0.0116952381, 0.3875095238, 
      0.8819571429], [0.0059571429, 0.4086142857, 0.8828428571], 
     [0.0165142857, 0.4266, 0.8786333333], [0.032852381, 0.4430428571, 
      0.8719571429], [0.0498142857, 0.4585714286, 0.8640571429], 
     [0.0629333333, 0.4736904762, 0.8554380952], [0.0722666667, 0.4886666667, 
      0.8467], [0.0779428571, 0.5039857143, 0.8383714286], 
     [0.079347619, 0.5200238095, 0.8311809524], [0.0749428571, 0.5375428571, 
      0.8262714286], [0.0640571429, 0.5569857143, 0.8239571429], 
     [0.0487714286, 0.5772238095, 0.8228285714], [0.0343428571, 0.5965809524, 
      0.819852381], [0.0265, 0.6137, 0.8135], [0.0238904762, 0.6286619048, 
      0.8037619048], [0.0230904762, 0.6417857143, 0.7912666667], 
     [0.0227714286, 0.6534857143, 0.7767571429], [0.0266619048, 0.6641952381, 
      0.7607190476], [0.0383714286, 0.6742714286, 0.743552381], 
     [0.0589714286, 0.6837571429, 0.7253857143], 
     [0.0843, 0.6928333333, 0.7061666667], [0.1132952381, 0.7015, 0.6858571429], 
     [0.1452714286, 0.7097571429, 0.6646285714], [0.1801333333, 0.7176571429, 
      0.6424333333], [0.2178285714, 0.7250428571, 0.6192619048], 
     [0.2586428571, 0.7317142857, 0.5954285714], [0.3021714286, 0.7376047619, 
      0.5711857143], [0.3481666667, 0.7424333333, 0.5472666667], 
     [0.3952571429, 0.7459, 0.5244428571], [0.4420095238, 0.7480809524, 
      0.5033142857], [0.4871238095, 0.7490619048, 0.4839761905], 
     [0.5300285714, 0.7491142857, 0.4661142857], [0.5708571429, 0.7485190476, 
      0.4493904762], [0.609852381, 0.7473142857, 0.4336857143], 
     [0.6473, 0.7456, 0.4188], [0.6834190476, 0.7434761905, 0.4044333333], 
     [0.7184095238, 0.7411333333, 0.3904761905], 
     [0.7524857143, 0.7384, 0.3768142857], [0.7858428571, 0.7355666667, 
      0.3632714286], [0.8185047619, 0.7327333333, 0.3497904762], 
     [0.8506571429, 0.7299, 0.3360285714], [0.8824333333, 0.7274333333, 0.3217], 
     [0.9139333333, 0.7257857143, 0.3062761905], [0.9449571429, 0.7261142857, 
      0.2886428571], [0.9738952381, 0.7313952381, 0.266647619], 
     [0.9937714286, 0.7454571429, 0.240347619], [0.9990428571, 0.7653142857, 
      0.2164142857], [0.9955333333, 0.7860571429, 0.196652381], 
     [0.988, 0.8066, 0.1793666667], [0.9788571429, 0.8271428571, 0.1633142857], 
     [0.9697, 0.8481380952, 0.147452381], [0.9625857143, 0.8705142857, 0.1309], 
     [0.9588714286, 0.8949, 0.1132428571], [0.9598238095, 0.9218333333, 
      0.0948380952], [0.9661, 0.9514428571, 0.0755333333], 
     [0.9763, 0.9831, 0.0538]]
                                               
    parula_map = LinearSegmentedColormap.from_list('parula', cm_data)
    return parula_map 

def derivative(x,y):
    '''
    Makes the derivative of a function y. 
    Make sure the x and y have the right dimensions of 2D arrays 
    '''
    dydx=[]
    for i in range(len(x[:,0])):
        dydx.append(np.diff(y[i]))
    return np.asarray(dydx)

def derivative_dvdi(x,y,z):
    '''
    Extract the Ic of a dvdi by doing the derivative of it 
    and taking the maximum values of the derivative as the Ic. 
    The Ic plus will be the maximum and the Ic minum will be the minimum
    Does not work for all the samples, depends on the shape of the dvdi
    x, y and z should be 3d arrays making a map
    dvdi is the derivative
    dvdipos and neg are the max and minimum
    Ic pos and neg will be the position in I of those minimum and maximum
    
    '''
    dydx=[]
    dydxpos=[]
    dydxneg=[]
    Icpos=[]
    Icneg=[]
    for i in range(len(x[:,0])):
        dydx.append(np.diff(z[i]))
        dydxpos.append(np.amax(dydx[i]))
        dydxneg.append(np.amin(dydx[i]))
    '''
    Extract the position in y of the derivative. 
    This works very bad, I think because of the few points in that direction
    '''
    for i in range(len(dydx)):
        for j in range(len(dydx[i])):
            if dydx[i][j]==dydxpos[i]:
                Icpos.append(y[i,j])
                continue
            elif dydx[i][j]==dydxneg[i]:
                Icneg.append(y[i,j])
                continue
    plot_derivative_map(x,y,z,dydxpos,dydxneg)
    
    return dydx,dydxpos,dydxneg

def plot_derivative_map(x,y,z,dvdipos,dvdineg,mycolor='RdBu_r'):
    '''
    Plot the extracted derivative values to see how it looks 
    '''
    
    plt.figure(figsize=(10,5))
    plt.pcolormesh(x,y,z,cmap=mycolor)
    plt.plot(x[:,0],dvdipos)
    plt.plot(x[:,0],dvdineg)
    plt.show()
    
def Ic_derivative(x,y,z,points_mult=2,window_smooth=3,smooth=True):
    
    '''
    Use the code of derivative_dvdi and then interpolates that data
    to get more points. The default is to get twice as many points. 
    This gives smoother data for fitting the model later. The interpolation
    also ensures that the data is equally spaced, which is needed
    to make any Fourier transform analysis. 
    Then I make a smoothing of the newly interpolated data with a running
    mean method to make it again easier to extrac a period
    with a Fourier transform method. 
    x has to be in the right direction in order to make the interpolation
    '''
    
    # This already makes a plot 
    dvdi1,dvdipos,dvdineg=derivative_dvdi(x,y,z)
    
    # Make the interpolated new data  
    if x[1,0]>x[0,0]:
        int_x=np.linspace(x[0,0],x[-1,0],points_mult*len(x))
        int_dvdipos=np.interp(int_x,x[:,0],dvdipos)
        int_dvdineg=np.interp(int_x,x[:,0],dvdineg)
        print('good order')
    else: 
        print('inverted order')
        int_x=np.linspace(x[-1,0],x[0,0],points_mult*len(x))
        int_dvdipos=np.interp(int_x,np.flip(x[:,0]),dvdipos)
        int_dvdineg=np.interp(int_x,np.flip(x[:,0]),dvdineg)
    
    # Plot the interpolated data to make sure that it fits 
    
    plt.figure(figsize=(10,5))
    plt.plot(x[:,0],dvdipos,'.-')
    plt.plot(int_x,int_dvdipos,'.-')
    plt.show()
    
    if smooth:
        # Now smooth the interpolated data 
        
        sm_dvdipos=runningMeanFast(int_dvdipos,window_smooth)
        sm_dvdineg=runningMeanFast(int_dvdineg,window_smooth)
        sm_x=int_x
    
        ## Correct the first and last values. If not corrected it will 
        ### average with the other side of the data and it will be incorrect
        
        sm_dvdipos[0]=dvdipos[0]
        sm_dvdipos[-1]=dvdipos[-1]
        sm_dvdipos[-2]=dvdipos[-2]
        sm_dvdineg[0]=dvdineg[0]
        sm_dvdineg[-1]=dvdineg[-1]
        sm_dvdineg[-2]=dvdineg[-2]
        
        # Plot it again to make sure it is not doing something wrong 
        
        plt.figure(figsize=(10,5))
        #plt.plot(xF1[num][:,0],dvdipos,'.-')
        #plt.plot(xF1[num][:,0],np.abs(dvdineg),'.-')
        plt.plot(sm_x,sm_dvdipos,'b.-',)
        plt.plot(int_x,int_dvdipos,'r.-')
        plt.show()
    
        return sm_x,sm_dvdipos,sm_dvdineg 
    else:
        return int_x,int_dvdipos,int_dvdineg 
    
def magnJJ_bulk_edge(B,w1,w2,w3,f1,f2,f3,
             phase1,phase2,phase3,phase4,phase5,edgeratio=1):
    '''
    Model to simulate the data combining a JJ and a SQUID
    Since the sample is 2D the SQUID will have to take into account
    the small JJs as well. 
    The model is basically: 
        Ic= w1 * JJ_bulk + w2* JJ_edge + w3 * SQUID
    To simulate the Magnetic behavior I add phase shifts at 
    certain positions, which match the observed jumps in the Fraunhofer
    maps. 
    w1, w2, w3: weigths corresponding to the different contributions.
                They should add up to 1. 
    f1, f2, f3: Frequency of the different effects. This frequency 
                is later matched to a size. This should make sense
    Phi1, Phi2: Global phase of the features 
    Phase1...Phase5: Phases added in the magnetic transitions. 
    '''
    ret_data=np.zeros_like((B))
    count=0
    ### Make the two edges different 
    wr=w2*edgeratio/2
    wl=w2*(1-edgeratio/2)
    print(wr)
    print(wl)
    for i in B:
        
        if i >=7.5:
            shift=phase1 #0.3
            den1=np.pi*2*f1*i
            fbulk=np.abs(np.sin(i*np.pi*2*f1+shift)/(den1+shift))
            
            sq1=np.abs(np.sin(i*np.pi*2*f3+shift))
            
            den2=np.pi*2*f2*i
            fedge=np.abs(np.sin(i*np.pi*2*f2+shift)/(den2+shift))
            
            ret_data[count]=(w1*fbulk+w3*sq1+wr*fedge+wl*fedge)
            
        elif i >=2.5 and i<7.5:
            shift=phase2 #np.pi/6
            #shift=shift
            den1=np.pi*2*f1*i
            fbulk=np.abs(np.sin(i*np.pi*2*f1+shift)/(den1+shift))
            
            sq1=np.abs(np.sin(i*np.pi*2*f3+shift))
            
            den2=np.pi*2*f2*i
            fedge=np.abs(np.sin(i*np.pi*2*f2+shift)/(den2+shift))
            
            ret_data[count]=(w1*fbulk+w3*sq1+wr*fedge+wl*fedge)

        elif i >=-2 and i<2.5:
            
            shift=phase3 # np.pi/2
            den1=np.pi*2*f1*i
            fbulk=np.abs(np.sin(i*np.pi*2*f1+shift)/(den1+shift))
            
            sq1=np.abs(np.sin(i*np.pi*2*f3+shift))
            
            den2=np.pi*2*f2*i
            fedge=np.abs(np.sin(i*np.pi*2*f2+shift)/(den2+shift))
            
            ret_data[count]=(w1*fbulk+w3*sq1+wr*fedge+wl*fedge)
            #if i>1.9:
            #    print('Total shift of '+str(0.4+np.pi/2) +' or '+str((0.4+np.pi/2)/np.pi)+' times pi')
        elif i >= -7.5 and i<-2: 
            ### Here it recovers the initial phase shift 
            shift=phase4 # -np.pi/2
            den1=np.pi*2*f1*i
            fbulk=np.abs(np.sin(i*np.pi*2*f1+shift)/(den1+shift))
                        
            sq1=np.abs(np.sin(i*np.pi*2*f3+shift))
            
            den2=np.pi*2*f2*i
            fedge=np.abs(np.sin(i*np.pi*2*f2+shift)/(den2+shift))
            
            ret_data[count]=(w1*fbulk+w3*sq1+wr*fedge+wl*fedge)
            
            #if i>-3:
            #    print('Total shift of '+str(0.4+np.pi/2-np.pi/2) +' or '+str((0.4+np.pi/2-np.pi/2)/np.pi)+' times pi')
                
        elif i<-7.5:
            shift=phase5 #2*np.pi
            den1=np.pi*2*f1*i
            fbulk=np.abs(np.sin(i*np.pi*2*f1+shift)/(den1+shift))

            sq1=np.abs(np.sin(i*np.pi*2*f3+shift))
                        
            den2=np.pi*2*f2*i
            fedge=np.abs(np.sin(i*np.pi*2*f2+shift)/(den2+shift))
            
            ret_data[count]=(w1*fbulk+w3*sq1+wr*fedge+wl*fedge)
        
        count=count+1
    total_shift=phase1+phase2+phase3+phase4+phase5
    print('Total phase shift is:' +str(total_shift/np.pi)+r'$\times \pi$')
    return ret_data

def magnJJ_bulkedge(B,w1,w2,w3,f1,f2,f3,
             phi1,phi2,phase1,phase2,phase3,phase4,phase5):
    '''
    SQUID formula is with a cosine
    All the formulas are multiplied by pi not 2 pi
    '''
    ret_data=np.zeros_like((B))
    count=0
    #print('Initial shift of '+str(0.4+np.pi/2-np.pi/2) +' or '+str((0.4+np.pi/2-np.pi/2)/np.pi)+' times pi')
    for i in B:
        if i >=7.5:
            shift=phase1 #0.3 
        elif i >=2.5 and i<7.5:
            shift=phase2 #np.pi/6
        elif i >=-2 and i<2.5:
            shift=phase3 # np.pi/2
        elif i >= -7.5 and i<-2: 
            ### Here it recovers the initial phase shift 
            shift=phase4 # -np.pi/2
        elif i<-7.5:
            shift=phase5 #2*np.pi

        flux1=np.pi*f1*i+shift
        fr1=np.abs(np.sin(flux1)/flux1)
        
        fluxsq=np.pi*f3*i+shift
        sq1=np.abs(np.cos(fluxsq))
        
        flux2=np.pi*f2*i+shift
        fr2=np.abs(np.sin(flux2)/flux2)
        
        ret_data[count]=(w1*fr1+w3*sq1+w2*fr2)
        count=count+1
    ## Fix for possible NaN values, usually at the zero field position 
    ## Once found just interpolate between previous and next value 
    nan_positions=np.where(np.isnan(ret_data))
    for n in nan_positions:
        print(n)
        ret_data[n]=(ret_data[n-1]+ret_data[n+1])/2+(ret_data[n-1]-ret_data[n-1])
    return ret_data
    
def magnJJ_bulkedge_up(B,w1,w2,w3,f1,f2,f3,
             phi1,phi2,phase1,phase2,phase3,phase4,phase5,hyst):

    ret_data=np.zeros_like((B))
    count=0
    #print('Initial shift of '+str(0.4+np.pi/2-np.pi/2) +' or '+str((0.4+np.pi/2-np.pi/2)/np.pi)+' times pi')
    for i in B:
        
        if i >=7.5:
            shift=phase1 #0.3
            den1=np.pi*2*f1*i
            fr1=np.abs(np.sin(i*np.pi*2*f1+phi1+shift)/(den1+phi1+shift))
            
            sq1=np.abs(np.sin(i*np.pi*2*f3+phi1+shift))
            
            den2=np.pi*2*f2*i
            fr2=np.abs(np.sin(i*np.pi*2*f2+phi2+shift)/(den2+phi2+shift))
            
            ret_data[count]=(w1*fr1+w3*sq1+w2*fr2)
            
        elif i >=2 and i<7.5:
            shift=phase2 #np.pi/6
            #shift=shift
            den1=np.pi*2*f1*i
            fr1=np.abs(np.sin(i*np.pi*2*f1+phi1+shift)/(den1+phi1+shift))
            
            sq1=np.abs(np.sin(i*np.pi*2*f3+phi1+shift))
            
            den2=np.pi*2*f2*i
            fr2=np.abs(np.sin(i*np.pi*2*f2+phi2+shift)/(den2+phi2+shift))
            
            ret_data[count]=(w1*fr1+w3*sq1+w2*fr2)
            
        elif i >=-3 and i<2:
            
            shift=phase3 # np.pi/2
            den1=np.pi*2*f1*i
            fr1=np.abs(np.sin(i*np.pi*2*f1+phi1+shift)/(den1+phi1+shift))
            
            sq1=np.abs(np.sin(i*np.pi*2*f3+phi1+shift))
            
            den2=np.pi*2*f2*i
            fr2=np.abs(np.sin(i*np.pi*2*f2+phi2+shift)/(den2+phi2+shift))
            
            ret_data[count]=(w1*fr1+w3*sq1+w2*fr2)
            #if i>1.9:
            #    print('Total shift of '+str(0.4+np.pi/2) +' or '+str((0.4+np.pi/2)/np.pi)+' times pi')
        elif i >= -6 and i<-3: 
            ### Here it recovers the initial phase shift 
            shift=phase4 # -np.pi/2
            den1=np.pi*2*f1*i
            fr1=np.abs(np.sin(i*np.pi*2*f1+phi1+shift)/(den1+phi1+shift))
                        
            sq1=np.abs(np.sin(i*np.pi*2*f3+phi1+shift))
            
            den2=np.pi*2*f2*i
            fr2=np.abs(np.sin(i*np.pi*2*f2+phi2+shift)/(den2+phi2+shift))
            
            ret_data[count]=(w1*fr1+w3*sq1+w2*fr2)
            
            #if i>-3:
            #    print('Total shift of '+str(0.4+np.pi/2-np.pi/2) +' or '+str((0.4+np.pi/2-np.pi/2)/np.pi)+' times pi')
                
        elif i<-6:
            shift=phase5 #2*np.pi
            den1=np.pi*2*f1*i
            fr1=np.abs(np.sin(i*np.pi*2*f1+phi1+shift)/(den1+phi1+shift))

            sq1=np.abs(np.sin(i*np.pi*2*f3+phi1+shift))
                        
            den2=np.pi*2*f2*i
            fr2=np.abs(np.sin(i*np.pi*2*f2+phi2+shift)/(den2+phi2+shift))
            
            ret_data[count]=(w1*fr1+w3*sq1+w2*fr2)

        count=count+1
    return ret_data
    
def magnJJ_bulkedge_vector(B,w1,w2,w3,f1,f2,f3,
              phase,jumps):
     '''
     SQUID formula is with a cosine
     All the formulas are multiplied by pi not 2 pi
     Phase and jumps have are a vector of numpy array shape, they have
     the same length
     '''
     ret_data=np.zeros_like((B))
     num_jumps=len(jumps)
     shift=0
     count=0
     #print('Initial shift of '+str(0.4+np.pi/2-np.pi/2) +' or '+str((0.4+np.pi/2-np.pi/2)/np.pi)+' times pi')
     for i in B:
         if i >=jumps[0]:
             shift=phase[0] #0.3 
         elif i <=jumps[1]:
             shift=phase[1] #0.3     
         else:
             shift=0
         flux1=np.pi*f1*i+shift
         fr1=np.abs(np.sin(flux1)/flux1)
         
         fluxsq=np.pi*f3*i+shift
         sq1=np.abs(np.cos(fluxsq))
         
         flux2=np.pi*f2*i+shift
         fr2=np.abs(np.sin(flux2)/flux2)
         
         ret_data[count]=(w1*fr1+w3*sq1+w2*fr2)
         count=count+1
     ## Fix for possible NaN values, usually at the zero field position 
     ## Once found just interpolate between previous and next value 
     nan_positions=np.where(np.isnan(ret_data))
     for n in nan_positions:
         print(n)
         ret_data[n]=(ret_data[n-1]+ret_data[n+1])/2+(ret_data[n-1]-ret_data[n-1])
     return ret_data      