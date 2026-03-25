24# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 10:56:20 2022

Code to calculate the distrubution of the current density in a JJ according to the
procedure they follow in Hart, S., Ren, H., et al. (Nature Physics 2014). "Induced superconductivity in the quantum spin Hall edge.
Which is based on the approach by Barone. 
The idea is to create an even and odd Current distributions and calculate the Fourier
transform of the sum of this. This complex Fourier transform should give you the 
current distribution in the sample. 
This code would be the opposite of IcB_SimulationFrom_Jcx.py, 
where you start from a current distribution and you calculate the Fraunhofer pattern. 
Here you start from the Fraunhofer pattern (real data) and you calculate the current distribution.
@author: Jaime Díez Mérida @ LDQM group at ICFO and LMU
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate 
from scipy.fft import fft, fftfreq, fftshift, ifft
from homemade_functions import plot_pre,plot_map,plot_line
from scipy.signal import savgol_filter

def subs_Fr(B,Fr,sm_window):
    ''' Make a background substraction of the Fraunhofer to make sure it gets
        the oscilations correctly '''
    if int(len(Fr)/sm_window)%2==0:
        window=int(len(Fr)/sm_window)+1
    else:
        window=int(len(Fr)/sm_window)
        
    backgrFr=savgol_filter(Fr,window,3)
    subsFr=Fr-backgrFr
    
    ax=plot_pre()
    plt.plot(B,backgrFr,'b')
    plt.plot(B,subsFr,'r')
    plot_line(ax,legx=0.4,legy=1.025,font_size=16,font_leg=2)
    plt.show()
    
    return subsFr

def Iceven_func(B,Fr,shift=0,sm_window=3,subs=True,center_w=5):
    '''
    Create the even part of the Fraunhofer pattern. It is basically the same shape
    as the original Fraunhofer but letting it change sign everytime there is an 
    oscillation. The code also gives back the flipping signal to create it. 
    Shift marks if the Fraunhofer is not centered at zero 
    '''
    flip=np.zeros_like(B)
    
    ## make a background substraction to make it easier to find the minimums 
    if subs:
        subsFr=subs_Fr(B,Fr,sm_window)
    else:
        subsFr=Fr
    
    # now create the flip function
    ## To be correct it has to start from the center always 
    
    
    count_flip=0
    for i in range(0,len(B)):
        # As it goes down is all +1,
        if subsFr[i]<subsFr[i-1] and count_flip%2==0:
            flip[i]=1
        # After it has reached the minimum it should change sign 
        elif subsFr[i]<subsFr[i-1] and count_flip%2!=0:
            flip[i]=-1
        elif subsFr[i]>subsFr[i-1] and count_flip%2==0: # New oscillation 
            flip[i]=-1
        elif subsFr[i]>subsFr[i-1] and count_flip%2!=0:
            flip[i]=+1
        try: ## Find the minimums
            if subsFr[i]>subsFr[i-1] and subsFr[i]>subsFr[i+1]:
                count_flip=count_flip+1
        except IndexError:
            continue

    ## Ensure there are no strange things in the central region 
    if center_w!=0:
        for i in range(int(len(B)/2)+shift-center_w,int(len(B)/2)+shift+center_w):
            flip[i]=+1
    else:
        flip[int(len(B)/2)]=flip[int(len(B)/2)+1]
    
    if flip[int(len(B)/2)]==-1:
        revert=-1
    else:
        revert=1
    
    ## Create the even function, changing sign everytime an oscillation finishes
    Ic_even=Fr*flip*revert
    
    ax=plot_pre()
    #ax.plot(B,Fr,'k',label='$I_{c,data}$')
    ax.plot(B,Ic_even,'b')
    ax.plot(B,flip,'r--')
    plot_line(ax,leg=False,ylabel='$I_E$ (nA)',font_size=16)
    plt.show()
    
    ax=plot_pre()
    #ax.plot(B,Fr,'k',label='$I_{c,data}$')
    ax.plot(B,Ic_even,'b')
    plot_line(ax,leg=False,ylabel='$I_E$ (nA)',font_size=16)
    plt.show()
    
    return Ic_even,flip

def Icodd_func(B,Fr):
    '''
    Create the odd part of the Fraunhofer pattern. This is the imaginary part
    of the complex current density.
    It finds the position where an oscillation finishes and it picks up its value
    with different sign each time. Then it interpolates in between these two
    values creating an odd funciton w.r.t. the field. This is the imaginary part
    of the complex current dentity as a function of field. 
    '''
       
    Ic_odd=np.zeros_like(B)
    # Find the minimums 
    flip_odd=1
    positions=[]
    for i in range(0,len(B)):
        try:
            if Fr[i]<Fr[i-1] and Fr[i]<Fr[i+1]: #it is a minimum
                Ic_odd[i]=Fr[i]*flip_odd
                flip_odd=flip_odd*(-1) ## Switch sign everytime there is a minimum
                positions.append(i)
        except IndexError:
            continue
    
    #print(positions)
    ## Interpolate in between these points  
    for i in range(len(positions)-1):
        Ic_odd[positions[i]:positions[i+1]]=np.linspace(Ic_odd[positions[i]],
                                                    Ic_odd[positions[i+1]],
                                               positions[i+1]-positions[i])
    #print(subsFr[295:305])

    ax=plot_pre()
    ax.plot(B,Ic_odd,'b')
    #ax.plot(B,Fr,'r',label='$I_{c,data}$')
    plot_line(ax,leg=False,ylabel='$I_O$ (nA)',font_size=16)
    plt.show()
   
    return Ic_odd

def Jc_func(B,Fr,shift=0,subs=True,center_w=5,sm_window=3):
    ''''
    The complex current distribution of the Fraunhofer. 
    '''
    Iceven,flip=Iceven_func(B,Fr,shift,sm_window,subs,center_w)
    Icodd=Icodd_func(B,Fr)
    
    ax=plot_pre()
    ax.plot(B,Icodd,'b',label='$I_{c,odd}$')
    ax.plot(B,Iceven,'r',label='$I_{c,even}$')
    plot_line(ax,legx=0.8,legy=1.025,font_size=16,font_leg=2)
    plt.show()
    
    
    return Iceven+(1j)*Icodd

def extract_Js(B,Fr,shift=0,subs=True,center_w=0,sm_window=3):
    '''
    Calculate the distribution of the current density in the sample and make 
    its Fourier transform to obtain the current distribution inside the sample
    '''
    ## Calculate the complex current density as a funcion of field 
    Jc=Jc_func(B,Fr,shift,subs,center_w,sm_window)
    ## Find if there are any infinites in the Fraunhfoer and convert them 
    
    ## Make the Fourier transform 
    ### number of signal points
    N = len(B)
    ### sample spacing
    dB=np.abs(B[1]-B[0])
    dx=2*np.pi/(dB*N)
    
    
    '''
    T = np.abs(B[1]-B[0])/N
    x = np.linspace(0.0, N*T, N, endpoint=False)
    '''
    
    y = Jc
    yf = fft(y)
    xf = fftfreq(N, dx)
    xf = fftshift(xf) ## fft shift reallocates to center around 0
    yplot = fftshift(yf)
    yplot=1/N * np.abs(yplot)
    
    ax=plot_pre()
    ax.plot(xf, yplot,'navy')
    plot_line(ax,legx=0.4,legy=1.025,font_size=16,font_leg=2,
              xlabel='$x$ ($\mu$m)',ylabel='$J_c$ [nA/$\mu$m]')
    plt.show()
    print(dx)
    return Jc,xf,yplot

