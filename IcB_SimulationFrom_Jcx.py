# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 10:48:18 2022


Code to simulate Fraunhofer patterns via the Dynes and Fulton 1971 PRB paper. This simulates a Fraunhofer pattern 
for a given current distribution, which can be homogenous or with edge states. 
The code also simulates the dVdI maps that you would obtain in a tunneling spectroscopy 
The code here was used for simulations published in https://doi.org/10.1038/s41467-023-38005-7. 

@author: Jaime Díez Mérida @ LDQM group at ICFO and LMU
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.integrate as integrate 
from scipy.fft import fft, fftfreq, fftshift, ifft
from homemade_functions import plot_pre,plot_map,plot_line
from scipy.signal import savgol_filter
from astropy.convolution import Gaussian1DKernel, convolve, Gaussian2DKernel
import cmath 
import LDQM_dataplot
import time 

"""
First you have to create the current distribution that you want
It can be homogenous or with edge states
"""

def exp_decay(x,tau):
    '''Function to plot a normalized exponential decay'''
    return np.exp(-x*tau)/np.amax(np.exp(-x*tau))

def square_exp(x,decay_l,decay_rate,jc=1,):
    '''
    Create a current distribution for a uniform distribution of curernt
    with an exponential decay on the edges of the sample
    x is the real space widht
    decay_l is the length of the decay
    decay_rate gives how fast you decay
    jc would be the critical current. it can be a.u. 
    '''
    ls=(x[-1]-x[0])/len(x) ## interval of points in the junction
    ly=int(decay_l/ls)-1  ## the number of points in the decaying region
    j=np.zeros_like(x)
    ''' The bulk will all have the same distribution '''
    for i in range(ly,len(x)-ly):
        j[i]=jc
    '''Now we do the decay'''
    j[0:ly]=exp_decay(-x[0:ly],decay_rate)
    j[len(x)-ly:]=exp_decay(x[len(x)-ly:],decay_rate)
    return j



def edges_exp(x,edge_w,ratio,edgeratio,decay_l,decay_rate,edgew2=0,jc=1,):
    '''
    Current distribution with edge states. The edges states have 
    an exponential decay both towards the end of the sample and towards
    the center of the sample. 
    '''
    ls=(x[-1]-x[0])/len(x) ## interval of points in the junction
    ly=int(decay_l/ls)-1  ## the number of points in the decaying region
    ledge=int(edge_w/ls)-1 ## number of points in the edge region
    edge_start=decay_l
    edge_end=decay_l+edge_w
    j=np.zeros_like(x)
    J_bulk=jc*ratio
    J_edge1=jc*(1-ratio)*edgeratio
    J_edge2=jc*(1-ratio)*(1-edgeratio)
    ''' The bulk will all have the same distribution '''
    for i in range(ly+ledge+ly,len(x)-ly-ledge-ly):
        j[i]=J_bulk
    '''Now make the edges'''
    for i in range(ly,ly+ledge):
        j[i]=J_edge1
    j[0:ly]=exp_decay(-x[0:ly],decay_rate)*J_edge1
    j[ly+ledge:ly+ledge+ly]=exp_decay(x[ly+ledge:ly+ledge+ly],decay_rate)*J_edge1+J_bulk
    for i in range(ly+ledge,ly+ledge+ly):
        if j[i]>J_edge1:
            j[i]=J_edge1
        
    '''Right edge'''
    for i in range(len(x)-ly-ledge,len(x)-ly):
        j[i]=J_edge2
    
    j[len(x)-ly-ledge-ly:len(x)-ledge-ly]=exp_decay(-x[len(x)-ly-ledge-ly:len(x)-ledge-ly],decay_rate)*J_edge2+J_bulk
    j[len(x)-ly:]=exp_decay(x[len(x)-ly:],decay_rate)*J_edge2
    
    for i in range(len(x)-ly-ledge-ly,len(x)-ledge-ly):
        if j[i]>J_edge2:
            j[i]=J_edge2
    return j

"""
Then you have to make the code for the DOS of the SC to simulate a dVdI


Tsme=np.zeros_like(V)
for i in range(len(V)):
    for j in range(len(integral2)):
        Tsme[i]=2*k*T
        if abs(V[i])>integral2[j]:
            DOSV[i,j]=V[i]
convoluted=[]
for i in range(len(DOSV)):
    convoluted.append(convolve(Ns(DOSV)[i],g,boundary='extend'))
"""  
e=1.602E-19
V=np.linspace(-0.5,0.5,401) # Voltage to give the energy range
E=e*V
B=np.linspace(-80,80,601) #Field range that you consider 
DOSV=np.zeros((len(V),len(B)))  # DOS, size of both because is a map
N0=1

def Ns(E,Delta,Is=1):
    return Is*(np.abs(E)/np.sqrt(np.abs((E)**2-Delta**2)))/N0

def Ns2D(E,Delta):
    """E is the voltage range times the electron charge 
    Delta would be your gap energy. In the simulation it would be the Ic
    """
    data=np.zeros_like((DOSV))
    #conv_data=np.zeros_like((data))
    for i in range(len(E)):
        for j in range(len(Delta)):
            data[i,j]=Is*(np.abs(E[i,j])/np.sqrt(np.abs((E[i,j])**2-(Delta[j])**2)))/N0
    return data

def Ns2D_IAsym(V,Icplus,Icminus, Is=35,N0=1):
    """
    Same as above but the positive and negative critical currents can be different
    It is already thought for dVdI simulation so instead of E you have Voltage
    And instead of gap we have the positive and negative critical current. 
    Is would be the critical current of the data, otherwise it is normalized
    """
    E=np.zeros((len(V),len(B)))
    data=np.zeros_like((E))
    Tsme=np.zeros_like(V)
    for i in range(len(V)):
        for j in range(len(Icplus)):
            Tsme[i]=2*k*T
            if V[i]>=0:
                if abs(V[i])>Icplus[j]:
                    E[i,j]=V[i]
            elif V[i]<0:
                if abs(V[i])>Icminus[j]:
                    E[i,j]=V[i]

    for i in range(len(V)):
        if V[i]>=0:
            for j in range(len(Icplus)):
                data[i,j]=Is*(np.abs(E[i,j])/np.sqrt(np.abs((E[i,j])**2-Icplus[j]**2)))/N0
        else:
            for j in range(len(Icplus)):
                data[i,j]=Is*(np.abs(E[i,j])/np.sqrt(np.abs((E[i,j])**2-Icminus[j]**2)))/N0
  
    return data

phi0=2.067E-15
def int_Js_brokenT(Js,B,x,edge_l,l,varphi1,varphi2,varphi3,t=0.25,Ic=1,edge_phase=0.5):
    '''
    Js would be the real space current density (obtained by edges_exp)
    This code will give a general equation to integrate to obtain the actual current density 
    edge_phase determines whether you take into account the decay with the phase or not 
    '''
    left_edge_st=x[0]+edge_phase*l
    left_edge_end=x[0]+edge_l+(1-edge_phase)*l
    right_edge_st=x[-1]-edge_l-(1-edge_phase)*l
    right_edge_end=x[-1]-edge_phase*l
    data=np.zeros_like((x),dtype='complex_')
    for i in range(len(x)):
        if x[i]<left_edge_end and x[i]>left_edge_st:
            beta=2*np.pi*t*(B+varphi1)/phi0   
        elif x[i]<right_edge_end and x[i]>right_edge_st:
            beta=2*np.pi*t*(B+varphi2)/phi0
        else:
            beta=2*np.pi*t*(B+varphi3)/phi0
            
        data[i]=Js[i]*np.exp(1j*beta*x[i])*Ic
        
    return data
savename='C:/Users/jdiez/LDQM/a. Data/Ju1/JJ_BF_2021/'
def Ic_calculate(x,B,edge_l,l,ratio,edgeratio,Ic,decay_rate,t,
                 phase_edge1,phase_edge2,phase_bulk, I_asym,
                 phase4=0,phase5=0,phase6=0,orig_cmap=matplotlib.cm.RdBu_r,
                 font=15,save=False,text_true=True,savename=savename):
    """
    x=np.linspace(-0.5,0.5,601), is an array of the width size of the Hall bar 
    B=np.linspace(-80,80,1001), is an array of the field you apply
    edge_l=0.05, is the width of the edges 
    l=0.05, is the decay length of the edges exponential 
    ratio=0.1, is the ratio of current through the bulk and the edges
    edgeratio=0.35, is the current ratio between the two edges
    Ic=19, is the critical current for the model
    decay_rate=120, is the decay rate of the exponential of the current distribution
    t=1/1.8, is the length of the junction. In the 2D formalism it is the width of the Hall bar 
    divided by 1.8
    phase_edge1=np.pi/4, is the phase of edge left, makes Fraunhfoer asymmetric
    phase_edge2=-np.pi/4, is the phase of edge right, makes Fraunhofer asymmetric
    phase_bulk=0, is the phase of the bulk. This just shift the Fraunhfoer left or right
    I_asym=True, makes the sign of the phases change sign for the Ic pos and negative 
    phase4,phase5 and phase6, could be used to make a complete different phase between the positive
    and negative critical currents
    The rest are all for plottings 
    If you make save False, then you can make the save outside of this function 
    """
    Icplus=[]
    Icminus=[]
    Jsx=edges_exp(x,edge_l+l,ratio,edgeratio,l,decay_rate,0,Ic)
    for i in B:
        ## Obtain the complex current density 
        int_trap_plus=np.trapz(int_Js_brokenT(Jsx,i,x,edge_l,l,phase_edge1,phase_edge2,phase_bulk,t,Ic),x)
        
        ## Calculate the modulus of the complex integral
        modulus=cmath.sqrt((int_trap_plus.real)**2+(int_trap_plus.imag)**2)

        Icplus.append(modulus )  
        
        if I_asym: # Do we add the curren phase asymmetry or not 
            int_trap_minus=np.trapz(int_Js_brokenT(Jsx,i,x,edge_l,l,-phase_edge1,
                                                   -phase_edge2,-phase_bulk,t,Ic),x)
            mod_minus=cmath.sqrt((int_trap_minus.real)**2+(int_trap_minus.imag)**2) 
            Icminus.append( mod_minus )
        else:
            int_trap_minus=np.trapz(int_Js_brokenT(Jsx,i,x,edge_l,l,
                                                   phase4,phase5,phase6,t,Ic),x)
            mod_minus=cmath.sqrt((int_trap_minus.real)**2+(int_trap_minus.imag)**2) 
            Icminus.append( mod_minus )
    """
    The final plotting will look as follows. I can 
    """
    ## Plot the results 
    ax=plot_pre(6,4)
    ax.plot(B,Icplus,'k')
    ax.plot(B,np.abs(Icminus),'b')
    ax.axvline(x=0,linestyle='dashed',color='silver',linewidth=1)
    plot_line(ax,leg=False)
    ax.set_ylim(-1,80)
    plt.show()

    ax=plot_pre(6,4)
    #ax.plot(x,Jsx,'navy')
    ax.plot(x,Jsx,'navy')
    plot_line(ax,xlabel='$x$ ($\mu$m)',ylabel='$J_c$ (nA/$\mu$m)',leg=False)
    ax.set_ylim(0,12)
    ax.set_xlim(-0.55,0.55)
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.25))
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.5))
    #if save:
        #plt.savefig(savename+'EdgeBulkAsymJs_TS_Fraunhofer_Real'+time.strftime("%H%M%d%m%y")+'.png',dpi=300) 
        #plt.savefig(savename+'EdgeBulkAsymJs_TS_Fraunhofer_Real'+time.strftime("%H%M%d%m%y")+'.pdf') 
        #savecount=savecount+1 

    plt.show()

    '''Code for the maps'''

    mycolor=LDQM_dataplot.shiftedColorMap(0.5,orig_cmap)


    fig=plt.figure(figsize=(6,4))
    ax=fig.add_subplot(111)
    im=plt.pcolormesh(B,V,Ns2D_IAsym(V,Icplus,Icminus)/15,cmap='RdBu_r',vmax=5,rasterized=True)

    ax.set_ylabel('$I$ (nA)',fontsize=font)
    ax.set_xlabel('$B$ (mT)',fontsize=font)  
    ### Add the colorbar 
    cbaxes = fig.add_axes([0.75, 0.90, 0.15,0.015]) #left,bottom,width,heigth
    cb = plt.colorbar(im,shrink=0.3,cax=cbaxes,orientation='horizontal') 

    cbaxes.xaxis.tick_top()
    cb.ax.tick_params(labelsize=font,direction='in')

    ax.tick_params(axis='both',which='major',labelsize=font,direction='in',pad=6,size=8,width=1.5)
    ax.tick_params(axis='both',which='minor',direction='in',size=4,width=1.5)
    #ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(20))
    #ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(50))
    if text_true:
        ax.text(20,40,r'$\psi_{left}(I>0)$'+'$={0:.2f}\pi$'.format(phase_edge1/np.pi),fontsize=font-3)
        ax.text(20,60,r'$\psi_{right}(I>0)$'+'$={0:.2f}\pi$'.format(phase_edge2/np.pi),fontsize=font-3)
        ax.text(20,80,r'$\psi_{bulk}(I>0)$'+'$={0:.2f}\pi$'.format(phase_bulk/np.pi),fontsize=font-3)
        if I_asym==True:
            ax.text(20,-40,r'$\psi_{left}(I>0)$'+'$={0:.2f}\pi$'.format(-phase_edge1/np.pi),fontsize=font-3)
            ax.text(20,-60,r'$\psi_{right}(I<0)$'+'$={0:.2f}\pi$'.format(-phase_edge2/np.pi),fontsize=font-3)
            ax.text(20,-80,r'$\psi_{bulk}(I>0)$'+'$={0:.2f}\pi$'.format(-phase_bulk/np.pi),fontsize=font-3)
        else:  
            ax.text(20,-40,r'$\psi_{left}(I>0)$'+'$={0:.2f}\pi$'.format(phase_edge1/np.pi),fontsize=font-3)
            ax.text(20,-60,r'$\psi_{right}(I<0)$'+'$={0:.2f}\pi$'.format(phase_edge2/np.pi),fontsize=font-3)
            ax.text(20,-80,r'$\psi_{bulk}(I>0)$'+'$={0:.2f}\pi$'.format(phase_bulk/np.pi),fontsize=font-3)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    #ax.axvline(x=0,linestyle='dashed',color='k',linewidth=1)
    if save:
        plt.savefig(savename+'FigS14_pi4_npi4_IB'+time.strftime("%H%M%d%m%y")+'.png',dpi=300) 
        plt.savefig(savename+'FigS14_pi4_npi4_IB'+time.strftime("%H%M%d%m%y")+'.pdf') 
        savecount=savecount+1    
        plt.show()

    return Icplus, Icminus, Jsx

