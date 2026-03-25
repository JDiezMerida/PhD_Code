# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 17:04:51 2021

@author: Jaime Díez Mérida @ LDQM group at ICFO and LMU
The code was based on the article "Solving the Generalized Poisson Equation Using the 
Finite-Difference Method (FDM)" by  James R. Nagel (February 15, 2012). The results used from this code
were used for simulations of TBG Josephson junctions published in https://doi.org/10.1038/s41467-023-38005-7

The code includes several versions of the same simulation, with different boundary conditions and solution approaches. 
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
#import ipympl
import time
import _pickle as pickle
import matplotlib 
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import basic_plotting
class Object():
    pass

# top hBN half the size as bottom hBN. So permittivity is double in the top than in the bottom 

# Now make a cut in the top one

# parallel plate capacitor

# top hBN half the size as bottom hBN. So permittivity is double in the top than in the bottom 

# Now make a cut in the top one

# parallel plate capacitor

# Set dimensions of the problem
L=1
Nx=81
Ny=61
ds=1

# Make the charge density matrix
rho0=-1
rho1=1
# One of the plates is set to 1 and the other one to -1. Lets put the plate in layer number 5
bot_plate=int(Ny/3)
top_plate=int(Ny-bot_plate-0.5*bot_plate)
# Also add the dielectric block inside the two capacitors
epsr1=5
epsr2=5
#x coordinate of the plates
plate_start=int(Nx/4)
cut_x=int(Nx/16)

orig_cmap = matplotlib.cm.jet
#save_folder='W:/a. Members/Jaime/f. Simulations/'

'''
Simulation with Dirichlet boundary Conditions, which makes it faster.
This is the final and best code
'''

def capacitors_Dirich(Nx,ds,dsy,rho_top_left,rho_bot,
               bot_hBN,top_hBN,bot_epsr,top_epsr,plate_start_x,
               cut_x,conv_error=10E-6,
               plot_names='', save=False,savename='',top_size=1,
               bot_size=1,hd=0,mid_split=False,
               rho_top_right=0,
               max_iter=1E4,cut_center=0,save_obj=False,onetop=False, notops=False,
               triplegate=False,triple_plate=0,rho_triple=0,
               ):
    
    print('Save settings is set to '+str(save))
    if cut_center==0:
        cut_center=int(Nx/2)
    
    if mid_split==False: 
        rho_top_right=rho_top_left
    
    ## The y-size for the proper simulations should be the exact size of the
    ## capacitor. Otherwise there are some artifact from the way it "grounds"
    ## the environment. 
    Ny=bot_hBN+top_hBN+top_size+bot_size
    if Ny%2==0:
        Ny=Ny+1
        bot_plate= 0+bot_size
        top_plate=Ny-top_size-1
    else:
        bot_plate= 0+bot_size
        top_plate=Ny-top_size-1
    
    x=np.arange(0,Nx,ds)
    y=np.arange(0,Ny,dsy)
    x_E=np.arange(0,Nx-ds,ds)
    y_E=np.arange(0,Ny-dsy,dsy)
    
    x_real=np.arange(-int(Nx/2),int(Nx/2)+ds,ds) #arange leaves out 1 values
    y_real=np.arange(-int(Ny/2),int(Ny/2)+dsy,dsy)
    x_Ereal=np.arange(-int(Nx/2),int(Nx/2),ds)
    y_Ereal=np.arange(-int(Ny/2),int(Ny/2),dsy)
    
    X,Y=np.meshgrid(x,y)
    X_E,Y_E=np.meshgrid(x_E,y_E)
    X_R,Y_R=np.meshgrid(x_real,y_real)
    X_ER,Y_ER=np.meshgrid(x_Ereal,y_Ereal)
    
    # Make the matrices to be filled. In np.zeros, the order is inverted
    rho=np.zeros((len(y),len(x)))
    epsilon=np.zeros((len(y),len(x)))+1
    # Make the initial guess for out matrix --> All 0s

    V=np.zeros((len(y),len(x)))
    R=np.zeros((len(y),len(x)))
    
    top_plate_size=np.arange(top_plate-top_size,top_plate+top_size+1)
    bot_plate_size=np.arange(bot_plate-bot_size,bot_plate+bot_size+1)
    plate_len=np.arange(plate_start_x,len(x)-plate_start_x)
    triple_plate_size=np.arange(triple_plate-bot_size,triple_plate+bot_size+1)

    ### Give the initial conditions. The charge distribution and the voltage of the electrodes 

    for i in plate_len:
        if cut_x==0.0:
            for j in top_plate_size:
                rho[j,i] = rho_top_left
                V[j,i]=rho_top_left    # Give the electrode the charge and not change it 
                
        else:
            if i <(cut_center-cut_x):
                for j in top_plate_size:
                    if notops:
                        rho[j,i] = 0
                        V[j,i]=0
                    else:
                        rho[j,i] = rho_top_left
                        V[j,i]=rho_top_left
            if i >(cut_center+cut_x):
                for j in top_plate_size:
                    if onetop==True or notops==True:
                        rho[j,i]=0
                        V[j,i]=0
                    else:
                        rho[j,i]=rho_top_right
                        V[j,i]=rho_top_right
                
        ### Same story as for the top plate, but now for the bottom plate
        
        for j in bot_plate_size:
            rho[j,i] = rho_bot
            V[j,i]=rho_bot
                
        if triplegate:
            for j in triple_plate_size:
                rho[j,i] = rho_triple
                V[j,i]=rho_triple


    ### Give the dielectric constant of the hBN inside the gates 
    for i in range(plate_start_x,len(x)-plate_start_x):
        for j in range(bot_plate,int(len(y)/2),1):
            epsilon[j,i] = bot_epsr  
        for j in range(int(len(y)/2),top_plate,1):
            epsilon[j,i] = top_epsr

    # Over relaxation constant. 0<w<2, if w is over 2 it will not converge

    t=np.cos(np.pi/Nx)+np.cos(np.pi/Ny)
    w=(8-np.sqrt(64-16*t**2))/t**2

    print('Over relaxation constant is: '+str(w))

    ### Solver ### 

    iterations = 0
    #eps = 10E-6 # Convergence threshold

    error = 1e4 # Large dummy error
    V_temp = np.copy(V)
    ### BCs are Dirichlet type, 0 for the boundaries and the charge for the electrodes
    
    while iterations < max_iter and error > conv_error:

        #V_temp = np.copy(V)
        error = 0 # we make this accumulate in the loop
        
        ### Left side of the electrodes and below the bottom electrode
        ### Keep the BCs always at zero voltage and the contacts at the given voltage 
  
        for i in range (1,len(x)-1):
            for j in range (1,len(y)-1):   
                # Calculate the permittivity matrices
                a0=epsilon[j,i]+epsilon[j-1,i]+epsilon[j,i-1]+epsilon[j-1,i-1]
                a1=0.5*(epsilon[j,i]+epsilon[j,i-1])
                a2=0.5*(epsilon[j-1,i]+epsilon[j,i])
                a3=0.5*(epsilon[j-1,i-1]+epsilon[j-1,i])
                a4=0.5*(epsilon[j,i-1]+epsilon[j-1,i-1])
                
                if j in bot_plate_size and i>plate_start_x and i < len(x)-plate_start_x:
                    R[j,i]=0
                    V[j,i]=rho_bot
                elif j in top_plate_size and i>plate_start_x and i < cut_center-cut_x and notops==False:
                    R[j,i]=0
                    V[j,i]=rho_top_left

                elif j in top_plate_size and i> cut_center+cut_x  and i < len(x)-plate_start_x and onetop==False and notops==False:
                    R[j,i]=0
                    V[j,i]=rho_top_right
                elif j in triple_plate_size and i>plate_start_x and i < len(x)-plate_start_x and triplegate==True:
                    R[j,i]=0
                    V[j,i]=rho_triple
                else: 
                    #Calculate the residual
                    R[j,i]=1/a0*(a1*V_temp[j+1,i] + a2*V_temp[j,i+1] +
                        a3*V_temp[j-1,i] + a4*V_temp[j,i-1] +rho[j,i])-V_temp[j,i]
                    
                    # Compute V using the overrelaxation method

                    V[j,i]=(V_temp[j,i]+w*R[j,i])
                    error += abs(V[j,i]-V_temp[j,i])

                    V_temp[j,i]=V[j,i] # Make it already work for the next iteration 

        iterations += 1
        error /= float(len(y))
    print (" iterations ="+str(iterations))
        
        
    ### Calculate the electric field
    ### it will always have one less component than the potentials, 
    ### just by the way it is defined
    Ex=np.zeros((len(y)-1,len(x)-1))
    Ey=np.zeros_like(Ex)
    Ey=np.zeros_like(Ex)
    Ex_c=np.zeros_like(Ex)
    Ey_c=np.zeros_like(Ex)
    E_total=np.zeros_like(Ex)

    ### This shifts the x and y with respect to each other
    
    for j in range (0,len(x) -1):
            for i in range (0,len(y) -1):

                Ex[i,j]=-(V[i,j+1]-V[i,j])/ds
                Ey[i,j]=-(V[i+1,j]-V[i,j])/dsy

                # Now we need to fix these shifts
    for j in range (0,len(x) -2):
            for i in range (0,len(y)-2):

                Ex_c[i,j]=0.5*(Ex[i+1,j]+Ex[i,j])
                Ey_c[i,j]=0.5*(Ey[i,j+1]+Ey[i,j])

                E_total[i,j]=np.sqrt(Ex_c[i,j]**2+Ey_c[i,j]**2)


    plt.figure(figsize=(8,5))
    plt.pcolormesh(X_R,Y_R,rho,cmap='jet')
    plt.colorbar()
    plt.title('Charges '+plot_names,fontsize=12)
    #savename='Charges'+plot_names
    
    plt.show()

    title_name='Voltage potentials '+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    basic_plotting.plot_map(X_R, Y_R, V, 8, 5, title_name, 
             x_label,y_label, z_label='', show=True, minmax=False,
             zmin=0, zmax=10, savename=save_name, save=save,hd=hd,
             cmap=orig_cmap, font=12, midpoint=0.5, start=0, stop=1,
             axis_lim=False,x_min=0,x_max=100,y_min=0,y_max=100,
             lineatzero=True,zerolineoffset=0)
     
    ###Make a quiver with smaller size

    new_XE=np.zeros((int(len(X_E)/3),int(len(X_E[0])/6)))
    new_YE=np.zeros_like(new_XE)
    new_Ex=np.zeros_like(new_XE)
    new_Ey=np.zeros_like(new_XE)
    count=0
    county=0
    for i in range(len(new_XE)):
        county=0
        for j in range(len(new_XE[0])):
            new_XE[i,j]=X_ER[count,county]
            new_YE[i,j]=Y_ER[count,county]
            new_Ex[i,j]=Ex[count,county]
            new_Ey[i,j]=Ey[count,county]
            county=county+6
        count=count+3
    
    title_name='E field with lines'+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    
    plt.figure(figsize=(8,5))
    plt.pcolormesh(X_ER,Y_ER,E_total,cmap='jet')
    plt.colorbar()
    plt.title('E field '+plot_names,fontsize=12)
    plt.streamplot(X_ER,Y_ER,Ex_c,Ey_c,density=1.5,color=Ey_c,cmap='hot_r')
    #plt.quiver(new_XE,new_YE,new_Ex,new_Ey)
    if save:
        print('Saving at '+save_name)
        plt.savefig(save_name+'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300) 

    plt.tight_layout()
        
    plt.show()
    
    plt.figure(figsize=(8,5))
    plt.plot(X_R[int(len(y)/2)],V[int(len(y)/2)],'b')
    plt.title('Linecut of V at the 0 position '+plot_names,fontsize=12)
    plt.xlabel('x coordinate')
    plt.ylabel('Potential')
    save_name=savename+'Linecut of V at the 0 position '+str(Y_R[int(len(y)/2),0])+plot_names
    plt.tight_layout()
    if save:
        plt.savefig(save_name+'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300)
        
    plt.show()
    
    
    ### Make an object to store all the data
    data_cap=Object()
    
    setattr(data_cap,'charge_top',rho_top_left)
    setattr(data_cap,'charge_bottom',rho_bot)
    setattr(data_cap,'charge_top_right',rho_top_right)
    setattr(data_cap,'charge',rho)
    
    setattr(data_cap,'hbn_top',top_plate)
    setattr(data_cap,'hbn_bottom',bot_plate)
    setattr(data_cap,'dielectric_1',bot_epsr)
    setattr(data_cap,'dielectric_2',top_epsr)
    setattr(data_cap,'dielectric',epsilon)
    setattr(data_cap,'plate_size', 2*plate_start)
    setattr(data_cap,'cut_size', 2*cut_x)
       
    
    setattr(data_cap,'x',X)
    setattr(data_cap,'y',Y)
    setattr(data_cap,'x_E',X_E)
    setattr(data_cap,'y_E',Y_E)
    setattr(data_cap,'x_real',X_R)
    setattr(data_cap,'y_real',Y_R)
    setattr(data_cap,'x_Ereal',X_ER)
    setattr(data_cap,'y_Ereal',Y_ER)
    
    setattr(data_cap,'V',V)
    setattr(data_cap,'E_tot',E_total)
    setattr(data_cap,'Ex',Ex_c)
    setattr(data_cap,'Ey',Ey_c)
    
    setattr(data_cap,'convergence_requirement',conv_error)
    setattr(data_cap,'name',plot_names)
    if save_obj:
        save_obj_filename=savename+plot_names+'dataobj.pkl'
        with open (save_obj_filename, 'wb') as output: # Overwrite any existing file
            pickle.dump(data_cap,output,-1)
    
    return data_cap





'''
Code incluying the environment area
'''
def capacitors_Dirich_env(Nx,Ny,ds,dsy,rho_top_left,rho_bot,bot_plate,
               top_plate,bot_epsr,top_epsr,plate_start_x,
               cut_x,conv_error=10E-6,
               plot_names='', save=False,savename='',top_size=1,
               bot_size=1,hd=0,mid_split=False,
               rho_top_right=0,
               max_iter=1E4,cut_center=0,save_obj=False,onetop=False, notops=False,
               triplegate=False,triple_plate=0,rho_triple=0,
               ):
    
    print('Save settings is set to '+str(save))
    if cut_center==0:
        cut_center=int(Nx/2)
    
    if mid_split==False: 
        rho_top_right=rho_top_left
    
    ## The y-size for the proper simulations should be the exact size of the
    ## capacitor. Otherwise there are some artifact from the way it "grounds"
    ## the environment. 
    Ny=top_plate+bot_plate+top_size+bot_size
    if Ny%2!=0:
        Ny=Ny+1
    
    x=np.arange(0,Nx,ds)
    y=np.arange(0,Ny,dsy)
    x_E=np.arange(0,Nx-ds,ds)
    y_E=np.arange(0,Ny-dsy,dsy)
    
    x_real=np.arange(-int(Nx/2),int(Nx/2)+ds,ds) #arange leaves out 1 values
    y_real=np.arange(-int(Ny/2),int(Ny/2)+dsy,dsy)
    x_Ereal=np.arange(-int(Nx/2),int(Nx/2),ds)
    y_Ereal=np.arange(-int(Ny/2),int(Ny/2),dsy)
    
    X,Y=np.meshgrid(x,y)
    X_E,Y_E=np.meshgrid(x_E,y_E)
    X_R,Y_R=np.meshgrid(x_real,y_real)
    X_ER,Y_ER=np.meshgrid(x_Ereal,y_Ereal)
    
    # Make the matrices to be filled. In np.zeros, the order is inverted
    rho=np.zeros((len(y),len(x)))
    epsilon=np.zeros((len(y),len(x)))+1
    # Make the initial guess for out matrix --> All 0s

    V=np.zeros((len(y),len(x)))
    R=np.zeros((len(y),len(x)))
    
    top_plate_size=np.arange(top_plate-top_size,top_plate+top_size+1)
    bot_plate_size=np.arange(bot_plate-bot_size,bot_plate+bot_size+1)
    plate_len=np.arange(plate_start_x,len(x)-plate_start_x)
    triple_plate_size=np.arange(triple_plate-bot_size,triple_plate+bot_size+1)

    ### Give the initial conditions. The charge distribution and the voltage of the electrodes 

    for i in plate_len:
        if cut_x==0.0:
            for j in top_plate_size:
                rho[j,i] = rho_top_left
                V[j,i]=rho_top_left    # Give the electrode the charge and not change it 
                
        else:
            if i <(cut_center-cut_x):
                for j in top_plate_size:
                    if notops:
                        rho[j,i] = 0
                        V[j,i]=0
                    else:
                        rho[j,i] = rho_top_left
                        V[j,i]=rho_top_left
            if i >(cut_center+cut_x):
                for j in top_plate_size:
                    if onetop==True or notops==True:
                        rho[j,i]=0
                        V[j,i]=0
                    else:
                        rho[j,i]=rho_top_right
                        V[j,i]=rho_top_right
                
        ### Same story as for the top plate, but now for the bottom plate
        
        for j in bot_plate_size:
            rho[j,i] = rho_bot
            V[j,i]=rho_bot
                
        if triplegate:
            for j in triple_plate_size:
                rho[j,i] = rho_triple
                V[j,i]=rho_triple


    ### Give the dielectric constant of the hBN inside the gates 
    for i in range(plate_start_x,len(x)-plate_start_x):
        for j in range(bot_plate,int(len(y)/2),1):
            epsilon[j,i] = bot_epsr  
        for j in range(int(len(y)/2),top_plate,1):
            epsilon[j,i] = top_epsr

    # Over relaxation constant. 0<w<2, if w is over 2 it will not converge

    t=np.cos(np.pi/Nx)+np.cos(np.pi/Ny)
    w=(8-np.sqrt(64-16*t**2))/t**2

    print('Over relaxation constant is: '+str(w))

    ### Solver ### 

    iterations = 0
    #eps = 10E-6 # Convergence threshold

    error = 1e4 # Large dummy error
    V_temp = np.copy(V)
    ### BCs are Dirichlet type, 0 for the boundaries and the charge for the electrodes
    
    while iterations < max_iter and error > conv_error:

        #V_temp = np.copy(V)
        error = 0 # we make this accumulate in the loop
        
        ### Left side of the electrodes and below the bottom electrode
        ### Keep the BCs always at zero voltage and the contacts at the given voltage 
  
        for i in range (1,len(x)-1):
            for j in range (1,len(y)-1):   
                # Calculate the permittivity matrices
                a0=epsilon[j,i]+epsilon[j-1,i]+epsilon[j,i-1]+epsilon[j-1,i-1]
                a1=0.5*(epsilon[j,i]+epsilon[j,i-1])
                a2=0.5*(epsilon[j-1,i]+epsilon[j,i])
                a3=0.5*(epsilon[j-1,i-1]+epsilon[j-1,i])
                a4=0.5*(epsilon[j,i-1]+epsilon[j-1,i-1])
                
                if j in bot_plate_size and i>plate_start_x and i < len(x)-plate_start_x:
                    R[j,i]=0
                    V[j,i]=rho_bot
                elif j in top_plate_size and i>plate_start_x and i < cut_center-cut_x and notops==False:
                    R[j,i]=0
                    V[j,i]=rho_top_left

                elif j in top_plate_size and i> cut_center+cut_x  and i < len(x)-plate_start_x and onetop==False and notops==False:
                    R[j,i]=0
                    V[j,i]=rho_top_right
                elif j in triple_plate_size and i>plate_start_x and i < len(x)-plate_start_x and triplegate==True:
                    R[j,i]=0
                    V[j,i]=rho_triple
                else: 
                    #Calculate the residual
                    R[j,i]=1/a0*(a1*V_temp[j+1,i] + a2*V_temp[j,i+1] +
                        a3*V_temp[j-1,i] + a4*V_temp[j,i-1] +rho[j,i])-V_temp[j,i]
                    
                    # Compute V using the overrelaxation method

                    V[j,i]=(V_temp[j,i]+w*R[j,i])
                    error += abs(V[j,i]-V_temp[j,i])

                    V_temp[j,i]=V[j,i] # Make it already work for the next iteration 

        iterations += 1
        error /= float(len(y))
    print (" iterations ="+str(iterations))
        
        
    ### Calculate the electric field
    ### it will always have one less component than the potentials, 
    ### just by the way it is defined
    Ex=np.zeros((len(y)-1,len(x)-1))
    Ey=np.zeros_like(Ex)
    Ey=np.zeros_like(Ex)
    Ex_c=np.zeros_like(Ex)
    Ey_c=np.zeros_like(Ex)
    E_total=np.zeros_like(Ex)

    ### This shifts the x and y with respect to each other
    
    for j in range (0,len(x) -1):
            for i in range (0,len(y) -1):

                Ex[i,j]=-(V[i,j+1]-V[i,j])/ds
                Ey[i,j]=-(V[i+1,j]-V[i,j])/dsy

                # Now we need to fix these shifts
    for j in range (0,len(x) -2):
            for i in range (0,len(y)-2):

                Ex_c[i,j]=0.5*(Ex[i+1,j]+Ex[i,j])
                Ey_c[i,j]=0.5*(Ey[i,j+1]+Ey[i,j])

                E_total[i,j]=np.sqrt(Ex_c[i,j]**2+Ey_c[i,j]**2)


    plt.figure(figsize=(8,5))
    plt.pcolormesh(X_R,Y_R,rho,cmap='jet')
    plt.colorbar()
    plt.title('Charges '+plot_names,fontsize=12)
    #savename='Charges'+plot_names
    
    plt.show()

    title_name='Voltage potentials '+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    basic_plotting.plot_map(X_R, Y_R, V, 8, 5, title_name, 
             x_label,y_label, z_label='', show=True, minmax=False,
             zmin=0, zmax=10, savename=save_name, save=save,hd=hd,
             cmap=orig_cmap, font=12, midpoint=0.5, start=0, stop=1,
             axis_lim=False,x_min=0,x_max=100,y_min=0,y_max=100,
             lineatzero=True,zerolineoffset=0)
     
    ###Make a quiver with smaller size

    new_XE=np.zeros((int(len(X_E)/3),int(len(X_E[0])/6)))
    new_YE=np.zeros_like(new_XE)
    new_Ex=np.zeros_like(new_XE)
    new_Ey=np.zeros_like(new_XE)
    count=0
    county=0
    for i in range(len(new_XE)):
        county=0
        for j in range(len(new_XE[0])):
            new_XE[i,j]=X_ER[count,county]
            new_YE[i,j]=Y_ER[count,county]
            new_Ex[i,j]=Ex[count,county]
            new_Ey[i,j]=Ey[count,county]
            county=county+6
        count=count+3
    
    title_name='E field with lines'+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    
    plt.figure(figsize=(8,5))
    plt.pcolormesh(X_ER,Y_ER,E_total,cmap='jet')
    plt.colorbar()
    plt.title('E field '+plot_names,fontsize=12)
    plt.streamplot(X_ER,Y_ER,Ex_c,Ey_c,density=1.5,color=Ey_c,cmap='hot_r')
    #plt.quiver(new_XE,new_YE,new_Ex,new_Ey)
    if save:
        print('Saving at '+save_name)
        plt.savefig(save_name+'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300) 

    plt.tight_layout()
        
    plt.show()
    
    plt.figure(figsize=(8,5))
    plt.plot(X_R[int(len(y)/2)],V[int(len(y)/2)],'b')
    plt.title('Linecut of V at the 0 position '+plot_names,fontsize=12)
    plt.xlabel('x coordinate')
    plt.ylabel('Potential')
    save_name=savename+'Linecut of V at the 0 position '+str(Y_R[int(len(y)/2),0])+plot_names
    plt.tight_layout()
    if save:
        plt.savefig(save_name+'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300)
        
    plt.show()
    
    
    ### Make an object to store all the data
    data_cap=Object()
    
    setattr(data_cap,'charge_top',rho_top_left)
    setattr(data_cap,'charge_bottom',rho_bot)
    setattr(data_cap,'charge_top_right',rho_top_right)
    setattr(data_cap,'charge',rho)
    
    setattr(data_cap,'hbn_top',top_plate)
    setattr(data_cap,'hbn_bottom',bot_plate)
    setattr(data_cap,'dielectric_1',bot_epsr)
    setattr(data_cap,'dielectric_2',top_epsr)
    setattr(data_cap,'dielectric',epsilon)
    setattr(data_cap,'plate_size', 2*plate_start)
    setattr(data_cap,'cut_size', 2*cut_x)
       
    
    setattr(data_cap,'x',X)
    setattr(data_cap,'y',Y)
    setattr(data_cap,'x_E',X_E)
    setattr(data_cap,'y_E',Y_E)
    setattr(data_cap,'x_real',X_R)
    setattr(data_cap,'y_real',Y_R)
    setattr(data_cap,'x_Ereal',X_ER)
    setattr(data_cap,'y_Ereal',Y_ER)
    
    setattr(data_cap,'V',V)
    setattr(data_cap,'E_tot',E_total)
    setattr(data_cap,'Ex',Ex_c)
    setattr(data_cap,'Ey',Ey_c)
    
    setattr(data_cap,'convergence_requirement',conv_error)
    setattr(data_cap,'name',plot_names)
    if save_obj:
        save_obj_filename=savename+plot_names+'dataobj.pkl'
        with open (save_obj_filename, 'wb') as output: # Overwrite any existing file
            pickle.dump(data_cap,output,-1)
    
    return data_cap
    

#### 
def capacitors(Nx,Ny,ds,dsy,rho_top_left,rho_bot,rho_mid,bot_plate,
               top_plate,bot_epsr,top_epsr,plate_start_x,
               cut_x,conv_error=10E-6,
               BCs=0, plot_names='', save=False,savename='',top_size=1,
               bot_size=1,hd=0,mid=False,mid_split=False,
               rho_top_right=0,
               max_iter=1E4,cut_center=0,save_obj=False,
               ):
    
    print('Save settings is set to '+str(save))
    if cut_center==0:
        cut_center=int(Nx/2)
    
    if mid_split==False: 
        rho_top_right=rho_top_left
    
    x=np.arange(0,Nx,ds)
    y=np.arange(0,Ny,dsy)
    x_E=np.arange(0,Nx-ds,ds)
    y_E=np.arange(0,Ny-dsy,dsy)
    
    x_real=np.arange(-int(Nx/2),int(Nx/2)+ds,ds) #arange leaves out 1 values
    y_real=np.arange(-int(Ny/2),int(Ny/2)+dsy,dsy)
    x_Ereal=np.arange(-int(Nx/2),int(Nx/2),ds)
    y_Ereal=np.arange(-int(Ny/2),int(Ny/2),dsy)
    
    X,Y=np.meshgrid(x,y)
    X_E,Y_E=np.meshgrid(x_E,y_E)
    X_R,Y_R=np.meshgrid(x_real,y_real)
    X_ER,Y_ER=np.meshgrid(x_Ereal,y_Ereal)
    
    # Make the matrices to be filled. In np.zeros, the order is inverted
    rho=np.zeros((len(y),len(x)))
    epsilon=np.zeros((len(y),len(x)))+1
    # Make the initial guess for out matrix --> All 0s

    V=np.zeros((len(y),len(x)))
    R=np.zeros((len(y),len(x)))

    for i in range(plate_start_x,len(x)-plate_start_x):

        ### Top plate charge is given. If we assume that it is floating
        ### then we should give it a charge of 0
        ### If we assume that it is GNDed, then it should react to the 
        ### bottom plate charge. 
        
        if cut_x==0.0:
            for j in [top_plate-top_size,top_plate,top_plate+top_size]:
                rho[j,i] = rho_top_left
        else:
            if i <(cut_center-cut_x):
                for j in [top_plate-top_size,top_plate,top_plate+top_size]:
                    rho[j,i] = rho_top_left
            elif i >(cut_center+cut_x):
                for j in [top_plate-top_size,top_plate,top_plate+top_size]:
                    rho[j,i]=rho_top_right
            else: 
                rho[j,i]=0
                
        ### Same story as for the top plate, but now for the bottom plate
        
        for j in [bot_plate-bot_size,bot_plate,bot_plate+bot_size]:
            rho[j,i] = rho_bot
        ### This would be the graphene and it is a bit more tricky 
        ### the charge might react to both plates, but this should be specified 
        ### in the running code, not here. 
        if mid==True:
            for j in [int(len(y)/2)]:
                ### The response is normalized by the bottom gate 
                rho[j,i]=-(rho[top_plate,i]*(top_epsr/bot_epsr)+rho[bot_plate,i])
            #for j in [int(len(y)/2)-1]:
            #    rho[j,i]=rho_midbot
            
        ## To make sure that there is no crazy fields outside of the capacitor, 
        ## which is actually not realistic
        '''
        for j in [bot_plate-1,top_plate+1]:
            rho[j,i]=0
        '''
    for i in range(plate_start_x,len(x)-plate_start_x):
        for j in range(bot_plate,int(len(y)/2),1):
            epsilon[j,i] = bot_epsr  
        for j in range(int(len(y)/2),top_plate,1):
            epsilon[j,i] = top_epsr

    # Over relaxation constant. 0<w<2, if w is over 2 it will not converge

    t=np.cos(np.pi/Nx)+np.cos(np.pi/Ny)
    w=(8-np.sqrt(64-16*t**2))/t**2

    print('Over relaxation constant is: '+str(w))

    # Solver

    iterations = 0
    #eps = 10E-6 # Convergence threshold

    error = 1e4 # Large dummy error

    ### BCs of value 0 means Dirichlet BCs 
    
    if BCs==0:
        while iterations < max_iter and error > conv_error:

            V_temp = np.copy(V)
            error = 0 # we make this accumulate in the loop
            for j in range (1,len(x)-1):
                for i in range (1,len(y)-1):
                    # Calculate the permittivity matrices
                    a0=epsilon[i,j]+epsilon[i-1,j]+epsilon[i,j-1]+epsilon[i-1,j-1]
                    a1=0.5*(epsilon[i,j]+epsilon[i,j-1])
                    a2=0.5*(epsilon[i-1,j]+epsilon[i,j])
                    a3=0.5*(epsilon[i-1,j-1]+epsilon[i-1,j])
                    a4=0.5*(epsilon[i,j-1]+epsilon[i-1,j-1])

                    #Calculate the residual
                    R[i,j]=1/a0*(a1*V_temp[i+1,j] + a2*V_temp[i,j+1] +
                        a3*V_temp[i-1,j] + a4*V_temp[i,j-1] +rho[i,j])-V_temp[i,j]

                    # Compute V using the overrelaxation method

                    V[i,j]=V_temp[i,j]+w*R[i,j]
                    error += abs(V[i,j]-V_temp[i,j])

                    V_temp[i,j]=V[i,j] # Make it already work for the next iteration 

            iterations += 1
            error /= float(len(y))
        print (" iterations ="+str(iterations))
        
        
     # BCs of 1 is a Neumann BC
        
    if BCs==1:  
        while iterations < max_iter and error > conv_error:

            V_temp = np.copy(V)
            error = 0 # we make this accumulate in the loop
            for j in range (0,len(x)-1):
                for i in range (0,len(y)-1):
                    
                    # Calculate the permittivity matrices.
                    ### Add he Neumann BCs
                    if j == 0 or j == len(x):
                        if i == 0 or i == len(y):
                            epsilon[i,j]=epsilon[i+1,j]
                            epsilon[i-1,j]=epsilon[i,j]
                        epsilon[i,j-1]=epsilon[i,j]
                        epsilon[i-1,j-1]=epsilon[i-1,j]

                    a0=epsilon[i,j]+epsilon[i-1,j]+epsilon[i,j-1]+epsilon[i-1,j-1]
                    a1=0.5*(epsilon[i,j]+epsilon[i,j-1])
                    a2=0.5*(epsilon[i-1,j]+epsilon[i,j])
                    a3=0.5*(epsilon[i-1,j-1]+epsilon[i-1,j])
                    a4=0.5*(epsilon[i,j-1]+epsilon[i-1,j-1])

                    #Calculate the residual
                    ### Add he Neumann BCs
                    if j == 0:
                        V_temp[i,j-1]=V_temp[i,j]
                    if i==0:
                        V_temp[i-1,j]=V_temp[i,j]
                    if j == Nx-2:
                        V_temp[i,j+1]=V_temp[i,j]
                    if i==Ny-2:
                        V_temp[i+1,j]=V_temp[i,j]
                            
                    R[i,j]=1/a0*(a1*V_temp[i+1,j] + a2*V_temp[i,j+1] +
                        a3*V_temp[i-1,j] + a4*V_temp[i,j-1] +rho[i,j])-V_temp[i,j]

                    # Compute V using the overrelaxation method

                    V[i,j]=V_temp[i,j]+w*R[i,j]
                    error += abs(V[i,j]-V_temp[i,j])

                    V_temp[i,j]=V[i,j] # Make it already work for the next iteration               
                    
       
            iterations += 1
            error /= float(len(y))
        print (" iterations ="+str(iterations))
        
    ### Calculate the electric field
    ### it will always have one less component than the potentials, 
    ### just by the way it is defined
    Ex=np.zeros((len(y)-1,len(x)-1))
    Ey=np.zeros_like(Ex)
    Ey=np.zeros_like(Ex)
    Ex_c=np.zeros_like(Ex)
    Ey_c=np.zeros_like(Ex)
    E_total=np.zeros_like(Ex)

    ### This shifts the x and y with respect to each other
    
    for j in range (0,len(x) -1):
            for i in range (0,len(y) -1):

                Ex[i,j]=-(V[i,j+1]-V[i,j])/ds
                Ey[i,j]=-(V[i+1,j]-V[i,j])/dsy

                # Now we need to fix these shifts
    for j in range (0,len(x) -2):
            for i in range (0,len(y)-2):

                Ex_c[i,j]=0.5*(Ex[i+1,j]+Ex[i,j])
                Ey_c[i,j]=0.5*(Ey[i,j+1]+Ey[i,j])

                E_total[i,j]=np.sqrt(Ex_c[i,j]**2+Ey_c[i,j]**2)

    plt.figure(figsize=(12,7))
    plt.pcolormesh(X_R,Y_R,rho,cmap='jet')
    plt.colorbar()
    plt.title('Charges '+plot_names,fontsize=12)
    savename='Charges'+plot_names
    
    plt.show()
    title_name='Voltage potentials '+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    basic_plotting.plot_map(X_R, Y_R, V, 12, 7, title_name, 
             x_label,y_label, z_label='', show=True, minmax=False,
             zmin=0, zmax=10, savename=save_name, save=save,hd=hd,
             cmap=orig_cmap, font=12, midpoint=0.5, start=0, stop=1,
             axis_lim=False,x_min=0,x_max=100,y_min=0,y_max=100,
             lineatzero=True,zerolineoffset=0)
    '''
    plt.figure(figsize=(12,7))
    plt.pcolormesh(X_R,Y_R,V,cmap='jet')
    plt.colorbar()
    plt.title('Voltage potentials'+plot_names,fontsize=12)
    savename='Voltage potentials'+plot_names
    
    plt.show()
    
    '''
    # Plot the total magnitude of the electric field
    '''
    title_name='E field '+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    basic_plotting.plot_map(X_ER, Y_ER, E_total, 12, 7, title_name, 
             x_label,y_label, z_label='', show=True, minmax=False,
             zmin=0, zmax=10, savename=save_folder+save_name, save=save,hd=hd,
             cmap=orig_cmap, font=12, midpoint=0.5, start=0, stop=1)
    
    '''
    ###Make a quiver with smaller size

    new_XE=np.zeros((int(len(X_E)/3),int(len(X_E[0])/6)))
    new_YE=np.zeros_like(new_XE)
    new_Ex=np.zeros_like(new_XE)
    new_Ey=np.zeros_like(new_XE)
    count=0
    county=0
    for i in range(len(new_XE)):
        county=0
        for j in range(len(new_XE[0])):
            new_XE[i,j]=X_ER[count,county]
            new_YE[i,j]=Y_ER[count,county]
            new_Ex[i,j]=Ex[count,county]
            new_Ey[i,j]=Ey[count,county]
            county=county+6
        count=count+3
    
    title_name='E field with lines'+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    
    plt.figure(figsize=(12,7))
    plt.pcolormesh(X_ER,Y_ER,E_total,cmap='jet')
    plt.colorbar()
    plt.title('E field '+plot_names,fontsize=12)
    plt.streamplot(X_ER,Y_ER,Ex_c,Ey_c,density=1.5,color=Ey_c,cmap='hot_r')
    #plt.quiver(new_XE,new_YE,new_Ex,new_Ey)
    if save:
        if hd==0:
            plt.savefig(save_name+'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300) 
        elif hd==1:
            plt.savefig(save_name +'_'+time.strftime("%H%M%d%m%y")+'.png', dpi=600)
        elif hd==2:
            plt.savefig(save_name +'_'+time.strftime("%H%M%d%m%y")+'.pdf')   

    plt.tight_layout()
        
    plt.show()
    
    
    # plot inverse of electric field, which is related to capacitance
    '''
    title_name='Capacitance '+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    basic_plotting.plot_map(X_ER, Y_ER, 1/E_total, 12, 7, title_name, 
             x_label,y_label, z_label='', show=True, minmax=False,
             zmin=0, zmax=10, savename=save_folder+save_name, save=save,hd=hd,
             cmap=orig_cmap, font=12, midpoint=0.5, start=0, stop=1)
    '''
    '''
    plt.figure(figsize=(12,7))
    plt.title('Capacitance model $C\propto 1/E$'+plot_names,fontsize=12)
    plt.pcolormesh(X_ER,Y_ER,1/E_total,cmap='jet',vmax=5)
    plt.colorbar()
    savename='Capacitance '+plot_names
    if save:
        if hd==0:
            plt.savefig(savename+'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300) 
        elif hd==1:
            plt.savefig(savename +'_'+time.strftime("%H%M%d%m%y")+'.png', dpi=600)
        elif hd==2:
            plt.savefig(savename +'_'+time.strftime("%H%M%d%m%y")+'.pdf')   

        np.savetxt(savename+time.strftime("%H%M%d%m%y")+'_metadata.txt',saveData, delimiter='\t', header = saveheader, comments='')
        
    plt.show()
    '''
    
    plt.figure(figsize=(12,7))
    plt.plot(X_R[int(len(y)/2)],V[int(len(y)/2)],'b')
    plt.title('Linecut of V at the 0 position '+plot_names,fontsize=12)
    plt.xlabel('x coordinate')
    plt.ylabel('Potential')
    save_name=savename+'Linecut of V at the 0 position '+str(Y_R[int(len(y)/2),0])+plot_names
    plt.tight_layout()
    if save:
        plt.savefig(save_name+'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300)
        
    plt.show()
    
    
    ### Make an object to store all the data
    data_cap=Object()
    

    setattr(data_cap,'charge_top',rho_top_left)
    setattr(data_cap,'charge_bottom',rho_bot)
    setattr(data_cap,'charge_top_right',rho_top_right)
    setattr(data_cap,'charge',rho)
    
    setattr(data_cap,'hbn_top',top_plate)
    setattr(data_cap,'hbn_bottom',bot_plate)
    setattr(data_cap,'dielectric_1',bot_epsr)
    setattr(data_cap,'dielectric_2',top_epsr)
    setattr(data_cap,'dielectric',epsilon)
    setattr(data_cap,'plate_size', 2*plate_start)
    setattr(data_cap,'cut_size', 2*cut_x)
       
    
    setattr(data_cap,'x',X)
    setattr(data_cap,'y',Y)
    setattr(data_cap,'x_E',X_E)
    setattr(data_cap,'y_E',Y_E)
    setattr(data_cap,'x_real',X_R)
    setattr(data_cap,'y_real',Y_R)
    setattr(data_cap,'x_Ereal',X_ER)
    setattr(data_cap,'y_Ereal',Y_ER)
    
    setattr(data_cap,'V',V)
    setattr(data_cap,'E_tot',E_total)
    setattr(data_cap,'Ex',Ex_c)
    setattr(data_cap,'Ey',Ey_c)
    
    setattr(data_cap,'convergence_requirement',conv_error)
    setattr(data_cap,'plot_names',plot_names)
    if save_obj:
        save_obj_filename=savename+plot_names+'dataobj.pkl'
        with open (save_obj_filename, 'wb') as output: # Overwrite any existing file
            pickle.dump(data_cap,output,-1)
    
    return data_cap
    

def capacitors_equipotential(Nx,Ny,ds,dsy,rho_top_left,rho_bot,rho_mid,bot_plate,
               top_plate,bot_epsr,top_epsr,plate_start_x,
               cut_x,conv_error=10E-6,
               BCs=0, plot_names='', save=False,savename='',top_size=1,
               bot_size=1,hd=0,mid=False,mid_split=False,
               rho_top_right=0,
               max_iter=1E4,cut_center=0,save_obj=False,
               ):
    
    print('Save settings is set to '+str(save))
    if cut_center==0:
        cut_center=int(Nx/2)
    
    if mid_split==False: 
        rho_top_right=rho_top_left
    
    x=np.arange(0,Nx,ds)
    y=np.arange(0,Ny,dsy)
    x_E=np.arange(0,Nx-ds,ds)
    y_E=np.arange(0,Ny-dsy,dsy)
    
    x_real=np.arange(-int(Nx/2),int(Nx/2)+ds,ds) #arange leaves out 1 values
    y_real=np.arange(-int(Ny/2),int(Ny/2)+dsy,dsy)
    x_Ereal=np.arange(-int(Nx/2),int(Nx/2),ds)
    y_Ereal=np.arange(-int(Ny/2),int(Ny/2),dsy)
    
    X,Y=np.meshgrid(x,y)
    X_E,Y_E=np.meshgrid(x_E,y_E)
    X_R,Y_R=np.meshgrid(x_real,y_real)
    X_ER,Y_ER=np.meshgrid(x_Ereal,y_Ereal)
    
    # Make the matrices to be filled. In np.zeros, the order is inverted
    rho=np.zeros((len(y),len(x)))
    epsilon=np.zeros((len(y),len(x)))+1
    # Make the initial guess for out matrix --> All 0s

    V=np.zeros((len(y),len(x)))
    R=np.zeros((len(y),len(x)))

    for i in range(plate_start_x,len(x)-plate_start_x):

        ### Top plate charge is given. If we assume that it is floating
        ### then we should give it a charge of 0
        ### If we assume that it is GNDed, then it should react to the 
        ### bottom plate charge. 
        
        if cut_x==0.0:
            for j in [top_plate-top_size,top_plate,top_plate+top_size]:
                rho[j,i] = rho_top_left
        else:
            if i <(cut_center-cut_x):
                for j in [top_plate-top_size,top_plate,top_plate+top_size]:
                    rho[j,i] = rho_top_left
            elif i >(cut_center+cut_x):
                for j in [top_plate-top_size,top_plate,top_plate+top_size]:
                    rho[j,i]=rho_top_right
            else: 
                rho[j,i]=0
                
        ### Same story as for the top plate, but now for the bottom plate
        
        for j in [bot_plate-bot_size,bot_plate,bot_plate+bot_size]:
            rho[j,i] = rho_bot
            
        ### This would be the graphene and it is a bit more tricky 
        ### the charge might react to both plates, but this should be specified 
        ### in the running code, not here. 
        if mid==True:
            for j in [int(len(y)/2)]:
                ### The response is normalized by the bottom gate 
                rho[j,i]=-(rho[top_plate,i]*(top_epsr/bot_epsr)+rho[bot_plate,i])
            #for j in [int(len(y)/2)-1]:
            #    rho[j,i]=rho_midbot
            
        ## To make sure that there is no crazy fields outside of the capacitor, 
        ## which is actually not realistic
        '''
        for j in [bot_plate-1,top_plate+1]:
            rho[j,i]=0
        '''
    for i in range(plate_start_x,len(x)-plate_start_x):
        for j in range(bot_plate,int(len(y)/2),1):
            epsilon[j,i] = bot_epsr  
        for j in range(int(len(y)/2),top_plate,1):
            epsilon[j,i] = top_epsr

    # Over relaxation constant. 0<w<2, if w is over 2 it will not converge

    t=np.cos(np.pi/Nx)+np.cos(np.pi/Ny)
    w=(8-np.sqrt(64-16*t**2))/t**2

    print('Over relaxation constant is: '+str(w))

    # Solver

    iterations = 0
    #eps = 10E-6 # Convergence threshold

    error = 1e4 # Large dummy error

    ### BCs of value 0 means Dirichlet BCs 
    
    if BCs==0:
        while iterations < max_iter and error > conv_error:

            V_temp = np.copy(V)
            error = 0 # we make this accumulate in the loop
            for j in range (1,len(x)-1):
                for i in range (1,len(y)-1):
                    # Calculate the permittivity matrices
                    a0=epsilon[i,j]+epsilon[i-1,j]+epsilon[i,j-1]+epsilon[i-1,j-1]
                    a1=0.5*(epsilon[i,j]+epsilon[i,j-1])
                    a2=0.5*(epsilon[i-1,j]+epsilon[i,j])
                    a3=0.5*(epsilon[i-1,j-1]+epsilon[i-1,j])
                    a4=0.5*(epsilon[i,j-1]+epsilon[i-1,j-1])

                    #Calculate the residual
                    R[i,j]=1/a0*(a1*V_temp[i+1,j] + a2*V_temp[i,j+1] +
                        a3*V_temp[i-1,j] + a4*V_temp[i,j-1] +rho[i,j])-V_temp[i,j]

                    # Compute V using the overrelaxation method

                    V[i,j]=V_temp[i,j]+w*R[i,j]
                    error += abs(V[i,j]-V_temp[i,j])

                    V_temp[i,j]=V[i,j] # Make it already work for the next iteration 
                    
                    current_pos=i
                    if i in [top_plate-top_size,top_plate,top_plate+top_size]:
                        for j in range(plate_start_x,cut_center-cut_x):
                            V_temp[i,j]=V[current_pos,j]
                        for j in range(cut_center+cut_x,len(x)-plate_start_x):
                            V_temp[i,j]=V[current_pos,j]
                            
                    elif i in [bot_plate-bot_size,bot_plate,bot_plate+bot_size]:
                        for j in range(plate_start_x,len(x)-plate_start_x):
                            V_temp[i,j]=V[current_pos,j]
 
            iterations += 1
            error /= float(len(y))
        print (" iterations ="+str(iterations))
        
        
     # BCs of 1 is a Neumann BC
        
    if BCs==1:  
        while iterations < max_iter and error > conv_error:

            V_temp = np.copy(V)
            error = 0 # we make this accumulate in the loop
            for j in range (0,len(x)-1):
                for i in range (0,len(y)-1):
                    
                    # Calculate the permittivity matrices.
                    ### Add he Neumann BCs
                    if j == 0 or j == len(x):
                        if i == 0 or i == len(y):
                            epsilon[i,j]=epsilon[i+1,j]
                            epsilon[i-1,j]=epsilon[i,j]
                        epsilon[i,j-1]=epsilon[i,j]
                        epsilon[i-1,j-1]=epsilon[i-1,j]

                    a0=epsilon[i,j]+epsilon[i-1,j]+epsilon[i,j-1]+epsilon[i-1,j-1]
                    a1=0.5*(epsilon[i,j]+epsilon[i,j-1])
                    a2=0.5*(epsilon[i-1,j]+epsilon[i,j])
                    a3=0.5*(epsilon[i-1,j-1]+epsilon[i-1,j])
                    a4=0.5*(epsilon[i,j-1]+epsilon[i-1,j-1])

                    #Calculate the residual
                    ### Add he Neumann BCs
                    if j == 0:
                        V_temp[i,j-1]=V_temp[i,j]
                    if i==0:
                        V_temp[i-1,j]=V_temp[i,j]
                    if j == Nx-2:
                        V_temp[i,j+1]=V_temp[i,j]
                    if i==Ny-2:
                        V_temp[i+1,j]=V_temp[i,j]
                            
                    R[i,j]=1/a0*(a1*V_temp[i+1,j] + a2*V_temp[i,j+1] +
                        a3*V_temp[i-1,j] + a4*V_temp[i,j-1] +rho[i,j])-V_temp[i,j]

                    # Compute V using the overrelaxation method

                    V[i,j]=V_temp[i,j]+w*R[i,j]
                    error += abs(V[i,j]-V_temp[i,j])

                    V_temp[i,j]=V[i,j] # Make it already work for the next iteration               
                    
       
            iterations += 1
            error /= float(len(y))
        print (" iterations ="+str(iterations))
        
    ### Calculate the electric field
    ### it will always have one less component than the potentials, 
    ### just by the way it is defined
    Ex=np.zeros((len(y)-1,len(x)-1))
    Ey=np.zeros_like(Ex)
    Ey=np.zeros_like(Ex)
    Ex_c=np.zeros_like(Ex)
    Ey_c=np.zeros_like(Ex)
    E_total=np.zeros_like(Ex)

    ### This shifts the x and y with respect to each other
    
    for j in range (0,len(x) -1):
            for i in range (0,len(y) -1):

                Ex[i,j]=-(V[i,j+1]-V[i,j])/ds
                Ey[i,j]=-(V[i+1,j]-V[i,j])/dsy

                # Now we need to fix these shifts
    for j in range (0,len(x) -2):
            for i in range (0,len(y)-2):

                Ex_c[i,j]=0.5*(Ex[i+1,j]+Ex[i,j])
                Ey_c[i,j]=0.5*(Ey[i,j+1]+Ey[i,j])

                E_total[i,j]=np.sqrt(Ex_c[i,j]**2+Ey_c[i,j]**2)
    '''
    plt.figure(figsize=(12,7))
    plt.title('Voltage potentials')
    plt.pcolormesh(X,Y,V,cmap='jet')
    plt.colorbar()
    plt.show()

    '''
    title_name='Voltage potentials '+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    basic_plotting.plot_map(X_R, Y_R, V, 12, 7, title_name, 
             x_label,y_label, z_label='', show=True, minmax=False,
             zmin=0, zmax=10, savename=save_name, save=save,hd=hd,
             cmap=orig_cmap, font=12, midpoint=0.5, start=0, stop=1,
             axis_lim=False,x_min=0,x_max=100,y_min=0,y_max=100,
             lineatzero=True,zerolineoffset=0)
    '''
    plt.figure(figsize=(12,7))
    plt.pcolormesh(X_R,Y_R,V,cmap='jet')
    plt.colorbar()
    plt.title('Voltage potentials'+plot_names,fontsize=12)
    savename='Voltage potentials'+plot_names
    
    plt.show()
    
    '''
    # Plot the total magnitude of the electric field
    '''
    title_name='E field '+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    basic_plotting.plot_map(X_ER, Y_ER, E_total, 12, 7, title_name, 
             x_label,y_label, z_label='', show=True, minmax=False,
             zmin=0, zmax=10, savename=save_folder+save_name, save=save,hd=hd,
             cmap=orig_cmap, font=12, midpoint=0.5, start=0, stop=1)
    
    '''
    ###Make a quiver with smaller size

    new_XE=np.zeros((int(len(X_E)/3),int(len(X_E[0])/6)))
    new_YE=np.zeros_like(new_XE)
    new_Ex=np.zeros_like(new_XE)
    new_Ey=np.zeros_like(new_XE)
    count=0
    county=0
    for i in range(len(new_XE)):
        county=0
        for j in range(len(new_XE[0])):
            new_XE[i,j]=X_ER[count,county]
            new_YE[i,j]=Y_ER[count,county]
            new_Ex[i,j]=Ex[count,county]
            new_Ey[i,j]=Ey[count,county]
            county=county+6
        count=count+3
    
    title_name='E field with lines'+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    
    plt.figure(figsize=(12,7))
    plt.pcolormesh(X_ER,Y_ER,E_total,cmap='jet')
    plt.colorbar()
    plt.title('E field '+plot_names,fontsize=12)
    plt.streamplot(X_ER,Y_ER,Ex_c,Ey_c,density=1.5,color=Ey_c,cmap='hot_r')
    #plt.quiver(new_XE,new_YE,new_Ex,new_Ey)
    if save:
        if hd==0:
            plt.savefig(save_name+'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300) 
        elif hd==1:
            plt.savefig(save_name +'_'+time.strftime("%H%M%d%m%y")+'.png', dpi=600)
        elif hd==2:
            plt.savefig(save_name +'_'+time.strftime("%H%M%d%m%y")+'.pdf')   

    plt.tight_layout()
        
    plt.show()
    
    
    # plot inverse of electric field, which is related to capacitance
    '''
    title_name='Capacitance '+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    basic_plotting.plot_map(X_ER, Y_ER, 1/E_total, 12, 7, title_name, 
             x_label,y_label, z_label='', show=True, minmax=False,
             zmin=0, zmax=10, savename=save_folder+save_name, save=save,hd=hd,
             cmap=orig_cmap, font=12, midpoint=0.5, start=0, stop=1)
    '''
    '''
    plt.figure(figsize=(12,7))
    plt.title('Capacitance model $C\propto 1/E$'+plot_names,fontsize=12)
    plt.pcolormesh(X_ER,Y_ER,1/E_total,cmap='jet',vmax=5)
    plt.colorbar()
    savename='Capacitance '+plot_names
    if save:
        if hd==0:
            plt.savefig(savename+'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300) 
        elif hd==1:
            plt.savefig(savename +'_'+time.strftime("%H%M%d%m%y")+'.png', dpi=600)
        elif hd==2:
            plt.savefig(savename +'_'+time.strftime("%H%M%d%m%y")+'.pdf')   

        np.savetxt(savename+time.strftime("%H%M%d%m%y")+'_metadata.txt',saveData, delimiter='\t', header = saveheader, comments='')
        
    plt.show()
    '''
    
    plt.figure(figsize=(12,7))
    plt.plot(X_R[int(len(y)/2)],V[int(len(y)/2)],'b')
    plt.title('Linecut of V at the 0 position '+plot_names,fontsize=12)
    plt.xlabel('x coordinate')
    plt.ylabel('Potential')
    save_name=savename+'Linecut of V at the 0 position '+str(Y_R[int(len(y)/2),0])+plot_names
    plt.tight_layout()
    if save:
        plt.savefig(save_name+'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300)
        
    plt.show()
    
    
    ### Make an object to store all the data
    data_cap=Object()
    

    setattr(data_cap,'charge_top',rho_top_left)
    setattr(data_cap,'charge_bottom',rho_bot)
    setattr(data_cap,'charge_top_right',rho_top_right)
    setattr(data_cap,'charge',rho)
    
    setattr(data_cap,'hbn_top',top_plate)
    setattr(data_cap,'hbn_bottom',bot_plate)
    setattr(data_cap,'dielectric_1',bot_epsr)
    setattr(data_cap,'dielectric_2',top_epsr)
    setattr(data_cap,'dielectric',epsilon)
    setattr(data_cap,'plate_size', 2*plate_start)
    setattr(data_cap,'cut_size', 2*cut_x)
       
    
    setattr(data_cap,'x',X)
    setattr(data_cap,'y',Y)
    setattr(data_cap,'x_E',X_E)
    setattr(data_cap,'y_E',Y_E)
    setattr(data_cap,'x_real',X_R)
    setattr(data_cap,'y_real',Y_R)
    setattr(data_cap,'x_Ereal',X_ER)
    setattr(data_cap,'y_Ereal',Y_ER)
    
    setattr(data_cap,'V',V)
    setattr(data_cap,'E_tot',E_total)
    setattr(data_cap,'Ex',Ex_c)
    setattr(data_cap,'Ey',Ey_c)
    
    setattr(data_cap,'convergence_requirement',conv_error)
    setattr(data_cap,'plot_names',plot_names)
    if save_obj:
        save_obj_filename=savename+plot_names+'dataobj.pkl'
        with open (save_obj_filename, 'wb') as output: # Overwrite any existing file
            pickle.dump(data_cap,output,-1)
    
    return data_cap


def capacitors_half(Nx,Ny,ds,dsy,rho_top,rho_bot,rho_midtop,rho_midbot,bot_plate,top_plate,epsr1,epsr2,plate_start_x,
               cut_x,eps=10E-6,
               BCs=0, plot_names='', save=False,savename='',top_size=1,bot_size=1,hd=0,mid=False,
               max_iter=1E4,cut_center=0
               ):
    
    if cut_center==0:
        cut_center=0
    
    x=np.arange(0,Nx,ds)
    y=np.arange(0,Ny,dsy)
    x_E=np.arange(0,Nx-ds,ds)
    y_E=np.arange(0,Ny-dsy,dsy)
    
    x_real=np.arange(0,Nx,ds)
    y_real=np.arange(0,Ny,dsy)
    x_Ereal=np.arange(0,Nx-ds,ds)
    y_Ereal=np.arange(0,Ny-dsy,dsy)
    
    X,Y=np.meshgrid(x,y)
    X_E,Y_E=np.meshgrid(x_E,y_E)
    X_R,Y_R=np.meshgrid(x_real,y_real)
    X_ER,Y_ER=np.meshgrid(x_Ereal,y_Ereal)
    
    # Make the matrices to be filled. In np.zeros, the order is inverted
    rho=np.zeros((len(y),len(x)))
    epsilon=np.zeros((len(y),len(x)))+1
    # Make the initial guess for out matrix --> All 0s

    V=np.zeros((len(y),len(x)))
    R=np.zeros((len(y),len(x)))

    for i in range(0,plate_start_x):

        ### Top plate charge is given. If we assume that it is floating
        ### then we should give it a charge of 0
        ### If we assume that it is GNDed, then it should react to the 
        ### bottom plate charge. 
        
        for j in [top_plate-top_size,top_plate,top_plate+top_size]:
            if i >(cut_center+cut_x):
                rho[j,i] = rho_top
            else: 
                rho[j,i]=0
                
        ### Same story as for the top plate, but now for the bottom plate
        
        for j in [bot_plate-bot_size,bot_plate,bot_plate+bot_size]:
            rho[j,i] = rho_bot
            
        ### This would be the graphene and it is a bit more tricky 
        ### the charge might react to both plates, but this should be specified 
        ### in the running code, not here. 
        if mid==True:
            for j in [int(len(y)/2)]:
                rho[j,i]=rho_midtop
            for j in [int(len(y)/2)-1]:
                rho[j,i]=rho_midbot
            
        ## To make sure that there is no crazy fields outside of the capacitor, 
        ## which is actually not realistic
        '''
        for j in [bot_plate-1,top_plate+1]:
            rho[j,i]=0
        '''
    for i in range(0,plate_start_x):
        
        for j in range(bot_plate,int(len(y)/2),1):
            epsilon[j,i] = epsr1   
            
        for j in range(int(len(y)/2),top_plate,1):
            epsilon[j,i] = epsr2 

    # Over relaxation constant. 0<w<2, if w is over 2 it will not converge

    t=np.cos(np.pi/Nx)+np.cos(np.pi/Ny)
    w=(8-np.sqrt(64-16*t**2))/t**2

    print(w)

    # Solver

    iterations = 0
    #eps = 10E-6 # Convergence threshold

    error = 1e4 # Large dummy error

    ### BCs of value 0 means Dirichlet BCs 
    
    if BCs==0:
        while iterations < max_iter and error > eps:

            V_temp = np.copy(V)
            error = 0 # we make this accumulate in the loop
            for j in range (0,len(x)-1):
                for i in range (0,len(y)-1):
                    # Calculate the permittivity matrices
                    a0=epsilon[i,j]+epsilon[i-1,j]+epsilon[i,j-1]+epsilon[i-1,j-1]
                    a1=0.5*(epsilon[i,j]+epsilon[i,j-1])
                    a2=0.5*(epsilon[i-1,j]+epsilon[i,j])
                    a3=0.5*(epsilon[i-1,j-1]+epsilon[i-1,j])
                    a4=0.5*(epsilon[i,j-1]+epsilon[i-1,j-1])

                    #Calculate the residual
                    R[i,j]=1/a0*(a1*V_temp[i+1,j] + a2*V_temp[i,j+1] +
                        a3*V_temp[i-1,j] + a4*V_temp[i,j-1] +rho[i,j])-V_temp[i,j]

                    # Compute V using the overrelaxation method

                    V[i,j]=V_temp[i,j]+w*R[i,j]
                    error += abs(V[i,j]-V_temp[i,j])

                    V_temp[i,j]=V[i,j] # Make it already work for the next iteration 

            iterations += 1
            error /= float(len(y))
        print (" iterations ="+str(iterations))
        
        
     # BCs of 1 is a Neumann BC
        
    if BCs==1:  
        while iterations < max_iter and error > eps:

            V_temp = np.copy(V)
            error = 0 # we make this accumulate in the loop
            for j in range (0,len(x)-1):
                for i in range (0,len(y)-1):
                    
                    # Calculate the permittivity matrices.
                    ### Add he Neumann BCs
                    if j == 0 or j == len(x):
                        if i == 0 or i == len(y):
                            epsilon[i,j]=epsilon[i+1,j]
                            epsilon[i-1,j]=epsilon[i,j]
                        epsilon[i,j-1]=epsilon[i,j]
                        epsilon[i-1,j-1]=epsilon[i-1,j]

                    a0=epsilon[i,j]+epsilon[i-1,j]+epsilon[i,j-1]+epsilon[i-1,j-1]
                    a1=0.5*(epsilon[i,j]+epsilon[i,j-1])
                    a2=0.5*(epsilon[i-1,j]+epsilon[i,j])
                    a3=0.5*(epsilon[i-1,j-1]+epsilon[i-1,j])
                    a4=0.5*(epsilon[i,j-1]+epsilon[i-1,j-1])

                    #Calculate the residual
                    ### Add he Neumann BCs
                    if j == 0:
                        V_temp[i,j-1]=V_temp[i,j]
                    if i==0:
                        V_temp[i-1,j]=V_temp[i,j]
                    if j == Nx-2:
                        V_temp[i,j+1]=V_temp[i,j]
                    if i==Ny-2:
                        V_temp[i+1,j]=V_temp[i,j]
                            
                    R[i,j]=1/a0*(a1*V_temp[i+1,j] + a2*V_temp[i,j+1] +
                        a3*V_temp[i-1,j] + a4*V_temp[i,j-1] +rho[i,j])-V_temp[i,j]

                    # Compute V using the overrelaxation method

                    V[i,j]=V_temp[i,j]+w*R[i,j]
                    error += abs(V[i,j]-V_temp[i,j])

                    V_temp[i,j]=V[i,j] # Make it already work for the next iteration               
                    
       
            iterations += 1
            error /= float(len(y))
        print (" iterations ="+str(iterations))
        
    ### Calculate the electric field
    ### it will always have one less component than the potentials, 
    ### just by the way it is defined
    Ex=np.zeros((len(y)-1,len(x)-1))
    Ey=np.zeros_like(Ex)
    Ey=np.zeros_like(Ex)
    Ex_c=np.zeros_like(Ex)
    Ey_c=np.zeros_like(Ex)
    E_total=np.zeros_like(Ex)

    ### This shifts the x and y with respect to each other
    
    for j in range (0,len(x) -1):
            for i in range (0,len(y) -1):

                Ex[i,j]=-(V[i,j+1]-V[i,j])/ds
                Ey[i,j]=-(V[i+1,j]-V[i,j])/ds

                # Now we need to fix these shifts
    for j in range (0,len(x) -2):
            for i in range (0,len(y)-2):

                Ex_c[i,j]=0.5*(Ex[i+1,j]+Ex[i,j])
                Ey_c[i,j]=0.5*(Ey[i,j+1]+Ey[i,j])

                E_total[i,j]=np.sqrt(Ex_c[i,j]**2+Ey_c[i,j]**2)
    '''
    plt.figure(figsize=(12,7))
    plt.title('Voltage potentials')
    plt.pcolormesh(X,Y,V,cmap='jet')
    plt.colorbar()
    plt.show()

    '''
    title_name='Voltage potentials '+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    basic_plotting.plot_map(X_R, Y_R, V, 12, 7, title_name, 
             x_label,y_label, z_label='', show=True, minmax=False,
             zmin=0, zmax=10, savename=save_name, save=save,hd=hd,
             cmap=orig_cmap, font=12, midpoint=0.5, start=0, stop=1)
    '''
    plt.figure(figsize=(12,7))
    plt.pcolormesh(X_R,Y_R,V,cmap='jet')
    plt.colorbar()
    plt.title('Voltage potentials'+plot_names,fontsize=12)
    savename='Voltage potentials'+plot_names
    
    plt.show()
    
    '''
    # Plot the total magnitude of the electric field
    
    title_name='E field '+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    basic_plotting.plot_map(X_ER, Y_ER, E_total, 12, 7, title_name, 
             x_label,y_label, z_label='', show=True, minmax=False,
             zmin=0, zmax=10, savename=save_name, save=save,hd=hd,
             cmap=orig_cmap, font=12, midpoint=0.5, start=0, stop=1)
    
    
    title_name='E field with lines'+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    
    plt.figure(figsize=(12,7))
    plt.pcolormesh(X_ER,Y_ER,E_total,cmap='jet')
    plt.colorbar()
    plt.title('E field '+plot_names,fontsize=12)
    #plt.quiver(X_ER,Y_ER,Ex_c,Ey_c)
    if save:
        if hd==0:
            plt.savefig(save_name+'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300) 
        elif hd==1:
            plt.savefig(save_name +'_'+time.strftime("%H%M%d%m%y")+'.png', dpi=600)
        elif hd==2:
            plt.savefig(save_name +'_'+time.strftime("%H%M%d%m%y")+'.pdf')   

       
        
    plt.show()
    
    
    # plot inverse of electric field, which is related to capacitance
    
    title_name='Capacitance '+plot_names
    x_label=''
    y_label=''
    save_name=savename+title_name
    basic_plotting.plot_map(X_ER, Y_ER, 1/E_total, 12, 7, title_name, 
             x_label,y_label, z_label='', show=True, minmax=False,
             zmin=0, zmax=10, savename=save_name, save=save,hd=hd,
             cmap=orig_cmap, font=12, midpoint=0.5, start=0, stop=1)
    '''
    plt.figure(figsize=(12,7))
    plt.title('Capacitance model $C\propto 1/E$'+plot_names,fontsize=12)
    plt.pcolormesh(X_ER,Y_ER,1/E_total,cmap='jet',vmax=5)
    plt.colorbar()
    savename='Capacitance '+plot_names
    if save:
        if hd==0:
            plt.savefig(savename+'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300) 
        elif hd==1:
            plt.savefig(savename +'_'+time.strftime("%H%M%d%m%y")+'.png', dpi=600)
        elif hd==2:
            plt.savefig(savename +'_'+time.strftime("%H%M%d%m%y")+'.pdf')   

        np.savetxt(savename+time.strftime("%H%M%d%m%y")+'_metadata.txt',saveData, delimiter='\t', header = saveheader, comments='')
        
    plt.show()
    '''
    
    plt.figure(figsize=(12,7))
    plt.plot(X_ER[int(len(y)/2)],1/E_total[int(len(y)/2)],'bo-',label='C at y position '+str(Y_ER[int(len(y)/2),0]))
    plt.plot(X_ER[int(len(y)/2)],1/E_total[int(len(y)/2)],'ro-',label='C at y position '+str(Y_ER[int(len(y)/2)-1,0]))
    plt.plot(X_ER[int(len(y)/2)],E_total[int(len(y)/2)],'ko-',label='E at y position '+str(Y_ER[int(len(y)/2),0]))
    plt.plot(X_ER[int(len(y)/2)],E_total[int(len(y)/2)],'go-',label='E at y position '+str(Y_ER[int(len(y)/2)-1,0]))
    plt.title('Linecut of C and E at the 0 position '+plot_names,fontsize=12)
    plt.legend()
    
    save_name=savename+'Linecut of C and E at the 0 position '+plot_names
    if save:
        plt.savefig(save_name+'_'+time.strftime("%H%M%d%m%y")+'.png',dpi=300)
    plt.show()
    
    
    ### Make an object to store all the data
    data_cap=Object()
    

    setattr(data_cap,'charge_top',rho_top)
    setattr(data_cap,'charge_bottom',rho_bot)
    setattr(data_cap,'charge_middle_top',rho_midtop)
    setattr(data_cap,'charge_middle_bot',rho_midbot)
    
    setattr(data_cap,'hbn_top',top_plate)
    setattr(data_cap,'hbn_bottom',bot_plate)
    setattr(data_cap,'dielectric_1',epsr1)
    setattr(data_cap,'dielectric_2',epsr2)
    setattr(data_cap,'plate_size', 2*plate_start)
    setattr(data_cap,'cut_size', 2*cut_x)
       
    
    setattr(data_cap,'x',X)
    setattr(data_cap,'y',Y)
    setattr(data_cap,'x_E',X_E)
    setattr(data_cap,'y_E',Y_E)
    setattr(data_cap,'x_real',X_R)
    setattr(data_cap,'y_real',Y_R)
    setattr(data_cap,'x_Ereal',X_ER)
    setattr(data_cap,'y_Ereal',Y_ER)
    
    setattr(data_cap,'V',V)
    setattr(data_cap,'E_tot',E_total)
    setattr(data_cap,'Ex',Ex)
    setattr(data_cap,'Ey',Ex)
    
    setattr(data_cap,'convergence_requirement',eps)
    
    if save:
        save_obj_filename=savename+plot_names+'dataobj.pkl'
        with open (save_obj_filename, 'wb') as output: # Overwrite any existing file
            pickle.dump(data_cap,output,-1)
    
    return data_cap
    

