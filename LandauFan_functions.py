# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 10:58:55 2023

This includes the typical functions to analyze Landau fan data of magic angle graphene transport devices. 
It includes functions to find the zero field, fix the gate voltage axis, find the CNP and CI positions, plot the Landau levels and extract the twist angle.
@author: Jaime Diez @LDQM-LS-Efetov group
"""
import matplotlib.pyplot as plt
from homemade_functions import plot_pre
import numpy as np

def zero_field(Bz):
    zerofield=np.amin(abs(Bz[:,0]))
    print('The lowest field value is: '+str(zerofield)+' T')
    if Bz[0,0]<0 and Bz[0,0]<Bz[1,0]:
        ascending=True
    elif Bz[0,0]<0 and Bz[0,0]>Bz[1,0]:
        ascending=False
    elif Bz[0,0]>0 and Bz[0,0]>Bz[1,0]:
        ascending=False
    elif Bz[0,0]>0 and Bz[0,0]<Bz[1,0]:
        ascending=True
    print(ascending)
    if ascending==True: 
        for i in range(len(Bz)):
            #print(Bz[i+1,0])
            if Bz[i,0]>=zerofield:
                print(i)
                
                zerofieldposition=i
                break
    else:
        for i in range(len(Bz)):
            if Bz[i,0]<=zerofield:
                print(i)
                zerofieldposition=i
                break
    return zerofield,zerofieldposition 



def fix_zero(Bx,By,Bz,zerofieldposition, CNP_st=30, CNP_end=20):
    # Range to find CNP, it should be close to the center, careful if there is 1/4 feature
    
    CNPgaterange=Bx[0,(int(len(Bz[0])/2)-CNP_st):(int(len(Bz[0])/2)+CNP_end)]
    
    ax=plot_pre(4,4)
    ax.plot(CNPgaterange,By[zerofieldposition,(int(len(Bz[0])/2)-CNP_st):(int(len(Bz[0])/2)+CNP_end)])
    plt.show()
    
    CNP_value=np.amax(By[zerofieldposition][(int(len(Bz[0])/2)-CNP_st):(int(len(Bz[0])/2)+CNP_end)])

    for i in range((int(len(Bz[0])/2)-CNP_st),(int(len(Bz[0])/2)+CNP_end)):
        if By[zerofieldposition,i]==CNP_value:
            CNP=i

    print('CNP value is: '+str(CNP_value)+' V')

    print('CNP gate position is: '+str(CNP))
    
    # Value of the gate at the CNP 
    
    zerogate=Bx[0][CNP]

    print('CNP gate value (Real zero) is: '+str(zerogate)+' V' )

    ## First fixe the back gate, then plot

    Bx_fixed=np.zeros_like(Bx)
    for i in range(len(Bx)):
        Bx_fixed[i]=Bx[i]-zerogate
    
    ax=plot_pre(4,4)
    ax.plot(Bx_fixed[0],By[zerofieldposition])
    ax.set_yscale('log')
    ax.axvline(x=0,linestyle='dashed',linewidth=1,color='silver')
    plt.show()
    return Bx_fixed

def find_CI(Bx,By,Bz,zerofieldposition,p1q_true=True,p1qst=10, pq1end=40, # DOes the sample ahve a 1/4 position 
             p1h_true=True,p1hst=25,p1hend=80):
    
    ## Find the -1/4 position 

    if p1q_true==True:
        p1qgaterange=Bx[0,(int(len(Bz[0])/2)+p1qst):(int(len(Bz[0])/2)+pq1end)]   
        ax=plot_pre(4,4)
        ax.plot(p1qgaterange,By[zerofieldposition,(int(len(Bz[0])/2)+p1qst):(int(len(Bz[0])/2)+pq1end)])
        ax.set_title('+1/4 position, Gate vs R at zerofield')
        plt.show()

        p1q_value=np.amax(By[zerofieldposition][(int(len(Bz[0])/2)+p1qst):(int(len(Bz[0])/2)+pq1end)])

        for i in range((int(len(Bz[0])/2)+p1qst),(int(len(Bz[0])/2)+pq1end)):
            if By[zerofieldposition,i]==p1q_value:
                p1q=i
        p1q_gate=Bx[0][p1q]


        print('+1/4 value is: '+str(p1q_value)+' V')

        print('+1/4 gate position is: '+str(p1q))

        print('+1/4 gate value (Real zero) is: '+str(p1q_gate)+' V' )
    else:
        print('This sample does not have 1/4 CI')
        p1q_gate=0
        
    # Find the +1/2
    if p1h_true==True:
        
        p1qgaterange=Bx[0,(int(len(Bz[0])/2)+p1hst):(int(len(Bz[0])/2)+p1hend)]   
        ax=plot_pre(4,4)
        ax.plot(p1qgaterange,By[zerofieldposition,(int(len(Bz[0])/2)+p1hst):(int(len(Bz[0])/2)+p1hend)])
        ax.set_title('+1/2 position, Gate vs R at zerofield')
        plt.show()

        p1h_value=np.amax(By[zerofieldposition][(int(len(Bz[0])/2)+p1hst):(int(len(Bz[0])/2)+p1hend)])

        for i in range((int(len(Bz[0])/2)+p1hst),(int(len(Bz[0])/2)+p1hend)):
            if By[zerofieldposition,i]==p1h_value:
                p1h=i

            #print(By67[zerofieldposition,i])

        print('p1h value is: '+str(p1h))

        p1h_gate=Bx[0][p1h]


        print('+1/2 value is: '+str(p1h_value)+' V')

        print('+1/2 gate position is: '+str(p1h))

        print('+1/2 gate value (Real zero) is: '+str(p1h_gate)+' V' )
    else:
        print('This sample does not have 1/2 CI')
        p1h_gate=0
    
    return p1q_gate,p1h_gate

def fix_moving_zero(Bx,By,Bz, CNP_st=30, CNP_end=20):
    # Range to find CNP, it should be close to the center, careful if there is 1/4 feature
    
    Bx_fixed=np.zeros_like(Bx)
    for j in range(len(Bz)):
        #ax=plot_pre(4,4)
        #ax.plot(CNPgaterange,By[j,(int(len(Bz[0])/2)-CNP_st):(int(len(Bz[0])/2)+CNP_end)])
        #ax.plot(CNPgaterange,By[-1,(int(len(Bz[0])/2)-CNP_st):(int(len(Bz[0])/2)+CNP_end)])
        #plt.show()
        CNP_value=np.amax(abs(By[j,(int(len(Bz[0])/2)-CNP_st):(int(len(Bz[0])/2)+CNP_end)]))
        for i in range((int(len(Bz[0])/2)-CNP_st),(int(len(Bz[0])/2)+CNP_end)):
            if By[j,i]==CNP_value:
                CNP=i
        zerogate=Bx[j,CNP]
        Bx_fixed[j]=Bx[j]-zerogate
    return Bx_fixed
def level_fit(x,slope,a):
    return a+slope*x

def Landau_plot(ax,levels,x_st,x_end,slope,degeneracy=4,CI_value=0):
    count=0
    B_start=[]
    for i in levels:
        if i==degeneracy:
            colors='tab:orange'
        else:
            colors='tab:red'
        if i%degeneracy==0:
            linestyles='-'
        else:
            linestyles='--'
        x=np.linspace(x_st[count],x_end[count],51)
        y=level_fit(x,slope/i,(-CI_value)*slope/i)
        ax.plot(x,y,color=colors,linestyle=linestyles,linewidth=1)   
        count=count+1
        B_start.append(np.amin(y))
    return B_start
        
def angle_extraction(BI_value,slope,):
    a=0.246E-9
    e=1.602*1E-19
    Cg_constant=3.87E-5
    
    Cg=Cg_constant*slope
    ns=Cg*BI_value/e
    theta=np.sqrt(np.sqrt(3)/8*abs(ns)*a**2)*180/np.pi
    print('C_g='+str(Cg))
    print('ns='+str(ns))
    print('Twist Angle={0:.4f}'.format(theta))
    return Cg, ns,theta

def span_magic(ax, ns=4,spansize=0.2):
    ## plus minus half filling
    ax.axvspan((ns/2)-spansize, (ns/2)+spansize, facecolor='salmon',alpha=0.1)
    ax.axvspan(-(ns/2)+spansize, -(ns/2)-spansize, facecolor='salmon',alpha=0.1)
    ## plus minus three quarters
    ax.axvspan(-(ns/4*3)+spansize, -(ns/4*3)-spansize, facecolor='b',alpha=0.1)
    ax.axvspan((ns/4*3)-spansize, (ns/4*3)+spansize, facecolor='b',alpha=0.1)
    ## plus one quarter
    ax.axvspan((ns/4)-spansize, (ns/4)+spansize, facecolor='skyblue',alpha=0.1)
    #ax.axvspan(-(ns/2)+0.1, -(ns/2)-0.1, facecolor='b',alpha=0.1)
    ## CNP
    ax.axvspan(-spansize, spansize, facecolor='silver',alpha=0.2)
    ## full_filling
    ax.axvspan((ns)-spansize, (ns)+spansize, facecolor='firebrick',alpha=0.1)
    ax.axvspan((-ns)-spansize, (-ns)+spansize, facecolor='firebrick',alpha=0.1)