# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 11:26:03 2021

@author: Jaime Díez Mérida @ LDQM group at ICFO and LMU

This code was used to extract the critical current of the Fraunhofer patterns using different criteria, 
such as the minimum value of the linecut, or a defined threshold. 
It can be used to extract the critical current from any map, as long as the data is given in the right format. 
The code also includes several parameters to define the limits of the plots, the font size, the color of the plots, etc.

"""
import numpy as np
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import time

savename='' # write the folder where you want to save the images 

'''
# x, y and z are the variables of the map 
# start is the point at which you start the extraction
# end is the point at which in ends
# middle is the middle point. So the zero current point. But if the variable, correct_middle is True
### then middle is not needed since it will find the real middle
# z_max, x_lim and y_lim are just the boundaaries of the plots +
# use min will use the minimum value of the linecut to start the threshold from there. Otherwise the threshold
# is just 0 + limit 
# limit2 was in case Ic- and Ic+ were very vdifferent. But I always use them as the same. 
# limit3 was an idea but is not used 
# normalize will plot divided by the maximum Ic 
# markers is the size of the markers 
# interp will interpolate the dvdi data so that it has many points to take and it will extract
# the exact value where the therehold is defined. Otherwise it will extract the closest
# data point to that value, but that will be wrong if you look at the actual values
'''

'''
Typical way of using it would be: 
Ic_pos,Ic_neg,D_Ic= get_Ic_minimum(x_data,y_data,z_data,limit=0.1,start=70,end=1000,
              middle=middle,z_max=5,x_lim=30,y_lim=100,save=False,font=18,window=7,use_min=True,limit2=0.1,
              correct_middle=False,normalize=False)

 It will extract the Ic and plot them, both in the map and also as plots. It will return 
 the Ic+, Ic- and the difference between them DIc
    
'''

def get_Ic_minimum(x,y,z,limit,start,end,middle,z_max=5,x_lim=50,y_lim=100,save=False,font=15,window=7,
                   use_min=False,limit2=0,limit3=0,plotting=True,correct_middle=False,normalize=False,
                   smooth=False,mycolor='RdBu_r',markers=5,interp=True,save_format='png'):
    
    #Extract it from a minimum value
    Ic_pos=np.zeros_like(x[:,0])
    Ic_neg=np.zeros_like(Ic_pos)
    middle_large=np.zeros_like(Ic_pos)
    
    if interp==True:
        new_x,new_y,new_z=interp_dvdi_data(x,y,z)
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
        if save:
            plt.savefig(savename+'JJBF_SCS_-Vbg1p12_Bdn_I7-18_Rxx1517'+time.strftime("%H%M%d%m%y")+'.'+save_format,dpi=300) 

        plt.show()


        plt.figure(figsize=(8,5))
        plt.plot(y[:,0],Ic_pos,'go--',label='$I_c^+(B)$',markersize=markers)
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
        plt.plot(y[:,0],Ic_pos_sm,'bo--',label='$I_c^+(B)$',markersize=markers)
        plt.plot(-y[:,0],-Ic_neg_sm,'ro--',label='$|I_c^-(-B)|$',markersize=markers)
        plt.title('Invert one of them',fontsize=font)
        plt.legend(loc='upper center', bbox_to_anchor=(0.2, 1), shadow=False, ncol=1,fontsize=font,frameon=False)
        plt.xlim(-x_lim,x_lim)
        plt.xlabel('$B$ (mT)',fontsize=font)
        plt.ylabel('$I_c$ (nA)',fontsize=font)
        plt.tick_params(axis='both',which='major',labelsize=font,direction='in')
        plt.show()
        
        
    return Ic_pos,Ic_neg,D_Ic


def get_Ic_minimum_xy1D(x,y,z,limit,start,end,middle,z_max=5,x_lim=50,y_lim=100,save=False,font=15,window=7,
                   use_min=False,limit2=0,limit3=0,plotting=True,correct_middle=False,normalize=False,
                   smooth=False,mycolor='RdBu_r',markers=5,interp=True,save_format='png'):
    
    #Extract it from a minimum value
    Ic_pos=np.zeros_like(x)
    Ic_neg=np.zeros_like(Ic_pos)
    middle_large=np.zeros_like(Ic_pos)
    
    if interp==True:
        new_x,new_y,new_z=interp_dvdi_data(x,y,z)
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
                    Ic_pos[count]=x[j]
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
        if save:
            plt.savefig(savename+'JJBF_SCS_-Vbg1p12_Bdn_I7-18_Rxx1517'+time.strftime("%H%M%d%m%y")+'.'+save_format,dpi=300) 

        plt.show()


        plt.figure(figsize=(8,5))
        plt.plot(y[:,0],Ic_pos,'go--',label='$I_c^+(B)$',markersize=markers)
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
        plt.plot(y[:,0],Ic_pos_sm,'bo--',label='$I_c^+(B)$',markersize=markers)
        plt.plot(-y[:,0],-Ic_neg_sm,'ro--',label='$|I_c^-(-B)|$',markersize=markers)
        plt.title('Invert one of them',fontsize=font)
        plt.legend(loc='upper center', bbox_to_anchor=(0.2, 1), shadow=False, ncol=1,fontsize=font,frameon=False)
        plt.xlim(-x_lim,x_lim)
        plt.xlabel('$B$ (mT)',fontsize=font)
        plt.ylabel('$I_c$ (nA)',fontsize=font)
        plt.tick_params(axis='both',which='major',labelsize=font,direction='in')
        plt.show()
        
        
    return Ic_pos,Ic_neg,D_Ic

def interp_dvdi_data(x,y,z,interpx_size=1001, interpy_size=601,interp_y=False):
    # Interpolate the data to get many more points 
    ### The interpolation has to go in ascending x order 
    y_interp=[]
    x_interp=[]
    new_x=np.linspace(x[0,0],x[0,-1],interpx_size)
    new_data=np.zeros((len(z),len(new_x))) 
    new_y_large=np.zeros((len(z),len(new_x)))
    new_x_large=np.zeros((len(z),len(new_x)))
    for i in range(0,len(z)):
      new_data[i]=np.interp(new_x,x[i],z[i])
      new_y_large[i]=y[i,0]
    for i in range(len(new_data)):
        new_x_large[i]=new_x  
        
    if interp_y: # this is usually set to False because it does not add anything
        ## Now interpolate the field 
        new_y=np.linspace(y[-1,0],y[0,0],interpy_size)
        
        new_dataB=np.zeros((len(new_y),len(new_x))) 
        new_y_largeB=np.zeros((len(new_y),len(new_x)))
        new_x_largeB=np.zeros((len(new_y),len(new_x)))
        for i in range(0,len(new_data[0])):
            new_dataB[:,i]=np.interp(new_y,np.flip(new_y_large[:,i]),np.flip(new_data[:,i]))
        
        for i in range(0,len(new_data[0])):
            #new_dataB[:,i]=np.flip(new_dataB[:,i])
            new_y_largeB[:,i]=new_y
        
        for i in range(len(new_dataB)):
            new_x_largeB[i]=new_x
        
        new_x_large=new_x_largeB
        new_y=new_y_largeB
        new_data=new_dataB
        
        
    return new_x_large,new_y_large,new_data
