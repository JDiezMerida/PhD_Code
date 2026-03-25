# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 11:33:41 2022

Code with homemade functions which can be used in general for plotting faster. The code includes
several functions which I use frequently, such as a function to make the line plots look nice, or to make the map plots look nice.

@author: Jaime Diez Mérida @LS-Efetov group at LMU
"""

import numpy as np
import matplotlib.pyplot as plt
import time
#from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib

font=16
orig_cmap = matplotlib.cm.RdBu_r

def runningMeanFast(x, N):
    '''
    Code to do an smoothing with a running average
    x : the data which you want to smooth 
    N : the number of values which you want to take as window to smooth
    '''
    return np.convolve(x, np.ones((N,))/N)[(N-1):]

def plot_pre(sizex=8,sizey=5):
    '''
    Make an empty figure and return ax, to add extra parameters 
    ax can be used to set the plot as one likes
    The default size is 8 x 5
    The way to use it combined with plot_line would be:
        ax=plot pre() # this creates the figure
        ax.plot(datax,datay,label=legend) #this plots the actual data
        plot_line(ax) #this makes the plot look nice
        
    '''
    fig=plt.figure(figsize=(sizex,sizey))
    ax=fig.add_subplot(111)
    return ax

def plot_premap(sizex=8,sizey=5):
    '''
    Make an empty figure and return ax, and fig, to add extra parameters 
    ax can be used to set the plot as one likes
    fig is needed for the colorbar of a map, a linecut does not need it
    The default size is 8 x 5
    The way to use it combined with plot_map would be:
        ax=plot pre() # this creates the figure
        im=ax.pcolormesh(datax,datay,dataz,cmap=mycolor) #this plots the actual data
        plot_map(ax,im,fig) #this makes the plot look nice
    '''
    fig=plt.figure(figsize=(sizex,sizey))
    ax=fig.add_subplot(111)
    return ax,fig

def plot_map(ax,im,fig,xlabel='$B$ (mT)',ylabel='$I$ (nA)',font_size=font):
    '''
    Adds the most common parameters to make the map plot look nice
    By using ax. extra parameters, everything here can be overwritten
    The plot is given by ax.
    Returns the colorbar, so that it can be modified 
    '''
    ax.set_ylabel(ylabel,fontsize=font_size)
    ax.set_xlabel(xlabel,fontsize=font_size)  
    ### Add the colorbar 
    cbaxes = fig.add_axes([0.75, 0.9, 0.15,0.015]) #left,bottom,width,heigth
    cb = plt.colorbar(im,shrink=0.3,cax=cbaxes,orientation='horizontal') 
    #cb.set_label('d$V/$d$I$ (k$\Omega$)',fontsize=font,position=(-0.4,1), labelpad=-15,rotation=0)
    cbaxes.xaxis.tick_top()
    cb.ax.tick_params(labelsize=font_size,direction='in')
    #cb.set_ticks([0,25])
    ax.tick_params(axis='both',which='major',labelsize=font_size,direction='in',pad=6,size=8,width=1.5)
    ax.tick_params(axis='both',which='minor',direction='in',size=4,width=1.5)
    #ax1.tick_params(axis='x',which='both',bottom=False,labelbottom=False)
    #ax1.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.1))
    #ax1.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.5))
    #ax1.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.5))
    #ax1.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.1))
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    return cb 

def plot_line(ax,ylabel='$I_c$ (nA)',xlabel='$B$ (mT)',leg=False,legx=0.4,legy=1,
              font_size=font,font_leg=2):
    '''
    Adds the most common parameters to make the line plot look nice
    By using ax. extra parameters, everything here can be overwritten
    The plot is given by ax. 
    '''
    ax.set_ylabel(ylabel,fontsize=font_size)
    ax.set_xlabel(xlabel,fontsize=font_size)
    if leg:
        ax.legend(loc='upper center', bbox_to_anchor=(legx, legy), shadow=False, ncol=1,fontsize=font_size-font_leg,frameon=False)
    ax.tick_params(axis='both',which='major',labelsize=font_size,direction='in',pad=6,size=8,width=1.5)
    ax.tick_params(axis='both',which='minor',direction='in',size=4,width=1.5)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
        
def saving_figure(savecount,name,savename,dpi_depth=300):
    '''
    Saves the figures both in png and pdf
    Savecount is used as a counter, to only save the images that you really 
    want. When you save the image it will not save again unless you define
    savecount=0, as the number would now be 1
    The way to use it would be:
        savecount=0
        savecount=saving_figure(savecount,name,savename) 
    
    '''
    if savecount==0:
        plt.savefig(savename+name+time.strftime("%H%M%d%m%y")+'.png',dpi=dpi_depth,bbox_inches='tight') 
        plt.savefig(savename+name+time.strftime("%H%M%d%m%y")+'.pdf',bbox_inches='tight')       
        savecount=savecount+1
    return savecount

def shiftedColorMap(midpoint=0.5,cmap=orig_cmap, start=0,  stop=1.0,
                    division=258,name='shiftedcmap'):
    '''
    I am not sure why, but this does not always works. Maybe depends on the 
    python version (?)
    
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


        
def saving_sourcedata(savedata,name,savename,header=''):
    '''
    Saves the sourcedata of a figure into a txt file 

    
    '''
