# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 15:48:44 2023

@author: J.Diez
"""

"""Loading a Dataset"""
class Object():
    pass
# Make everything into an object 
general_directory='C:/jdiez_local/LDQM/a. Data/2023 April/'
folder=Object()
directory=Object()
files=Object()
dataobj=Object()
rawdata=Object()
variables=Object()

# Create the object where the data is stored  and plot all the datasets 
datatype='DT'
setattr(folder,datatype,'/TBGTB19/RxxTvsD/')
setattr(directory,datatype,general_directory+getattr(folder,datatype))

setattr(files,datatype,['J-TB19_RTvsD_40nA_I-4-16_Rxy1-8-3_Rxx2-9-11_Rxx3-8-9_VSiG-n10V_Vtg(0p535xMan+1p87(Vbg+Man))_mVbgT_T1p5.txt',
        'J-TB19_RTvsD_40nA_I-4-16_Rxy1-8-3_Rxx2-9-11_Rxx3-8-9_VSiG-n10V_Vtg(0p535xMan+1p87(Vbg+Man))_mVbgT_T2.txt',
         'J-TB19_RTvsD_40nA_I-4-16_Rxy1-8-3_Rxx2-9-11_Rxx3-8-9_VSiG-n10V_Vtg(0p535xMan+1p87(Vbg+Man))_mVbgT_T2p5.txt',     
         'J-TB19_RTvsD_40nA_I-4-16_Rxy1-8-3_Rxx2-9-11_Rxx3-8-9_VSiG-n10V_Vtg(0p535xMan+1p87(Vbg+Man))_mVbgT_T3.txt',
         'J-TB19_RTvsD_40nA_I-4-16_Rxy1-8-3_Rxx2-9-11_Rxx3-8-9_VSiG-n10V_Vtg(0p535xMan+1p87(Vbg+Man))_mVbgT_T3p5.txt',  
        'J-TB19_RTvsD_40nA_I-4-16_Rxy1-8-3_Rxx2-9-11_Rxx3-8-9_VSiG-n10V_Vtg(0p535xMan+1p87(Vbg+Man))_mVbgT_T4.txt',
         'J-TB19_RTvsD_40nA_I-4-16_Rxy1-8-3_Rxx2-9-11_Rxx3-8-9_VSiG-n10V_Vtg(0p535xMan+1p87(Vbg+Man))_mVbgT_T4p5.txt',  
          'J-TB19_RTvsD_40nA_I-4-16_Rxy1-8-3_Rxx2-9-11_Rxx3-8-9_VSiG-n10V_Vtg(0p535xMan+1p87(Vbg+Man))_mVbgT_T5.txt',               
                       ])

setattr(dataobj,datatype,{})  #stores all the data
setattr(rawdata,datatype,{}) #store the xaxis data for plotting
setattr(variables,datatype,{}) #store the data to plot in the right units

# All IV curves
# scan type = 0
for i in range(len(getattr(files,datatype))):
    file=getattr(directory,datatype)+getattr(files,datatype)[i]
    getattr(dataobj,datatype)['data'+str(i)], getattr(rawdata,datatype)['data'+str(i)],getattr(variables,datatype)['data'+str(i)]=LDQM_dataplot.plot_map_fromfiles(file)
    
    
### 
### TEMPERATURE FIGURES IN LINEPLOT 

ax=plot_pre()
num=5
for i in range(0,len(Rxy_temp)):
    ax.plot(n[str(num)],Rxy[num][i],color=colors[i])
ax.set_yscale('log')
plot_line(ax,xlabel='$n$',ylabel='$R_{xx}$ (k$\Omega$)')
ax.text(-2.88,100,'$D/\epsilon_0$={0:.2f} V/nm'.format(D[str(num)][100]),fontsize=font,color='k')
ax.set_xlim(-4.4,4.4)
ax.set_ylim(0.8,200)

nu_y=30
#ax.text(2.7/cap*1E16,nu_y,r'$\nu=4$',fontsize=font,rotation=-15)
ax.text(+2,nu_y,r'$\nu=3$',fontsize=font,rotation=-15)
ax.text(+1.15,nu_y,r'$\nu=2$',fontsize=font,rotation=-15)
ax.text(+0.35,nu_y,r'$\nu=1$',fontsize=font,rotation=-15)
ax.text(-1,nu_y,r'$\nu=0$',fontsize=font,rotation=-15)
ax.text(-3.4,nu_y,r'$\nu=-4$',fontsize=font,rotation=-15)
ax.text(-4.1,100,r'$hBN$',fontsize=font,rotation=-15)
ax.text(3.4,100,'$hBN$',fontsize=font,rotation=-15)

Vns=3.
VhBN=3.64
span_w=0.15
#ax.axvspan(-(Vns/2)+span_w, -(Vns/2)-span_w, facecolor='salmon',alpha=0.1)
#ax.axvspan(-(Vns/4*3)+span_w, -(Vns/4*3)-span_w, facecolor='b',alpha=0.1)
#ax.axvspan(-(ns/2)+0.1, -(ns/2)-0.1, facecolor='b',alpha=0.1)
ax.axvspan((Vns/2)-span_w, (Vns/2)+span_w, facecolor='salmon',alpha=0.1)
ax.axvspan((Vns/4)-span_w, (Vns/4)+span_w, facecolor='skyblue',alpha=0.1)
ax.axvspan((Vns/4*3)-span_w, (Vns/4*3)+span_w, facecolor='b',alpha=0.1)
ax.axvspan(-span_w, span_w, facecolor='silver',alpha=0.2)
ax.axvspan((Vns)-span_w, (Vns)+span_w, facecolor='firebrick',alpha=0.1)
ax.axvspan((-Vns)-span_w, (-Vns)+span_w, facecolor='firebrick',alpha=0.1)

ax.axvspan((VhBN)-span_w, (VhBN)+span_w, facecolor='royalblue',alpha=0.1)
ax.axvspan((-VhBN)-span_w, (-VhBN)+span_w, facecolor='royalblue',alpha=0.1)

plt.show()
