#%%
#%%
#import scipy
import numpy as np 
import matplotlib.pyplot as plt

#3D JJ
B=np.arange(-50,50,0.001)
w=0.8E-6   # Width of the Hall bar 
L=100E-9   #Length of the JJ 
lam=0
A=(w-2*lam)*(L+2*lam)  # Area of the JJ 
phi=B*A
I_c=1

phi_0=2.067833848E-12
ws=250E-9   # width of the SQUID arms (in our case it might be the wideness of the edge currents )

As=1E-6*300E-9  # Area inside of the SQUID
#As=(w-ws)*L*3  # Area inside of the SQUID
#As= (0.6E-6)*(0.6E-6)
area_totalSQUID=4*ws**2+As
a0=0.43 
font=20
def Is(B,A):
    return I_c*np.abs(np.sin(np.pi*B*A/phi_0)/(np.pi*B*A/phi_0))

def I3DSQUID(B,A):
    return I_c*np.sin((np.pi*B*A/phi_0))
# 2D Jj correct with Bessel function

def IB(B,w):
    h=w**2*B/phi_0   # Reduced field. Eq. 9 of Moshe Paper 
    return abs(scipy.special.jv(0, 0.43*h))
# 2D JJ
def Im(B,w):
    #return 0.61*np.sqrt(phi_0/B)*np.abs(np.cos(1.72*(B*L^2/phi_0)-(np.pi/4)))
    #return 0.61*np.sqrt(phi_0/np.abs(B))*np.abs(np.cos(1.72*(B*w**2/phi_0)-(np.pi/4)))
    return I_c*np.abs(np.sin(np.pi*B*w**2/(1.8*phi_0))/(np.pi*B*w**2/(1.8*phi_0)))

# 2D SQUID Bessel 
def ISquid(B,ws,As):
    #phi=B*((w-2*ws)**2)
    phi=B*area_totalSQUID
    J=scipy.special.jv(0, 4*a0*(ws**2/As)*phi/phi_0)
    #J=scipy.special.jv(0, 4*a0*(ws**2/As)*B)
    #phi=B*(w-2*ws)**2/1.8
    return I_c*abs(J*np.cos(np.pi*phi/phi_0))
    #return Ic*abs(J*np.cos(np.pi*B))

plt.plot(B,Im(B,w),'k')
plt.plot(B,ISquid(B,ws,As))
# %%
