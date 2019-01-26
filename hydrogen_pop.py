import matplotlib.pyplot as plt
import numpy as np
import math


num_levels = 2
max_levels_to_plot = 2
keV = 8.61734e-5      #boltzman constant eV/K
k = 1.38065e-23          #boltzman constant j/K
ion_energy = 13.6
m_e = 9.10938e-31
h = 6.62607e-34
H_density = 2.1e20



def boltzman_pop(temp,j):  #n_j/N ratio of electrons in level j to all particles in ionisation state
    g_j = 2.0 * float(j*j)
    norm = keV * float(temp)
    E_j = -ion_energy / float(j*j)
    relnum = ((g_j / particion_func(temp)) * np.exp(-(E_j - (-ion_energy)) / norm))
    return relnum

def delta_E(i,j):
    return ((- ion_energy / float(j*j)) - (- ion_energy / float(i*i)))


def particion_func(temp):
    Z = 0.0
    norm = keV * temp
    for i in range(1,num_levels+1):
        E_i = -ion_energy / float(i*i)
        g_i = 2.0 * float(i*i)
        Z += g_i * np.exp(-(E_i - (-ion_energy)) / norm)
    return Z


def saha_quad(temp,density):
    Z = particion_func(temp)

    const = (((2.0 * np.pi * m_e * k * temp) / (h*h))**(3.0/2.0)) * 2.0 * np.exp(-(ion_energy) / (keV *temp))
    const = const / (Z * density)

    relnum = (-const + math.sqrt(const**2 + (4*const))) * 0.5
    if(relnum < 0.0):
        print("ERROR: ratio negative")
        exit()
    return relnum





count = 0
t_min = 5000
t_max = 85000
t_step = 100

len=int((t_max-t_min)/t_step)
ints_levels = np.zeros((max_levels_to_plot,len))

temps = np.arange(t_min,t_max,t_step)
plt.figure(1)
plt.subplot(212)
for j in range(1,max_levels_to_plot+1):
    t_count = 0
    for t in range(t_min,t_max,t_step):
        ints_levels[count,t_count] = boltzman_pop(t,j) * (1.0-saha_quad(t,H_density))

        t_count += 1
    plt.plot(temps,ints_levels[count,:],label="state = {}".format(j))
    count += 1
plt.yscale("log")
plt.xscale("log")
plt.xlabel("Temperature [K]")
plt.ylabel("Fraction of H in state")
plt.axhline(1.0,color='k',linestyle='--',linewidth=1,alpha=0.5)
plt.legend()


len = int((t_max-t_min) / t_step)
data = np.zeros((len,1))
data1 = np.zeros((len,3))
count = 0
for t in range(t_min,t_max,t_step):
    data[count,0]=t
    data1[count,0]=boltzman_pop(t,1)
    data1[count,2]=boltzman_pop(t,2)
    data1[count,1]= saha_quad(t,H_density)

    count += 1


plt.subplot(222)
plt.plot(data[:,0],data1[:,0],'g',label="Ratio of state {} and all states".format(1))
plt.plot(data[:,0],data1[:,2],'b',label="Ratio of state {} and all states".format(2))
plt.axhline(1.0,color='k',linestyle='--',linewidth=1,alpha=0.5)
plt.yscale("log")
plt.ylabel("Fraction of H in state")
plt.xlabel("Temperature [K]")
plt.legend()
plt.subplot(221)
plt.plot(data[:,0],data1[:,1],'r',label="population of HII")
plt.axhline(1.0,color='k',linestyle='--',linewidth=1,alpha=0.5)
plt.ylabel("Fraction of H in state")
plt.legend()
plt.xlabel("Temperature [K]")

plt.show()






###DEBUGING### ###DEBUGING### ###DEBUGING### ###DEBUGING###

"""
t = 10000
print("Partition function for temperature {} is {:.5f}".format(t,particion_func(t)))
print("Relative pop between levels {} and {} is {:.5f}%".format(i,j,100*boltzman_pop_state(t,i,j)))
print("Ionised atoms to non ionised atoms is {}".format(saha_pop(t,20.0)))
print("Ionised atoms to non ionised atoms is {}".format(saha_quad(t,H_density)))
"""
"""
states = np.arange(1,num_levels+1,1)
fig=plt.figure()
plt.plot(states,ints_levels,'r.')
plt.vlines(states, [0],ints_levels)
plt.yscale("log")
plt.show()
"""
