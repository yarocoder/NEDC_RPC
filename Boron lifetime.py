#Author: Zirui Jiang
#Contributers: Tyler Bailey, Yaro Kaminskiy

#README
#This is a simple script that takes into account the high neutron capture economy of Cd-113
#We are assuming that the loss of Cd-112 is pretty much small negiligble 
#Main assumptions: 
#we are only assuming natural Cd
#we are only taking account of Cd-112 and Cd-113
#Cd-112 has an abundance of 24.11 at%
#Cd-113 has an abundance of 12.23 at%
#For now we want the metal to have > 10.00 at% Cd-113; This is arbitrary

#Updated the code on calculation for B-10
#B-10 natural abundance: 19.8%
#Assuming natural B-10
#Thermal energy: 2.6 meV=0.026 eV = 2.6*10**-6 MeV

import math
import numpy as np
import matplotlib.pyplot as plt

#The majority of this data was extracted from https://www-nds.iaea.org/exfor/endf.htm 
#Units are CGS
density = 2.08
moles_per_volume = 2.08/10.81
atoms_per_volume = moles_per_volume*6.02*10**23
B_10 = atoms_per_volume*.198 #intital amount of B-10
sigma_10 = 4.5 * 10**3*10**-24
flux = 10**9/(4*math.pi*1)      # 1 cm away from a 10^9 [neutron/s] thermal neutron source

#This only computes and prints our start amount of B-10 and our end amount of B-10 
print('starting amount of B 10')
print(B_10)
print('10% of B 10')
print(atoms_per_volume*.1)

def main(t):
    '''
    Main function
    This computes the the neutron+target reaction rate for both isotopes
    negative rate is a loss while postive rate is a gain
    This is discretized with time intervals of .1 gigaseconds
    x is the amount of B-10 atoms used up after an interval
    y is the amount of B-10 atoms that remain after an interval
    this only returns z and y
    
    '''
    #rate_10 = sigma_10*B_10*flux
    #net_rate = 0-rate_10
    #x = net_rate * t
    y = B_10*math.exp(-sigma_10*flux*t)
    return y

#graph preparation. Time from t=0s to t=.5*10^11s with time interval 10**8s
B_10_lst = []
time = []
time.append(0)
B_10_lst.append(100)
i = 10**8
while i <= 10**14:
    time.append(i)
    i+=10**8

#more graph preparation
i = 10**8
while i <= 10**14:
    B_10_lst.append(main(i)*100/B_10)
    i+=10**8

#graph. The intersection is an estimated lifetime
horizontal = []
horizontal.append(10)

'''
i = 10**8
while i <= 10**14:
    horizontal.append(10**1)
    i+=10*8
'''

plt.semilogx(time, B_10_lst)
plt.ylim([0,100])
plt.xlabel('Time elapsed (s)', fontsize=18)
plt.ylabel('Percentage of B-10 Isotope Left (%)', fontsize=18)
plt.title('Boron Lifetime Calculation', fontsize=18)

plt.show()
