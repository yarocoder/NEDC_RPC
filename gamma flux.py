
# coding: utf-8

# In[64]:

# Author: Yaro Kaminskiy

# This code calculates the gamma current coming out of a sheet of boron reacting with incident neutrons to produce Li-7 in
# an excited state. That excited state has a 42 femtosecond half-life, so the decay can be assumed instant after the
# formation of Li-7.

# Import the number pi and the functions exp and plot from their respective libraries.
from math import pi, exp
import matplotlib.pyplot as plt

# Import the numpy library.
import numpy as np

# Definitions                # Units, Brief Descriptor

# Info given by Dr. Lee Bernstein.
I = 10e9                     # [# neutrons from source per second]=[no 'real' units], approx. top neutron current (HFNG)
r = 0.5                      # [m], assumed distance from target boron slab given general size of HFNG assembly


# Obtained using WolframAlpha.com
Na = 6.022141e23             # [# of B-10 nuclei/mol B-10]=[1/mol], Avogadro's number
rho = 2.46e6                 # [mass natural boron/volume it occupies]=[g/m^3], density of natural boron
x = 0.198                    # [mass B-10/mass natural boron]=[no 'real' units], mass fraction of B-10 in natural boron
M = 10.012936992             # [mass B-10/mol B-10]=[g/mol]=[u], molar mass of B-10

# Obtained using ENDF website: http://www.nndc.bnl.gov/exfor/endf00.jsp
sigma = 3600.69e-28          # [m^2], cross-section for thermal energy peak at approximately 0.0253 eV

# Define a function to calculate the gammas current produced per unit volume of the boron slab (gammas produced/(m^3*s)) as a
# function of time.

def gamma(t):
    # Calculate the initial B-10 nuclei density (# of initial B-10 nucei/volume of material).
    # Assume the slab is entirely natural, elemental boron.
    rho_knot = Na*rho*x / M
    
    # Calculate the neutron flux at the target using the assumed neutron current from the source and distance from the
    # source.
    omega_knot = I / (4 * pi * r * r)
    
    # Return the gammas produced per unit volume of the material at the given input time.
    return rho_knot * sigma * omega_knot * exp(-sigma * omega_knot * t)

print("The gamma production density at the start of the experiment is", str(gamma(0)), "gammas/s/cm^3.")
print("The gamma production density in 415 million years is", str(gamma(int(1.31e16))/1e6), "gammas/s/cm^3.")

# Initialize the array for gamma ray production density for plotting purposes, as well as a
# time array and blank line array.
gamma_array = []
time_array = range(1,int(1.5e16),int(1e13))

# Create an array with a constant gamma count of 10/s/cm^2 for the lower bound of activity/cm^3.
null_rad = []

# Fill the gamma production density value array.
for i in range(1,int(1.5e16),int(1e13)):
    gamma_array.append(gamma(i)/1e6)
    
# Fill the lower bound array.    
for i in range(1,int(1.5e16),int(1e13)):
    null_rad.append(int(10))
    
# Plot the gamma production density as a function of time.

plt.figure()
plt.loglog(time_array, gamma_array, 'b-', label='Gamma activity density (1/s/cm^3)')
plt.loglog(time_array, null_rad, 'r--', label='Lower bound activity')
plt.title('Gamma Production Current Density in Natural Boron Over Time', fontsize=40)
plt.ylabel('Gamma Production Current Density (1/s/cm^3)', fontsize=30)
plt.xlabel('Time(s)', fontsize=30)
plt.legend(loc=(0.2,0.2))

plt.show()


# In[ ]:



