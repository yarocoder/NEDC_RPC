# Author: Yaro Kaminskiy

# This code calculates the gamma flux density coming out of a sheet of boron reacting with incident neutrons to produce Li-7 in
# an excited state. That excited state has a 42 femtosecond half-life, so the decay can be assumed instant after the formation of Li-7.

# Import the number pi and the functions exp and plot from their respective libraries.
from math import pi, exp
import matplotlib.pyplot as plt

# Import the numpy library.
import numpy as np

# Definitions                # Units, Brief Descriptor. All units are in the CGS (centimeters-gram-second) unit system.

# Info given by Dr. Lee Bernstein.
I = 10e9                     # [# neutrons from source per second]=[no 'real' units], approx. top neutron current (HFNG)
r = 0.005                    # [cm], assumed distance from target boron slab given general size of HFNG assembly


# Obtained using WolframAlpha.com
Na = 6.022141e23             # [# of B-10 nuclei/mol B-10]=[1/mol], Avogadro's number
rho = 2.46                   # [mass natural boron/volume it occupies]=[g/cm^3], density of natural boron
x = 0.198                    # [mass B-10/mass natural boron]=[no 'real' units], mass fraction of B-10 in natural boron
M = 10.012936992             # [mass B-10/mol B-10]=[g/mol]=[u], molar mass of B-10

# Obtained using ENDF website: http://www.nndc.bnl.gov/exfor/endf00.jsp
sigma = 3.60069e-21          # [cm^2], cross-section for B-10(n,a)Li-7 reaction at the thermal energy peak at approximately 0.0253 eV.

# Define a function to calculate the gammas current produced per unit volume of the boron slab (gammas produced/(cm^3*s)) as a
# function of time.
def gamma(t):
    '''
    Calculate the initial B-10 nuclei density (# of initial B-10 nucei/volume of material).
    '''
    
    # The formula used falls out of dimensional analysis ( [# B-10 nuclei/mol B-10] * [mass nat. boron/volume] * [mass B-10/mass nat. boron]
    # / [mass B-10/mol B-10] = [# B-10 nuclei/volume]). The equation below calls values from the definitions section of the code above.
    # Assume the slab is entirely natural, elemental boron.
    rho_knot = Na*rho*x / M
    
    '''
    Calculate the neutron flux at the target using the assumed neutron current (neutrons per second NOT neutrons per second per area)
    from the source and distance from the source. 
    '''
    
    # The formula used falls out of the physical reasoning that flux is the current per unit area, and the total area we are evaluating
    # for is the surface area of a sphere with a radius equal to the seperation from the source and target. The equation below calls
    # values from the definitions section of the code above.
    omega_knot = I / (4 * pi * r * r)
    
    '''
    Return the gammas produced per unit volume of the material at the given input time, or the gamma flux density.
    '''
    
    # The formula below is derived from the equation for reaction rate, R = rho * I * sigma, by recognizing that R = -dN/dt, in which N
    # is the number of nuclei undergoing the reaction in the sample and -dN/dt is the negative of the time derivative of the number of
    # nuclei (assuming 1 nucleus is consumed per nuclear reaction), and rho is N/volume. We divide both sides by volume and solve the
    # resulting differential equation. The result is the number of gammas produced per unit time per unit volume of the boron slab
    # (in other words the gamma flux density).
    return rho_knot * sigma * omega_knot * exp(-sigma * omega_knot * t)


# Initialize the array for gamma flux density for plotting purposes, as well as a time array.
gamma_array = []
time_array = range(1,int(1.5e16),int(1e13))

# Create an array with a constant gamma count of 10 gammas/s/cm^3 for the lower bound of activity/cm^3.
null_rad = []

# Fill the gamma flux density value array.
for i in range(1,int(1.5e16),int(1e13)):
    # Append each new gamma flux density value to the end of the array, gradually building it.
    gamma_array.append(gamma(i))
    
# Fill the lower bound array.    
for i in range(1,int(1.5e16),int(1e13)):
    # Append each new 'null' gamma flux density value to the end of the array, gradually building it. Note that the value is the same
    # for every entry in the array.
    null_rad.append(int(10))

# Print relevant info on th start of the experiment.
print("The gamma production density at the start of the experiment is", str(gamma(0)), "gammas/s/cm^3.")

# The particular time indicated is used because it is the time at which the gamma flux nearly reaches the value of 10/s/cm^3, considered
# the 'null' value of gamma flux because it is so low. 
print("The gamma production density in 415 million years is", str( gamma( int(1.31e16) ) ), "gammas/s/cm^3.")

'''
Plot the gamma production density as a function of time.
'''

# Create a new figure window for plotting the graph.
plt.figure()

# Plot the gamma flux density as a function of time on a log-log plot. The label will appear in the plot legend.
plt.loglog(time_array, gamma_array, 'b-', label='Gamma flux density (1/s/cm^3)')

# Plot the 'null' gamma flux density as a function of time on the same log-log plot. The label will appear in the plot legend.
plt.loglog(time_array, null_rad, 'r--', label='Lower bound activity')

# Add the title of the plot in a reasonably sized font.
plt.title('Gamma Production Current Density in Natural Boron Over Time', fontsize=40)

# Add the y-axis and x-axis labels in a reasonably sized font to the plot.
plt.ylabel('Gamma Production Current Density (1/s/cm^3)', fontsize=30)
plt.xlabel('Time(s)', fontsize=30)

# Add the legend to the plot with a convenient location specified so the both gamma flux density curves can be seen in their entirety.
plt.legend(loc=(0.2,0.2))

# Show the created plot.
plt.show()
