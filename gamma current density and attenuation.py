# Author: Yaro Kaminskiy

# This code calculates the gamma fluence density coming out of a sheet of boron reacting with incident neutrons to produce Li-7 in
# an excited state. That excited state has a 42 femtosecond half-life, so the decay can be assumed instant after the formation of Li-7.

# Import the number pi and the functions exp and plot from their respective libraries.
from math import pi, exp
import matplotlib.pyplot as plt

# Import the numpy library.
import numpy as np

'''
Definitions            		 # Units, Brief Descriptor. All units are in the CGS (centimeters-gram-second) unit system.
'''

# Info given by Dr. Lee Bernstein.
I = 10e9                    # [# neutrons from source per second]=[no 'real' units], approx. top neutron current (HFNG)
r = 1                       # [cm], assumed distance from target boron slab given general size of HFNG assembly


# Obtained using WolframAlpha.com
Na = 6.022141e23             # [# of B-10 nuclei/mol B-10]=[1/mol], Avogadro's number
rho = 2.46                   # [mass natural boron/volume it occupies]=[g/cm^3], density of natural boron
x = 0.198                    # [mass B-10/mass natural boron]=[no 'real' units], mass fraction of B-10 in natural boron
M = 10.012936992             # [mass B-10/mol B-10]=[g/mol]=[u], molar mass of B-10

# Obtained using ENDF website: http://www.nndc.bnl.gov/exfor/endf00.jsp
sigma = 3.60069e-21          # [cm^2], cross-section for B-10(n,a)Li-7 reaction at the thermal energy peak at approximately 0.0253 eV.

# Obtained using WolframAlpha.com
ironDen = 7.874				 # [g/cm^3], density of natural iron
polyDen = 0.95				 # [g/cm^3], density of polyethylene
leadDen = 11.34				 # [g/cm^3], density of lead
ferricOxideDen = 5.26		 # [g/cm^3], density of iron (III) oxide (ferric oxide)
ferrousOxideDen = 5.7		 # [g/cm^3], density of iron (II) oxide (ferrous oxide)
limestoneDen = 2.93			 # [g/cm^3], density of limestone (calcium carbonate)
stainlessSteelDen = 7.9		 # [g/cm^3], density of stainless steel (assume mean value)

# Derived from molar and atomic masses by taking the atomic masses of each individual element, multiplying by the number of atoms of that
# element in a molecule of that compound, and dividing by the molar mass of the compound available on WolframAlpha.com
wtFracArrayFerricOxide = [0.6994176, 0.3005711]
wtFracArrayFerrousOxide = [0.7773092, 0.2226964]
wtFracArrayLimestone = [0.4004196, 0.119999, 0.4795504]

# Values for the weight fractions of steel assumed as 89.5% iron, 10.5% chromium, and an insignificant amount of carbon ( <1% ) based off of
# the definition of stainless steel and carbon content information available on https://en.wikipedia.org/wiki/Stainless_steel
wtFracArrayStainlessSteel = [0.895, 0.105]

# Obtained using the NIST X-Ray Mass Attenuation Coefficients database: https://www.nist.gov/pml/x-ray-mass-attenuation-coefficients
attenCoeffPoly = 9.947e-02		 # [cm^2/g], mass attenuation coefficient for polyethylene
attenCoeffLead = 1.614e-01		 # [cm^s/g], mass attenuation coefficient for natural lead
attenCoeffIron = 8.414e-02		 # [cm^2/g], mass attenuation coefficient for natural iron
attenCoeffOxygen = 8.729e-02	 # [cm^s/g], mass attenuation coefficient for natural oxygen
attenCoeffCalcium = 8.851e-02	 # [cm^s/g], mass attenuation coefficient for natural calcium
attenCoeffCarbon = 8.715e-02	 # [cm^s/g], mass attenuation coefficient for natural carbon
attenCoeffChromium = 8.281e-02	 # [cm^s/g], mass attenuation coefficient for natural chromium

# Define mass attenuation arrays for the compounds of interest, to be inputed into an overall mass attenuation coefficient function later in
# the script.
attenArrayFerricOxide = [attenCoeffIron, attenCoeffOxygen]
attenArrayFerrousOxide = [attenCoeffIron, attenCoeffOxygen]
attenArrayLimestone = [attenCoeffCalcium, attenCoeffCarbon, attenCoeffOxygen]
attenArrayStainlessSteel = [attenCoeffIron, attenCoeffChromium]

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
    Calculate the neutron fluence at the target using the assumed neutron current (neutrons per second NOT neutrons per second per area)
    from the source and distance from the source. 
    '''
    
    # The formula used falls out of the physical reasoning that fluence is the current per unit area, and the total area we are evaluating
    # for is the surface area of a sphere with a radius equal to the seperation from the source and target. The equation below calls
    # values from the definitions section of the code above.
    omega_knot = I / (4 * pi * r * r)
    
    '''
    Return the gammas produced per unit volume of the material at the given input time, or the gamma fluence density.
    '''
    
    # The formula below is derived from the equation for reaction rate, R = rho * I * sigma, by recognizing that R = -dN/dt, in which N
    # is the number of nuclei undergoing the reaction in the sample and -dN/dt is the negative of the time derivative of the number of
    # nuclei (assuming 1 nucleus is consumed per nuclear reaction), and rho is N/volume. We divide both sides by volume and solve the
    # resulting differential equation. The result is the number of gammas produced per unit time per unit volume of the boron slab
    # (in other words the gamma fluence density).
    return rho_knot * sigma * omega_knot * exp(-sigma * omega_knot * t)

# Define a function to calculate the secondary gamma ray fluence per unit volume of the boron slab attenuated in a gamma absorbing material slab
# placed to be just directly touching the surface of the boric acid.
def gammaAtten(attenCoeff, dens, x):

	# Return the value of the gamma fluence. Formula developed from one given from http://physics.nist.gov/PhysRefData/XrayMassCoef/chap2.html
	return gamma(0) * ( 1 - exp(-attenCoeff * dens * x) ) # x in cm

# Calculate the percentage of gamma ray fluence attenuated in a material as a function of the mass attenuation coefficient, depth in the
# material and the time since the start of the experiment.
def gammaAttenPercent(attenCoeff, dens, x, t):

	# Divide the gamma fluence attenuation by the original gamma fluence.
	return 100 * gammaAtten(attenCoeff, dens, x) / gamma(t)	# x in cm, t in seconds


# Iterate through the arrays of mass attenuation coefficient and weight fraction material property arrays to calculate the mass attenuation
# coefficient (for the secondary gamma rays emitted for the B-10 isotope) for compounds and mixtures at a certain depth in the material.
def attenCoeffOverall(attenArray, wtFracArray):

	# Initialize a value for the mass attenuation coefficient of a compound or mixture.
	massAttenCoeff = 0

	# Iterate through each material property array to calculate the mass attenuation coefficient for the compound or mixture of interest.
	for i in np.arange(0, len(attenArray)):

		# Add the weighted mass attenuation coefficient of an element or component in the compound or mixture to the previous overall
		# mass attenuation coefficient. Formula obtained from http://physics.nist.gov/PhysRefData/XrayMassCoef/chap2.html
		massAttenCoeff = massAttenCoeff + wtFracArray[i] * attenArray[i]

	# Return the overall mass attenuation coefficient.
	return massAttenCoeff

# Set the values to end at and the incremental value for the gamma fluence density, time, and lower bound arrays.
final_val = int(8e12)
increm_val = int(1e8)

# Initialize the array for gamma fluence density for plotting purposes, as well as a time array.
gamma_array = []
time_array = np.arange(1, final_val, increm_val)

# Create an array with a constant gamma count of 10 gammas/s/cm^3 for the lower bound of activity/cm^3.
null_rad = []

# Fill the gamma fluence density value and lower bound array.
for i in np.arange(1, final_val, increm_val):
    
    # Append each new gamma fluence density value to the end of the array, gradually building it.
    gamma_array.append(gamma(i))

    # Append each new 'null' gamma fluence density value to the end of the array, gradually building it. Note that the value is the same
    # for every entry in the array.
    null_rad.append(int(10))

# Create arrays for the gamma fluence percentages of the relevant materials to be analyzed.
poly_array = []
lead_array = []
ferric_oxide_array = []
ferrous_oxide_array = []
limestone_array = []
stainless_steel_array = []

# Create an array for the depth values in a slab of a particular material.
x_array = np.arange(0, 1.22, 0.06)

# Calculate the mass attenuation coefficients for the relevant compounds.
attenCoeffFerricOxide = attenCoeffOverall(attenArrayFerricOxide, wtFracArrayFerricOxide)
attenCoeffFerrousOxide = attenCoeffOverall(attenArrayFerrousOxide, wtFracArrayFerrousOxide)
attenCoeffLimestone = attenCoeffOverall(attenArrayLimestone, wtFracArrayLimestone)
attenCoeffStainlessSteel = attenCoeffOverall(attenArrayStainlessSteel, wtFracArrayStainlessSteel)

# Create arrays for the various materials for plotting purposes.
for i in np.arange(0, 1.22, 0.06):
	poly_array.append( gammaAttenPercent(attenCoeffPoly, polyDen, i, 0) )
	lead_array.append( gammaAttenPercent(attenCoeffLead, leadDen, i, 0) )
	ferric_oxide_array.append( gammaAttenPercent(attenCoeffFerricOxide, ferricOxideDen, i, 0) )
	ferrous_oxide_array.append( gammaAttenPercent(attenCoeffFerrousOxide, ferrousOxideDen, i, 0) )
	limestone_array.append( gammaAttenPercent(attenCoeffLimestone, limestoneDen, i, 0) )
	stainless_steel_array.append( gammaAttenPercent(attenCoeffStainlessSteel, stainlessSteelDen, i, 0) )


# Print relevant info on the start of the experiment.
print("The gamma activity density at the start of the experiment is", str(gamma(0)), "Bq/cm^3.")

# The particular time indicated is used because it is the time at which the gamma fluence nearly reaches the value of 10 Bq/cm^3, considered
# the 'null' value of gamma fluence because it is so low. 
print("The gamma activity density in 9.25 years is", str( gamma( int(2.9182e8) ) ), "Bq/cm^3.")

# Print the gamma fluence percentage of the original gamma fluence remaining at 5 mm in the gamma absorbing material at the start of the experiment.
print( 'The percentage of the original gamma fluence remaining at 5 mm in natural iron is', gammaAttenPercent(attenCoeffIron, ironDen, 0.5, 0) )
print( 'The percentage of the original gamma fluence remaining at 5 mm in polyethylene is', gammaAttenPercent(attenCoeffPoly, polyDen, 0.5, 0) )
print( 'The percentage of the original gamma fluence remaining at 5 mm in lead is', gammaAttenPercent(attenCoeffLead, leadDen, 0.5, 0) )

'''
Plot the gamma production density as a function of time.
'''

# Create a new figure window for plotting the graph.
plt.figure()

# Plot the gamma fluence density as a function of time on a log-log plot. The label will appear in the plot legend.
plt.loglog(time_array, gamma_array, 'b-', label='Gamma current density (Bq/cm^3)')

# Plot the 'null' gamma fluence density as a function of time on the same log-log plot. The label will appear in the plot legend.
plt.loglog(time_array, null_rad, 'r--', label='Lower bound activity')

# Add the title of the plot in a reasonably sized font.
plt.title('Gamma Current Density from B-10(n,a + g)L-7 Reaction')

# Add the y-axis and x-axis labels in a reasonably sized font to the plot.
plt.xlabel('Time(s)')
plt.ylabel('Gamma Current Density (Bq/cm^3)')

# Add the legend to the plot with a convenient location specified so the both gamma fluence density curves can be seen in their entirety.
plt.legend(loc=(0.2,0.2))

# Show the created plot.
plt.show()

'''
Plot the gamma fluence attenuation percentage for different materials at the start of the experiment as functions of the depth within the material.
'''

plt.figure()

plt.plot(x_array, lead_array, 'r--', label = 'Lead')
plt.plot(x_array, stainless_steel_array, 'b', label = 'Stainless steel')
plt.plot(x_array, ferrous_oxide_array, 'm:', label = 'Ferrous oxide')
plt.plot(x_array, ferric_oxide_array, 'c+', label = 'Ferric oxide')
plt.plot(x_array, limestone_array, 'k.', label = 'Limestone')
plt.plot(x_array, poly_array, 'gx', label = 'Polyethylene')

plt.title('478 keV Gamma Fluence Attenuation')

plt.xlabel('Depth in material [cm]')
plt.ylabel('Gamma fluence attenuation [%]')

plt.legend(loc = (0,0.6))

plt.show()
