# Author: Yaro Kaminskiy

'''
Calculates the percentage of thermal neutrons absorbed with different thicknesses and compositions
of the composite material.
'''

# Import the relevant functionality from the correct libraries.
from math import pi, exp
import matplotlib.pyplot as plt
import numpy as np

'''
Definitions                         Units, Brief Descriptor. All units are in the CGS unit system.
'''

# Obtained using WolframAlpha.com, Wikipedia.org, item description (for steel and epoxy components)
Na = 6.022141e23                    # [# of B-10 nuclei/mol B-10]=[1/mol], Avogadro's number

iso_frac = 0.198                    # [mass B-10/mass boron]=[no 'real' units], mass fraction of B-10 in natural boron

boron_molar_mass = 10.811           # [mass boron/mol boron]=[g/mol]=[u], molar mass of natural boron
B_10_molar_mass = 10.012936992      # [mass B-10/mol B-10]=[g/mol]=[u], molar mass of B-10
boric_acid_molar_mass = 61.833      # [mass boric acid/mol boric acid]=[g/mol]=[u], molar mass of boric acid

boric_acid_density = 1.435          # [mass boric acid/volume boric acid]=[g/cm^3], density of boric acid
steel_density = 7.9                 # [mass steel/volume steel]=[g/cm^3], density of steel
resin_density = 1.1628              # [mass epoxy resin/volume epoxy resin]=[g/cm^3], density of epoxy resin
hardener_density = 0.922            # [mass epoxy hardener/volume epoxy hardener]=[g/cm^3], density of epoxy hardener

# Obtained using ENDF website: http://www.nndc.bnl.gov/exfor/endf00.jsp
sigma = 3.8468e-21                  # [cm^2], cross-section for B-10 nonelastic neutron absorption reaction at the thermal energy peak at approximately 0.0253 eV.

# Define the relative mass fractions of the composite components for varying compostions of the composite material.
# Order of mass fractions: [boric acid, steel, epoxy resin, epoxy hardener]
mass_frac_array_1 = [0.500, 0, 0.306, 0.194]
mass_frac_array_2 = [0.333, 0.333, 0.206, 0.128]
mass_frac_array_3 = [0.250, 0, 0.458, 0.292]
mass_frac_array_4 = [0.170, 0.330, 0.306, 0.194]
mass_frac_array_5 = [0.200, 0, 0.489, 0.311]
mass_frac_array_6 = [0.200, 0.300, 0.306, 0.194]
mass_frac_array_7 = [0.333, 0.167, 0.306, 0.194]

'''
Calculate the percentage of thermal neutrons attenuated in the composite material.
'''

def neut_atten_percent(x, boric_acid_frac, steel_frac, resin_frac, hardener_frac):


    # Calculate the density of the composite material.
    rho = 1 / ( boric_acid_frac/boric_acid_density + steel_frac/steel_density + resin_frac/resin_density + hardener_frac/hardener_density )

    # Calculate the number density of B-10 nuclei.
    N = Na * iso_frac * boric_acid_frac * rho * boron_molar_mass / B_10_molar_mass / boric_acid_molar_mass

    # Calculate the percentage of neutrons attenuated assuming a narrow beam.
    # Also assume only the B-10 nonelastic neutron absoption is a significant
    # contributer to reduction of the neutron fluence.
    return 100 * (1 - exp(-N * sigma * x))

# Set the final and incremental values for the composite thickness and neutron percentage arrays.
final_val = 1.4    # [cm]
increm_val = 0.05  # [cm]

# Initialize the thickness and neutron percentage arrays.
thickness_array = np.arange(0, final_val, increm_val)
comp_array_1 = []
comp_array_2 = []
comp_array_3 = []
comp_array_4 = []
comp_array_5 = []
comp_array_6 = []
comp_array_7 = []

# Create arrays for the various composite compositions for plotting purposes.
for i in np.arange(0, final_val, increm_val):


    comp_array_1.append( neut_atten_percent(i, mass_frac_array_1[0],
                                             mass_frac_array_1[1],
                                             mass_frac_array_1[2],
                                             mass_frac_array_1[3]) )
    

    comp_array_2.append( neut_atten_percent(i, mass_frac_array_2[0],
                                             mass_frac_array_2[1],
                                             mass_frac_array_2[2],
                                             mass_frac_array_2[3]) )


    comp_array_3.append( neut_atten_percent(i, mass_frac_array_3[0],
                                             mass_frac_array_3[1],
                                             mass_frac_array_3[2],
                                             mass_frac_array_3[3]) )


    comp_array_4.append( neut_atten_percent(i, mass_frac_array_4[0],
                                             mass_frac_array_4[1],
                                             mass_frac_array_4[2],
                                             mass_frac_array_4[3]) )


    comp_array_5.append( neut_atten_percent(i, mass_frac_array_5[0],
                                             mass_frac_array_5[1],
                                             mass_frac_array_5[2],
                                             mass_frac_array_5[3]) )


    comp_array_6.append( neut_atten_percent(i, mass_frac_array_6[0],
                                             mass_frac_array_6[1],
                                             mass_frac_array_6[2],
                                             mass_frac_array_6[3]) )


    comp_array_7.append( neut_atten_percent(i, mass_frac_array_7[0],
                                             mass_frac_array_7[1],
                                             mass_frac_array_7[2],
                                             mass_frac_array_7[3]) )
'''
Plot the percentages of neutron current remaining at a particular thickness
in the absorber for the four compositions of interest.
'''
plt.figure()

plt.plot(thickness_array, comp_array_1, 'g+', label = "50.0% B, 0% S, 30.6% R, and 19.4% H")
plt.plot(thickness_array, comp_array_2, 'r--', label = "33.3% B, 33.3% S, 20.6% R, and 12.8% H")
plt.plot(thickness_array, comp_array_7, 'k*', label = "33.3% B, 16.7% S, 30.6% R, 19.4% H")
plt.plot(thickness_array, comp_array_6, 'y', label = "20.0% B, 30.0% S, 30.6% R, 19.4% H")
plt.plot(thickness_array, comp_array_3, 'c^', label = "25.0% B, 0% S, 45.8% R, and 29.2% H")
plt.plot(thickness_array, comp_array_4, 'bx', label = "17.0% B, 33.0% S, 30.6% R, and 19.4% H")
plt.plot(thickness_array, comp_array_5, 'm.', label = "20.0% B, 0% S, 48.9% R, 31.1% H")

plt.title('Thermal Neutron Fluence Attenuation', fontsize=24)

plt.xlabel('Thickness of material [cm]', fontsize=24)
plt.ylabel('Neutron fluence attenuated [%]', fontsize=24)

plt.legend(loc = (0.25,0), fontsize=18)

plt.show()
