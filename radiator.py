# Spacecraft Radiator thermal simulation
# (c) John Nguyen 2025
# MIT License

from vpython import *
import numpy as np
import matplotlib.pyplot as plt

# UNIVERSAL CONSTANTS

# Stefan-Boltzmann constant
stefboltz = 5.67e-8 # Watt / meter^2 Kelvin^4

# TEMPERATURES

# space
T_inf = 0 # Kelvin

# fluid inlet
T_inlet = 320 # Kelvin

# initial surface temperature
T_surf = 200 # Kelvin

# THERMAL RESISTANCES

# conductive thermal resistance of the radiator panel
conduct_rt = 8.98e-6 # Kelvin meter^2 / Watt

# convective thermal resistance of working fluid inside panel
convect_rt = 0.0511 # Kelvin meter^2 / Watt

# MATERIAL PROPERTIES OF RADIATOR PANEL

# coefficient of thermal expansion of panel material
ctexpand = 2.36e-5

# emissivity of the radiator panel
emiss = 0.8

# specific heat capacity of panel
cp_panel = 896 # Joule / kilogram Kelvin

# DIMENSIONS OF RADIATOR PANEL

width = 0.58 # meter
length = 1.73 # meter
thickness = 0.0035 # meter

# area density of panel
dens_panel = 6.78 # kilogram / meter^2

# tube diameter
diameter = 0.0025 # meter

# tube spacing
spacing = 0.005 # meter

# number of tubes
n_tubes = 112

# path length
path_length = 201.81 # meter

# MATERIAL PROPERTIES OF WORKING FLUID

density = 696 # kilogram / meter^3

viscosity = 0.00013 # kilogram / meter second (alternately Pascal second)

cp_fluid = 4744 # Joule / kilogram Kelvin

velocity = 1.46 # meter / second

# SIMULATION OVERVIEW

# This program assumes a rectangular radiator panel
# with many coolant tubes running parallel in the 
# same direction through it, assumed to be served
# by inlet and outlet manifolds at either end.
# The properties of this panel, detailed above, are
# calculated in a spreadsheet also in this folder.

# This program aims to calculate the heat transfer
# from this radiator at steady state in darkness.

# For now, no heat will bleed between adjacent tubes.

# The thermal circuit is composed of 3 stages:
# 1. Convective transfer between working fluid
#    and the radiator tube
# 2. Conductive transfer between the tube wall
#    and the exterior radiator surface
# 3. Radiative transfer from the exterior surface
#    into space.

# Which are modeled as such:
# 1.    The Swamee-Jain approximation is used to find
#    the Darcy friction factor.  This is used in
#    conjunction with the Prandtl number and Reynolds
#    number to find the average Nusselt number, through
#    Gnielinski's correlation.  Then the heat transfer
#    coefficient is calculated from the Nusselt number.
#       This is all done in the spreadsheet.
# 2.    Fourier's Law of conduction is used assuming that
#    the average effective thickness of the panel is
#    approximately equal to the average of the tubing
#    wall thickness and the total wall thickness.
# 3.    Stefan-Boltzmann law of radiation is used,
#    assuming that the area is the entire half-spacing
#    on either side of the tube element, and that the 
#    background temperature is absolute zero.

# SIMULATION VARIABLES

# time
# Note that time step dt also directly affects number of elements!
# This simulation takes O(1/dt^2) time to run.
t = 0 # seconds
dt = 0.1 # seconds
tmax = 1 # seconds

# number of elements
# sized for 1 element representing length fluid flows in 0.01 seconds
n = ceil(path_length / velocity * 100)

# sim tick counter
m = 0
m_max = ceil(tmax/dt)

# ELEMENT PROPERTIES

# element length
len_ea = path_length/n # meter

# element area
area_fluid = len_ea * diameter * np.pi # meter^2, area of inside of tube
area_radiate = len_ea * spacing * 2 # meter^2, area for radiation and conduction (both sides included)

# Adjust the thermal resistances for the area of the element
conductive_rt = conduct_rt / area_radiate # Kelvin / Watt
convective_rt = convect_rt / area_fluid # Kelvin / Watt

# heat capacity of the fluid in element
cp_ea_fluid = cp_fluid * density * diameter**2 / 4 * np.pi * len_ea # Joule / Kelvin

# heat capacity of the solid panel in element
cp_ea_panel = cp_panel * dens_panel * len_ea * spacing # Joule / Kelvin

# ELEMENT PROPERTY MATRIX DEFINITION

# 0-n y z = Element position
# x 0 z = Fluid temperature
# x 1 z = Surface temperature
# x 2 z = Radiative heat transfer
# x y 0-m = Simulation tick

shape = (n, 3, m_max)

matrix = np.zeros(shape)

nn = 0
while nn < n:
    matrix[nn,0,0] = T_surf
    matrix[nn,2,0] = T_surf
    nn += 1
# SIMULATION LOOP
while m < ceil(tmax/dt):
    nn = 0
    while nn < n:
        # if at fluid inlet, clamp to fluid inlet temperature
        if nn == 0:
            fluid_T_inlet = T_inlet
            matrix[0,0,m] = T_inlet
        # otherwise, take the temperature propagated from previous position
        else:
            fluid_T_inlet = matrix[nn,0,m]

        # Take the current surface temperature (this is written to by the previous timestep)
        surf_T_initial = matrix[nn,1,m]

        # Thermal resistance of radiation, Kelvin / Watt
        radiative_rt = 1 / (stefboltz * emiss * area_radiate * (surf_T_initial**2 + T_inf**2) * (surf_T_initial + T_inf))

        print((radiative_rt+conductive_rt+convective_rt)*area_radiate)
        # Total heat transfer out of this element, Joule
        heat_transfer = 1/(radiative_rt + conductive_rt + convective_rt) * (surf_T_initial - T_inf) * dt

        # Fluid outlet temperature, Kelvin
        fluid_T_outlet = fluid_T_inlet - heat_transfer / cp_ea_fluid

        # Heat transfer rate out of this element, Watt
        matrix[nn,2,m] = heat_transfer / dt

        # Propagate to next fluid element unless at outlet
        if nn+1 < n:
            matrix[nn+1,0,m] = fluid_T_outlet

        # Final surface temperature, Kelvin
        # Temperature difference = Effective thermal resistance * Heat transfer rate
        surf_T_final = fluid_T_inlet - (conductive_rt + convective_rt)* heat_transfer / dt
        matrix[nn,1,m] = surf_T_final

        # Propagate to next time step unless at final time
        if m+1 < m_max:
            matrix[nn,1,m+1] = surf_T_final

        nn += 1
    m += 1

# this function displays a single snapshot in time of the radiator.
# choose your property by proper array indexing as shown in the example below
# matrix[0:n-1:1,1,m_max] grabs the surface temperature snapshot from the final time step
def displayreshaped_straight(mat, length, len_ea, width, figure,title):
    parallel = mat
    ng = 0
    # duplicates single tube across full width to achieve to-scale display
    while ng < n_tubes:
        parallel = np.vstack([parallel,mat])
        ng += 1

    # plot final surface temperature
    fig = plt.figure(figure)
    plot = plt.imshow(parallel,cmap='hot',interpolation='nearest',aspect=spacing/len_ea)
    fig.add_artist(plot)
    plt.title(title)
    plt.colorbar()

# this function displays a single snapshot in time of the radiator.
# choose your property by proper array indexing as shown in the example below
# matrix[0:n-1:1,1,m_max] grabs the surface temperature snapshot from the final time step
def displayreshaped(mat, length, len_ea, width, figure,title):
    # number of elements per straight length of tube
    n_per_straight = ceil(length / len_ea)
    # total length of array
    final_slice_length = mat.shape[0]
    # empty spaces to pad array with to enable reshaping into rectangle for display
    padding_number = n_per_straight - final_slice_length % n_per_straight
    # pad final array with values copied from the edge to avoid zeros messing up the colormap
    padded = np.pad(mat,(0,padding_number),'edge')
    # reshape final surface temperature slice into rectangle
    final_grid = np.reshape(padded,[-1,n_per_straight])

    # correct for zigzag tubing by flipping every other row around
    max_col = final_grid.shape[0]
    col = 1
    while col < max_col:
        final_grid[col,:] = np.flip(final_grid[col,:],0)
        col += 2

    # plot final surface temperature
    fig = plt.figure(figure)
    plot = plt.imshow(final_grid,cmap='hot',interpolation='nearest')
    fig.add_artist(plot)
    plt.title(title)
    plt.colorbar()

# surface temperature
surftemp = matrix[0:n-1:1,1,m_max-1]

# fluid temperature
fluidtemp = matrix[0:n-1:1,0,m_max-1] # Kelvin

# radiative heat transfer rate
radiation = matrix[0:n-1:1,2,m_max-1]/area_radiate # Watt / meter^2
print(surftemp.shape)
print(f"Total heat transfer: {np.sum(radiation*area_radiate):.2f} Watts")
print(f"Total area: {length*width:.2f} meter^2")
print(f"Total dry mass: {dens_panel*length*width:.2f} kg")
print(f"Total coolant mass: {path_length*diameter**2/4*np.pi*density:.2f} kg")

displayreshaped(surftemp,length,len_ea,width,1,"Final surface temperature (K)")
displayreshaped(fluidtemp,length,len_ea,width,2,"Final fluid temperature (K)")
displayreshaped(radiation,length,len_ea,width,3,"Final radiative heat transfer rate (Watt/meter^2)")
plt.show()
