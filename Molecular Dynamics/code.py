import numpy as np
from matplotlib import *

# set some basic stuff of the simulation
numParticles = 64
dim = 2   # dimensionlaity of system, 2 for ease of visualisation
kBT = 1.0 # temperature
mass = 40 * 1.6 * (10**-27)  # mass of argon gas

# set interaction parameters: sigma and epsilon 
sigma = 1.0 
epsilon = 1.0

# set integrator parameters
dt = 0.01

# setup arrays for storing stuff
positons = numpy.zeros( (numParticles, dim) )
momentas = numpy.zeros( (numParticles, dim) )
forces   = numpy.zeros( (numParticles, dim) )



# initialise_config
def initialize_positions( boxLength_in ):
    global positions, numParticles, dim, boxLength
    boxLength = boxLength_in
    # set positions: lattice start inside a box of length 'boxLength'
    positions = np.random.uniform(0,boxLength,(numParticles,dim))
    
    return


# initialization of momentas
def initialize_momentas(temp):
    global positions, momentas, forces, numParticles, dim, kBT, mass
    # set momenta, as Maxwell-Boltzmann distribution at temperaute 'temp'
    kBT = temp    
    momentas = np.random.normal(0, np.sqrt(2*mass*1.38*(10**-23)*kBT), (numParticles,dim))

    return



def calc_forces():
    global positions, sigma, epsilon, forces
    for idx, vec in enumerate(positions):
        r = np.linalg.norm(vec - positions, axis = 1)
        r = np.delete(r, idx)
        position_temp = np.delete(vec - positions, idx, 0)
        force_mag = 24*epsilon*( 2*((sigma/r)**12) - ((sigma/r)**6) )/(r*r)
        forces[idx]= np.sum(position_temp*force_mag [:, None], axis = 0)
    return     


def potential_curve():
    global positions, sigma, epsilon, mass, momentas
    Pot_Energy = 0
    for idx, vec in enumerate(positions):
        r = np.linalg.norm(vec-positions, axis = 1)
        r = np.delete(r, idx)
        Pot_Energy += 4*epsilon*( ((sigma/r)**12) - ((sigma/r)**6) )
    Pot_Energy = np.sum(Pot_Energy)
    
    Kin_Ener = np.sum(np.einsum('ij,ij->i', momentas, momentas))/(2*mass)
    
    Energy = Pot_Energy + Kin_Ener
    
    return Pot_Energy, Kin_Ener, Energy



# integrator
def euler_integrate(dt):
    global positions, momentas, forces, boxLength, mass
    calc_forces()
    positions += dt * (momentas / mass)
    # positions = np.mod(positions, boxLength)
    momentas += dt * forces
    return



def configDraw():
    global positions
    figure(figsize=(8,8)) # setup figure of 8 inches x 8 inches
    axis = gca()         # get current axes
    axis.set_xlim(0,(30))  # xlim (-15,15)
    axis.set_ylim(0,(30))  # ylim (-15,15)
    circles = []
    for i in range(numParticles):
        circles.append( axis.add_patch( Circle( positions[i], radius=0.5,
                                                linewidth=2, edgecolor='black') ) )
    grid()
    show()
    
    return




#### 
initialize_positions(2)
initialize_momentas(273)
dt = 1/np.power(10, 9)
tot_steps = 1000
Energy_plot = []
Kinetic = []
Tot = []
for steps in range(tot_steps):
    euler_integrate(dt)
    Energy_plot.append(potential_curve()[0])
    Kinetic.append(potential_curve()[1])
    Tot.append(potential_curve()[2])
    if (steps % 100==0):
        configDraw()