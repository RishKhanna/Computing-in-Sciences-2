
# set some basic stuff of the simulation
numParticles = 64
dim = 2   # dimensionlaity of system, 2 for ease of visualisation
kBT = 1.0 # temperature

# set interaction parameters: sigma and epsilon 
sigma = 1.0
epsilon = 1.0

# set integrator parameters
dt = 0.01

# setup arrays for storing stuff
positons = numpy.zeros( (numParticles, dim) )
momentas = numpy.zeros( (numParticles, dim) )
forces   = numpy.zeros( (numParticles, dim) )

### initailistaion functions

# initialization of position
def initialize_positions(boxLength):
    global positions, momentas, forces, numParticles, dim
    # set positions: lattice start inside a box of length 'boxLength'
    
    return 

# initialization of momentas
def initialize_momentas(temp):
    global positions, momentas, forces, numParticles, dim
    # set momenta, as Maxwell-Boltzmann distribution at temperaute 'temp'

    return

# funciton to compute forces on all particles in system
def calc_forces():
    global numParticles, sigma, epsilon
    
    for i in range(....):
        for j in range(....):
            forces[i] += compute_pairwiseLJ(i, j, sigma, epsilon)

    return

# function to compute pairwise forces
def compute_pairwiseLJ(i,j,sigma,epsilon):
    
    return

# function to compute total Energy of system
def compute_totalEnergy():
    global positions, momentas
    # calculate kinetic energy from momenta
    # calculate potential energy  using calc_energy_pairwiseLJ(i,j,sigma,epsilon)

    return

# integrator
def euler_integrate(dt):
    global positions, momemta, forces
    calc_forces()
    positions += dt * momenta / mass
    momenta += dt * forces
    return

def configDraw():
    global positions
    figure(figsize=(5,5)) # setup figure of 5 inches x 5 inches
    axis = gca()         # get current axes
    axis.set_xlim(0,10)  # xlim (0,10)
    axis.set_ylim(0,10)  # ylim (0,10)
    circles = []
    for i in range(self.numParticles):
        circles.append( axis.add_patch( Circle( positions[i], radius=0.5,
                                                linewidth=2, edgecolor='black') ) )
    show()
    return


#### 
initialize_positions(numParticles, dim)
initialize_momentas(numParticles, dim)

for steps in range(tot_steps):
    euler_integrate(dt)
    if (steps % 100 ==0) configDraw()

