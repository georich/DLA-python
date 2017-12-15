# DLA Square Border Growth, George Richards 4228068
import numpy as np
import matplotlib.pyplot as plt 
import time

# Initial values and movements possible in 2D (a plane)
t0 = time.time()
movements = np.array([[1,0], [-1,0], [0,1], [0,-1]])
length = 500
particleCap = 10001

def initialise(length, movements):
    # Marking all seed particles as the entire border of the square lattice/array
    lattice = np.zeros(shape=(length+1,length+1))
    lattice[0,:] = 1
    lattice[-1,:] = 1
    lattice[:,0] = 1
    lattice[:,-1] = 1
    stickyLattice = np.zeros(shape=(length+1,length+1))
    stickyLattice[1:-1,1] = 1
    stickyLattice[-2,1:-1] = 1
    stickyLattice[1:-1,-2] = 1
    stickyLattice[1,1:-1] = 1
    return lattice, stickyLattice

def generate_path(steps, startingPosition):
    # Path of 1000 steps generated via randmly picking one of four directions
    particlePosition = startingPosition
    direction = np.random.choice(4,steps)
    path = np.zeros((steps,2))
    path[0] = startingPosition
    for n in range(steps-1):
        path[n+1] = particlePosition + movements[direction[n]]
        particlePosition = path[n+1]
    return path

def particle_collision(logicCheck, path, lattice, stickyLattice, particleNumber):
    # Mark in the lattice which site the particle would have colliede at with
    # successively higher numbers in order to generate the coloured plot
    logicCheck2 = np.nonzero(logicCheck)
    collision = path[logicCheck2[0]]
    lattice[int(collision[0,0]),int(collision[0,1])] = particleNumber + 1
    stickySitesTest = []
    for n in range(4):
        stickySitesTest.append((collision[0,0],collision[0,1]) + movements[n])
    for n in range(4):
        stickyLattice[int(stickySitesTest[n][0]),int(stickySitesTest[n][1])] = 1
    # Remove sites where particles reside from the sticky lattice to remove
    # a bug where particles could deposit on same site
    latticeCheck = np.nonzero(lattice)
    latticeCheckLength = len(latticeCheck[0])
    for n in range(latticeCheckLength):
        stickyLattice[latticeCheck[0][n],latticeCheck[1][n]] = 0
    return lattice, stickyLattice

# More initial values and setting up arrays
lattice, stickyLattice = initialise(length,movements)
particleNumber = 1
steps = 2000

while particleNumber < particleCap:
    moving = 1
    # In this case starting position will always be set as the centre of the lattice
    startingPosition = np.array([length/2,length/2])
    while moving:
        path = generate_path(steps, startingPosition)
        logicCheck = np.zeros(steps)
        for n in range(steps):
            try:
                logicCheck[n] = stickyLattice[int(path[n,0]),int(path[n,1])] == 1
            except IndexError:
                continue
            if logicCheck[n] == 1:
                # Call functions to update both lattice arrays upon collision
                lattice, stickyLattice = particle_collision(logicCheck, path, lattice, stickyLattice, particleNumber)
                moving = 0
                particleNumber += 1
                print('Particle', particleNumber-1, 'has collided')
                break
        if moving:
            startingPosition = path[-1]
# Plots a coloured plot of the resultant aggregate
plt.figure()
cmap = plt.cm.plasma
cmap.set_under(color='black')
plt.pcolormesh(lattice, cmap=cmap, vmin = 0.0001)
plt.title('Fractal growth with a square border')
t1 = time.time()
t2 = round(t1-t0,2)
print('Finished in',t2/60,'minutes')         