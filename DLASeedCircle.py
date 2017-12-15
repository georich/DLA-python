# DLA Including a circular seed, George Richards 4228068
import numpy as np
import matplotlib.pyplot as plt
import time

# Initial values and movements possible in 2D (a plane)
length = 1000
movements = np.array([[1,0], [-1,0], [0,1], [0,-1]])
particleCap = 4001
steps = 1000

def initialise(length, movements):
    # Seed particle is built from a circle of radius 10, 'sticky' sites will
    # as usual be all sites adjacent
    radiusCircle = 10
    circleX, circleY = np.ogrid[-radiusCircle:radiusCircle+1, -radiusCircle:radiusCircle+1]
    circleSelection = circleX*circleX + circleY*circleY <= radiusCircle*radiusCircle
    coordsToCheck = np.where(circleSelection == True)
    coords = np.zeros(shape=(len(coordsToCheck[0]),2))
    for n in range(len(coordsToCheck[0])):
        coords[n,0], coords[n,1] =  coordsToCheck[0][n], coordsToCheck[1][n]
    lattice = np.zeros((length+1,length+1))
    for n in range(len(coordsToCheck[0])):
        lattice[int(length/2-radiusCircle+coords[n,0]),int(length/2-radiusCircle+coords[n,1])] = 1
    stickyInitial = np.nonzero(lattice)
    stickySites = []
    stickyInitialLength = len(stickyInitial[0])
    for site in range(stickyInitialLength):
        for n in range(4):
            stickySites.append((stickyInitial[0][site],stickyInitial[1][site]) + movements[n])
    stickyLattice = np.zeros((length+1,length+1))
    for n in range(4*stickyInitialLength):
        stickyLattice[stickySites[n][0],stickySites[n][1]] = 1
    return lattice, stickyLattice

def pick_starting_position(radius):
    # Starting position generated via parametric equations and a random angle
    startingPosition = np.zeros(2)
    randomAngle = np.random.randint(0,359)
    angleRadians = np.deg2rad(randomAngle)
    x = radius*np.cos(angleRadians)
    y = radius*np.sin(angleRadians)
    startingPosition[0], startingPosition[1] = int((length/2)+x), int((length/2)+y)
    return startingPosition

def generate_path(steps, startingPosition):
    # Path of 1000 steps generated via randomly picking one of four directions
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
    stickyInitialTest = np.nonzero(lattice)
    stickySitesTest = []
    for n in range(4):
        toDo = len(stickyInitialTest[0]) # How many sites to build around
        for o in range(toDo):
            stickySitesTest.append((stickyInitialTest[0][o],stickyInitialTest[1][o]) + movements[n])
    for n in range(4*toDo):
        stickyLattice[stickySitesTest[n][0],stickySitesTest[n][1]] = 1
    # Remove sites where particles reside from the sticky lattice to remove
    # a bug where particles could deposit on same site
    latticeCheck = np.nonzero(lattice)
    latticeCheckLength = len(latticeCheck[0])
    for n in range(latticeCheckLength):
        stickyLattice[latticeCheck[0][n],latticeCheck[1][n]] = 0
    return lattice, stickyLattice

def lattice_radius_check(lattice, length, radius, killRadius):
    # Calculate distances of all particles from the seed particle, record
    # highest value as rmax for later recording
    latticeRadiusInitial = np.nonzero(lattice)
    latticeRadiusY = length/2 - latticeRadiusInitial[0]
    latticeRadiusX = latticeRadiusInitial[1] - length/2
    latticeRadiusValues = np.sqrt(latticeRadiusX**2 + latticeRadiusY**2)
    latticeRadiusMax = max(latticeRadiusValues)
    # Increase generation radius of new particles if rmax of the aggregate
    # is encroaching on it
    if latticeRadiusMax > radius:
        radius = latticeRadiusMax*1.5
        killRadius = latticeRadiusMax*2
        return latticeRadiusMax, radius, killRadius
    else:
        return latticeRadiusMax, radius, killRadius

def kill_check(steps, path, length, startingPosition, moving, killRadius):
    # Since particle did not collide with another, remove it if it walked
    # outside of the kill radius, otherwise continue with same particle
    killCheck = np.zeros(steps)
    killCheck = np.sqrt(((path[:,0]-length/2)**2) + ((path[:,1]-length/2)**2))
    toKill = np.any(killCheck > killRadius)
    if toKill: # Kill particle
        kill = 1
        moving = 0
        return kill,moving,startingPosition
    else: # Keep going with current particle
        kill = 0
        startingPosition = path[-1] # Last position as new starting position
        return kill, moving, startingPosition

# More initial values and setting up both arrays
t0 = time.time()
lattice, stickyLattice = initialise(length, movements)
particleNumber = 1
kill = 1 # Set as 1 to generate a starting position on first iteration
radius = 25
killRadius = radius*2
moving = 1

while particleNumber < particleCap:
    
    # Generate starting position again if needed
    if kill or not moving:
        startingPosition = pick_starting_position(radius)
            
    moving = 1
    kill = 0
    
    while moving:
        path = generate_path(steps, startingPosition)
        logicCheck = np.zeros(steps)
        for n in range(steps):
            try:
                # Check if any point on the path is equal to a 'sticky' site
                logicCheck[n] = stickyLattice[int(path[n,0]),int(path[n,1])] == 1
            except IndexError:
                # This allows the code to continue if particles 'exit' the
                # lattice, particles will be killed so has no adverse
                # effect
                continue
            if logicCheck[n] == 1:
                # Call functions to update both lattice arrays upon collision
                lattice, stickyLattice = particle_collision(logicCheck, path, lattice, stickyLattice, particleNumber)
                latticeRadiusMax, radius, killRadius = lattice_radius_check(lattice,length,radius,killRadius)
                
                moving = 0
                particleNumber += 1
                #print('Particle', particleNumber-1, 'collided, generation radius is', radius)
                break
                
        if moving:
            # Kill particle or carry on with same one
            kill, moving, startingPosition = kill_check(steps, path, length, startingPosition, moving, killRadius)

# Plot a coloured plot of the resultant aggregate
plt.figure()
cmap = plt.cm.plasma
cmap.set_under(color='black')
plt.pcolormesh(lattice, cmap=cmap, vmin = 0.1)
plt.title('Fractal: Circular Seed')
plt.xlim((length/2)-(latticeRadiusMax*1.25),(length/2)+(latticeRadiusMax*1.25))
plt.ylim((length/2)-(latticeRadiusMax*1.25),(length/2)+(latticeRadiusMax*1.25))
t1 = time.time()
t2 = round((t1-t0)/60,1)
print('Finished in', t2, 'minutes.')