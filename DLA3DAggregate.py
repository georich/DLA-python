# DLA 3D Aggregate, George Richards 4228068
import numpy as np
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D

# Initial values and movements possible in 3D
length = 200
particleCap = 4001
steps = 1000
attempts = 10
# Imagining a cube, can move towards each face, corner or midpoint of each edge
# listed in that order respectively
movements = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1],
                      [1,1,1],[1,1,-1],[1,-1,1],[1,-1,-1],[-1,1,1],[-1,-1,1],
                      [-1,-1,-1],[-1,1,-1],[-1,0,1],[0,-1,1],[1,0,1],[0,1,1],
                      [-1,-1,0],[-1,1,0],[1,1,0],[1,-1,0],[1,0,-1],[-1,0,-1],
                      [0,1,-1],[0,-1,-1]])
# Sticky sites as each site connected to a face of the deposited particle
stickyBuilding = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
dataStore = np.zeros(shape=(particleCap-1,attempts))

def initialise():
    # Seed particle needs to be marked in the middle of the array/lattice
    # 'Sticky' site are the ones adjacent to the seed, marked in a seperate
    # array to not interfere with colouring of the later plot
    lattice = np.zeros(shape=(length+1,length+1,length+1))
    initialSeed = [(length/2),(length/2),(length/2)]
    lattice[int(initialSeed[0]),int(initialSeed[1]),int(initialSeed[2])] = 1
    stickyInitial = np.nonzero(lattice)
    stickySites = []
    for n in range(6):
        stickySites.append((stickyInitial[0][0],stickyInitial[1][0],stickyInitial[2][0]) + stickyBuilding[n])
    stickyLattice = np.zeros(shape=(length+1,length+1,length+1))
    for n in range(6):
        stickyLattice[stickySites[n][0],stickySites[n][1],stickySites[n][2]] = 1
    return lattice, stickyLattice

def pick_starting_position(radius):
    # Starting position generated via spherical polar coordinates
    startingPosition = np.zeros(3)
    randomAnglePhi = np.random.randint(0,359)
    angleRadiansPhi = np.deg2rad(randomAnglePhi)
    randomAngleTheta = np.random.randint(0,179)
    angleRadiansTheta = np.deg2rad(randomAngleTheta)
    x = radius*np.cos(angleRadiansPhi)*np.sin(angleRadiansTheta)
    y = radius*np.sin(angleRadiansPhi)*np.sin(angleRadiansTheta)
    z = radius*np.cos(angleRadiansTheta)
    startingPosition[0], startingPosition[1], startingPosition[2] = int((length/2)+y), int((length/2)+z), int((length/2)+x)
    return startingPosition

def generate_path(steps, startingPosition):
    # Path of 1000 steps generated via randomly picking one of the directions
    particlePosition = startingPosition
    direction = np.random.choice(len(movements),steps)
    path = np.zeros(shape=(steps,3))
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
    lattice[int(collision[0,0]), int(collision[0,1]), int(collision[0,2])] = particleNumber + 1
    stickyInitialTest = np.nonzero(lattice)
    stickySitesTest = []
    for n in range(6):
        toDo = len(stickyInitialTest[0]) # How many sites to build around
        for o in range(toDo):
            stickySitesTest.append((stickyInitialTest[0][o], stickyInitialTest[1][o], stickyInitialTest[2][o]) + stickyBuilding[n])
    for n in range(6*toDo):
        stickyLattice[stickySitesTest[n][0],stickySitesTest[n][1],stickySitesTest[n][2]] = 1
    # Remove sites where particles reside from the sticky lattice to remove
    # a bug where particles could deposit on same site
    latticeCheck = np.nonzero(lattice)
    latticeCheckLength = len(latticeCheck[0])
    for n in range(latticeCheckLength):
        stickyLattice[latticeCheck[0][n], latticeCheck[1][n], latticeCheck[2][n]] = 0
    return lattice, stickyLattice

def lattice_radius_check(lattice, length, radius, killRadius):
    # Calculate distances of all particles from the seed particle, record
    # highest value as rmax for later recording
    latticeRadiusInitial = np.nonzero(lattice)
    latticeRadiusY = length/2 - latticeRadiusInitial[0]
    latticeRadiusX = latticeRadiusInitial[2] - length/2
    latticeRadiusZ = latticeRadiusInitial[1] - length/2
    latticeRadiusValues = np.sqrt((latticeRadiusX**2) + (latticeRadiusY**2) + (latticeRadiusZ**2))
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
    killCheck = np.sqrt(((path[:,0]-length/2)**2) + ((path[:,1]-length/2)**2) + ((path[:,2]-length/2)**2))
    toKill = np.any(killCheck > killRadius)
    if toKill: # Kill particle
        kill = 1
        moving = 0
        return kill, moving, startingPosition
    else: # Keep going with current particle
        kill = 0
        startingPosition = path[-1] # Last position as new starting position
        return kill, moving, startingPosition

for run in range(attempts):
    t0 = time.time()
    
    # More initial values
    lattice, stickyLattice = initialise()
    particleNumber = 1
    kill = 1 # Set as 1 to generate a starting position on first iteration
    radius = 25
    killRadius = radius*2
    moving = 1
    latticeRadiusData = []
    
    while particleNumber < particleCap:
        
        # Generate starting position again if needed
        if kill or not moving:
            startingPosition = pick_starting_position(radius)
        
        moving = 1
        kill = 0
        
        while moving:
            path = generate_path(steps,startingPosition)
            logicCheck = np.zeros(steps)
            for n in range(steps):
                try:
                    # Check if any point on the path is equal to a 'sticky' site
                    logicCheck[n] = stickyLattice[int(path[n,0]),int(path[n,1]),int(path[n,2])] == 1
                except IndexError:
                    # This allows the code to continue if particles 'exit' the
                    # lattice, particles will be killed so has no adverse
                    # effect
                    continue
                if logicCheck[n] == 1:
                    
                    # Call functions to update both lattice arrays upon collision
                    lattice, stickyLattice = particle_collision(logicCheck, path, lattice, stickyLattice, particleNumber)
                    latticeRadiusMax, radius, killRadius = lattice_radius_check(lattice, length, radius, killRadius)
                    
                    latticeRadiusData.append([particleNumber, latticeRadiusMax])
                    moving = 0
                    #print('Particle', particleNumber, 'collided, generation radius is', radius)
                    particleNumber += 1
                    
                    break
            if moving:
                # Kill particle or carry on with same one
                kill, moving, startingPosition = kill_check(steps, path, length, startingPosition, moving, killRadius)

    N, r = zip(*latticeRadiusData)
    dataStore[:,run] = r
    
    # Record positions of the particles within the lattice
    latticePlot = np.nonzero(lattice)
    ValuesLength = len(latticePlot[0])
    latticeValues = []
    for n in range(ValuesLength):
        latticeValues.append(lattice[(latticePlot[0][n],latticePlot[1][n],latticePlot[2][n])])

    # Produce an array of colours based on when the particle was deposited
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    cmap = plt.cm.plasma
    colourArrangement = np.zeros(particleCap)
    for n in range(particleCap):
        colourArrangement[n] = lattice[latticePlot[0][n],latticePlot[1][n],latticePlot[2][n]]
    ax.scatter(latticePlot[2],latticePlot[1],latticePlot[0],c=colourArrangement,cmap=cmap,marker='s')
    plt.title('FCM 3D Fractal no=%i' %(run+1))
    t1 = time.time()
    t2 = round((t1-t0)/60,1)
    print('Finished run', run+1, 'in', t2, 'minutes.')

# Plotting the log-log plot of N vs r
dataRadiusAverage = np.zeros(particleCap-1)
for n in range(particleCap-1):
    dataRadiusAverage[n] = np.mean(dataStore[n,:])
plt.figure()
plt.plot(np.log(dataRadiusAverage),np.log(N),label='Averaged Data')
plt.xlabel('ln(r)')
plt.ylabel('ln(N)')
plt.title('FCM 3D fractal growth over time')
# Plot the gradient and print its value
gradient, intercept = np.polyfit(np.log(dataRadiusAverage),np.log(N),1)
plt.plot(np.log(r),gradient*np.log(r)+intercept, 'r--', label='Regression Line')
plt.legend()
print('Value of the fractal dimension is', gradient)