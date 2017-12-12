import numpy as np
import matplotlib.pyplot as plt
import time
from matplotlib.animation import FuncAnimation

length = 150
movements = np.array([[1,0],[-1,0],[0,1],[0,-1]])
particleCap = 1001
steps = 1000
#collisionStore = np.zeros(shape=(particleCap,2))
#collisionStore[0]=initialSeed

def initialise(length,movements):
    lattice = np.zeros((length+1,length+1))
    initialSeed = [(length/2),(length/2)]
    lattice[int(initialSeed[0]),int(initialSeed[1])] = 1
    collisionStore = np.zeros(shape=(particleCap,2))
    collisionStore[0] = initialSeed
    #movements = np.array([[1,0],[-1,0],[0,1],[0,-1]])
    stickyInitial = np.nonzero(lattice)
    stickySites = []
    for n in range(4):
        stickySites.append((stickyInitial[0][0],stickyInitial[1][0]) + movements[n])
    stickyLattice = np.zeros((length+1,length+1))
    for n in range(4):
        stickyLattice[stickySites[n][0],stickySites[n][1]] = 1
    return lattice, stickyLattice, collisionStore

def pick_starting_position(radius):
    startingPosition = np.zeros(2)
    randomAngle = np.random.randint(0,359)
    angleRadians = np.deg2rad(randomAngle)
    x = radius*np.cos(angleRadians)
    y = radius*np.sin(angleRadians)
    startingPosition[0], startingPosition[1] = int((length/2)+x), int((length/2)+y)
    return startingPosition

def generate_path(steps,startingPosition):
    particlePosition = startingPosition
    direction = np.random.choice(4,steps)
    path = np.zeros((steps,2))
    path[0] = startingPosition
    for n in range(steps-1):
        path[n+1] = particlePosition + movements[direction[n]]
        particlePosition = path[n+1]
    return path

def particle_collision(logicCheck,path,lattice,stickyLattice,particleNumber,collisionStore):
    logicCheck2 = np.nonzero(logicCheck)
    collision = path[logicCheck2[0]]
    collisionStore[particleNumber] = collision
    lattice[int(collision[0,0]),int(collision[0,1])] = particleNumber + 1
    stickyInitialTest = np.nonzero(lattice)
    stickySitesTest = []
    for n in range(4):
        toDo = len(stickyInitialTest[0]) # how many sites to build around
        for o in range(toDo):
            stickySitesTest.append((stickyInitialTest[0][o],stickyInitialTest[1][o]) + movements[n])
    #stickyLattice = np.zeros((length+1,length+1))
    for n in range(4*toDo):
        stickyLattice[stickySitesTest[n][0],stickySitesTest[n][1]] = 1
    latticeCheck = np.nonzero(lattice) #######################
    latticeCheckLength = len(latticeCheck[0])
    for n in range(latticeCheckLength):
        stickyLattice[latticeCheck[0][n],latticeCheck[1][n]] = 0 ######################
    return lattice, stickyLattice, collisionStore

def lattice_radius_check(lattice,length,radius,killRadius):
    latticeRadiusInitial = np.nonzero(lattice)
    latticeRadiusY = length/2 - latticeRadiusInitial[0]
    latticeRadiusX = latticeRadiusInitial[1] - length/2
    latticeRadiusValues = np.sqrt(latticeRadiusX**2 + latticeRadiusY**2)
    latticeRadiusMax = max(latticeRadiusValues)
    
    if latticeRadiusMax > radius:
        radius = latticeRadiusMax*1.5
        killRadius = latticeRadiusMax*2
        #steps = steps
        return latticeRadiusMax, radius, killRadius
    else:
        return latticeRadiusMax, radius, killRadius

def kill_check(steps,path,length,startingPosition,moving,killRadius):
    killCheck = np.zeros(steps)
    killCheck = np.sqrt(((path[:,0]-length/2)**2) + ((path[:,1]-length/2)**2))
    toKill = np.any(killCheck > killRadius)
    if toKill: # Kill particle
        kill = 1
        moving = 0
        return kill,moving,startingPosition
    else: # Keep going with current particle
        kill = 0
        startingPosition = path[-1] # Last position as new starting position for next n steps
        return kill,moving,startingPosition

t0 = time.time()
lattice, stickyLattice, collisionStore = initialise(length, movements)
particleNumber = 1
kill = 1
radius = 15
killRadius = radius*2
moving = 1
#steps = 1000

while particleNumber < particleCap:
    
    if kill or not moving:
        startingPosition = pick_starting_position(radius)
    
    moving = 1
    kill = 0
    #particlePosition = startingPosition
    while moving:
        path = generate_path(steps,startingPosition)
        logicCheck = np.zeros(steps)
        for n in range(steps):
            try:
                logicCheck[n] = stickyLattice[int(path[n,0]),int(path[n,1])] == 1
            except IndexError:
                continue
            if logicCheck[n] == 1:
                
                lattice, stickyLattice, collisionStore = particle_collision(logicCheck, path, lattice, stickyLattice, particleNumber, collisionStore)
                latticeRadiusMax, radius, killRadius = lattice_radius_check(lattice,length,radius,killRadius) #,steps)
                
                moving = 0
                particleNumber += 1
                print('Particle', particleNumber-1, 'collided, generation radius is', radius)
                break
        if moving:
            kill,moving,startingPosition = kill_check(steps,path,length,startingPosition,moving,killRadius)

#plt.figure()
#cmap = plt.cm.plasma
#cmap.set_under(color='black')
#plt.pcolormesh(lattice, cmap=cmap, vmin = 0.0001)
#plt.title('2D Fractal')
#plt.xlim((length/2)-(latticeRadiusMax*1.25),(length/2)+(latticeRadiusMax*1.25))
#plt.ylim((length/2)-(latticeRadiusMax*1.25),(length/2)+(latticeRadiusMax*1.25))

fig, ax = plt.subplots()
cmap = plt.cm.plasma
cmap.set_under(color='black')
cax = ax.pcolormesh(lattice, cmap=cmap, vmin=0.1)

animationLattice = np.zeros(shape=(particleCap,length+1,length+1))
for n in range(particleCap):
    particle = 1
    for k in range(n+1):
        animationLattice[n,int(collisionStore[k,0]),int(collisionStore[k,1])] = particle
        particle += 1
 
t1 = time.time()
t2 = round((t1-t0)/60,1)
print('Finished in', t2, 'minutes.')
    
def animate(i):
    cax.set_array(animationLattice[i,:,:].flatten())
    plt.title('Animation of 2D fractal growth')

anim = FuncAnimation(fig, animate, interval=30, frames=particleCap)
plt.draw
plt.show

#t1 = time.time()
#t2 = round((t1-t0)/60,1)
#print('Finished in', t2, 'minutes.')
