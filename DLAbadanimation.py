import numpy as np
import matplotlib.pyplot as plt
import time
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
t0 = time.time()
# Setting up array for lattice, initial seed, starting position and movements
length = 100
particles = 50
lattice = np.zeros((length+1,length+1))
initialSeed = [(length/2)+1,(length/2)+1]
lattice[int(initialSeed[0]),int(initialSeed[1])] = -1
collisionStore = np.zeros(shape=(particles+1,2))
collisionStore[0]=initialSeed

cax = ax.pcolormesh(lattice)

movements = np.array([[1,0],[-1,0],[0,1],[0,-1]])

# Loop which simulates model
for n in range(particles):
    # Starting at boundary
    startingPosition = np.zeros(2)
    startingSide = np.random.choice(4) # 0 = top edge, 1 = right edge, 2 = bottom edge, 3 = left edge
    if startingSide == 0:
        startingPosition[0] = 0
        startingPosition[1] = np.random.randint(0,length)
    elif startingSide == 1:
        startingPosition[0] = np.random.randint(0,length)
        startingPosition[1] = length
    elif startingSide == 2:
        startingPosition[0] = length
        startingPosition[1] = np.random.randint(0,length)
    else:
        startingPosition[0] = np.random.randint(0,length)
        startingPosition[1] = 0
    # Movement of particle
    particlePosition = startingPosition # y position = [0], x position = [1]
    moving = 1
    #particlePath = [] # Debugging
    while moving == 1:
        inspectPosition = particlePosition
        # Moving particle
        inspectPosition = inspectPosition + movements[np.random.choice(4)]
        # Checking periodic boundaries of y values
        if inspectPosition[0] == -1:
            inspectPosition[0] = length
        if inspectPosition[0] == length + 1:
            inspectPosition[0] = 0
        # Checking periodic boundaries of x values
        if inspectPosition[1] == -1:
            inspectPosition[1] = length
        if inspectPosition[1] == length + 1:
            inspectPosition[1] = 0
        # Checking if it collides with seed
        if lattice[int(inspectPosition[0]),int(inspectPosition[1])] == -1:
            moving = 0
            #print('Particle', n+1, 'has collided.')
            lattice[int(particlePosition[0]),int(particlePosition[1])] = -1
            collisionStore[n+1] = particlePosition
            continue
        particlePosition = inspectPosition
        #particlePath.append(particlePosition)

animationLattice = np.zeros(shape=(particles+1,length+1,length+1))
for n in range(particles):
    for k in range(n+1):
        animationLattice[n,int(collisionStore[k,0]),int(collisionStore[k,1])] = -1

def animate(i):
    cax.set_array(animationLattice[i,:,:].flatten())
    plt.title('Animation of fractal growth')
    
anim = FuncAnimation(fig, animate, interval=100, frames=particles+1)
plt.draw
plt.show

t1 = time.time()
tFinal = t1-t0
print('Time taken %.2f' % tFinal)