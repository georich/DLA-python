import numpy as np
import matplotlib.pyplot as plt
import time

length = 2000
movements = np.array([[1,0],[-1,0],[0,1],[0,-1]])
particleCap = 10001
steps = 1000
attempts = 10
dataStore = np.zeros(shape=(particleCap-1,attempts))
filledDensityMeanStore = np.zeros(shape=(10,attempts))
radiusDensity = np.linspace(2,100,15)
latticeDataStore = np.zeros(shape=(attempts,length+1,length+1))

def initialise(length,movements):
    lattice = np.zeros((length+1,length+1))
    initialSeed = [(length/2),(length/2)]
    lattice[int(initialSeed[0]),int(initialSeed[1])] = 1
    #movements = np.array([[1,0],[-1,0],[0,1],[0,-1]])
    stickyInitial = np.nonzero(lattice)
    stickySites = []
    for n in range(4):
        stickySites.append((stickyInitial[0][0],stickyInitial[1][0]) + movements[n])
    stickyLattice = np.zeros((length+1,length+1))
    for n in range(4):
        stickyLattice[stickySites[n][0],stickySites[n][1]] = 1
    return lattice, stickyLattice

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

def particle_collision(logicCheck,path,lattice,stickyLattice,particleNumber):
    logicCheck2 = np.nonzero(logicCheck)
    collision = path[logicCheck2[0]]
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
    return lattice, stickyLattice

def lattice_radius_check(lattice,length,radius,killRadius):
    latticeRadiusInitial = np.nonzero(lattice)
    latticeRadiusY = length/2 - latticeRadiusInitial[0]
    latticeRadiusX = latticeRadiusInitial[1] - length/2
    latticeRadiusValues = np.sqrt(latticeRadiusX**2 + latticeRadiusY**2)
    latticeRadiusMax = max(latticeRadiusValues)
    #latticeRadiusData = []
    #latticeRadiusData.append([particleNumber,latticeRadiusMax])
    
    if latticeRadiusMax > radius:
        radius = latticeRadiusMax*1.5
        killRadius = latticeRadiusMax*2
        #steps = steps
        return latticeRadiusMax, radius, killRadius,# steps
    else:
        return latticeRadiusMax, radius, killRadius,# steps

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

for run in range(attempts):
    
    t0 = time.time()
    lattice, stickyLattice = initialise(length, movements)
    latticeRadiusData = []
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
                    
                    lattice, stickyLattice = particle_collision(logicCheck, path, lattice, stickyLattice, particleNumber)
                    latticeRadiusMax, radius, killRadius = lattice_radius_check(lattice,length,radius,killRadius) #,steps)
                    
                    latticeRadiusData.append([particleNumber,latticeRadiusMax])    
                    moving = 0
                    particleNumber += 1
                    #print('Particle', particleNumber-1, 'collided, generation radius is', radius)
                    break
            if moving:
                kill,moving,startingPosition = kill_check(steps,path,length,startingPosition,moving,killRadius)
    
    latticeDataStore[run,:,:] = lattice[:,:]
    plt.figure()
    cmap = plt.cm.plasma
    cmap.set_under(color='black')
    plt.pcolormesh(lattice, cmap=cmap, vmin = 0.0001)
    plt.title('Fractal no=%i' %(run+1))
    plt.xlim((length/2)-(latticeRadiusMax*1.25),(length/2)+(latticeRadiusMax*1.25))
    plt.ylim((length/2)-(latticeRadiusMax*1.25),(length/2)+(latticeRadiusMax*1.25))
    plt.pause(0.1)
    print('The radius of the fractal', run+1, 'is', latticeRadiusData[-1][1])
    N, r = zip(*latticeRadiusData)
    #dataStore = np.zeros(shape=(particleCap-1,10))
    dataStore[:,run] = r
    #dataStore[:,1+(run*2)] = N
    
    t1 = time.time()
    t2 = round((t1-t0)/60,1)
    print('Finished run', run+1, 'in', t2, 'minutes.')
    
    tCorrelation1 = time.time()
    particlePositions = np.nonzero(lattice)
    particlePositionsRows = particlePositions[0][:]
    particlePositionsColumns = particlePositions[1][:]
    particleCoords = np.zeros(shape=(particleCap,2))
    particleCoords[:,0],particleCoords[:,1] = particlePositionsRows, particlePositionsColumns
    filledDensityStore = np.zeros(shape=(10,particleCap))
    #filledDensityMeanStore = np.zeros(shape(10,10))
    for particle in range(particleCap-1):
        referenceParticleRow = particleCoords[particle,0]
        referenceParticleColumn = particleCoords[particle,1]
        
        #radiusDensity = np.linspace(2,20,10)
        filledDensity = np.zeros(10)
        #filledDensityStore = np.zeros(shape=(10,particleCap-1))
        for width in range(10):
            radiusCircle = int(radiusDensity[width])
    #position = 5
            repeats = 1
            emptyStore = np.zeros(repeats)
            filledStore = np.zeros(repeats)
            correlationLattice = np.zeros(shape=((2*radiusCircle)+1,(2*radiusCircle)+1))
            for row in range((2*radiusCircle)+1):
                correlationLattice[row] = lattice[int(referenceParticleRow+radiusCircle-row),int(referenceParticleColumn-radiusCircle):int(referenceParticleColumn+radiusCircle+1)]
            circleX, circleY = np.ogrid[-radiusCircle:radiusCircle+1, -radiusCircle:radiusCircle+1]
            circleSelection = circleX*circleX + circleY*circleY <= radiusCircle*radiusCircle
            coordsToCheck = np.where(circleSelection == True)
            coords = np.zeros(shape=(len(coordsToCheck[0]),2))
            for n in range(len(coordsToCheck[0])):
                coords[n,0], coords[n,1] =  coordsToCheck[0][n], coordsToCheck[1][n]
            filled = 0
            empty = 0
            for n in range(len(coords)):
                row = int(coords[n,0])
                column = int(coords[n,1])
                check = correlationLattice[row,column]
                if check == 0:
                    empty += 1
                else:
                    filled += 1
            emptyStore = empty
            filledStore = filled
            filledDensity[width] = filled/(empty+filled)
        filledDensityStore[:,particle] = filledDensity
        #print('Finished particle', particle+1)
    filledDensityMean = np.zeros(10)
    for n in range(10):
        filledDensityMean[n] = np.mean(filledDensityStore[n,:])
    filledDensityMeanStore[:,run] = filledDensityMean[:]

    tCorrelation2 = time.time()
    tCorrelationFinish = round(tCorrelation2 - tCorrelation1,2)
    print('Finished correlation', run+1, 'in', tCorrelationFinish, 'seconds')

N = np.linspace(1,particleCap-1,particleCap-1)
dataRadiusAverage = np.zeros(particleCap-1)
for n in range(particleCap-1):
    #dataRadiusAverage[n] = np.average(dataStore[n,0]+dataStore[n,2]+dataStore[n,4]+dataStore[n,6]+dataStore[n,8])
    dataRadiusAverage[n] = np.mean(dataStore[n,:])
plt.figure()
plt.plot(np.log(dataRadiusAverage),np.log(N), label='Averaged Data')
plt.xlabel('ln(r)')
plt.ylabel('ln(N)')
plt.title('Fractal growth over time')
gradient, intercept = np.polyfit(np.log(dataRadiusAverage),np.log(N),1)
plt.plot(np.log(dataRadiusAverage),gradient*np.log(dataRadiusAverage)+intercept, 'r--', label='Regression Line')
plt.legend()
print('Value of the fractal dimension is', gradient)
correlationFunctionValues = np.zeros(attempts)
for n in range(attempts):
    correlationFunctionValues[n] = np.mean(filledDensityMeanStore[n,:])
logFilledDensityMeans = np.log(filledDensityMeanStore)
errorRanges = np.zeros(attempts)
errorDensityAverages = np.zeros(attempts)
for n in range(attempts):
    errorRanges[n] = (np.max(logFilledDensityMeans[n,:]) - np.min(logFilledDensityMeans[n,:])) / 2
    errorDensityAverages[n] = errorRanges[n]/np.sqrt(attempts) 
plt.figure()
#plt.plot(np.log(radiusDensity),np.log(correlationFunctionValues),'x',label='Averaged Data')
plt.errorbar(np.log(radiusDensity),np.log(correlationFunctionValues),
             errorDensityAverages,marker='x',ls='none',label='Averaged Data',capsize=5)
plt.title('Radial density function for all particles')
plt.xlabel('ln(r)')
plt.ylabel('ln(C(r))')
gradient2, intercept2 = np.polyfit(np.log(radiusDensity),np.log(correlationFunctionValues),1)
plt.plot(np.log(radiusDensity),gradient2*np.log(radiusDensity)+intercept2,'r', label='Regression Line')
plt.legend()
print('Power law relation between C(r) and r is', gradient2)
powerlawDimension = 2 - abs(gradient2)
print('Fractal dimension according to power law relation is:', powerlawDimension)
print('Finished')