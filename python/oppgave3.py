

import numpy as np

class VVSolver():
    
    def __init__(self, position, velocity, dt, t0, t):
        
        self.dt = dt
        self.t0 = t0
        self.t = t
        self.timepoints = np.linspace(t0, t, int((t-t0)/dt))
        
        
        # generate position vectors
        self.N = len(position)
        
        self.positions = np.zeros((len(timepoints),N, 3))

        self.positions[0] = position
        self.velocities = np.zeros((len(timepoints),N, 3))
        self.velocities[0] = velocity
        
    
    
    def a(self, ri, rj):
        if (ri == rj).all() or np.linalg.norm(ri-rj)>3:
            acc = np.zeros(len(ri))
        else:
            
            accfacc = 24* (2*(np.linalg.norm(ri-rj))**(-12) - (np.linalg.norm(ri-rj))**(-6))
            acc = accfacc * (ri-rj) /((np.linalg.norm(ri-rj))**(2))
            
        return acc
        
        
    def getaccelerations(self, position):
        
        N = len(position)
        

        
        accelerations = np.zeros((N, N, 3))
        calculated = []
        for i in range(N):
            for j in range(N):
                if f"{j}{i}" in calculated:
                    accelerations[i,j] = -1 * accelerations[j,i]
                else:
                    accelerations[i,j] = self.a(position[i], position[j])
                    calculated.append(f"{i}{j}")
                    
    
        return accelerations
    
    
    def run(self):
        start = time.time()
        
        l = len(self.timepoints)
        positions = self.positions
        velocities = self.velocities
        for t_, val in enumerate(positions):
            if t_ != l-1:
            
                prevacc = np.array([sum(ac) for ac in self.getaccelerations(val)])
        
                for ind, pos in enumerate(val):
                    
                    positions[t_+1][ind] = positions[t_][ind] + velocities[t_][ind]*dt  + 0.5*prevacc[ind]*dt**2

                    
                    nextacc = np.array([sum(ac) for ac in self.getaccelerations(positions[t_+1])])
                    
                    velocities[t_+1][ind] = velocities[t_][ind] + 0.5*(prevacc[ind]+nextacc[ind])*dt
                    
                
                    prevacc = nextacc
                break
            
        self.positions = positions
        self.velocities = velocities
        print(f"Runtime: {int(time.time() - start)} s")
        return positions, velocities
    


dt =0.01
t0 = 0
t = 5
timepoints = np.linspace(t0, t, int((t-t0)/dt))


# generate position vectors
N = 10
import random
import time


init_positions = np.zeros((N, 3))
for i in range(N):
    
    position = np.array(random.sample(range(-100, 100), 3))/100
    init_positions[i] = position


velocities = np.zeros((N, 3))

machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, t0=t0, t=t)
pos, vel = machine.run()

with open("C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/results2.txt", "w") as outfile:
    
    for ind, timeval in enumerate(pos):
        if ind%5 ==0:
            outfile.write(f"{N}\n")
            outfile.write("type x y z\n")
            for i, cord in enumerate(timeval):
                outfile.write(f"Ar {cord[0]} {cord[1]} {cord[2]}\n")
            
outfile.close()

# No the change does not impact



