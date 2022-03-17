# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 13:13:52 2022

@author: Anders sin PC
"""



import numpy as np
import time
import matplotlib.pyplot as plt

class VVSolver():
    
    def __init__(self, position, velocity, dt, t0, t):
        
        
        self.dt = dt
        self.t0 = t0
        self.t = t
        self.timepoints = np.linspace(t0, t, int((t-t0)/dt))
        
        
        # generate position vectors
        self.N = len(position)
        
        self.positions = np.zeros((len(self.timepoints),self.N, 3))

        self.positions[0] = position
        self.velocities = np.zeros_like(self.positions)
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
            for j in range(i+1, N):
                accl = self.a(position[i], position[j])
                accelerations[i,j] = accl
                accelerations[j,i] = -accl
                    
        #print(accelerations)
        return accelerations
    
    
    def run(self, prog=False):
        start = time.time()
        
        l = len(self.timepoints)
        positions = self.positions
        velocities = self.velocities
        dt = self.dt
        prevacc = np.array([sum(ac) for ac in self.getaccelerations(positions[0])])
        
        
        if prog:
            for t_, val in enumerate(positions):
                print(t_)
                if t_ != l-1:
                    positions[t_+1] = positions[t_] + velocities[t_]*dt  + 0.5*prevacc*dt**2
                    nextacc = np.array([sum(ac) for ac in self.getaccelerations(positions[t_+1])])
                    velocities[t_+1] = velocities[t_] + 0.5*(prevacc+nextacc)*dt
                    prevacc = nextacc
    
        
            self.positions = positions
            self.velocities = velocities
            print(f"Runtime: {int(time.time() - start)} s")
        else:
            
            for t_, val in enumerate(positions):
                # if t_ == 50:
                #     break
                if t_ != l-1:
                    positions[t_+1] = positions[t_] + velocities[t_]*dt  + 0.5*prevacc*dt**2
                    nextacc = np.array([sum(ac) for ac in self.getaccelerations(positions[t_+1])])
                    velocities[t_+1] = velocities[t_] + 0.5*(prevacc+nextacc)*dt
                    prevacc = nextacc
    
        
            self.positions = positions
            self.velocities = velocities
            print(f"Runtime: {int(time.time() - start)} s")
        return positions, velocities
    
    
    
def Uscaled_shifted(r):
    
    constant_shift = 4 * (3**(-12) - 3**(-6))
    
    if r >3:
        return 0
    else:
        
        res = 4 * (r**(-12) - r**(-6)) - constant_shift
    
    return res
    
    
    
def bi():
    
    dt =0.01
    t0 = 0
    t = 5
    timepoints = np.linspace(t0, t, int((t-t0)/dt))
    
    
    # generate position vectors
    
    N = 4
    
    init_positions = np.zeros((N, 3))
        
    
    init_positions[0] = [1.5, 0, 0]
    init_positions[1] = [0, 0, 0]
    init_positions[2] = [0.95, 4, 0]
    init_positions[3] = [0, 4, 0]
    
    velocities = np.zeros((N, 3))
    
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, t0=t0, t=t)
    pos, vel = machine.run()
    
    with open("C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/results2.txt", "w") as outfile:
        
        for ind, timeval in enumerate(pos):
            outfile.write(f"{N}\n")
            outfile.write("type x y z\n")
            for i, cord in enumerate(timeval):
                outfile.write(f"Ar {cord[0]} {cord[1]} {cord[2]}\n")
            
    outfile.close()




def bii():
    # generate position vectors
    # 3b ii
    
    dt =0.01
    t0 = 0
    t = 5
    timepoints = np.linspace(t0, t, int((t-t0)/dt))
    
    N = 4
    
    init_positions = np.zeros((N, 3))
        
    init_positions[0] = [1, 0, 0]
    init_positions[1] = [0, 1, 0]
    init_positions[2] = [-1, 0, 0]
    init_positions[3] = [0, -1, 0]
    
    
    velocities = np.zeros((N, 3))
    
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, t0=t0, t=t)
    pos, vel = machine.run()
    
    with open("C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/results.txt", "w") as outfile:
        
        for ind, timeval in enumerate(pos):
            outfile.write(f"{N}\n")
            outfile.write("type x y z\n")
            for i, cord in enumerate(timeval):
                outfile.write(f"Ar {cord[0]} {cord[1]} {cord[2]}\n")
            
    outfile.close()
    
def biv():
    dt =0.01
    t0 = 0
    t = 5
    timepoints = np.linspace(t0, t, int((t-t0)/dt))
    
    N = 4
    
    init_positions = np.zeros((N, 3))
        
    init_positions[0] = [1, 0, 0]
    init_positions[1] = [0, 1, 0]
    init_positions[2] = [-1, 0, 0]
    init_positions[3] = [0, -1, 0]
    
    
    velocities = np.zeros((N, 3))
    
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, t0=t0, t=t)
    pos, vel = machine.run()
    
    posit = pos[0]
    
    
    pot_all = []
    for posit in pos:
        calculated =  []
        potential = np.zeros((len(posit), len(posit)))
        for i in range(len(posit)):
            for j in range(len(posit)):
                if i != j:
                    if f"{i}{j}" in calculated:
                        potential[i,j] = potential[j,i]
                    else:
        
                        potential[i,j] = Uscaled_shifted(np.linalg.norm(posit[i]-posit[j]))
                        calculated.append(f"{i}{j}")
        
        pot_all.append(sum([sum(vec) for vec in potential]))
        
    kin_all = []
    for veloc in vel:
            kin_all.append(sum([np.dot(v, v) for v in veloc]))
    
    
    plt.clf()
    plt.close()
    pot_all = np.array(pot_all)
    kin_all = np.array(kin_all)
    plt.plot(timepoints, pot_all, label = "Potential Energy")
    plt.plot(timepoints, kin_all, label = "Kinetic Energy")
    plt.plot(timepoints, pot_all+kin_all, label = "Total Energy")
    plt.grid()
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Energy")
    
    # plt.show()



def bvi():

    dt =0.01
    t0 = 0
    t = 5
    timepoints = np.linspace(t0, t, int((t-t0)/dt))
    
    N = 4
    
    init_positions = np.zeros((N, 3))
        
    init_positions[0] = [1, 0.1, 0]
    init_positions[1] = [0, 1, 0]
    init_positions[2] = [-1, 0, 0]
    init_positions[3] = [0, -1, 0]
    
    
    velocities = np.zeros_like(init_positions)
    
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, t0=t0, t=t)
    pos, vel = machine.run()
    
    with open("C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/results3d.txt", "w") as outfile:
        
        for ind, timeval in enumerate(pos):
            outfile.write(f"{N}\n")
            outfile.write("type x y z\n")
            for i, cord in enumerate(timeval):
                outfile.write(f"Ar {cord[0]} {cord[1]} {cord[2]}\n")
            
    outfile.close()

# 3b iv
def biv():
    
    dt =0.01
    t0 = 0
    t = 5
    timepoints = np.linspace(t0, t, int((t-t0)/dt))
    
    N = 4
    
    init_positions = np.zeros((N, 3))
        
    init_positions[0] = [1, 0, 0]
    init_positions[1] = [0, 1, 0]
    init_positions[2] = [-1, 0, 0]
    init_positions[3] = [0, -1, 0]
    
    
    velocities = np.zeros((N, 3))
    
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, t0=t0, t=t)
    pos, vel = machine.run()
    
    
    
    posit = pos[0]
    
    pot_all = []
    for posit in pos:
        calculated =  []
        potential = np.zeros((len(posit), len(posit)))
        for i in range(len(posit)):
            for j in range(len(posit)):
                if i != j:
                    if f"{i}{j}" in calculated:
                        potential[i,j] = potential[j,i]
                    else:
        
                        potential[i,j] = Uscaled_shifted(np.linalg.norm(posit[i]-posit[j]))
                        calculated.append(f"{i}{j}")
    
        pot_all.append(sum([sum(vec) for vec in potential]))
        
    
    kin_all = []
    for veloc in vel:
            kin_all.append(sum([np.dot(v, v) for v in veloc]))
    
    
    plt.clf()
    plt.close()
    pot_all = np.array(pot_all)
    kin_all = np.array(kin_all)
    plt.plot(timepoints, pot_all, label = "Potential Energy")
    plt.plot(timepoints, kin_all, label = "Kinetic Energy")
    plt.plot(timepoints, pot_all+kin_all, label = "Total Energy")
    plt.grid()
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Energy")
    
    plt.show()



def lulz():
    # Just for fun
    N = 200
    import random
    import time
    
    
    init_positions = np.zeros((N, 3))
    for i in range(N):
        
        position = np.array(random.sample(range(-500, 500), 3))/100
        init_positions[i] = position
    
    dt =0.01
    t0 = 0
    t = 7
    timepoints = np.linspace(t0, t, int((t-t0)/dt))
    
    velocities = np.zeros_like(init_positions)
    
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, t0=t0, t=t)
    pos, vel = machine.run(prog=True)
    
    with open(f"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/resultsN{N}s{t}.txt", "w") as outfile:
        
        for ind, timeval in enumerate(pos):
            outfile.write(f"{N}\n")
            outfile.write("type x y z\n")
            for i, cord in enumerate(timeval):
                outfile.write(f"Ar {cord[0]} {cord[1]} {cord[2]}\n")
            
    outfile.close()
    
    
    
    def Uscaled_shifted(r):
        constant_shift = 4 * (3**(-12) - 3**(-6))
        if r >3:
            return 0
        else:
            
            res = 4 * (r**(-12) - r**(-6)) - constant_shift
        
        return res
    
    
    posit = pos[0]
    
    pot_all = []
    for posit in pos:
        calculated =  []
        potential = np.zeros((len(posit), len(posit)))
        for i in range(len(posit)):
            for j in range(len(posit)):
                if i != j:
                    if f"{i}{j}" in calculated:
                        potential[i,j] = potential[j,i]
                    else:
        
                        potential[i,j] = Uscaled_shifted(np.linalg.norm(posit[i]-posit[j]))
                        calculated.append(f"{i}{j}")
    
        pot_all.append(sum([sum(vec) for vec in potential]))
        
    
    kin_all = []
    for veloc in vel:
            kin_all.append(sum([np.dot(v, v) for v in veloc]))
    
    
    plt.clf()
    plt.close()
    pot_all = np.array(pot_all)
    kin_all = np.array(kin_all)
    plt.plot(timepoints, pot_all, label = "Potential Energy")
    plt.plot(timepoints, kin_all, label = "Kinetic Energy")
    plt.plot(timepoints, pot_all+kin_all, label = "Total Energy")
    plt.grid()
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Energy")
    
    plt.show()
    
    

def optimize():
    
    N = 50
    import random
    import time
    
    
    init_positions = np.zeros((N, 3))
    for i in range(N):
        
        position = np.array(random.sample(range(-500, 500), 3))/100
        init_positions[i] = position
    
    dt =0.01
    t0 = 0
    t = 3
    timepoints = np.linspace(t0, t, int((t-t0)/dt))
    
    velocities = np.zeros_like(init_positions)
    
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, t0=t0, t=t)
    pos, vel = machine.run(#prog=True
                           )
    
    
    
    # posit = pos[0]
    
    # pot_all = []
    # for posit in pos:
    #     calculated =  []
    #     potential = np.zeros((len(posit), len(posit)))
    #     for i in range(len(posit)):
    #         for j in range(len(posit)):
    #             if i != j:
    #                 if f"{i}{j}" in calculated:
    #                     potential[i,j] = potential[j,i]
    #                 else:
        
    #                     potential[i,j] = Uscaled_shifted(np.linalg.norm(posit[i]-posit[j]))
    #                     calculated.append(f"{i}{j}")
    
    #     pot_all.append(sum([sum(vec) for vec in potential]))
        
    
    # kin_all = []
    # for veloc in vel:
    #         kin_all.append(sum([np.dot(v, v) for v in veloc]))
    
    
    # plt.clf()
    # plt.close()
    # pot_all = np.array(pot_all)
    # kin_all = np.array(kin_all)
    # plt.plot(timepoints, pot_all, label = "Potential Energy")
    # plt.plot(timepoints, kin_all, label = "Kinetic Energy")
    # plt.plot(timepoints, pot_all+kin_all, label = "Total Energy")
    # plt.grid()
    # plt.legend()
    # plt.xlabel("Time")
    # plt.ylabel("Energy")
    
    # plt.show()




if __name__ == "__main__":
    optimize()








