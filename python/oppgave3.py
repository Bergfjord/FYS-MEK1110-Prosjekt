

import numpy as np
import matplotlib.pyplot as plt
import time
from tqdm import trange


class VVSolver():
    
    def __init__(self, position, velocity, dt, t0, t):
        
        self.dt = dt
        self.t0 = t0
        self.t = t
        timepoints = np.linspace(t0, t, int((t-t0)/dt))
        self.timepoints = timepoints
        
        # generate position vectors
        N = len(position)
        self.N = N
        
        self.positions = np.zeros((len(timepoints),N, 3))

        self.positions[0] = position
        self.velocities = np.zeros((len(timepoints),N, 3))
        self.velocities[0] = velocity
        
    
    
    def a(self, ri, rj):
        if (ri == rj).all() or np.linalg.norm(ri-rj)>3: # 3*(a) iii
            acc = np.zeros(len(ri))
        else:
            
            accfacc = 24* (2*(np.linalg.norm(ri-rj))**(-12) - (np.linalg.norm(ri-rj))**(-6))
            acc = accfacc * (ri-rj) /((np.linalg.norm(ri-rj))**(2))
            
        return acc
        
        
    def getaccelerations(self, pos):
        
        dr_matr = np.zeros((self.N, self.N , 3))
        
        for j in range(self.N):
        
            drs = pos[j+1:] - pos[j]
            
            dr_matr[j, j+1:] = drs
            dr_matr[j+1:, j] = -drs # Newtons third Law is applied here 
    
        #print(dr_matr)
        
        r = np.linalg.norm(dr_matr, axis = 2)
        
        r = np.where(r>3, 0, r) #Sets distance over 3 to 0, to be ignored
        #print(r)
        
        r1 = np.where(r!=0, r**-2, 0) # Ignores all 0's in acceleration calc
        #print(r1)
        
        
        r2 = 24*(2*r1**(6) - r1**(3))*r1 
        #print(r2)
        
        a = np.einsum('ij,ijk->jk', r2, dr_matr)
        
        
        return a
    
    
    def run(self):
        start = time.time()
        
        l = len(self.timepoints)
        positions = self.positions
        velocities = self.velocities
        dt = self.dt
        prevacc = self.getaccelerations(positions[0])
        print("Running simulation...")
        for i in trange(len(positions)-1):
            new_pos = positions[i] + velocities[i]*dt  + 0.5*prevacc*dt**2
            
            
            nextacc = self.getaccelerations(new_pos)
            velocities[i+1] = velocities[i] + 0.5*(prevacc+nextacc)*dt
            positions[i+1] = new_pos
            
            prevacc = nextacc
            
        
        self.positions = positions
        self.velocities = velocities
        print(f"Runtime: {int(time.time() - start)} s")
        return positions, velocities
    
    
    
    
    
    def ovito(self, name="delete_me"):
        pos = self.positions
        N = self.N
        with open(f"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/{name}.txt", "w") as outfile:
        
            for ind, timeval in enumerate(pos):
                outfile.write(f"{N}\n")
                outfile.write("type x y z\n")
                for i, cord in enumerate(timeval):
                    outfile.write(f"Ar {cord[0]} {cord[1]} {cord[2]}\n")



def Uscaled_shifted(r):
    
    
    constant_shift = 4 * (3**(-12) - 3**(-6))
    
    if r >3 or r == 0:
        return 0
    else:
        
        res = 4 * (r**(-12) - r**(-6)) - constant_shift
    
    return res


def Uscaled(r):    
    res = 4 * (r**(-12) - r**(-6))
    
    return res

# 3(c)i
def lattice(n, L):


    iterat = range(n)
    
    arr = []
    
    for i in iterat:
        for j in iterat:
            for k in iterat:
                arr.append([i,j,k])
                arr.append([i, 0.5+j,0.5+k])
                arr.append([0.5+i,j,0.5+k])
                arr.append([0.5+i,0.5+j,k])
    
    arr = np.array(arr)*(L/n)
    return arr

    

def task_3a():

    print("Running task_3a")
    def Uscaled(r):    
        res = 4 * (r**(-12) - r**(-6))
        
        return res
    
    
    points = np.linspace(2.5, 3.5, 100)
    org_pot = np.array([Uscaled(r) for r in points])
    new_pot = []
    for r in points:
        
        if r >= 3:
            new_pot.append(0)
        else:
            new_pot.append(Uscaled(r)-Uscaled(3))
        
    new_pot = np.array(new_pot)
    
    plt.plot(points, org_pot, label="Original potential")
    plt.plot(points, new_pot, label="New potential")
    plt.title("Shifted Potential")
    plt.xlabel("Sigma Distance")
    plt.ylabel("Potential")
    plt.legend()
    plt.grid()
    #plt.savefig()
    plt.show()
    
def task_3bi():
    print("Running task_3bi")
    dt =0.01
    t0 = 0
    t = 5
    timepoints = np.linspace(t0, t, int((t-t0)/dt))
    
    
    # generate position vectors
    
    N = 4
    
    init_positions = np.zeros((N, 3))
        
    # The two sets of atoms are place too far apart to interact (d>3) so
    # can simulate both cases simulatneously
    init_positions[0] = [1.5, 0, 0]
    init_positions[1] = [0, 0, 0]
    init_positions[2] = [0.95, 4, 0]
    init_positions[3] = [0, 4, 0]
    
    velocities = np.zeros((N, 3))
    
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, t0=t0, t=t)
    pos, vel = machine.run()
    
    machine.ovito("n4sigma095and150")
    
def task_3bii():
    print("Running task_3bii")
    
    dt =0.01
    t0 = 0
    t = 5
    
    N = 4
    
    init_positions = np.zeros((N, 3))
        
    init_positions[0] = [1, 0, 0]
    init_positions[1] = [0, 1, 0]
    init_positions[2] = [-1, 0, 0]
    init_positions[3] = [0, -1, 0]
    
    
    velocities = np.zeros((N, 3))
    
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, t0=t0, t=t)
    pos, vel = machine.run()
    
    machine.ovito("n4leafcloverformation")

def task_3biv():
    
    print("Running task_3biv")
    
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
        
        pot_all.append(np.sum(potential))
        
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

def task_3bv():
    
    print("Running task_3bv")
    
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
    
    machine.ovito("n4perturbation")
    
    
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
    
def task_3cii():

    print("Running task_3cii")
    
    n = 3
    L = 20
    
    init_positions = lattice(n=n, L=L)
    
    N = len(init_positions)
    
    # Indeed 108 atoms
    print(f"Box with sidelength {L} with {n} cells contain {N} atoms in total")
    
    
    dt =0.01
    t0 = 0
    t = 5

    velocities = np.zeros_like(init_positions)
    
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, t0=t0, t=t)
    pos, vel = machine.run()
    
    machine.ovito("latticestructure3cellsL20")

def task_3d():
    
    print("Running task_3d")
    
    
    n = 4
    L = n*1.7
    
    init_positions = lattice(n=n, L=L)
    
    N = len(init_positions)
    
    # Indeed 108 atoms
    print(f"Box with sidelength {L} with {n} cells contain {N} atoms in total")
    
    
    dt =0.01
    t0 = 0
    t = 5

    velocities = np.zeros_like(init_positions)
    
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, t0=t0, t=t)
    pos, vel = machine.run()
    
    machine.ovito("latticestructure3cellsL20")
    
    timepoints = machine.timepoints
    
    posit = pos[0]
    
    def calc_potential(pos, N):
    
        dr_matr = np.zeros((N, N , 3))
        #print(pos)
        for j in range(N):
        
            drs = pos[j+1:] - pos[j]
            dr_matr[j, j+1:] = drs
            dr_matr[j+1:, j] = -drs
    
        #print(dr_matr)
        
        r = np.linalg.norm(dr_matr, axis = 2)
        #print(r)
        r1 = np.where(r>3, 0, r)
        r1 = np.where(r!=0, r**-6, 0)
        
        #print(r1)
        
        r2 = np.nan_to_num(4*(r1**2 - r1))
        #print(r2)
        
        pot = np.sum(np.triu(r2))
        #print(pot)
        return pot  
    
    pot_all = []
    for posit in pos:
        pot_all.append(calc_potential(posit, N=N))
    
        
    kin_all = []
    for veloc in vel:
            kin_all.append(0.5*np.sum(veloc**2))
    
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
    
def task_3e():
    
    print("Running task_3e")
    
    # Here i will be defining a new solver, which i will be 
    # using (inititally at least) in task 4
    
    class VVSolver():
        
        def __init__(self, position, velocity, L, dt, t0, t):
            
            self.temp = []
            self.auto_cor = []
            self.L = L
            self.dt = dt
            self.t0 = t0
            self.t = t
            self.timepoints = np.linspace(t0, t, int((t-t0)/dt))
            
            
            # generate position vectors
            self.N = len(position)
            
            self.positions = np.zeros((len(self.timepoints),    self.N, 3))
            
            self.positions[0] = position
            self.velocities = np.zeros_like(self.positions)
            self.velocities[0] = velocity
            self.initial_velocity = velocity
        
        def getaccelerations(self, pos):
            
            dr_matr = np.zeros((self.N, self.N , 3))
            
            for j in range(self.N):
            
                drs = pos[j+1:] - pos[j]
                
                # Implementation of the "trick" to sense an atom on the "other side"
                # of the box
                drs = drs - np.round(drs/self.L)*self.L
    
                dr_matr[j, j+1:] = drs
                dr_matr[j+1:, j] = -drs # Newtons third Law is applied here 
        
            #print(dr_matr)
            
            r = np.linalg.norm(dr_matr, axis = 2)
            
            r = np.where(r>3, 0, r) #Sets distance over 3 to 0, to be ignored
            #print(r)
            
            r1 = np.where(r!=0, r**-2, 0) # Ignores all 0's in acceleration calc
            #print(r1)
            
            
            r2 = 24*(2*r1**(6) - r1**(3))*r1 
            #print(r2)
            
            a = np.einsum('ij,ijk->jk', r2, dr_matr)
            
            
            return a
        
        
        def run(self):
            start = time.time()
            
            l = len(self.timepoints)
            positions = self.positions
            velocities = self.velocities
            dt = self.dt
            prevacc = self.getaccelerations(positions[0])
            L = self.L
            print("Running simulation...")
            for i in trange(len(positions)-1):
                new_pos = positions[i] + velocities[i]*dt  + 0.5*prevacc*dt**2
                
                
                # Boundary conditions
                new_pos = new_pos - np.floor(new_pos/L)*L
                
                nextacc = self.getaccelerations(new_pos)
                velocities[i+1] = velocities[i] + 0.5*(prevacc+nextacc)*dt
                positions[i+1] = new_pos
                
                prevacc = nextacc
                
            
            self.positions = positions
            self.velocities = velocities
            print(f"Runtime: {int(time.time() - start)} s")
            return positions, velocities, self.timepoints
        
        def ovito(self, name="delete_me"):
            pos = self.positions
            N = self.N
            with open(f"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/{name}.txt", "w") as outfile:
            
                for ind, timeval in enumerate(pos):
                    outfile.write(f"{N}\n")
                    outfile.write("type x y z\n")
                    for i, cord in enumerate(timeval):
                        outfile.write(f"Ar {cord[0]} {cord[1]} {cord[2]}\n")
    
    
    dt =0.01
    t0 = 0
    t = 5
    
    N = 1
    
    L = 1.7
    
    
    init_positions = np.zeros((N, 3))
    init_positions[0] = [1, 0, 0]

    
    
    velocities = np.zeros_like(init_positions)
    velocities[0] = [1, 0, 0]
    
    
    machine = VVSolver(position=init_positions, velocity=velocities, L=L, dt=dt, t0=t0, t=t)
    pos, vel, timepoints = machine.run()
    
    machine.ovito("boundaryteleportation")
    
    
if __name__ == "__main__":
    # task_3a()
    # task_3bi()
    # task_3bii()
    # task_3biv()
    # task_3bv()
    # task_3cii()
    # task_3d()
    # task_3e()
    
    pass

