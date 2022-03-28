import numpy as np
import matplotlib.pyplot as plt
import random
import time
from tqdm import trange


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
    
    def plot_energy(self, mode="show"):
        def calc_potential(pos, N):
    
            dr_matr = np.zeros((N, N , 3))
            #print(pos)
            for j in range(N):
            
                drs = pos[j+1:] - pos[j]
                drs = drs - np.around(drs/self.L)*self.L
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


        pos = self.positions
        vel = self.velocities
        timepoints = self.timepoints
        pot_all = []
        for posit in pos:
            pot_all.append(calc_potential(posit, N=self.N))
            
        kin_all = []
        for veloc in vel:
                kin_all.append(0.5*np.sum(veloc**2))
        
        

        pot_all = np.array(pot_all)
        kin_all = np.array(kin_all)
        plt.plot(timepoints, pot_all, label = "Potential Energy")
        plt.plot(timepoints, kin_all, label = "Kinetic Energy")
        plt.plot(timepoints, pot_all+kin_all, label = "Total Energy")
        plt.grid()
        plt.legend()
        plt.title(f"Energy | N={self.N}")
        plt.xlabel("Time")
        plt.ylabel("Energy")
        
        if mode == "show":
            plt.show()
        else:
            plt.savefig(fr"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/figures/Energy_{mode}_{self.N}atoms.pdf")
        
        
        plt.clf()
        
    def plot_temp(self, mode="show"):
        velocities = self.velocities
        N = self.N
        timepoints = self.timepoints
        temp = []
        for veloc in velocities:
                temp.append(np.sum(veloc**2)/(3*N))
        

        temp = np.array(temp)
        temp = temp * 119.7
        self.temp = temp
        plt.plot(timepoints, temp, label = "Temperature")
        plt.grid()
        #plt.legend()
        plt.title(f"Avg. Temperature | N={self.N}")
        plt.xlabel("Time")
        plt.ylabel("Kelvin")
        
        if mode == "show":
            plt.show()
        else:
            plt.savefig(fr"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/figures/Temp_{mode}_{self.N}atoms.pdf")
        
        plt.clf()
        
    def ovito(self, name="delete_me"):
        pos = self.positions
        N = self.N
        with open(f"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/{name}.txt", "w") as outfile:
        
            for ind, timeval in enumerate(pos):
                outfile.write(f"{N}\n")
                outfile.write("type x y z\n")
                for i, cord in enumerate(timeval):
                    outfile.write(f"Ar {cord[0]} {cord[1]} {cord[2]}\n")
        
    def plot_vel_autocor(self, mode="show", runs=1):
        
        init_vel = self.initial_velocity
        
        init_dot = []
        
        for vel in init_vel:
            init_dot.append(np.dot(vel,vel))
        
        init_dot=np.sum(init_dot)
        
        #init_dot = np.einsum('ij->i',init_vel**2)
        velocities = self.velocities
        N = self.N
        timepoints = self.timepoints
        # velocities = velocities[start:]
        # timepoints = timepoints[start:]
        
        auto_cor = []
        for veloc in velocities:
            
            numerator_dot = []
            for ind, vel in enumerate(veloc):
                numerator_dot.append(np.dot(vel, init_vel[ind]))
            
            res = np.array(numerator_dot)
            res = np.sum(res)/init_dot
            #res = np.sum(np.einsum('ij->i',veloc*init_vel) / init_dot)
            
            auto_cor.append(res)
        
        auto_cor = np.array(auto_cor)/N
        
        
        
        for i in range(runs-1):
            self.initial_velocity = self.velocities[-1]
            self.velocities = np.zeros_like(self.velocities)
            self.velocities[0] = self.initial_velocity
            
            init_vel = self.initial_velocity
            init_dot = np.einsum('ij->i',init_vel**2)
            
            start_pos = self.positions[-1]
            
            self.positions = np.zeros_like(self.positions)
            self.positions[0] = start_pos
            
            self.positions, velocities, self.timepoints = self.run()
            


            auto_cor1 = []
            for veloc in velocities:
                
                res = np.sum(np.einsum('ij->i',veloc*init_vel) / init_dot)
                
                auto_cor1.append(res)
                
            auto_cor += np.array(auto_cor1)/N            


        auto_cor = auto_cor / runs

        plt.plot(timepoints, auto_cor, label = "Auto correleation")
        plt.grid()
        #plt.legend()
        plt.title(f"Auto correlation | N={self.N}")
        plt.xlabel("Time")
        plt.ylabel("A")
        self.auto_cor = auto_cor
        
        if mode == "show":
            plt.show()
        else:
            plt.savefig(fr"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/figures/Temp_{mode}_{self.N}atoms.pdf")
        
        plt.clf()     

    def diffusion(self, mode ="show"):
        auto_cor = self.auto_cor
        timepoints = self.timepoints
        diff = []
        fig = 0 
        for aut in auto_cor:
            fig +=aut
            diff.append(fig)
            
            
        diff =  np.array(diff)*self.dt/3
        
        
        plt.plot(timepoints, diff, label = "Diffusion")
        plt.grid()
        #plt.legend()
        plt.title(f"Diffusion | N={self.N}")
        plt.xlabel("Time")
        plt.ylabel("Displacement")
        self.diff = diff
        
        if mode == "show":
            plt.show()
        else:
            plt.savefig(fr"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/figures/diffusion_{mode}_{self.N}atoms.pdf")
        
        plt.clf()
            
    def msd(self, mode="show"):
            
        
        positions = self.positions
        init_pos = positions[0]
        N = self.N
        L = self.L
        timepoints = self.timepoints
        
        jumps = np.zeros_like(positions)
        prev_pos = positions[0]
        for ind, pos in enumerate(positions):
            cur_jump = np.where(pos-prev_pos<L*-0.6, -1, 0)
            cur_jump += np.where(pos-prev_pos>L*0.6, 1, 0)
            
            jumps[ind] = jumps[ind-1] + cur_jump

            prev_pos = pos
            

        msd = []
        for ind, pos in enumerate(positions):

            res = np.sum(np.linalg.norm(pos-L*jumps[ind] - init_pos, axis=1)**2)
            # print(pos)
            # print(L*jumps[ind])
            # print(init_pos)
            # print((pos-L*jumps[ind]) - init_pos)
            # print(res)

            
            msd.append(res)
        
        msd =  np.array(msd)/N
        
        
        plt.plot(timepoints, msd, label = "Mean Squared Displacement")
        plt.grid()
        #plt.legend()
        plt.title(f"Mean Squared Displacement | N={self.N}")
        plt.xlabel("Time")
        plt.ylabel("Mean Squared Displacement")
        
        
        if mode == "show":
            plt.show()
        else:
            plt.savefig(fr"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/figures/msd_{mode}_{self.N}atoms.pdf")
        
        plt.clf()

    def rdf(self, bin_edges, mode="show"):
        
        r = self.positions[-1]
        V = self.L**3
        """
        bin_edges = edges of bins. Typically np.linspace(0, rc, num_bins+1)
        for some cut-off rc.
        r = Nx3-array of positions of atoms at a given timestep.
        V = volume of system.
        """
        N = self.N
        bin_centres = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        bin_sizes = bin_edges[1:] - bin_edges[:-1]
        n = np.zeros_like(bin_sizes)
        for i in range(N):
            dr = r - r[i]
            dr = dr - np.round(dr/self.L)*self.L
            dr = np.linalg.norm(dr, axis=1) # Distances from atom i.
            
            n += np.histogram(dr, bins=bin_edges)[0] # Count atoms within each distance interval.
        n[0] = 0
        # Equation (7) on the preceding page:
        rdf = V / N**2 * n / (4 * np.pi * bin_centres**2 * bin_sizes)
        
        plt.plot(bin_centres, rdf, label = "Diffusion")
        #plt.legend()
        plt.title(f"Radial Distribution | N={self.N}")
        plt.xlabel("r")
        plt.ylabel("Radial Distribution")
        
        
        if mode == "show":
            plt.show()
        else:
            plt.savefig(fr"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/figures/rdf_{mode}.pdf")
        
        plt.clf()

    
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
    

def task_4_a():
    print("Running task_4_a")
    n=4
    L=1.7*n # Corresponds to d = 1.7
    
    init_positions = lattice(n=n, L=L)
    print((np.max(init_positions)))
    dt =0.01
    t0 = 0
    t = 5
    N = len(init_positions)
    print(N)
    # 256s
    T = 300/119.7
    velocities = np.random.normal(0, np.sqrt(T), size=(N,3))

    
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, L=L, t0=t0, t=t)
    pos, vel, timepoints = machine.run()
    machine.plot_temp()# Method in VVSolver class gives temperature over time 
    machine.ovito("temperature_measurement")
    print(np.mean(machine.temp)) 

    # Plotting with K=180 start temp to get 94.4 (ish) equilibrium temp
    
    T = 180/119.7
    velocities = np.random.normal(0, np.sqrt(T), size=(N,3))

    
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, L=L, t0=t0, t=t)
    pos, vel, timepoints = machine.run()
    machine.plot_temp()# Method in VVSolver class gives temperature over time 
    machine.ovito("temperature_measurement")
    print(np.mean(machine.temp)) 

def task_4_b():
    print("Running task_4_b")
    n=4
    L=1.7*n # Corresponds to d = 1.7
    
    init_positions = lattice(n=n, L=L)
    dt =0.01
    t0 = 0
    t = 10
    N = len(init_positions)
    print(N)
    T = 180/119.7
    velocities = np.random.normal(0, np.sqrt(T), size=(N,3))
    #velocities[50] = [1,1,1]
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, L=L, t0=t0, t=t)
    pos, vel, timepoints = machine.run()
    
    runs = 100 # Averaging over this many runs
    machine.plot_vel_autocor(mode=f"autcorr_N{N}_runs{runs}", runs=runs)
    
    machine.diffusion(mode=f"diffusionrandomstart")

def task_4_c():
    print("Running task_4_c")
    n=9
    L=1.7*n # Corresponds to d = 1.7
    
    init_positions = lattice(n=n, L=L)
    dt =0.01
    t0 = 0
    t = 15
    N = len(init_positions)
    print(f"{N} atoms are in play")
    T = 180/119.7
    velocities = np.random.normal(0, np.sqrt(T), size=(N,3))
    #velocities[50] = [1,1,1]
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, L=L, t0=t0, t=t)
    pos, vel, timepoints = machine.run()
    
    machine.msd(mode="fig")
    
def task_4_d():
    print("Running task_4_d")
    n=6
    L=1.7*n # Corresponds to d = 1.7
    
    init_positions = lattice(n=n, L=L)
    dt =0.01
    t0 = 0
    t = 5
    N = len(init_positions)
    # 256s
    T = 180/119.7
    velocities = np.random.normal(0, np.sqrt(T), size=(N,3))
    #velocities[50] = [1,1,1]
    machine = VVSolver(position=init_positions, velocity=velocities, dt=dt, L=L, t0=t0, t=t)
    pos, vel, timepoints = machine.run()
    runs=3
    machine.plot_vel_autocor(runs=runs)
    bin_edges = np.linspace(0, L*0.5, 200)
    machine.rdf(bin_edges = bin_edges, mode=f"rdf_binedge{round(bin_edges[0])}")
    
    

    
if __name__ == "__main__":
    # task_4_a()
    # task_4_b()
    # task_4_c()
    # task_4_d()
    pass
    
    
    
    
    
    
    
    
    