# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 15:15:55 2022

@author: Anders sin PC
"""



import numpy as np



class Solver():
    def __init__(self, ri0, rj0, t0, t, dt):
        
        self.ri0, self.rj0, self.t0, self.t, self.dt = ri0, rj0, t0, t, dt
    
    
    def a(self, ri, rj):
        
        accfacc = 24* (2*(np.linalg.norm(ri-rj))**(-12) - (np.linalg.norm(ri-rj))**(-6))
        acc = accfacc * (ri-rj) /((np.linalg.norm(ri-rj))**(2))
        
        return acc
        
    def Euler(self):
        print(f"Euler inititated with values ri0={self.ri0}, rj0={self.rj0}, t0={self.t0}, t={self.t}, dt={self.dt}")
        dt = self.dt
        t0 = self.t0
        t = self.t
        timepoints = np.linspace(t0, t, int((t-t0)/dt))
        
        vi = np.zeros((len(timepoints), 3))
        vj = np.zeros((len(timepoints), 3))
        
        ri = np.zeros((len(timepoints), 3))
        ri[0] = self.ri0
        rj = np.zeros((len(timepoints), 3))
        rj[0] = self.rj0
        
        a = np.zeros((len(timepoints), 3))
        
        dist = np.zeros(len(timepoints))
        dist[0] = np.linalg.norm(ri[0]-rj[0])
        
        
        for i in range(len(timepoints)-1):
            a[i] = self.a(ri[i], rj[i])
            # print(an)
            # print(vi[i] + an * dt)
            # print(vj[i] - an * dt)
            
            vi[i+1] = vi[i] + a[i] * dt
            vj[i+1] = vj[i] - a[i] * dt
            
            ri[i+1] = ri[i] + dt*vi[i]
            rj[i+1] = rj[i] + dt*vj[i]
            
            dist[i+1] = np.linalg.norm(ri[i+1]-rj[i+1])
            # print(dist[i+1])
        return (timepoints, ri, rj, vi, vj, a, dist)
    
    def EulerCromer(self):
        print(f"EulerCromer inititated with values ri0={self.ri0}, rj0={self.rj0}, t0={self.t0}, t={self.t}, dt={self.dt}")
        dt = self.dt
        t0 = self.t0
        t = self.t
        timepoints = np.linspace(t0, t, int((t-t0)/dt))
        
        vi = np.zeros((len(timepoints), 3))
        vj = np.zeros((len(timepoints), 3))
        
        ri = np.zeros((len(timepoints), 3))
        ri[0] = self.ri0
        rj = np.zeros((len(timepoints), 3))
        rj[0] = self.rj0
        
        a = np.zeros((len(timepoints), 3))
        
        dist = np.zeros(len(timepoints))
        dist[0] = np.linalg.norm(ri[0]-rj[0])
        
        
        
        for i in range(len(timepoints)-1):
            a[i] = self.a(ri[i], rj[i])
            # print(an)
            # print(vi[i] + an * dt)
            # print(vj[i] - an * dt)
            
            vi[i+1] = vi[i] + a[i] * dt
            vj[i+1] = vj[i] - a[i] * dt
            
            ri[i+1] = ri[i] + dt*vi[i+1]
            rj[i+1] = rj[i] + dt*vj[i+1]
            
            dist[i+1] = np.linalg.norm(ri[i+1]-rj[i+1])
            # print(dist[i+1])
        return (timepoints, ri, rj, vi, vj, a, dist)
    
    
    def VelocityVerlet(self):
        print(f"VelocityVerlet inititated with values ri0={self.ri0}, rj0={self.rj0}, t0={self.t0}, t={self.t}, dt={self.dt}")
        dt = self.dt
        t0 = self.t0
        t = self.t
        timepoints = np.linspace(t0, t, int((t-t0)/dt))
        
        vi = np.zeros((len(timepoints), 3))
        vj = np.zeros((len(timepoints), 3))
        
        ri = np.zeros((len(timepoints), 3))
        ri[0] = self.ri0
        rj = np.zeros((len(timepoints), 3))
        rj[0] = self.rj0
        
        rj = np.zeros((len(timepoints), 3))
        
        a = np.zeros((len(timepoints), 3))
        a[0] = self.a(ri[0], rj[0])
        
        dist = np.zeros(len(timepoints))
        dist[0] = np.linalg.norm(ri[0]-rj[0])
        
        
        for i in range(len(timepoints)-1):
            
            
            
            ri[i+1] = ri[i] + dt*vi[i] + 0.5 * a[i] * dt**2
            rj[i+1] = rj[i] + dt*vj[i] - 0.5 * a[i] * dt**2
            
            a[i+1] = self.a(ri[i+1], rj[i+1])
            
            vi[i+1] = vi[i] + 0.5*(a[i] + a[i+1])*dt
            vj[i+1] = vj[i] - 0.5*(a[i] + a[i+1])*dt
            
            
            
            dist[i+1] = np.linalg.norm(ri[i+1]-rj[i+1])
            # print(dist[i+1])
        return (timepoints, ri, rj, vi, vj, a, dist)
    







# Testrun
import matplotlib.pyplot as plt
ri0 = np.array([2**(1/6)-0.1,0,0])
rj0 = np.array([0,0,0])
Machine = Solver(ri0, rj0, 0, 2, 0.001)

(timepoints, ri, rj, vi, vj, a, dist) = Machine.Euler()


plt.plot(timepoints, dist)
plt.show()


(timepoints, ri, rj, vi, vj, a, dist) = Machine.EulerCromer()

plt.plot(timepoints, dist)
plt.show()


(timepoints, ri, rj, vi, vj, a, dist) = Machine.VelocityVerlet()


plt.plot(timepoints, dist)
plt.show()

# 2b i
# Running along x-axis for simplicity
ri0 = np.array([1.5,0,0])
rj0 = np.array([0,0,0])
Machine = Solver(ri0, rj0, t0=0, t=5, dt=0.01)
(timepoints, ri, rj, vi, vj, a, dist) = Machine.EulerCromer()

# ii
plt.plot(timepoints, dist)
plt.xlabel("Time")
plt.ylabel("Distance")
plt.title("Atom motion | Sigma = 1.5")
#plt.show()
plt.savefig(r"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/motion150sigma.pdf")

# iii
# Fits preciciely as postulated from task 1 a. See pdf for more thourough explenation

# iv
ri0 = np.array([0.95,0,0])
rj0 = np.array([0,0,0])
Machine = Solver(ri0, rj0, t0=0, t=5, dt=0.01)
(timepoints, ri, rj, vi, vj, a, dist) = Machine.EulerCromer()

plt.plot(timepoints, dist)
plt.xlabel("Time")
plt.ylabel("Distance")
plt.title("Atom motion | 0.95 Sigma")
plt.show()
plt.savefig(r"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/motion095sigma.pdf")
# We see that the particle start oo close, with positive potential. 
# This causes the particle to accelerate to velocities which 
# cannot be counteracted on the other side of the lennard jones potential .
# Thus, we see the particle move almost lenearly wrt time after te initial inteaction.



    
    



ri0 = np.array([1.5,0,0])
rj0 = np.array([0,0,0])
Machine = Solver(ri0, rj0, t0=0, t=5, dt=0.01)
(timepoints, ri, rj, vi, vj, a, dist) = Machine.EulerCromer()

# ii
plt.plot(timepoints, dist, label="1.5 Sigma")
ri0 = np.array([0.95,0,0])
rj0 = np.array([0,0,0])
Machine = Solver(ri0, rj0, t0=0, t=5, dt=0.01)
(timepoints, ri, rj, vi, vj, a, dist) = Machine.EulerCromer()

plt.plot(timepoints, dist, label="0.95 Sigma")
plt.xlabel("Time")
plt.ylabel("Distance")
plt.title("Atom motion")
plt.legend()
plt.show()
plt.savefig(r"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/motioncomparison.pdf")



# 2c


def Uscaled(r):    
    res = 4 * (r**(-12) - r**(-6))
    
    return res


ri0 = np.array([1.5,0,0])
rj0 = np.array([0,0,0])
Machine = Solver(ri0, rj0, t0=0, t=5, dt=0.01)
(timepoints, ri, rj, vi, vj, a, dist) = Machine.EulerCromer()

Ek = 0.5*np.array([np.linalg.norm(v)**2 for v in vi]) # Standard Kinetic energy formula

Ep = 0.5*np.array([Uscaled(d) for d in dist]) # Since dist is traversed half by atom i and half by atom j

Etot = Ek + Ep




# ii
plt.plot(timepoints, dist, label="Position")
plt.plot(timepoints, np.array([np.linalg.norm(v) for v in vi]), label="Speed")
plt.plot(timepoints, np.array([np.linalg.norm(v) for v in a]), label="Acceleration")
plt.legend()
plt.grid()
plt.show()










