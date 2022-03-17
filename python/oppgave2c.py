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
            rj[i+1] = rj[i] + dt*vj[i] + 0.5 - a[i] * dt**2
            
            a[i+1] = self.a(ri[i+1], rj[i+1])
            
            vi[i+1] = vi[i] + 0.5*(a[i] + a[i+1])*dt
            vj[i+1] = vj[i] - 0.5*(a[i] + a[i+1])*dt
            
            
            
            dist[i+1] = np.linalg.norm(ri[i+1]-rj[i+1])
            # print(dist[i+1])
        return (timepoints, ri, rj, vi, vj, a, dist)
    


# 2c

import matplotlib.pyplot as plt
def Uscaled(r):    
    res = 4 * (r**(-12) - r**(-6))
    
    return res


def getplots(sigma, dt):
    ri0 = np.array([sigma,0,0])
    rj0 = np.array([0,0,0])
    Machine = Solver(ri0, rj0, t0=0, t=5, dt=dt)
    (timepoints, ri, rj, vi, vj, a, dist) = Machine.EulerCromer()
    
    Ek = 0.5*np.array([np.linalg.norm(v)**2 for v in vi]) # Standard Kinetic energy formula
    
    Ep = 0.5*np.array([Uscaled(d) for d in dist]) # Since dist is traversed half by atom i and half by atom j
    
    Etot = Ek + Ep
    
    # ii
    plt.plot(timepoints, Ek, label="Kinetic energy")
    plt.plot(timepoints, Ep, label="Potential Energy")
    plt.plot(timepoints, Etot, label="Total energy")
    plt.grid()
    plt.legend()
    plt.title(f"Sigma {sigma} | dt = {dt}")
    plt.xlabel("Time")
    plt.ylabel("Energy")
    
    plt.savefig(fr"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/Energy_sig{round(sigma,2)}_dt{dt}.pdf")
    plt.clf()
    plt.close()
    

dts = [0.01, 0.0001]
sigmas = [1.5, 0.95]

for dt in dts:
    for sigma in sigmas:
        getplots(sigma=sigma, dt=dt)

