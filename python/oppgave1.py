# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 14:42:50 2022

@author: Anders sin PC
"""
# Oppgave 1
# i. 
#
import matplotlib.pyplot as plt
import numpy as np


def U(r):
    eps = 1
    sigma = 1
    
    res = 4*eps * ((sigma/r)**12 - (sigma/r)**6)
    
    return res

n = 101
r = np.linspace(0.9, 3, n)
potential = U(r)


plt.plot(r, potential, label="U(r)")
plt.title("Potential Graph U(r)")
plt.xlabel("Distance")
plt.ylabel("Potential")
plt.legend()
plt.show()
#plt.savefig(r"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/Potential.pdf")


# ii

# ic 



# opgave 3
plt.clf()
plt.close()

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
plt.show()
#plt.savefig(r"C:/Users/Anders sin PC/Documents/UIO/V22/FYS-MEK1110/Project_molekylær_dynamikk/Potential_shift.pdf")




