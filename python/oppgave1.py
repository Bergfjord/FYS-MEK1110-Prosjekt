# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 14:42:50 2022

@author: Anders sin PC
"""
# Oppgave 1
# i. 
#
import numpy as np
import matplotlib.pyplot as plt

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
#plt.savefig()



