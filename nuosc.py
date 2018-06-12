
"""
Anisha Kadri 2018
anishakadri1@gmail.com

A script for plotting probabilites for three- neutrino oscillation
using the PMNS model.

----------------------------------
References

 [1] Esteban, Ivan; Gonzalez-Garcia, M.C.; 
 Maltoni, Michele; Martinez Soler, Ivan; Schwetz, Thomas (2018).
 "Updated fit to three neutrino mixing: exploring the accelerator-reactor 
 complementarity".
 Retrieved 2018-05-01.

"""
import parameters as par
import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
import os

#parameter values

theta12 = 33.62
theta13 = 8.54
theta23 = 47.2
deltacp = 0 #4.154

pm12 = [0.78,-0.76]
pm23 = [1.9,-3.9]
pm13 = [0.15, -0.15]
pmcp = [.75,-0.54]

#arguments

print PMNS
#mass eigenvalues
#flavour eigenvalues
#create neutrino in muon flavour state

# P_mu_tau
# P_tau_mu
# P_e_mu
# P_mu_e
# m1
# m2
# m3

# psi_1 = np.exp((-1j)*)
# psi_2 = np.exp((-1j)*deltacp)
# psi_3 = np.exp((-1j)*deltacp)


# flavours = np.array([n])



