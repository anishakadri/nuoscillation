import numpy as np
#oscillation parameters
theta12 = 0.58677969
theta13 = 0.1490511
theta23 = 0.8237954
deltacp = 4.084


#uncertainty intervals
pm12 = [0.0136136,-0.0132645]
pm23 = [0.0331613,-0.0680678]
pm13 = [0.002617994, -0.002617994]
pmcp = [.75,-0.54]

#pmns terms
c12 = np.cos(theta12)
s12 = np.sin(theta12)
c13 = np.cos(theta13)
s13 = np.sin(theta13)
c23 = np.cos(theta23)
s23 = np.cos(theta23)
cp = np.exp((1j)*deltacp)
cpm = np.exp((-1j)*deltacp)

time = 1.51975683891e15 #1s in eV^-1
length = 5076142.13198 #1m in eV^-1

def natural_t(seconds):
	return(seconds*time)

def natural_l(meters):
	return(meters*length)

#PMNS Matrix

PMNS = np.array([[c12*c13,s12*c13,s13*cpm],[(-s12*c23)-(c12*s23*s13*cp), (c12*c23)-(s12*s23*s13*cp), s23*c13],[(s12*s23)-(c12*c23*s13*cp), (-c12*s23)-(s12*c23*s13*cp), c23*c13]])

#neutrino masses: (since only values of mass-squared 
# differences are relevant, can choose m1=0.12eV and 
#compute the rest accordingly)

mass = np.array([0.12 , np.sqrt(7.37e-5 + (0.12**2)), np.sqrt(2.56e-3 + (0.12**2))]) #coverted to MeV

print np.linalg.norm(PMNS)






