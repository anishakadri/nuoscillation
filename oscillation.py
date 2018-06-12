import parameters as par
import numpy as np 
import scipy.optimize as spo

#N.B. - mass states only defined for non-sterline states

im = -1j

electron = 0
muon = 1
tau = 2

m1 = 0
m2 = 1
m3 = 2

def prob(mass, flavour):
		return par.PMNS[flavour,mass]

class MassState(object):
	"""
	A class definining mass eigenstates;
    Where it can be initialised with the following variables
    	MassState(num, time=0., energy=0.6e9, amp=1.):
    		num = index of mass state (0,1,2) or predefined variable (m1, m2, m3)
    		time = propagation time of the state, in units of eV^-1
    		energy = energy of state in units of eV
    		amp = amplitude (sqrt(probability)) of the state
    """
	def __init__(self, m_num, time=0., energy=0.6e9, amp=1.):
		super(MassState, self).__init__()
		if m_num !=0 and m_num!=1 and m_num!=2:
			raise Exception("no defined mass state, choose m=1,2,3") 
		self.index = m_num
		self.mass = par.mass[self.index]

		if amp >1.:
			raise Exception("amplitude cannot be greater than one")
		self.decomp = amp*(np.conj(par.PMNS[:,self.index])) #spectral decomposition into flavour states
		self.time = time
		self.energy = energy
		self.amp = np.absolute(self.decomp)


	def __repr__(self):
		return "%s(m %r, time=%r,energy=%r, amplitude=%r)" % ("MassState", 
			self.index+1, self.time, self.energy, self.amp)

	def probability(self):
		k = self.index
		return np.absolute(self.decomp)		

	def propagate(self, L):
		self.time = self.time + L
		E = self.energy
		if self.energy ==0.:
			raise Exception("energy is 0") 

		propagator = np.exp((im*(self.mass**2)*L)/(2*E))
		print self.amp
		self.decomp = np.array([propagator*self.decomp[0], propagator*self.decomp[1], propagator*self.decomp[2]])
		self.amp = np.absolute(self.decomp)
		print self.amp
		return self

#---------------------------#
class FlavourState(object):
	"""
	A class definining flavour eigenstates;
    Where it can be initialised with the following variables
    	MassState(num, time=0., energy=0.6e9, amp=1.):
    		num = index of flavour state (0,1,2) or predefined variable (electron,muon, tau)
    		time = propagation time of the state, in units of eV^-1
    		energy = energy of state in units of eV
    		amp = amplitude (sqrt(probability)) of the state
    """

	def __init__(self, flavour, time=0., energy=0.6e9, amp=1.):
		super(FlavourState, self).__init__()
		self.index = flavour
		if self.index == 0:
			self.flavour = "electron"
		elif self.index == 1:
			self.flavour = "muon"
		elif self.index == 2:
			self.flavour = "tau"
		else:
			raise Exception("no defined flavour state, choose muon, electron or tau")

		if amp >1.:
			raise Exception("amplitude cannot be greater than one")
		self.decomp = amp*par.PMNS[self.index,:]
		self.time = 0
		self.energy = energy
		self.amp = np.sqrt(np.absolute(self.decomp))


	def __repr__(self):
		return "%s(%s, time=%r, energy=%r, amplitude=%r)" % ("FlavourState", self.flavour, self.time, self.energy, self.amp)

#---------------------------#

class Neutrino(object):
	"""
	N.B- The neutrino must be initialised in a definite flavour state

	"""
	def __init__(self, flavour, E = 0.6e9):
		super(Neutrino, self).__init__()
		self.time = 0.
		self.energy = E
		self.initial = FlavourState(flavour, time=0., energy=self.energy, amp=1.)

		self.flavour_spectrum = np.array([FlavourState(electron, time=0., energy=self.energy, amp=0.),FlavourState(muon, time=0., energy=self.energy, amp=0.),FlavourState(tau, time=0., energy=self.energy, amp=0.)])
		
		self.flavour_spectrum[flavour] = self.initial

		self.nu_e = self.flavour_spectrum[electron]
		self.nu_mu = self.flavour_spectrum[muon]
		self.nu_tau = self.flavour_spectrum[tau]

		self.mass1 = MassState(m1, time=0., energy=self.energy, amp=self.initial.decomp[m1])
		self.mass2 = MassState(m2, time=0., energy=self.energy, amp=self.initial.decomp[m2])
		self.mass3 = MassState(m3, time=0., energy=self.energy, amp=self.initial.decomp[m3])

		self.mass_spectrum = np.array([self.mass1,self.mass2,self.mass3])

	def __repr__(self):
		return "%s(time=%r,energy=%r, [e: %r, mu: %r, tau: %r] )" % ("Neutrino", 
			self.time, self.energy, self.nu_e.amp, self.nu_mu.amp, self.nu_tau.amp)

	def flav_amps(self):
		return np.array([self.nu_e.amp, self.nu_mu.amp, self.nu_tau.amp])

	def mass_amps(self):
		return np.array([self.mass1.amp, self.mass2.amp, self.mass3.amp])

	def propagate(self, L):
	 	self.time = self.time + L

	 	#print self.mass3.amp

	  	self.mass1 = self.mass1.propagate(L)
	  	self.mass2 = self.mass2.propagate(L)
	  	self.mass3 = self.mass3.propagate(L)

	  	#print self.mass3.amp

	  	self.mass_spectrum = np.array([self.mass1,self.mass2,self.mass3])

	  	for i in np.array([0,1,2]):
	  		amp = self.mass1.decomp[i]+ self.mass2.decomp[i] + self.mass3.decomp[i]
	  		self.flavour_spectrum[i] = FlavourState(i, self.time, self.energy, amp)

	  	self.nu_e = self.flavour_spectrum[electron]
	  	self.nu_mu = self.flavour_spectrum[muon]
	  	self.nu_tau = self.flavour_spectrum[tau]

	  	return self

new = Neutrino(muon)
new = new.propagate(295e3*par.length)


	# def observe(self):

	# 	self = 

		









		



		