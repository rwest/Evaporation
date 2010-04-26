#!/usr/bin/env python
# encoding: utf-8
"""
evaporation.py

Created by Richard West on 2010-04-26.
Copyright (c) 2010 MIT. All rights reserved.
"""

import sys
import os

import numpy

class Compound:
	"""A chemical compound"""
	def __init__(self, A, B, C, molar_density, MW):
		self.A = A
		self.B = B
		self.C = C
		self.molar_density = molar_density
		self.molecular_weight = MW
	
	def getPureComponentVaporPressure(self,Temperature):
		"""Use Antoine Equation to get saturated vapor pressure at Temperature.
		
		Note the units in Antoine Eqn is mmHg and C
		P = 10^(A-B/(C+T))
		"""
		A = self.A
		B = self.B
		C = self.C
		# transfer from mmHg to Pa
		# A = A0 + log10(101325/760)
		A = A + 2.124903
		# transfer from C to K
		C = C-273.15
		# Antoine's Equation
		return 10**(A-B / (C + Temperature))
		
class Layer:
	"""
	A uniform layer of several compounds, all mixed, at the same temperature.
	Unit surface area.
	"""
	def __init__(self, components, initial_amounts = None):
		self.components = components
		self.number_of_components = len(components)
		
		if initial_amounts is not None: # number of moles (per unit surface area)
			self.amounts = initial_amounts
		else:
			self.amounts = numpy.zeros(self.number_of_components)
			
		self.molecular_weights = numpy.array([ c.molecular_weight for c in self.components ])
	
	def getMoleFractions(self):
		"""return an array of the mole fractions"""
		return self.amounts / self.amounts.sum()
		
	def getPureVaporPressures(self,T):
		"""return an array of the pure component vapor pressures"""
		answer = list()
		for c in self.components:
			answer.append( c.getPureComponentVaporPressure(T) )
		return numpy.array(answer)
	
	def getVaporPressures(self,T):
		"""
		return an array of the vapor pressures, given the current mole fractions
		Assumes Henry's law.
		"""
		return self.getPureVaporPressures(T) * self.getMoleFractions()
		
	def getMolarFluxes(self,T):
		"""
		Return the molar flux in moles per second (per unit area) at a given T
		"""
		vp = self.getVaporPressures(T)
		mw = self.molecular_weights
		kb = 1.38E-23 # m^2 kg s^-2 K^-1 !!
		Je = vp / numpy.sqrt(2*numpy.pi * kb * T)
		return Je
		
	def rightSideOfODE(self, Y, t):
		"""
		The right hand side of the ODE that will be solved.
		
		Y is the vector of values at time t.
		Temperature is the last element of Y, the others are amounts.
		
		returns dY/dt
		"""
		temperature = Y[-1]
		amounts = Y[:-1]
		self.amounts = amounts
		# first get the species amounts
		dYdt = -1 * self.getMolarFluxes(temperature)
		# then get the temperature change
		dTdt = 
		dYdt.append(dTdt)

def main():
	pass

if __name__ == '__main__':
	main()


undecane = Compound(A=6.9722, B=1569.57, C=187.7, molar_density=4945, MW = 156.0)
c21 = Compound(A=7.0842, B=2054, C=120.1, molar_density =2729, MW=310.0)

print "Vapor pressure of pure undecane at 400K is ", undecane.getPureComponentVaporPressure(400)

list_of_compounds = [undecane, c21]

initial_amounts = numpy.array([1.0, 2.0]) # moles (per unit area)

layer = Layer(list_of_compounds, initial_amounts)

print "The initial mole fractions are ", layer.getMoleFractions()

print "The vapor pressures at 400K are", layer.getVaporPressures(400)

print "The molar fluxes at 400K are", layer.getMolarFluxes(400)

