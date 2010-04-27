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

from scipy.integrate import odeint

class Compound:
	"""
	A chemical compound
	
	A,B,C are Antoine parameters in Bar and Kelvin
	"""
	def __init__(self, A, B, C, molar_density, MW, Hvap, Cp):
		self.A = A
		self.B = B
		self.C = C
		self.molar_density = molar_density
		self.molar_mass = MW
		self.enthalpy_of_vaporization = Hvap
		self.molar_heat_capacity = Cp
	
	def getPureComponentVaporPressure(self,Temperature):
		"""
		Use Antoine Equation to get saturated vapor pressure at Temperature.
		returns P in Pa
		
		Assumes A,B,C are for bar and K (as provided on webbook.nist.gov)
		
		P = 10^(A-B/(C+T))
		"""
		A = self.A
		B = self.B
		C = self.C
		
		# Antoine's Equation
		Pbar =  10**(A - B / (C + Temperature))
		return Pbar * 1E5 # to get Pa
		
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
			
		self.molar_masses = numpy.array([ c.molar_mass for c in self.components ])
		self.enthalpies_of_vaporization  = numpy.array([ c.enthalpy_of_vaporization for c in self.components ])
		self.molar_heat_capacities = numpy.array([c.molar_heat_capacity for c in self.components])
		self.molar_densities = numpy.array([c.molar_density for c in self.components])
		
		# it would make more sense to take these values as variables (like components, etc.)
		# but for now we just hard-code them here: 
		self.heat_transfer_coefficient = 40  # J / K / m^2 / s    (i.e. W/m2K)
		# 40 W/m2K = ballpark figure eyeballed from a picture in doi:10.1016/0017-9310(80)90153-2
		self.T_wall = 500 # Kelvin. 
	
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
		sigma = 0.1 # Hertz-Knudsen condensation/evaporation coefficient
		# see doi:10.1103/PhysRevE.59.429
		# http://www.mie.utoronto.ca/labs/tkl/publications/WardFangExpress.pdf
		vp = self.getVaporPressures(T)
		mw = self.molar_masses*0.001 / 6.02E23  # mass of molecule in kg
		kb = 1.38E-23 # m^2 kg s^-2 K^-1 !!
		Je = vp / numpy.sqrt(2*numpy.pi * mw * kb * T) 
		return Je/ 6.02E23
	
	def getHeatFlux(self, T):
		"""
		Get the heat flux (Q) into the layer from the hot wall (per unit area)
		when the layer is at temperature T.
		"""
		Q = self.heat_transfer_coefficient * (self.T_wall - T)
		return Q
		
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
		# first get the species amount changes (as an array)
		dNdt = -1 * self.getMolarFluxes(temperature)
		# then get the temperature change, and append it to the array
		dUdt = sum( dNdt * self.enthalpies_of_vaporization )
		dUdt = dUdt + self.getHeatFlux(temperature)
		dTdt = dUdt / sum( self.amounts * self.molar_heat_capacities )
		return numpy.append(dNdt,dTdt)

def main():
	pass

if __name__ == '__main__':
	main()

undecane = Compound(A=4.101, B=1572.477, C=-85.128, molar_density=4945, MW=156.0, Hvap=56.4, Cp=341.1)
c21      = Compound(A=5.921, B=3571.218, C=-19.953, molar_density=2729, MW=310.0, Hvap=142,  Cp=666.4)

print "Vapor pressure of pure undecane at 400K is ", undecane.getPureComponentVaporPressure(400)
print "and its Enthlapy of vaporization is",undecane.enthalpy_of_vaporization

list_of_compounds = [undecane, c21]

layer = Layer(list_of_compounds)

initial_amounts = numpy.array([1.0, 2.0])   # relative molar amounts - will be scaled to get correct initial_thickness
current_thickness = sum( initial_amounts / layer.molar_densities )

initial_thickness = 0.5e-3 # m
initial_amounts = initial_amounts * initial_thickness/current_thickness # moles (per unit area)

# set the amounts
layer.amounts = initial_amounts

print "The initial mole fractions are ", layer.getMoleFractions()

print "The vapor pressures at 400K are", layer.getVaporPressures(400)

print "The molar fluxes at 400K are", layer.getMolarFluxes(400)

# add the temperature on as the last variable
vector_of_variables = numpy.append(layer.amounts, 400)

print "The right hand side of the ODE at t=0 is ", layer.rightSideOfODE(vector_of_variables, 0)

timepoints = numpy.linspace(0,5,501) # 501 points spread linearly between 0 and 5 seconds
results = odeint(layer.rightSideOfODE, vector_of_variables, timepoints)
print "After %g seconds the vector of results is %s"%(timepoints[-1], results[-1])

# slice off the last column of results, which are temperatures
amounts = results[:,:-1]
temperatures = results[:,-1]



#plot the results
import pylab
pylab.figure(1)
pylab.plot(timepoints,amounts)
pylab.show()

pylab.figure(2)
pylab.plot(timepoints,temperatures)
pylab.show()


