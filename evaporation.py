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
import csv

class Compound:
	"""
	A chemical compound
	
	A,B,C are Antoine parameters in mmHg and Kelvin
	"""
	def __init__(self, name, Antoine_params, mass_density, MW, Hvap, Cp):
		self.name = name
		self.Antoine_params = Antoine_params # a tuple or list: [A,B,C]
		self.mass_density = mass_density # kg/m^3
		self.molar_mass = MW # g/mol
		self.enthalpy_of_vaporization = Hvap
		self.molar_heat_capacity = Cp
		# derived properties
		self.molar_density = mass_density / (0.001*MW) # kg/m^3 / kg/mol = mol/m^3
		
	def __repr__(self):
		"""This function returns how the compound will look in the console."""
		return "<Compound %r>"%self.name
		
	def getPureComponentVaporPressure(self,Temperature):
		"""
		Use Antoine Equation to get saturated vapor pressure at Temperature.
		returns P in Pa
		
		A,B,C are for mmHg and K (as provided in Epstein 2009 Supp Info) and must subtract 273.15 from Temperature to convert to degrees C as required by Antoine Eqn.
		
		P = 10^(A-B/(C+T))
		"""
		A = self.Antoine_params[0]
		B = self.Antoine_params[1]
		C = self.Antoine_params[2]
		
		# Antoine's Equation
		PmmHg =  10**(A - B / (C + Temperature - 273.15))
		return PmmHg * 133.322 # to get Pa
		
class CompoundsDatabase(dict):
	"""
	A collection of information about compounds.
	Behaves like a dict (dictionary)
	"""
	
	def __init__(self,file_path):
		"""
		Reads in an excel-formatted CSV file with column headings in the first row.
		"""
		data_reader = csv.DictReader(file(file_path,'rU'))
		for row in data_reader:
			# we have to turn the strings into floating point numbers.
			c = Compound( name = row['Name'],
			              Antoine_params = [float(row['Antoine A']),float(row['Antoine B']),float(row['Antoine C'])],
			              mass_density = float(row['Mass Density']),
			              MW = float(row['Molecular Weight']),
			              Hvap = float(row['Enthalpy of Vaporization']),
			              Cp = float(row['Molar Heat Capacity']) )
			# place it in the dictionary
			#print "Have just read in ",c
			self[c.name] = c
		
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
		self.T_wall = 873.15 # Kelvin. 
	
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

compounds = CompoundsDatabase('compounds.csv')
undecane=compounds['undecane']

print "Vapor pressure of pure undecane at 300K is ", undecane.getPureComponentVaporPressure(300)
print "and its Enthlapy of vaporization is",undecane.enthalpy_of_vaporization

# this would be the full list in a random order:
list_of_compounds = compounds.values()
# this will be just the eight we specify, in the specified order:
list_of_compounds = [ compounds[name] for name in ['undecane', 'c21', 'acetophenone', 'diethyl phthalate', 'diethylene glycol', 'epsilon-caprolactone', 'gamma-butyrolactone', 'glycerol', 'triethylene glycol'] ]

layer = Layer(list_of_compounds)

initial_amounts = numpy.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])   # relative molar amounts - will be scaled to get correct initial_thickness
current_thickness = sum( initial_amounts / layer.molar_densities )

initial_thickness = 1.0e-5 # m
initial_amounts = initial_amounts * initial_thickness/current_thickness # moles (per unit area)

# set the amounts
layer.amounts = initial_amounts

print "The initial mole fractions are ", layer.getMoleFractions()

print "The vapor pressures at 300K are", layer.getVaporPressures(300)

print "The molar fluxes at 300K are", layer.getMolarFluxes(300)

# add the temperature on as the last variable
vector_of_variables = numpy.append(layer.amounts, 300)

print "The right hand side of the ODE at t=0 is ", layer.rightSideOfODE(vector_of_variables, 0)

timepoints = numpy.linspace(0,0.2,201) # 201 points spread linearly between 0 and 0.2 seconds
results = odeint(layer.rightSideOfODE, vector_of_variables, timepoints)
print "After %g seconds the vector of results is %s"%(timepoints[-1], results[-1])

# slice off the last column of results, which are temperatures
amounts = results[:,:-1]
temperatures = results[:,-1]

thicknesses = amounts / layer.molar_densities
thicknesses = thicknesses.sum(axis=1) # add the components up, at each time point


#plot the results
import pylab
pylab.figure(1)
pylab.plot(timepoints,amounts)
pylab.xlabel("time (s)")
pylab.ylabel("amount (mol/m2)")
pylab.show()

pylab.figure(2)
pylab.plot(timepoints,temperatures)
pylab.xlabel("time (s)")
pylab.ylabel("temperature (K)")
pylab.show()

pylab.figure(3)
pylab.plot(timepoints,thicknesses)
pylab.xlabel("time (s)")
pylab.ylabel("thickness (m)")
pylab.show()

