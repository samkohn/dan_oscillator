#!/usr/bin/env python
#
# Three-flavor Neutrino oscillation, using GLoBES Interface
#
# Created by: dadwyer@lbl.gov 2013/05/01
# Last Modified: dadwyer@lbl.gov 2013/05/01
# 


# Useful imports
from Oscillator import NeutrinoState
from Units import *
from math import sqrt

# Load GLoBES Interface
from ROOT import gROOT
gROOT.ProcessLine(".x loadInterfaces.C")
from ROOT import GLoBES_Interface

class Oscillator_GLoBES(object):
    def __init__(self, sinTh12, sinTh13, sinTh23, 
                 dmSq21, dmSq32, deltaCP,
                 rho_e, energy):
        "Constructor"
        self._globesInterface = GLoBES_Interface()
        self._globesInterface.setOscillation(sinTh12, sinTh13, sinTh23, 
                                             dmSq21, dmSq32, deltaCP)
        self._rho_e = rho_e
        self._energy = energy

    def evolve(self, initialState, distance ):
        "Evolve neutrino quantum state, and return final state vector"
        # FIXME: GLoBES cannot handle coherent oscillation.
        #        Only valid for pure flavor initial state
        finalAmps = [complex(0), complex(0), complex(0)]
        initFlavor = 0
        activeFlavors = 0
        initAmps = initialState.amplitudes()
        for flavorIdx in range(3):
            if initAmps[flavorIdx] != 0:
                activeFlavors += 1
                initFlavor = flavorIdx
        if activeFlavors != 1:
            print "Error: Oscillator_GLoBES can't do real oscillation.  Use pure flavor states only."
            return NeutrinoState(complex(1),
                                 complex(0),
                                 complex(0),
                                 initialState.isMatter())
        # Calculate Oscillation
        probabilities = []
        for finalFlavor in range(3):
            prob = self._globesInterface.probability(initFlavor+1,
                                                     finalFlavor+1,
                                                     initialState.isMatter(),
                                                     self._energy/GeV,
                                                     distance/km,
                                                     self._rho_e/(1/(m**3)))
            probabilities.append( prob )
        finalState = NeutrinoState(complex(sqrt(probabilities[0])),
                                   complex(sqrt(probabilities[1])),
                                   complex(sqrt(probabilities[2])),
                                   initialState.isMatter())
        return finalState
