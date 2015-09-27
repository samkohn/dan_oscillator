#!/usr/bin/env python
#
# Test script for neutrino oscillation tools
#
# Created by: dadwyer@lbl.gov 2013/03/12
# Last Modified: dadwyer@lbl.gov 2013/03/15
#

from Oscillator import Oscillator, NeutrinoState
from NeutrinoParameters import neutrinoParams_PDG2012
from Units import *

if "__main__" == __name__:
    # Run tests

    # Daya Bay Oscillation
    print ""
    print "Daya Bay:"
    oscillator_vacuum = Oscillator(neutrinoParams_PDG2012['sinTheta12'],
                                   neutrinoParams_PDG2012['sinTheta13'],
                                   neutrinoParams_PDG2012['sinTheta23'],
                                   neutrinoParams_PDG2012['deltaMSq21'],
                                   neutrinoParams_PDG2012['deltaMSq32'],
                                   neutrinoParams_PDG2012['deltaCP'],
                                   0, # Electron density
                                   4 * MeV) # Neutrino energy
    antiNuE = NeutrinoState(1.0, 0, 0)
    finalState = oscillator_vacuum.evolve(antiNuE, 1.6*km)
    amplitudes = finalState.amplitudes()
    probabilities = finalState.probabilities()
    print " Evolve: 1.6 km:"
    print "  Oscillated amplitudes:"
    print "    antiNuE   = ",amplitudes[0]
    print "    antiNuMu  = ",amplitudes[1]
    print "    antiNuTau = ",amplitudes[2]
    print "  Oscillated probabilities:"
    print "    antiNuE   = ",probabilities[0]
    print "    antiNuMu  = ",probabilities[1]
    print "    antiNuTau = ",probabilities[2]

    midState = oscillator_vacuum.evolve(antiNuE, 0.8*km)
    finalState = oscillator_vacuum.evolve(midState, 0.8*km)
    amplitudes = finalState.amplitudes()
    probabilities = finalState.probabilities()
    print " Evolve: 2 x 0.8 km:"
    print "  Oscillated amplitudes:"
    print "    antiNuE   = ",amplitudes[0]
    print "    antiNuMu  = ",amplitudes[1]
    print "    antiNuTau = ",amplitudes[2]
    print "  Oscillated probabilities:"
    print "    antiNuE   = ",probabilities[0]
    print "    antiNuMu  = ",probabilities[1]
    print "    antiNuTau = ",probabilities[2]


    # KamLAND Oscillation
    print ""
    print "KamLAND:"
    finalState = oscillator_vacuum.evolve(antiNuE, 180.0*km)
    amplitudes = finalState.amplitudes()
    probabilities = finalState.probabilities()
    print "  Oscillated amplitudes:"
    print "    antiNuE   = ",amplitudes[0]
    print "    antiNuMu  = ",amplitudes[1]
    print "    antiNuTau = ",amplitudes[2]
    print "  Oscillated probabilities:"
    print "    antiNuE   = ",probabilities[0]
    print "    antiNuMu  = ",probabilities[1]
    print "    antiNuTau = ",probabilities[2]


    # Earth density
    print ""
    print "Earth mantle density:"
    electronDensity_mantle = (4.5 * g / cm**3) * (0.5 / amu)
    oscillator_mantle = Oscillator(neutrinoParams_PDG2012['sinTheta12'],
                                   neutrinoParams_PDG2012['sinTheta13'],
                                   neutrinoParams_PDG2012['sinTheta23'],
                                   neutrinoParams_PDG2012['deltaMSq21'],
                                   neutrinoParams_PDG2012['deltaMSq32'],
                                   neutrinoParams_PDG2012['deltaCP'],
                                   electronDensity_mantle, # Electron density
                                   5 * GeV) # Neutrino energy
    nuMu = NeutrinoState(0, 1.0, 0)
    print " Oscillate 10,000 km:"
    finalState = oscillator_mantle.evolve(nuMu, 10000.*km)
    amplitudes = finalState.amplitudes()
    probabilities = finalState.probabilities()
    print "  Oscillated amplitudes:"
    print "    nuE   = ",amplitudes[0]
    print "    nuMu  = ",amplitudes[1]
    print "    nuTau = ",amplitudes[2]
    print "  Oscillated probabilities:"
    print "    nuE   = ",probabilities[0]
    print "    nuMu  = ",probabilities[1]
    print "    nuTau = ",probabilities[2]

    print " Oscillate 2 x 5,000 km:"
    midState = oscillator_mantle.evolve(nuMu, 5000.*km)
    finalState = oscillator_mantle.evolve(midState, 5000.*km)
    amplitudes = finalState.amplitudes()
    probabilities = finalState.probabilities()
    print "  Oscillated amplitudes:"
    print "    nuE   = ",amplitudes[0]
    print "    nuMu  = ",amplitudes[1]
    print "    nuTau = ",amplitudes[2]
    print "  Oscillated probabilities:"
    print "    nuE   = ",probabilities[0]
    print "    nuMu  = ",probabilities[1]
    print "    nuTau = ",probabilities[2]
