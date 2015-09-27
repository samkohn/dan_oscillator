#!/usr/bin/env python
#
# Test script: compare 3-flavor oscillation calculation with GLoBES
#
# Created by: dadwyer@lbl.gov 2013/05/01
# Last Modified: dadwyer@lbl.gov 2013/05/01
#

from Plotter import Plotter, drawGraphs
from Oscillator import Oscillator, NeutrinoState
from Oscillator_GLoBES import Oscillator_GLoBES
from NeutrinoParameters import neutrinoParams_PDG2012 as neutrinoParams
from math import sin, cos, asin, acos, sqrt, atan, tan, atan2, pi
from Units import *

#deltaCP = neutrinoParams['deltaCP']
#deltaCP = 0
#deltaCP = pi/2
deltaCP = 3*pi/2

if "__main__" == __name__:
    # Run tests
    plotter = Plotter()
    energyAtmo = 6 * GeV
    Ne_mantle = (4.5*g/cm**3)*0.5/amu
    oscillator_mantleMSW = Oscillator(neutrinoParams['sinTheta12'],
                                      neutrinoParams['sinTheta13'],
                                      neutrinoParams['sinTheta23'],
                                      neutrinoParams['deltaMSq21'],
                                      neutrinoParams['deltaMSq32'],
                                      #neutrinoParams['deltaCP'],
                                      deltaCP,
                                      Ne_mantle, # Electron density
                                      energyAtmo) # Neutrino energy
    oscillator_mantleMSW_Alt = Oscillator_GLoBES(neutrinoParams['sinTheta12'],
                                                 neutrinoParams['sinTheta13'],
                                                 neutrinoParams['sinTheta23'],
                                                 neutrinoParams['deltaMSq21'],
                                                 neutrinoParams['deltaMSq32'],
                                                 #neutrinoParams['deltaCP'],
                                                 deltaCP,
                                                 Ne_mantle, # Electron density
                                                 energyAtmo) # Neutrino energy
    nuMu = NeutrinoState(0, 1.0, 0, True)

    
    # Earth Mantle Oscillation
    oscNumericMSW = []
    oscNumericMSW_Alt = []
    baselines = []
    nSamples = 1000
    for i in range(nSamples):
        l = (15000./nSamples)*i*km
        baselines.append( l )
        finalStateMSW = oscillator_mantleMSW.evolve(nuMu, l)
        oscNumericMSW.append( finalStateMSW.probabilities()[0] )
        finalStateMSW_Alt = oscillator_mantleMSW_Alt.evolve(nuMu, l)
        oscNumericMSW_Alt.append( finalStateMSW_Alt.probabilities()[0] )
    
    g1 = plotter.makeGraph(baselines, oscNumericMSW)
    g2 = plotter.makeGraph(baselines, oscNumericMSW_Alt)

    g1.SetTitle("Atmospheric Neutrino Oscillation (E_{#nu}=6 GeV)")
    g1.GetXaxis().SetTitle("distance [m]")
    if nuMu.isMatter():
        g1.GetYaxis().SetTitle("P(#nu_{#mu} #rightarrow #nu_{e})")
    else:
        g1.GetYaxis().SetTitle("P(#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e})")
    for graph in [g1,g2]: graph.SetLineStyle(2)
    canvas = plotter.makeCanvas()
    labels = ["Numeric (3-#nu, MSW)", "GLoBES, mod. (3-#nu, MSW)"]
    drawGraphs([g1,g2],canvas,None,labels)
    g1.SetMarkerColor(4)
    g1.SetLineColor(4)
    g2.SetMarkerColor(2)
