#!/usr/bin/env python
#
# Neutrino Oscillation Physical Parameters
#
# Created by: dadwyer@lbl.gov 2013/03/06
# Last Modified: dadwyer@lbl.gov 2013/03/15
# 

from Units import *

# PDG2012 best estimates of neutrino oscillation parameters
neutrinoParams_PDG2012 = {'sinTheta12':0.553,
                          'sinTheta13':0.151,
                          'sinTheta23':0.648,
                          'deltaMSq21':7.58e-5 * eV**2,
                          'deltaMSq32':2.35e-3 * eV**2,
                          'deltaCP':0,
                          }

# Neutrino oscillation parameters: arXiv:1205.7071v5
neutrinoParams_arXiv12057071 = {'sinTheta12':0.5585696,
                                'sinTheta13':0.1581139,
                                'sinTheta23':0.6480741,
                                'deltaMSq21':7.60e-5 * eV**2,
                                'deltaMSq32':2.35e-3 * eV**2,
                                'deltaCP':0,
                                }

# Latest estimates of neutrino oscillation parameters
neutrinoParams_best = {'sinTheta12':0.553,
                       'sinTheta13':0.151,
                       'sinTheta23':0.648,
                       'deltaMSq21':7.5e-5 * eV**2,
                       'deltaMSq32':2.5e-3 * eV**2,
                       'deltaCP':0,
                       }
