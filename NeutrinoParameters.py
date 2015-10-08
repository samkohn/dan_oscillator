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

# NuFit 2.0 (2014)
nufit_NO = {'best': {'sinTheta12': 0.55136195,
                 'sinTheta13': 0.147648231,
                 'sinTheta23': 0.672309453,
                 'deltaMSq21': 7.5e-5 * eV**2,
                 'deltaMSq32': 2.382e-3 * eV**2,
                 'deltaCP': 0,
                 },
             '+3sigma': {'sinTheta12': 0.586515132,
                         'sinTheta13': 0.158113883,
                         'sinTheta23': 0.801872808,
                         'deltaMSq21': 8.09e-5 * eV**2,
                         'deltaMSq32': 2.5368e-3 * eV**2,
                         'deltaCP': 0,
                         },
             '-3sigma': {'sinTheta12': 0.519615242,
                         'sinTheta13': 0.136381817,
                         'sinTheta23': 0.618061486,
                         'deltaMSq21': 7.02e-5 * eV**2,
                         'deltaMSq32': 2.2361e-3 * eV**2,
                         'deltaCP': 0,
                         },
             }

nufit_IO = {'best':
                {'sinTheta12': 0.55136195,
                 'sinTheta13': 0.147986486,
                 'sinTheta23': 0.760920495,
                 'deltaMSq21': 7.5e-5 * eV**2,
                 'deltaMSq32': -2.449e-3 * eV**2,
                 'deltaCP':0,
                 },
            '+3sigma':
                {'sinTheta12': 0.586515132,
                 'sinTheta13': 0.158429795,
                 'sinTheta23': 0.802496106,
                 'deltaMSq21': 8.09e-5 * eV**2,
                 'deltaMSq32': -2.307e-3 * eV**2,
                 'deltaCP': 0,
                 },
            '-3sigma':
                {'sinTheta12': 0.519615242,
                 'sinTheta13': 0.137113092,
                 'sinTheta23': 0.623698645,
                 'deltaMSq21': 8.09e-5 * eV**2,
                 'deltaMSq32', -2.590e-3 * eV**2,
                 'deltaCP': 0,
                 },
            }
