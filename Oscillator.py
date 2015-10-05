#!/usr/bin/env python
#
# Three-flavor Neutrino oscillation tools
#
# Created by: dadwyer@lbl.gov 2013/03/12
# Last Modified: dadwyer@lbl.gov 2013/03/14
# 


# Useful imports
from sys import float_info
from cmath import exp, sqrt
from math import sin, cos, atan, atan2
from math import sqrt as sqrt_real
from LinearAlgebra import array, dot
from Units import *

class NeutrinoState(object):
    _normTolerance = 1.0e-4 # Tolerance to maintain state normalization
    def __init__(self, amplitude_e, amplitude_mu, amplitude_tau, isMatter=True):
        "Constructor"
        self._amplitude_e   = complex(amplitude_e)
        self._amplitude_mu  = complex(amplitude_mu)
        self._amplitude_tau = complex(amplitude_tau)
        self._isMatter = isMatter
        # Ensure proper state normalization
        #assert abs(self.normalization() - 1.0) < self._normTolerance, "State is not normalized, N = %f" % (self.normalization())
        if abs(self.normalization() - 1.0) > NeutrinoState._normTolerance:
            print "Warning: State is not normalized, N = %f" % (self.normalization())
        return

    def amplitudes(self):
        "Return the complex state vector (Psi_e, Psi_mu, Psi_tau)"
        return [self._amplitude_e, self._amplitude_mu, self._amplitude_tau]

    def probabilities(self):
        "Return the probabilities for detection of each flavor"
        return [(self._amplitude_e*self._amplitude_e.conjugate()).real,
                (self._amplitude_mu*self._amplitude_mu.conjugate()).real,
                (self._amplitude_tau*self._amplitude_tau.conjugate()).real]
    
    def normalization(self):
        "Return the state normalization"
        return sum(self.probabilities())

    def isMatter(self):
        "Return true for matter, false for antimatter"
        return self._isMatter

class OscillatorBase(object):
    def __init__(self, sinTh12, sinTh13, sinTh23,
                 dmSq21, dmSq32, deltaCP,
                 rho_e, energy):
        pass

    def evolve(self, initialState, distance):
        pass

    @classmethod
    def fromParameterSet(cls, parameterSet, rho_e, energy):
        "Get an oscillator using the parameters saved in parameterSet."
        oscillator = cls(parameterSet['sinTheta12'],
                parameterSet['sinTheta13'], parameterSet['sinTheta23'],
                parameterSet['deltaMSq21'], parameterSet['deltaMSq32'],
                parameterSet['deltaCP'], rho_e, energy)
        return oscillator

class Oscillator(OscillatorBase):
    def __init__(self, sinTh12, sinTh13, sinTh23, 
                 dmSq21, dmSq32, deltaCP,
                 rho_e, energy):
        "Constructor"
        self._sinTh12 = sinTh12
        self._sinTh13 = sinTh13
        self._sinTh23 = sinTh23
        self._dmSq21 = dmSq21
        self._dmSq32 = dmSq32
        self._deltaCP = deltaCP
        self._rho_e = rho_e
        self._energy = energy
        self._Ufm = None
        self._Umf = None
        self._updatePMNSMatrix()
        [lambdaEig, amplFactor] = self._updateAmplitudeFactors(1)
        self._lambdaMatter = lambdaEig
        self._amplFactorMatter = amplFactor
        [lambdaEigAnti, amplFactorAnti] = self._updateAmplitudeFactors(-1)
        self._lambdaAntiMatter = lambdaEigAnti
        self._amplFactorAntiMatter = amplFactorAnti

    def _updatePMNSMatrix(self):
        "Update PMNS neutrino mixing matrix"
        s12 = self._sinTh12
        s13 = self._sinTh13
        s23 = self._sinTh23
        c12 = sqrt( 1 - s12**2 )
        c13 = sqrt( 1 - s13**2 )
        c23 = sqrt( 1 - s23**2 )
        expCP = exp( 1J * self._deltaCP )
        expMCP = exp( -1J * self._deltaCP )
        U_e1 = c12*c13
        U_e2 = s12*c13
        U_e3 = s13 * expMCP
        U_m1 = -s12*c23 - c12*s23*s13 * expCP
        U_m2 =  c12*c23 - s12*s23*s13 * expCP
        U_m3 =  s23*c13
        U_t1 =  s12*s23 - c12*c23*s13 * expCP
        U_t2 = -c12*s23 - s12*c23*s13 * expCP
        U_t3 =  c23*c13
        Ufm = [[U_e1, U_e2, U_e3],
               [U_m1, U_m2, U_m3],
               [U_t1, U_t2, U_t3]]
        self._Ufm = array( Ufm )
        self._Umf = self._Ufm.conj().transpose()
        return

    def _updateAmplitudeFactors(self, matterSign):
        "Pre-calculate amplitude factors needed for oscillation"
        amplFactor = array( [[[complex(0)]*3]*3]*3 )
        # Consider only relative energy differences in relativistic limit
        E1 = 0.0
        E2 = self._dmSq21 / ( 2 * self._energy )
        E3 = (self._dmSq32 + self._dmSq21) / ( 2 * self._energy )
        Ufm = None
        Umf = None
        if matterSign > 0:
            # Neutrino mixing matrix
            Ufm = self._Ufm
            Umf = self._Umf
        else:
            # Antineutrino mixing matrix
            Ufm = self._Ufm.conj()
            Umf = self._Umf.conj()
        if self._rho_e == 0:
            # Vacuum oscillations.  Simple calculation
            lambdaEig = [E1, E2, E3]
            for mi in range(3): # mass eigenstate
                for ffl in range(3): # Final flavor
                    for ifl in range(3): # Initial flavor
                        amplFactor[ifl,ffl,mi] = (Umf[mi,ifl]
                                                  *Ufm[ffl,mi])
            return [lambdaEig, amplFactor]
        # Calculate matter-enhanced oscillation factors
        identity3x3 = array( [[1, 0, 0],
                              [0, 1, 0],
                              [0, 0, 1]] )
        Ae = matterSign * sqrt_real(2.) * Gf * self._rho_e
        energyScale = (E1 + E2 + E3 + Ae) / 4.
        E1 /= energyScale
        E2 /= energyScale
        E3 /= energyScale
        Ae /= energyScale
        Hm_vac = array( [[E1, 0, 0],
                         [0, E2, 0],
                         [0, 0, E3]] )
        Vf = array(  [[Ae, 0, 0],
                      [0, 0, 0],
                      [0, 0, 0]] )
        Vm = dot(Umf, dot(Vf, Ufm) )  
        Hm = Hm_vac + Vm
        trHm = E1 + E2 + E3 + Ae
        Tm = Hm - ((trHm/3.)*identity3x3)
        lambdaEig = self._eigenvalues( Tm ) 
        Tf = dot( Ufm, dot(Tm, Umf) )
        TfSq = dot( Tf, Tf )
        c1 = (Tm[0,0]*Tm[1,1] + Tm[0,0]*Tm[2,2] + Tm[1,1]*Tm[2,2]
              - (Tm[0,1]*Tm[1,0] + Tm[0,2]*Tm[2,0] + Tm[1,2]*Tm[2,1]) )
        for mi in range(3): # matter-enhanced energy eigenstate
            lambda_a = lambdaEig[mi]
            denom = 3*lambda_a*lambda_a + c1
            amplFactorM = ((lambda_a*lambda_a + c1) * identity3x3 
                           + lambda_a * Tf
                           + TfSq) * (1 / denom)
            for ffl in range(3): # Final flavor
                for ifl in range(3): # Initial flavor
                    amplFactor[ifl,ffl,mi] = amplFactorM[ffl,ifl]
        # Unscale matter-enhanced energy eigenvalues
        lambdaEig = [eigval*energyScale for eigval in lambdaEig]
        return [lambdaEig, amplFactor]

    def _eigenvalues(self, Tm):
        "Find energy eigenvalues of 3x3 complex Hamiltonian"
        # Calculate eigenvalues
        c0 = -(self._determinant(Tm))
        c1 = (Tm[0,0]*Tm[1,1] + Tm[0,0]*Tm[2,2] + Tm[1,1]*Tm[2,2]
              - (Tm[0,1]*Tm[1,0] + Tm[0,2]*Tm[2,0] + Tm[1,2]*Tm[2,1]) )
        if abs(c0.imag) > NeutrinoState._normTolerance:
            raise ArithmeticError("c0 = %s"% c0)
        if abs(c1.imag) > NeutrinoState._normTolerance:
            raise ArithmeticError("c1 = %s"% c1)

        c0 = c0.real
        c1 = c1.real
        arcTanArg = sqrt_real(-c0*c0 - (4/27.)*c1*c1*c1) / c0
        arcTanTerm = (1/3.)*atan(arcTanArg)
        sqrtMC1 = sqrt_real(-c1)
        sqrtMC1D3 = sqrtMC1/sqrt_real(3.)
        cosATT = cos( arcTanTerm )
        sinATT = sin( arcTanTerm )
        # Calculate final eigenvalues
        eigenvals = [(-sqrtMC1D3*cosATT + sqrtMC1*sinATT),
                     (-sqrtMC1D3*cosATT - sqrtMC1*sinATT),
                     (2*sqrtMC1D3*cosATT)]
        if c0 > 0:
            eigenvals = [-1*eigenval for eigenval in eigenvals]
        #self._checkEigenvalues(eigenvals, Tm)
        #eg = eigenvals
        #if abs(eg[0]*eg[1]*eg[2] + c0) > 1.0e-5:
        #    #if True:
        #    TrTm = Tm[0,0]+Tm[1,1]+Tm[2,2]
        #    print "   Bad zeros: %f %f %f %f %f" % (
        #        TrTm.real,
        #        TrTm.imag,
        #        sum(eigenvals), 
        #        eg[0]*eg[1]+eg[0]*eg[2]+eg[1]*eg[2] - c1,
        #        eg[0]*eg[1]*eg[2] + c0)
        #    #print "   eigen:",eigenvals
        return eigenvals


    def _determinant(self, T):
        "Determinant of 3x3 matrix"
        det = (T[0,0]*(T[1,1]*T[2,2] - T[1,2]*T[2,1])
               -T[0,1]*(T[1,0]*T[2,2] - T[1,2]*T[2,0])
               +T[0,2]*(T[1,0]*T[2,1] - T[1,1]*T[2,0]))
        return det

    def evolve(self, initialState, distance ):
        "Evolve neutrino quantum state, and return final state vector"
        finalAmps = [complex(0), complex(0), complex(0)]
        amplFactor = self._amplFactorMatter
        energyEigvals = self._lambdaMatter
        if not initialState.isMatter():
            amplFactor = self._amplFactorAntiMatter
            energyEigvals = self._lambdaAntiMatter
        initialAmplitudes = initialState.amplitudes()
        #traceHm = sum(energyEigvals)
        #totalPhase = exp( -1j*distance*traceHm/3.0 )
        for finalFlavorIdx in range(3):
            for initFlavorIdx in range(3):
                for massIdx in range(3):
                    finalAmps[finalFlavorIdx] += (
                        exp( -1j * distance * energyEigvals[massIdx] / hbarc)
                        * amplFactor[initFlavorIdx, 
                                     finalFlavorIdx,
                                     massIdx]
                        * initialAmplitudes[initFlavorIdx] )
            # Neglect overall phase, since it cancels in any prob. calculation
            #finalAmps[finalFlavorIdx] *= totalPhase
        finalState = NeutrinoState(finalAmps[0],
                                   finalAmps[1],
                                   finalAmps[2],
                                   initialState.isMatter())
        return finalState

    def _checkEigenvalues(self, eig, Tm):
        "Use characteristic eq to confirm eigenvalue calculation"
        identity3x3 = array( [[1, 0, 0],
                              [0, 1, 0],
                              [0, 0, 1]] )
        printResult = False
        result = []
        for eVal in eig:
            res = self._determinant(Tm-(eVal*identity3x3))
            result.append( res )
            if abs(res)>1.0e-14: 
                printResult=True
        if printResult:
            print "  eig:",eig
            print "  res:",result
        return
