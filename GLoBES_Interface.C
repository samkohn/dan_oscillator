//
// Rough Interface to interact with GLoBES oscillation calculator
// 
//
// Created by: dadwyer@lbl.gov 2013/05/01
// Last Modified: dadwyer@lbl.gov 2013/05/01
//

#include <globes/globes.h>

#ifndef GLOBES_INTERFACE_H
#define GLOBES_INTERFACE_H

class GLoBES_Interface {
 public:
  GLoBES_Interface();
  virtual ~GLoBES_Interface();
  void setOscillation(double sinTh12, double sinTh13, double sinTh23,
		      double dmSq21, double dmSq32, double deltaCP);
  double probability(int flavorInit, int flavorFinal, bool isMatter, 
		     double energy, double distance, double elecDensity=0);
};

#endif //GLOBES_INTERFACE_H

GLoBES_Interface::GLoBES_Interface()
{
  char procName[17] = "GLoBES_Interface";
  glbInit(procName);
}

GLoBES_Interface::~GLoBES_Interface()
{
  ;
}

void GLoBES_Interface::setOscillation(double sinTh12, double sinTh13, 
				      double sinTh23, double dmSq21, 
				      double dmSq32, double deltaCP)
{
  // Set the oscillation parameters for GLoBES
  //   Note GLoBES units:
  //     radians for angles and CP-phase 
  //     eV^2 for mass squared differences
  double theta12 = asin(sinTh12);
  double theta13 = asin(sinTh13);
  double theta23 = asin(sinTh23);
  double dmSq31 = dmSq32 + dmSq21;
  glb_params true_values = glbAllocParams();
  glbDefineParams(true_values,theta12,theta13,theta23,deltaCP,dmSq21,dmSq31);
  glbSetOscillationParameters(true_values);
  return;
}

double GLoBES_Interface::probability(int flavorInit, int flavorFinal, 
				     bool isMatter, double energy, 
				     double distance, double elecDensity/*=0*/)
{
  // Calculate oscillation probability for vacuum or constant density
  //   Note GLoBES units:
  //     flavor indices: e=1, mu=2, tau=3
  //     energy: GeV
  //     distance: km
  //     elecdensity: electrons/m^3
  int matterFlag=1;
  if( !isMatter ) matterFlag = -1;
  if(elecDensity){
    // Matter-enhanced oscillation
    //   Convert electron density to globes 'density'
    double globesVfactor = 7.5e-14;
    double globesEfactor = 0.5;
    double sqrt2 = 1.4142135623730951;
    double Gf = 8.961808701567238e-44;
    double density = elecDensity*sqrt2*Gf/(globesVfactor*globesEfactor);
    return glbConstantDensityProbability(flavorInit,flavorFinal,
					 matterFlag,energy,distance,density);
  }
  // Vacuum oscillation
  return glbVacuumProbability(flavorInit,flavorFinal,
			      matterFlag,energy,distance);
}

