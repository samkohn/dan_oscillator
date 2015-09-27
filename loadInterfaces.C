// Compile and load interface classes
{

  TString globesLibPrefix = "/Users/dandwyer/Tools/installs/globes-3.1.11/lib";  
  TString globesIncPrefix = "/Users/dandwyer/Tools/installs/globes-3.1.11/include";  

  gSystem->Load((globesLibPrefix + "/libglobes.dylib").Data());
  
  TString include = ".include ";
  gROOT->ProcessLine( (include + globesIncPrefix).Data() ); 

  TString load = ".L ";
  TString interfacePrefix = ".";  
  gROOT->ProcessLine( (load + interfacePrefix + "/GLoBES_Interface.C+").Data() );

}
