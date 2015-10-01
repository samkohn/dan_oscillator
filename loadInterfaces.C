// Compile and load interface classes
{

  TString globesLibPrefix = "/Users/skohn/Documents/DUNE/software/install/lib";
  TString globesIncPrefix = "/Users/skohn/Documents/DUNE/software/install/include";

  gSystem->Load((globesLibPrefix + "/libglobes.dylib").Data());
  
  TString include = ".include ";
  gROOT->ProcessLine( (include + globesIncPrefix).Data() ); 

  TString load = ".L ";
  TString interfacePrefix = ".";  
  gROOT->ProcessLine( (load + interfacePrefix + "/GLoBES_Interface.C+").Data() );

}
