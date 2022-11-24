#include <iostream>

#include <TROOT.h>
#include <TStopwatch.h>
#include <TFile.h>
#include <TString.h>

#include "FemtoDstAnalyzer_New.C"
//#include "functions.C"

int main(int argc, char** argv)
{
  TString inFileName, inFileNameRec, inFileNameFlat, outFileName, WorkMode, Energy;
  TStopwatch timer;
  if (argc < 12)
  {
    std::cerr << "./FemtoDstAnalyzer_New -i INPUTFILE -o OUTPUTFILE -r INPUTFILEREC -f INPUTFILEFLAT -m WORKMODE -g ENERGY" << std::endl;
    return 13;
  }

  for (int i=1;i<argc;i++)
  {
    if (std::string(argv[i]) != "-i" &&
        std::string(argv[i]) != "-o" &&
        std::string(argv[i]) != "-m" &&
        std::string(argv[i]) != "-r" &&
        std::string(argv[i]) != "-f" &&
        std::string(argv[i]) != "-g" )
    {
      std::cerr << "\nUnknown parameter: " << argv[i] << std::endl;
      return 14;
    }
    else
    {
      if (std::string(argv[i]) == "-i" && i != argc-1)
      {
        inFileName = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-i" && i == argc-1)
      {
        std::cerr << "\nInput file was not specified!" << std::endl;
	return 15;
      }
      if (std::string(argv[i]) == "-o" && i != argc-1)
      {
        outFileName = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-o" && i == argc-1)
      {
        std::cerr << "\nOutput file was not specified!" << std::endl;
	return 16;
      }
      if (std::string(argv[i]) == "-r" && i != argc-1)
      {
        inFileNameRec = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-r" && i == argc-1)
      {
        std::cerr << "\ninFile_recentrening file was not specified!" << std::endl;
  return 17;
      }
      if (std::string(argv[i]) == "-f" && i != argc-1)
      {
        inFileNameFlat = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-f" && i == argc-1)
      {
        std::cerr << "\nInput file (flattening) was not specified!" << std::endl;
  return 18;
      }
      if (std::string(argv[i]) == "-m" && i != argc-1)
      {
        WorkMode = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-m" && i == argc-1)
      {
        std::cerr << "\nmode not specified!" << std::endl;
  return 19;
      }
      if (std::string(argv[i]) == "-g" && i != argc-1)
      {
        Energy = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-g" && i == argc-1)
      {
        std::cerr << "\nenergy not specified!" << std::endl;
  return 20;
      }
    }
  }
  if (inFileName == "" || outFileName == "" || WorkMode == "" || Energy == "")
  {
    std::cerr << "\nInput/Output file has not been set properly!" << std::endl;
    return 21;
  }

  FemtoDstAnalyzer_New(inFileName.Data(),outFileName.Data(),inFileNameRec.Data(),inFileNameFlat.Data(),WorkMode.Data(),Energy.Atoi());

  return 0;
}

