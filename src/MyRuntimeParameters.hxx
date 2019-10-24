#ifndef MyRuntimeParameters_hxx_seen
#define MyRuntimeParameters_hxx_seen

#include "TString.h"
#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <set>
#include "MyUnitsTableParser.hxx"
#include "Log.hxx"
#include "ErrorCore.hxx"

// General DummyDatabase exceptions
MYEXCEPTION(EOARuntimeParameters,COMET::EoaCore);
// Exception for missing database parameter
MYEXCEPTION(ENonexistantDatabaseParameter,EOARuntimeParameters);
// Exception for reading parameter file
MYEXCEPTION(EBadParameterFile,EOARuntimeParameters);
MYEXCEPTION(EBadParameterConversion,EOARuntimeParameters);
MYEXCEPTION(EDuplicateParameterInCurrentFile,EOARuntimeParameters);

/// This class is meant to provide access to a set of parameters
/// that are defined in a set of text files; in particular, the text files
/// are used to store 'software runtime parameters' in calibration and 
/// reconstruction algorithms.
class MyRuntimeParameters {
 public:
  ~MyRuntimeParameters();

  ///  Get a reference to the singleton instance of dummy 
  ///  database information.  If this is first attempt at reference
  ///  then singleton is instantiated and parameters are 
  ///  read from text files.
  static MyRuntimeParameters& Get(void) {
    if (!fMyRuntimeParameters)
      fMyRuntimeParameters = new MyRuntimeParameters();
    return *fMyRuntimeParameters;
  }

  /// Check if Parameter is stored in database
  bool HasParameter(std::string);

  /// Check if parameter with prefix is stored in database
  bool HasPrefix(std::string);

  /// Get parameter.  Value is returned as bool. 
  bool GetParameterB(std::string);

  /// Get parameter.  Value is returned as integer. 
  int GetParameterI(std::string);

  /// Get parameter.  Value is returned as double. 
  double GetParameterD(std::string);

  /// Get parameter.  Value is returned as string. 
  std::string GetParameterS(std::string);

  /// Get parameter. Value type is defined by template
  /// throw exception when parameter is not found
  template<typename T> void GetParameter(std::string, T&);

  void ClearMapOfMyRuntimeParameters();

  MyUnitsTableParser* GetUnitsTableParser(){
    return fUnitsTableParser;
  }

    /// This command allows the user to set parameters from the
    /// command line; the command is different from the standard
    /// file reading, in that the parameters that are set are 'fixed'.
    /// Ie, they are immutable and cannot be changed, even if they
    /// exist in some other parameters file that is read in later.
    void ReadParamOverrideFile(TString filename);


  /// Reads parameters from input files.
  /// Function can be used to read in extra parameter files.
  void ReadInputFile(TString filename, TString dirName="", bool tryFile = false, bool fixParameters = false);  

  /// Prints list of saved parameters
  void PrintListOfParameters();
 private:

  MyRuntimeParameters();

  // Add a non-functional private copy constructor.
  MyRuntimeParameters(const MyRuntimeParameters &){
      MyLog("MyRuntimeParameters copy constructor is private.  Shouldn't be called!!!");
  }


  /// Helper method to attempt to read parameters file based on parameters name.
  /// The return value indicates whether we can now find this parameter after 
  /// trying to open this file (return 1 means parameter found).
  int TryLoadingParametersFile(std::string parameterName);

  /// map containing list of parameters and their values.
  std::map<std::string, std::string, std::less<std::string> > mapOfMyRuntimeParameters;
  typedef std::map<std::string, std::string, std::less<std::string> >::iterator mapIterator;

  /// A set of all the 'fixed'.  'Fixed' parameters cannot be changed, even if they reappear
  /// in a different parameters file.  This is used for the parameter override files.
  std::set<std::string> fixedParameters;

  MyUnitsTableParser *fUnitsTableParser;

  /// The static pointer to the singleton instance.
  static MyRuntimeParameters* fMyRuntimeParameters;
};

#endif
