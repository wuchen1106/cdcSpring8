#include "MyRuntimeParameters.hxx"
#include <cstdlib>

#include "Log.hxx"
#include "HEPUnits.hxx"
#include "HEPConstants.hxx"

using namespace std;

MyRuntimeParameters* MyRuntimeParameters::fMyRuntimeParameters = NULL;

MyRuntimeParameters::MyRuntimeParameters() {

  MyLog("*** Initializing MyRuntimeParameters");

  // Build the set of units for Geant4.
  fUnitsTableParser = new MyUnitsTableParser();

  // We do not open any parameter files by default in the constructor.  
  // Instead, parameter files are open automatically when an entry
  // in the file is requested by the user.

}


MyRuntimeParameters::~MyRuntimeParameters() {

}

/// Read in each parameter file.
void MyRuntimeParameters::ReadInputFile(TString fileName, TString dirName, bool tryFile, bool fixParameters) {

  TString fullFileName;

  // Make full file name.
  if(!dirName.CompareTo(""))
      fullFileName = fileName.Data();    
  else
      fullFileName = dirName.Data() + TString("/") + fileName.Data();
  
  MyInfo("Opening parameter file: " << fullFileName.Data());
  
  ifstream inputFile(fullFileName.Data(), ios::in);
  
  // If input file doesn't exist.
  if (!inputFile){
    // Just return if 'tryFile' is specified.
    if(tryFile){
      MyInfo("\n***** MyRuntimeParameters::ReadInputFile *****\n"
                 << "Cannot open input file '" << fullFileName.Data() << "'.");
      return;
    }else{ // Otherwise, throw exception.
      MyError("\n***** MyRuntimeParameters::ReadInputFile *****\n"
                 << "Cannot open input file '" << fullFileName.Data() << "'.");
      throw EBadParameterFile();
    }
  }
  
  int inputState = 0;
  string inputString;
  string parameterName;
  string parameterValueUnit;
  string parameterValue;

  while (inputFile >> inputString) {
    if (inputState == 0) {
      if (inputString == "<")
        inputState = 1;
    }
    else if (inputState == 1) {
      parameterName = inputString;
      inputState = 2;
    }
    else if (inputState == 2) {
      if (inputString == "=")
        inputState = 3;
      else {
        
        MyError("\n***** MyRuntimeParameters::ReadInputFile *****\n"
                   "Input file '" << fileName << "'. Last parameter '"
                   << parameterName << "'.\n" 
                   << "Cannot find symbol '='.\n"
                   << "Badly formatted parameters file.");
        throw EBadParameterFile();
      }
    }
    else if (inputState == 3) {	
      //        parameterValue = atof(inputString.c_str());
        parameterValue = inputString.c_str();
	parameterValueUnit = inputString;
	inputState = 4;
    }
    else if (inputState == 4) {	
      if (inputString == ">"){
          // Finished reading. Save parameter; but only if the parameter
          // isn't already 'fixed'
          if(fixedParameters.find(parameterName) == fixedParameters.end()){
              mapOfMyRuntimeParameters[parameterName] = parameterValue;
              // If fixParameters bool is set, fix this parameter now.
              if(fixParameters)
                  fixedParameters.insert(parameterName);
          }
          inputState = 0;
      }
      else if (inputString == "<") {
        MyError("\n***** MyRuntimeParameters::ReadInputFile *****\n"
                   << "Input file '" << fileName << "'. Last parameter '"
                   << parameterName << "'.\n" 
                   << "Unexpected symbol '<'.\n"
                   << "Badly formatted parameters file.");
        throw EBadParameterFile();
      }
      else {

	// The parameter must have a unit.  Resave the value with the correct unit.
	parameterValueUnit.append(" ");
	parameterValueUnit.append(inputString);

	// Use MyUnitsTableParser class to convert string of value+unit to 
	// a double with value+unit.
	parameterValue = fUnitsTableParser->Convert2DoubleWithUnit(parameterValueUnit);
      }
    }
  }
  
  if (inputState != 0) {    
    MyError("\n***** MyRuntimeParameters::ReadInputFile *****\n"
               << "Input file '" << fileName << "'. Last parameter '"
               << parameterName << "'.\n"
               << "Cannot find symbol '>' at the end of file.\n"
               << "Badly formatted parameters file."); 
    throw EBadParameterFile();
  }  
  inputFile.close();
}


void MyRuntimeParameters::PrintListOfParameters() {

  MyInfo("***** MyRuntimeParameters::PrintListOfParameters *****");
  MyInfo("List of parameters:");
  
  for (mapIterator i = mapOfMyRuntimeParameters.begin();
       i != mapOfMyRuntimeParameters.end(); i++)
    MyInfo("  " << (*i).first << " = " << (*i).second);
  MyInfo("");

}

int MyRuntimeParameters::TryLoadingParametersFile(std::string parameterName){

  MyInfo("Trying to load parameters file for parameter: " << parameterName);

  // Figure out the name of the package
  int pos = parameterName.find(".");
  if (pos<0){
      MyWarn("Cannot find the packageName from parameter "<<parameterName);
      return 0;
  }
  TString packageName(parameterName.c_str(), pos);

  // and the file name for parameters file
  TString fileName = packageName + ".parameters.dat";

  // and the directory of parameters file.
  packageName.ToUpper();
  TString packageROOT = packageName + "ROOT";
  TString dirName =  getenv(packageROOT.Data()) 
    + TString("/parameters/");
  if (dirName == "/parameters/"){
      MyWarn("Cannot find enviroment variable for package "<<packageName);
      dirName = "parameters/";
  }
  // Now try reading in this file.  Last input variable is set to true,
  // indicating that we don't want to throw exception if a file is not found.
  ReadInputFile(fileName,dirName,true);
 
  // Now try to find this parameter again
  mapIterator i = mapOfMyRuntimeParameters.find(parameterName);  
  if(i != mapOfMyRuntimeParameters.end())
    return 1;
  else{
    return 0;
  }
}


void MyRuntimeParameters::ReadParamOverrideFile(TString filename){

    MyLog("Using MyRuntimeParameters override file = " << filename);

    // Setting final input variable to true forces the parameters
    // that are loaded to be 'fixed'; ie immutable.
    ReadInputFile(filename,"",false,true);

}


bool MyRuntimeParameters::HasParameter(string parameterName) {

  mapIterator i = mapOfMyRuntimeParameters.find(parameterName);

  if(i != mapOfMyRuntimeParameters.end())
    return 1;
  else{
    
    // Okay, we didn't find this parameter on the first try.  
    // Let's see if we can try loading the parameters file.
    // The function will return whether this parameter 
    // was found afterwards.
    return TryLoadingParametersFile(parameterName);

  }
}


bool MyRuntimeParameters::HasPrefix(string paramPrefix){
    // Search for element whose key would not go before the prefex, i.e. in 
    // greater than or equal to
    mapIterator i = mapOfMyRuntimeParameters.lower_bound(paramPrefix);
    if (i != mapOfMyRuntimeParameters.end()) {
        const string& key = i->first;
        // Check if the prefix is a match to the element
        if (key.compare(0, paramPrefix.size(), paramPrefix) == 0) return true;
    };
    // Otherwise return false
    return false;
}

int MyRuntimeParameters::GetParameterI(string parameterName) {
  
  if(HasParameter(parameterName)) {
    return atoi(mapOfMyRuntimeParameters[parameterName].c_str());
  }else{
    MyError("\n***** MyRuntimeParameters::GetParameterAsInteger *****\n"
                << "Cannot find parameter '" << parameterName << "'.");
    throw ENonexistantDatabaseParameter();
  }
  return -1;

}


double MyRuntimeParameters::GetParameterD(string parameterName) {

  if(HasParameter(parameterName)) {
    return atof(mapOfMyRuntimeParameters[parameterName].c_str());
  }else{
    MyError("\n***** MyRuntimeParameters::GetParameterAsDouble *****\n"
                << "Cannot find parameter '" << parameterName << "'.");
    throw ENonexistantDatabaseParameter();
  }
  return -1;

}


std::string MyRuntimeParameters::GetParameterS(string parameterName) {

  if(HasParameter(parameterName)) {
    return mapOfMyRuntimeParameters[parameterName];
  }else{
    MyError("\n***** MyRuntimeParameters::GetParameterAsString *****\n"
                << "Cannot find parameter '" << parameterName << "'.");
    throw ENonexistantDatabaseParameter();
  }
  return std::string();

}


void MyRuntimeParameters::ClearMapOfMyRuntimeParameters() {
  mapOfMyRuntimeParameters.clear();
}
