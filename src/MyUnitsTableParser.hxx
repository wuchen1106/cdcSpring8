#ifndef MyUnitsTableParser_h
#define MyUnitsTableParser_h 
#include "TString.h"
#include <iostream>
#include <string>
#include <fstream>
#include <map>

/// This class provides a method for converting a 
/// string like "1.5 cm" into a double with the 
/// appropriate unit.  To do so it defines a set 
/// of units, using the same base units as in 
/// oaEvent/src/HEPUnits.hxx: ie mm, ns, MeV...
/// Only a fairly limited set of units is defined.

class MyUnitsTableParser {

public:

  /// Constructor.  Creates list of units.
  MyUnitsTableParser();
  ~MyUnitsTableParser();
  
  /// Converts a string like "1.5 cm" into a double
  /// with the appropriate units.
  //  double Convert2DoubleWithUnit(std::string line);
  std::string Convert2DoubleWithUnit(std::string line);

  /// Prints all the defined units.
  void PrintListOfUnits();

private:

  std::map<std::string, double> units;

};

#endif
