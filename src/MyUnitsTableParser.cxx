#include "MyUnitsTableParser.hxx"
#include <iostream>
#include <sstream>
#include <string>

#include <math.h>

MyUnitsTableParser::MyUnitsTableParser() {

  // Define units.  Add to this list as needed.

  // Length 
  units["pc"] = 3.08568e+19 ;
  units["km"] = 1e+06 ;
  units["m"] = 1000 ;
  units["cm"] = 10 ;
  units["mm"] = 1 ;
  units["um"] = 0.001 ;
  units["nm"] = 1e-06 ;
  units["Ang"] = 1e-07  ;
  units["fm"] = 1e-12 ;

  // Area
  units["km2"] = 1e+12 ;
  units["m2"] = 1e+06 ;
  units["cm2"] = 100 ;
  units["mm2"] = 1 ;
  units["barn"] = 1e-22 ;
  units["mbarn"] = 1e-25 ;
  units["mubarn"] = 1e-28  ;
  units["nbarn"] = 1e-31 ;
  units["pbarn"] = 1e-34 ;

  // Volume
  units["km3"] = 1e+18  ;
  units["m3"] = 1e+09 ;
  units["cm3"] = 1000  ;
  units["mm3"] = 1  ;

  // Degree
  units["rad"] = 1 ;
  units["mrad"] = 0.001  ;
  units["sr"] = 1 ;
  units["deg"] = 0.0174533 ;

  // Time
  units["s"] = 1e+09 ;
  units["ms"] = 1e+06 ;
  units["mus"] = 1000  ;
  units["ns"] = 1 ;
  units["ps"] = 0.001 ;

  // Frequency
  units["Hz"] = 1e-09 ;
  units["kHz"] = 1e-06  ;
  units["MHz"] = 0.001  ;

  // Electric Charge
  units["e+"] = 1  ;
  units["C"] = 6.24151e+18 ;

  // Energy
  units["eV"] = 1e-06 ;
  units["keV"] = 0.001  ;
  units["MeV"] = 1  ;
  units["GeV"] = 1000  ;
  units["TeV"] = 1e+06  ;
  units["PeV"] = 1e+09  ;
  units["J"] = 6.24151e+12 ;

  // Energy/Length
  units["GeV/cm"] = 100  ;
  units["MeV/cm"] = 0.1  ;
  units["keV/cm"] = 0.0001  ;
  units["eV/cm"] = 1e-07 ;

  // Inverse energy  
  units["1/eV"] = 1.0/units["eV"];  ;
  units["1/keV"] = 1.0/units["keV"];  ;
  units["1/MeV"] = 1.0/units["MeV"];  ;
  units["1/GeV"] = 1.0/units["GeV"];  ;

  // Mass
  units["mg"] = 6.24151e+18  ;
  units["g"] = 6.24151e+21 ;
  units["kg"] = 6.24151e+24  ;

  // Volumic Mass
  units["g/cm3"] = 6.24151e+18 ;
  units["mg/cm3"] = 6.24151e+15  ;
  units["kg/m3"] = 6.24151e+15 ;

  // Mass/Surface
  units["g/cm2"] = 6.24151e+19 ;
  units["mg/cm2"] = 6.24151e+16  ;
  units["kg/cm2"] = 6.24151e+22  ;

  // Surface/Mass
  units["cm2/g"] = 1.60218e-20  ;

  // Energy*Surface/Mass
  units["eV*cm2/g"] = 1.60218e-26 ;
  units["keV*cm2/g"] = 1.60218e-23  ;
  units["MeV*cm2/g"] = 1.60218e-20  ;
  units["GeV*cm2/g"] = 1.60218e-17  ;

  // Power
  units["W"] = 6241.51  ;

  // Force
  units["N"] = 6.24151e+09  ;

  // Pressure
  units["Pa"] = 6241.51 ;
  units["bar"] = 6.24151e+08  ;
  units["atm"] = 6.32421e+08  ;

  // Electric current
  units["A"] = 6.24151e+09 ;
  units["mA"] = 6.24151e+06 ;
  units["muA"] = 6241.51  ;
  units["nA"] = 6.24151 ;

  // Electric potential
  units["V"] = 1e-06 ;
  units["kV"] = 0.001  ;
  units["MV"] = 1  ;

  // Magnetic flux
  units["Wb"] = 1000  ;

  // Magnetic flux density
  units["T"] = 0.001 ;
  units["kG"] = 0.0001  ;
  units["G"] = 1e-07 ;

  // Speed
  units["cm/us"] = units["cm"]/units["mus"]  ;
  units["cm/ns"] = units["cm"]/units["ns"]  ;
  units["mm/ns"] = units["mm"]/units["ns"]  ;

  // Length/Energy
  units["mm/MeV"] = units["mm"]/units["MeV"];
  units["mm/keV"] = units["mm"]/units["keV"];
  units["cm/MeV"] = units["cm"]/units["MeV"];
  units["cm/keV"] = units["cm"]/units["keV"];

  // Dummy units for diffusion coefficient
  units["cm/sqrt(cm)"] = units["cm"]/sqrt(units["cm"]);
  units["mm/sqrt(cm)"] = units["mm"]/sqrt(units["cm"]);
  units["um/sqrt(cm)"] = units["um"]/sqrt(units["cm"]);

  // Dummy units for electron mobility
  units["cm2/(Vs)"] = units["cm2"]/(units["V"]*units["s"]);
}


MyUnitsTableParser::~MyUnitsTableParser() {

}


std::string MyUnitsTableParser::Convert2DoubleWithUnit(std::string input){

  double value;
  std::string unit;
  
  std::istringstream line(input);
  if(!(line >> value >> unit)){
    std::cerr << "MyUnitsTableParser: badly formatted input string. Returning 0."<<std::endl;
    return "0.0";
  }
  
  // Check if requested unit is in map.
  if(units.find(unit) == units.end()){
    std::cerr << "MyUnitsTableParser: requested unit '"
              << unit << "' not found. Returning 0."<<std::endl;
    return "0.0";
  }

  
  value = value * units[unit];

  char s[256];
  sprintf(s,"%f",value);


  return s;

}


void MyUnitsTableParser::PrintListOfUnits() {

  std::cout << std::endl;
  std::cout << "***** List of available units *****" << std::endl <<std::endl;  
  for (std::map<std::string, double>::iterator unit = units.begin();
       unit != units.end(); unit++)
    std::cout << (*unit).first <<  std::endl;
  std::cout << std::endl;

}
