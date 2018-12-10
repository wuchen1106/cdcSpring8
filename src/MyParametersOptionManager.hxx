#ifndef TParametersOptionManager_hxx_seen
#define TParametersOptionManager_hxx_seen

#include "MyRuntimeParameters.hxx"

#include <string>

/// Class to parse options for runtime parameters.
/// Specifically, the class provides a method for overriding 
/// certain parameters using a command-line specified text file.
class MyParametersOptionManager {
public:
  
  ///Default Constructor
  MyParametersOptionManager(){};

  ///Virtual D'tor
  virtual ~MyParametersOptionManager(){};

  /// This function checks if the option is relevant for parameter 
  /// overriding.
  bool IsRelevantOption(std::string option); 

  /// This function actually does the overriding.
  void UseRelevantOption(std::string value); 

  ///Prints the usage information.  Should be called by the Usage()
  ///method of any event loop that incorporates these options
  void Usage();


private:



};

#endif //TParametersOptionManager_hxx_seen
