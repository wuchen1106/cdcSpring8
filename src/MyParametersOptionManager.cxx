#include "MyParametersOptionManager.hxx"
#include "Log.hxx"

bool MyParametersOptionManager::IsRelevantOption(std::string option){

  if(option == "par_override")
    return true;
  else
    return false;

  return false;
}

void MyParametersOptionManager::UseRelevantOption(std::string value){

  MyLog("Using the runtime parameter override file: " << value);
  MyRuntimeParameters::Get().ReadParamOverrideFile(value);

}


void MyParametersOptionManager::Usage(){

  std::cout << "Options provided by MyParametersOptionManager:" << std::endl;

  std::cout << "    -O par_override=<parameter file> : sets an override parameter file." << std::endl;
  std::cout << "        This allows the user to override the default " 
            << "sets of parameters using the specified file"<< std::endl;                                         


}

