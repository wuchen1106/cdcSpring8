#include "MyRuntimeParamBlock.hxx"

MyRuntimeParamBlock::MyRuntimeParamBlock(
        const std::string& section, const std::string& package):
    fPackage(package),fSection(section){
}

std::string MyRuntimeParamBlock::MakeName(const std::string& name)const{
    std::string full_name;
    if(!fPackage.empty())full_name+=fPackage+".";
    if(!fSection.empty())full_name+=fSection+".";
    full_name+=name;
    return full_name;
}
