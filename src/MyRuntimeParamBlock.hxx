#ifndef MyRUNTIMEPARAMETERS_MyRUNTIMEPARAMBLOCK_HXX
#define MyRUNTIMEPARAMETERS_MyRUNTIMEPARAMBLOCK_HXX


#include <MyRuntimeParameters.hxx>
#include <sstream>

class MyRuntimeParamBlock{
    public:
        MyRuntimeParamBlock(const std::string& section, const std::string& package="");
        ~MyRuntimeParamBlock(){}

        int GetInt(const std::string& name)const{
            return MyRuntimeParameters::Get().GetParameterI(MakeName(name));
        }
        double GetDouble(const std::string& name)const{
            return MyRuntimeParameters::Get().GetParameterD(MakeName(name));
        }
        std::string GetString(const std::string& name)const{
            return MyRuntimeParameters::Get().GetParameterS(MakeName(name));
        }

        int HasParameter(const std::string& name)const{
            std::string full_name=MakeName(name);
            return MyRuntimeParameters::Get().HasParameter(full_name);
        }
        int GetInt(const std::string& name, int def_val)const{
            std::string full_name=MakeName(name);
            if(!MyRuntimeParameters::Get().HasParameter(full_name)) return def_val;
            return MyRuntimeParameters::Get().GetParameterI(full_name);
        }
        double GetDouble(const std::string& name, double def_val)const{
            std::string full_name=MakeName(name);
            if(!MyRuntimeParameters::Get().HasParameter(full_name)) return def_val;
            return MyRuntimeParameters::Get().GetParameterD(full_name);
        }
        std::string GetString(const std::string& name, const std::string& def_val)const{
            std::string full_name=MakeName(name);
            if(!MyRuntimeParameters::Get().HasParameter(full_name)) return def_val;
            return MyRuntimeParameters::Get().GetParameterS(full_name);
        }

        const std::string& GetPackage()const{return fPackage;}
        const std::string& GetSection()const{return fSection;}
        std::string PrintParam(const std::string& param, const std::string& value) const{
            return " < " + MakeName(param) + " = " + value + " > ";}
        std::string PrintParam(const std::string& param, 
                               const int value) const{
            std::stringstream ss;
            ss << value;
            return PrintParam(param, ss.str());}
        std::string PrintParam(const std::string& param, 
                               const double value) const{
            std::stringstream ss;
            ss << value;
            return PrintParam(param, ss.str());}



    private:
        std::string MakeName(const std::string&)const;
        std::string fPackage, fSection;
};

#endif // MyRUNTIMEPARAMETERS_MyRUNTIMEPARAMBLOCK_HXX
