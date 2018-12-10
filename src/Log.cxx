#include <vector>
#include <sstream>
#include <fstream>
#include <time.h>

#include "Log.hxx"

Log::ErrorPriority Log::fErrorPriority = Log::ErrorLevel;
Log::LogPriority Log::fLogPriority = Log::LogLevel;
std::ostream* Log::fDebugStream = NULL;
std::ostream* Log::fLogStream = NULL;
std::map<std::string,Log::ErrorPriority> Log::fErrorTraces;
std::map<std::string,Log::LogPriority> Log::fLogTraces;
int Log::fIndentation = 0;

Log::Log() { }
Log::~Log() { }

void Log::SetDebugLevel(const char* trace, 
                              Log::ErrorPriority level) {
    fErrorTraces[trace] = level;
}

Log::ErrorPriority Log::GetDebugLevel(const char* trace) {
    std::map<std::string,ErrorPriority>::iterator elem = fErrorTraces.find(trace);
    if (elem == fErrorTraces.end()) return fErrorPriority;
    return elem->second;
}

namespace {
    std::string MakeTimeStamp() {
        std::string stamp = "Unknown Time";
        time_t t = time(NULL);
        struct tm *local = localtime(&t);
        if (!local) return stamp;
        char localTime[80];
        if (!strftime(localTime,sizeof(localTime),"%c",local)) return stamp;
        struct tm *utc = gmtime(&t);
        if (!utc) return stamp;
        char utcTime[80];
        if (!strftime(utcTime,sizeof(utcTime),"%Y-%m-%d %H:%M (UTC)",utc)) 
            return stamp;
        stamp = localTime;
        stamp += " [";
        stamp += utcTime;
        stamp += "]";
        return stamp;
    }
}

void Log::SetDebugStream(std::ostream* err) {
    Log::fDebugStream = err;
    if (!fDebugStream) return;
    std::ofstream* ofile = dynamic_cast<std::ofstream*>(err);
    if (ofile && !(ofile->is_open())) {
        fDebugStream = NULL;
        MySevere("Debug stream is not open.");
    }
    *fDebugStream << std::endl
                  << "##################################################" 
                  << std::endl
                  << "# ERROR LOG STARTS AT: " << MakeTimeStamp()
                  << std::endl
                  << "##################################################" 
                  << std::endl
                  << std::endl;
}

std::ostream& Log::GetDebugStream() {
    if (!Log::fDebugStream ) return GetLogStream();
    return *Log::fDebugStream;
}
    
void Log::SetLogLevel(const char* trace, 
                            Log::LogPriority level) {
    fLogTraces[trace] = level;
}

Log::LogPriority Log::GetLogLevel(const char* trace) {
    std::map<std::string,LogPriority>::iterator elem = fLogTraces.find(trace);
    if (elem == fLogTraces.end()) return fLogPriority;
    return elem->second;
}

void Log::SetLogStream(std::ostream* log) {
    Log::fLogStream = log;
    if (!fLogStream) return;
    std::ofstream* ofile = dynamic_cast<std::ofstream*>(log);
    if (ofile && !(ofile->is_open())) {
        fLogStream = NULL;
        MySevere("Log stream is not open.");
    }
    *fLogStream << std::endl
                << "##################################################" 
                << std::endl
                << "# LOG STARTS AT: " << MakeTimeStamp()
                << std::endl
                << "##################################################" 
                << std::endl
                << std::endl;
}

std::ostream& Log::GetLogStream() {
    if (!Log::fLogStream) return std::cout;
    return *Log::fLogStream;
}

void Log::SetIndentation(int i) {
    Log::fIndentation = std::max(i,0);
}

void Log::IncreaseIndentation() {
    ++Log::fIndentation;
}

void Log::DecreaseIndentation() {
    if (Log::fIndentation>0) --Log::fIndentation;
}
    
void Log::ResetIndentation() {
    Log::fIndentation = 0;
}


std::string Log::MakeIndent() {
    if (fIndentation<1) return "";
    std::string indent = "";
    for (int i=0; i<fIndentation; ++i) {
        indent += "..";
    }
    indent += " ";
    return indent;
}

namespace {
    bool TranslateLogLevel(const std::string& name, 
                           Log::LogPriority& level) {
        if (name == "QuietLevel") {
            level = Log::QuietLevel;
            return true;
        }
        if (name == "LogLevel") {
            level = Log::LogLevel;
            return true;
        }
        if (name == "InfoLevel") {
            level = Log::InfoLevel;
            return true;
        }
        if (name == "VerboseLevel") {
            level = Log::VerboseLevel;
            return true;
        }
        return false;
    }

    bool TranslateErrorLevel(const std::string& name, 
                             Log::ErrorPriority& level) {
        if (name == "SilentLevel") {
            level = Log::SilentLevel;
            return true;
        }
        if (name == "ErrorLevel") {
            level = Log::ErrorLevel;
            return true;
        }
        if (name == "SevereLevel") {
            level = Log::SevereLevel;
            return true;
        }
        if (name == "WarnLevel") {
            level = Log::WarnLevel;
            return true;
        }
        if (name == "DebugLevel") {
            level = Log::DebugLevel;
            return true;
        }
        if (name == "TraceLevel") {
            level = Log::TraceLevel;
            return true;
        }
        return false;
    }

    std::ostream* StreamPointer(const std::string& name) {
        if (name == "STDCOUT") return &std::cout;
        if (name == "STDCERR") return &std::cerr;
        if (name[0] != '"') return NULL;
        if (name[name.size()-1] != '"') return NULL;
        std::string file = name.substr(1,name.size()-2);
        std::ofstream* output = new std::ofstream(file.c_str(),
                                                  std::ios::out|std::ios::app);
        if (output->is_open()) return output;
        return NULL;
    }

    bool ReadConfigurationFile(const char* config) {
        std::ifstream input(config);
        if (!input.is_open()) return false;

        int inputLine = 0;
        for (;;) {
            std::string line;
            std::getline(input,line);
            if (input.eof()) break;

            // Save the current line number and cache the value so error
            // messages can be printed later.
            std::string cache(line);
            ++inputLine;
            
            // Strip the comments out of the file.
            std::string::size_type position = line.find("#");
            if (position != std::string::npos) line.erase(position);

            // Strip the white space at the beginning of the line.
            line.erase(0,line.find_first_not_of("\t "));

            // Skip lines that are too short.
            if (line.size()==0) continue;

            // Split the line into fields and a value.
            position = line.find("=");
            if (position == std::string::npos) {
                // Houston, we have a problem... There isn't a value.
                std::cerr << "WARNING: " << config << ":" << inputLine << ": "
                          << "Configuration line missing an '='"
                          << std::endl;
                std::cerr << "  Line: <" << cache << ">"
                          << std::endl;
                std::cerr << "  Configuration line has been skip" << std::endl;
                continue;
            }

            // Split the value off the end of the line.
            std::string value = line.substr(position+1);
            line.erase(position);

            // Strip the white space at the beginning of the value.
            value.erase(0,value.find_first_not_of("\t "));

            // Strip the white space at the end of the value.
            position = value.find_last_not_of("\t ");
            if (position != std::string::npos) value.erase(position+1);
            
            // Strip the white space at the end of the fields.
            position = line.find_last_not_of("\t ");
            if (position != std::string::npos) line.erase(position+1);
            
            // Split the remaining line in to fields.
            std::vector<std::string> fields;
            for (;;) {
                position = line.find(".");
                if (position == std::string::npos) {
                    fields.push_back(line);
                    break;
                }
                fields.push_back(line.substr(0,position));
                line.erase(0,position+1);
            } 
            
            // Process the fields and value.
            if (fields.size() == 2
                && fields[0] == "log"
                && fields[1] == "file") {
                // Set the log file name.
                std::ostream* str = StreamPointer(value);
                if (!str) {
                    std::cerr << "WARNING: " << config << ":" 
                              << inputLine << ": "
                              << "Cannot open log stream."
                              << std::endl;
                    std::cerr << "  Line: <" << cache << ">"
                              << std::endl;
                    std::cerr << "  Configuration line has been skip" 
                              << std::endl;
                    continue;
                }
                Log::SetLogStream(str);
            }
            else if (fields.size() == 2
                     && fields[0] == "error"
                     && fields[1] == "file") {
                // Set the error file name.
                std::ostream* str = StreamPointer(value);
                if (!str) {
                    std::cerr << "WARNING: " << config << ":" 
                              << inputLine << ": "
                              << "Cannot open error stream."
                              << std::endl;
                    std::cerr << "  Line: <" << cache << ">"
                              << std::endl;
                    std::cerr << "  Configuration line has been skip" 
                              << std::endl;
                    continue;
                }
                Log::SetDebugStream(str);
            }
            else if (fields.size() == 3
                     && fields[0] == "log"
                     && fields[1] == "default"
                     && fields[2] == "level") {
                // Set the default log level.
                Log::LogPriority level;
                if (!TranslateLogLevel(value,level)) {
                    std::cerr << "WARNING: " << config << ":" 
                              << inputLine << ": "
                              << "Unknown log level name."
                              << std::endl;
                    std::cerr << "  Line: <" << cache << ">"
                              << std::endl;
                    std::cerr << "  Configuration line has been skip" 
                              << std::endl;
                    continue;
                }
                Log::SetLogLevel(level);
            }
            else if (fields.size() == 3
                     && fields[0] == "error"
                     && fields[1] == "default"
                     && fields[2] == "level") {
                // Set the default error level.
                Log::ErrorPriority level;
                if (!TranslateErrorLevel(value,level)) {
                    std::cerr << "WARNING: " << config << ":" 
                              << inputLine << ": "
                              << "Unknown error level name."
                              << std::endl;
                    std::cerr << "  Line: <" << cache << ">"
                              << std::endl;
                    std::cerr << "  Configuration line has been skip" 
                              << std::endl;
                    continue;
                }
                Log::SetDebugLevel(level);
            }
            else if (fields.size() == 3
                     && fields[0] == "log"
                     && fields[2] == "level") {
                // Set the log level.
                Log::LogPriority level;
                if (!TranslateLogLevel(value,level)) {
                    std::cerr << "WARNING: " << config << ":" 
                              << inputLine << ": "
                              << "Unknown log level name."
                              << std::endl;
                    std::cerr << "  Line: <" << cache << ">"
                              << std::endl;
                    std::cerr << "  Configuration line has been skip" 
                              << std::endl;
                    continue;
                }
                Log::SetLogLevel(fields[1].c_str(),level);
            }
            else if (fields.size() == 3
                     && fields[0] == "error"
                     && fields[2] == "level") {
                // Set the error level.
                Log::ErrorPriority level;
                if (!TranslateErrorLevel(value,level)) {
                    std::cerr << "WARNING: " << config << ":" 
                              << inputLine << ": "
                              << "Unknown error level name."
                              << std::endl;
                    std::cerr << "  Line: <" << cache << ">"
                              << std::endl;
                    std::cerr << "  Configuration line has been skip" 
                              << std::endl;
                    continue;
                }
                Log::SetDebugLevel(fields[1].c_str(),level);
            }
            else {
                std::cerr << "WARNING: " << config << ":" << inputLine << ": "
                          << "Unknown command."
                          << std::endl;
                std::cerr << "  Line: <" << cache << ">"
                          << std::endl;
                std::cerr << "  Configuration line has been skip" << std::endl;
            }                
        }

        return true;
    }
}

void Log::Configure(const char* conf) {
    // Try to read a local configuration file.  
    ReadConfigurationFile("./cometlog.config");
    if (conf) {
        bool success = ReadConfigurationFile(conf);
        if (!success) MyLog("MyLog configuration file was not read.");
    }
}
