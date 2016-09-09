//
//  Common.h
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 8/31/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#ifndef Common_h
#define Common_h

class ConfigManager {

private:
    
    boost::program_options::options_description commandLine;
    
    boost::program_options::options_description fileOptions;
    
    boost::program_options::variables_map optionMap;
    
public:
    
    ConfigManager();
    
    void ParseCommandLine(int argc, const char* argv[]);
    
    void ParseOptionsFile();
    
    void HelpMessage();
    
    template <typename Type>
    Type Get(std::string path) {
        try {
            return optionMap[path].as<Type>();
        }
        catch (...) {
            throw std::exception();
        }
    }
};

extern ConfigManager globalConfig;

#endif
