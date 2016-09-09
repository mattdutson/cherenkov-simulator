// common.h
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#ifndef Common_h
#define Common_h

namespace cherenkov_simulator {
    
    class ConfigManager {
        
    private:
        
        boost::program_options::options_description command_line_;
        
        boost::program_options::options_description file_options_;
        
        boost::program_options::variables_map option_map_;
        
    public:
        
        ConfigManager();
        
        void ParseCommandLine(int argc, const char* argv[]);
        
        void ParseOptionsFile();
        
        void HelpMessage();
        
        template <typename Type>
        Type Get(std::string path) {
            try {
                return option_map_[path].as<Type>();
            }
            catch (...) {
                throw std::exception();
            }
        }
    };

    extern ConfigManager globalConfig;
}


#endif
