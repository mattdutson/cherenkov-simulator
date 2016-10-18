// common.cpp
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include "Common.h"

using namespace boost::program_options;
using namespace std;

namespace cherenkov_simulator
{
    
    ConfigManager::ConfigManager(std::string name): allowed_options(name) {}

    std::string ConfigManager::HelpMessage()
    {
        stringstream output;
        output << allowed_options;
        return output.str();
    }
    
    CommandOptions::CommandOptions(): ConfigManager("Command line options") {
        allowed_options.add_options()
        ("config_file", value<std::string>()->default_value("config.txt"), "The path to the configuration file")
        ("help", "Show help message");
    }
    
    void CommandOptions::ParseCommand(int argc, const char **argv) {
        store(parse_command_line(argc, argv, allowed_options), option_map);
    }
    
    FileOptions::FileOptions(): ConfigManager("File options") {
        allowed_options.add_options()
        ("mirror_radius", value<double>()->required(), "The radius of curvature of the telescope mirror");
        notify(option_map);
    }
    
    void FileOptions::ParseFile(std::string filename) {
        ifstream file(filename.c_str());
        store(parse_config_file(file, allowed_options), option_map);
        notify(option_map);
    }
}
