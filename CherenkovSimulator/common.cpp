// common.cpp
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
//

#include "Common.h"

using namespace boost::program_options;
using namespace std;

namespace cherenkov_simulator {
    
    ConfigManager::ConfigManager(): command_line_("Command Line Options"), file_options_("Configuration File Options") {
        command_line_.add_options()
        ("config_file", value<std::string>()->default_value("config.txt"), "The path to the configuration file")
        ("help", "Show help message");
        
        file_options_.add_options()
        ("mirror_radius", value<double>(), "The radius of curvature of the telescope mirror");
        
        options_description sharedOptions;
        sharedOptions.add_options()
        ("test", value<std::string>(), "The test to be performed");
        file_options_.add(sharedOptions);
        command_line_.add(sharedOptions);
    }

    void ConfigManager::ParseCommandLine(int argc, const char* argv[]) {
        store(parse_command_line(argc, argv, command_line_), option_map_);
        notify(option_map_);
    }

    void ConfigManager::ParseOptionsFile() {
        string filename = Get<string>("config_file");
        ifstream file(filename.c_str());
        store(parse_config_file(file, file_options_), option_map_);
        notify(option_map_);
    }

    void ConfigManager::HelpMessage() {
        cout << command_line_ << endl;
    }
        
    ConfigManager globalConfig;
}