//
//  Common.cpp
//  CherenkovSimulator
//
//  Created by Matthew Dutson on 9/8/16.
//  Copyright © 2016 Matthew Dutson. All rights reserved.
//

#include "Common.h"

using namespace boost::program_options;
using namespace std;

ConfigManager::ConfigManager(): commandLine("Command Line Options"), fileOptions("Configuration File Options") {
    commandLine.add_options()
        ("config_file", value<std::string>()->default_value("config.txt"), "The path to the configuration file")
        ("help", "Show help message");
    
    fileOptions.add_options()
        ("mirror_radius", value<double>(), "The radius of curvature of the telescope mirror");
    
    options_description sharedOptions;
    sharedOptions.add_options()
    ("test", value<std::string>(), "The test to be performed");
    fileOptions.add(sharedOptions);
    commandLine.add(sharedOptions);
}

void ConfigManager::ParseCommandLine(int argc, const char* argv[]) {
    store(parse_command_line(argc, argv, commandLine), optionMap);
    notify(optionMap);
}

void ConfigManager::ParseOptionsFile() {
    string filename = Get<string>("config_file");
    ifstream file(filename.c_str());
    store(parse_config_file(file, fileOptions), optionMap);
    notify(optionMap);
}

void ConfigManager::HelpMessage() {
    cout << commandLine << endl;
}

ConfigManager globalConfig;