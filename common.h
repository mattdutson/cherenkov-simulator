// common.h
// cherenkov_simulator
//
// Created by Matthew Dutson on 9/8/16.
//
// Contains a configuration manager and globally useful methods

#ifndef Common_h
#define Common_h

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <iostream>
#include <fstream>
#include <string>

namespace cherenkov_simulator
{
    constexpr double light_speed = 3e8;

    // A template for a class containing the simulator configuration.
    class ConfigManager
    {
        
    protected:
        
        boost::program_options::options_description allowed_options;
        
        boost::program_options::variables_map option_map;
        
        ConfigManager(std::string name);
        
    public:
        
        std::string HelpMessage();
        
        template <typename Type>
        Type Contains(std::string path)
        {
            return option_map.count(path) > 0;
        }
        
        template <typename Type>
        Type Get(std::string path)
        {
            try
            {
                return option_map[path].as<Type>();
            }
            catch (...)
            {
                throw std::exception();
            }
        }
    };

    // A class containing options parsed from the command line.
    class CommandOptions: public ConfigManager
    {
        
    public:
        
        CommandOptions();
        
        void ParseCommand(int argc, const char** argv);
        
    };
    
    class FileOptions: public ConfigManager
    {
        
    public:
        
        FileOptions();
        
        void ParseFile(std::string filename);
    };
}


#endif
