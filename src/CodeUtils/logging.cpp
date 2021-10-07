#include <ctime>
#include <fstream>
#include <iostream>
#include "includes/CodeUtils/logging.hpp"

void gmml::log(int line, std::string file_path, LogLevel level, std::string msg, std::string out_file_name)
{
	std::ofstream file;
	if(out_file_name == "")
	{
		char* gemshome_env_var = std::getenv("GEMSHOME");
		// Check if the environment variables exist.
		if(!gemshome_env_var)
		{
			std::cerr << "\nMust set GEMSHOME environment variable.\n\n    BASH:   export GEMSHOME=/path/to/gems\n    SH:     setenv GEMSHOME /path/to/gems\n" << std::endl;
			return;
		}
		std::string GEMSHOME(gemshome_env_var);
		out_file_name = GEMSHOME + "/gmml/GMML_Log.txt";
	}
	file.open(out_file_name.c_str(), std::ios_base::app);
	time_t t = time(0);
	std::string time_str = std::asctime(std::localtime(&t));
	file << time_str.substr(0, time_str.size() - 1) << " >>> " << file_path << ":" << line << " >>>";
	switch(level)
	{
	case INF:
		file << " [INFO]: ";
		break;
	case ERR:
		file << " [ERROR]: ";
		break;
	case WAR:
		file << " [WARNING]: ";
		break;
	}
	file << msg << std::endl;
	file.close();
}
