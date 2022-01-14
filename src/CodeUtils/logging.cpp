#include <ctime>
#include <fstream>
#include <iostream>
#include "includes/CodeUtils/logging.hpp"

void gmml::log(int line, std::string file_path, LogLevel level, std::string msg, std::string out_file_name)
{
	std::ofstream file;
	if(out_file_name == "")
	{
		char* gmmlLogEnvVar = std::getenv("GMMLLOG");
		// Check if the environment variables exist.
		if(!gmmlLogEnvVar)
		{
			return; // Log nothing unless GMMLLOG is set in the user's environment.
		}
		std::string gmmlLog(gmmlLogEnvVar);
		out_file_name = gmmlLog;
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
