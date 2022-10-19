#include <fstream>
#include <iostream>
#include <iomanip>
#include <ios>

#include "includes/common.hpp"
#include "includes/utils.hpp"
#include "includes/ParameterSet/PrepFile/prepFile.hpp"
#include "includes/ParameterSet/PrepFile/prepResidue.hpp"
#include "includes/ParameterSet/PrepFile/prepAtom.hpp"
#include "includes/CodeUtils/files.hpp" // ensureFileExists
#include "includes/CodeUtils/strings.hpp" // split

using prep::PrepFile;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
PrepFile::PrepFile(const std::string& prep_file)
{
	std::cout << "Constructor entered" << std::endl;
	codeUtils::ensureFileExists(prep_file);
	std::ifstream in_file(prep_file.c_str());
	if(in_file.is_open())
	{
		this->ReadAllResidues(in_file);
		in_file.close();            /// Close the prep file
	}
	else
	{
		throw std::runtime_error("Prep file exists but couldn't be opened.");
	}
}

PrepFile::PrepFile(const std::string& prep_file, const std::vector<std::string> queryNames)
{
	codeUtils::ensureFileExists(prep_file);
	std::ifstream in_file(prep_file.c_str());
	if(in_file.is_open())
	{
		this->ReadQueryResidues(in_file, queryNames);
		in_file.close();            /// Close the prep files
	}
	else
	{
		throw std::runtime_error( "Prep file exists but couldn't be opened.");
	}
}
//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
//std::vector<std::string> PrepFile::GetAllResidueNames()
//{
//	std::vector<std::string> residue_names;
//	for(PrepFile::ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++){
//		std::string residue_name = (*it).first;
//		residue_names.push_back(residue_name);
//	}
//	return residue_names;
//}

//gmml::ResidueNameMap PrepFile::GetAllResidueNamesMap()
//{
//	gmml::ResidueNameMap residue_names = gmml::ResidueNameMap();
//	for(PrepFile::ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++){
//		std::string residue_name = (*it).first;
//		residue_names[residue_name] = residue_name;
//	}
//	return residue_names;
//}

//std::vector<std::string> PrepFile::GetAllAtomNamesOfResidue(std::string residue_name)
//{
//    std::vector<std::string> atom_names_of_residue;
//    ResidueMap residue_map = GetResidues();
//    prep::PrepResidue* prep_file_residue = residue_map[residue_name];
//    std::vector<PrepFileAtom*> atoms = prep_file_residue->GetAtoms();
//    for(std::vector<PrepFileAtom*>::iterator it = atoms.begin(); it != atoms.end(); it++)
//    {
//        PrepFileAtom* atom = (*it);
//        atom_names_of_residue.push_back(atom->GetName());
//    }
//    return atom_names_of_residue;
//}

//////////////////////////////////////////////////////////
//                         MUTATORS                     //
//////////////////////////////////////////////////////////
void PrepFile::SetAtomConnectivities()
{
    for ( auto &residue : this->getResidues() )
    {
        residue->SetConnectivities();
    }
    return;
}

void PrepFile::Generate3dStructures()
{
    for ( auto &residue : this->getResidues() )
    {
        residue->Generate3dStructure();
    }
    return;
}
//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
void PrepFile::ReadAllResidues(std::ifstream &in_file)
{
	std::string line;
	getline(in_file, line);
	getline(in_file, line); // first two lines are always blank apparently. smh.
	getline(in_file, line); // This should be first line of residue entry.
    while (gmml::Trim(line).find("STOP") == std::string::npos)           /// End of file
    {
        this->addResidue(std::make_unique<PrepResidue>(in_file, line));
        getline(in_file, line); // This should be first line of next residue entry or STOP.
        //std::cout << "Back out and line is: " << line << std::endl;
    }
	//std::cout << "Ok this is done with line as:\n " << line << std::endl;
}

// Reads each line of the file. If it finds one of the query residues it reads it in. Won't read in query repeats twice.
void PrepFile::ReadQueryResidues(std::ifstream &in_file, const std::vector<std::string>& queryNames)
{
	std::string line;
	getline(in_file, line);
	getline(in_file, line); // first two lines are always blank apparently. smh.
	getline(in_file, line); // This should be first line of residue entry.
    while (gmml::Trim(line).find("STOP") == std::string::npos) // While not at end of file
    {
    	std::streampos firstResidueLinePosition = in_file.tellg(); // save correct position to start reading residue
    	// Need to move to line with residue name on it.
    	getline(in_file, line); // title
    	getline(in_file, line); // residue name appears here
    	std::vector<std::string> residueNameLine = codeUtils::split(line, ' '); // front() string will be name
    	//std::cout << "Current residue name is: " << residueNameLine.front() << std::endl;
    	if (std::find(queryNames.begin(), queryNames.end(), residueNameLine.front()) != queryNames.end() )
    	{
    		std::cout << "Found query residue: " << residueNameLine.front() << "\n";
    		in_file.seekg(firstResidueLinePosition);  //go back here so the residue constructor works
    		this->addResidue(std::make_unique<PrepResidue>(in_file, line));
    	}
    	else
    	{ // need to flush the lines until we find a residue we want.
    		while(gmml::Trim(line).find("DONE") == std::string::npos)
    		{
    			getline(in_file, line);
    		}
    	}
        getline(in_file, line); // This should be first line of next residue entry or STOP.
        //std::cout << "Back out and line is: " << line << std::endl;
    }
	//std::cout << "Ok this is done with line as:\n " << line << std::endl;
}

void PrepFile::Write(const std::string& prep_file)
{
	std::ofstream out_file;
	try
	{
		out_file.open(prep_file.c_str());
	}
	catch(...)
	{
		throw std::runtime_error("PrepFile could not be created for writing");
	}
	try
	{
		this->Write(out_file);
	}
	catch(...)
	{
		out_file.close();
	}
}

void PrepFile::Write(std::ofstream &stream)
{
	stream << "\n" << "\n";
	for(auto &residue: this->getResidues())
	{
		residue->Write(stream);
	}
	stream << "STOP\n";
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
std::string PrepFile::Print() const
{
    std::string out;
	for(auto &residue : this->getResidues() )
	{
		out += "**********************************************************************************\n";
		out += residue->Print();
	}
	return out;
}
