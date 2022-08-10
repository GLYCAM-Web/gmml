#include <fstream>
#include <iostream>
#include <iomanip>
#include <ios>

#include "includes/common.hpp"
#include "includes/utils.hpp"
#include "includes/ParameterSet/PrepFile/prepFile.hpp"
#include "includes/ParameterSet/PrepFile/PrepResidue.hpp"
#include "includes/ParameterSet/PrepFile/prepAtom.hpp"
#include "includes/CodeUtils/files.hpp" // ensureFileExists

using prep::PrepFile;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
PrepFile::PrepFile(const std::string& prep_file)
{
	codeUtils::ensureFileExists(prep_file);
	std::ifstream in_file(prep_file.c_str());
	if(in_file.is_open())
	{
		ReadAllResidues(in_file);
		in_file.close();            /// Close the prep file
	}
	else
	{
		throw std::runtime_error("Prep file exists but couldn't be opened.");
	}
}

PrepFile::PrepFile(const std::string& prep_file, std::vector<std::string>& query_residue_names)
{
	codeUtils::ensureFileExists(prep_file);
	std::ifstream in_file(prep_file.c_str());
	if(in_file.is_open())
	{
		ReadOnlyQueryResidues(in_file, query_residue_names);
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


//void PrepFile::SetResidues(ResidueMap residues)
//{
//    residues_.clear();
//    for(ResidueMap::iterator it = residues.begin(); it != residues.end(); it++)
//    {
//        prep::PrepResidue* residue = (*it).second;
//        std::string residue_name = (*it).first;
//        residues_[residue_name] = residue;
//    }
//}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
void PrepFile::ReadAllResidues(std::ifstream &in_file)
{
	std::string header1, header2;
	getline(in_file, header1);
	getline(in_file, header2);
    this->addResidue(std::make_unique<PrepResidue>(in_file));
}

void PrepFile::ReadOnlyQueryResidues(std::ifstream &in_file, std::vector<std::string>& query_residue_names)
{
	std::string query_names_str;
	std::string header1, header2;
	getline(in_file, header1);
	getline(in_file, header2);
	//Have to reach the line containing residue name to tell if this residue should be read. However, residue processing starts from the residue title line, which is two lines above
	//Sadly, for getline() there is no easy way back.
	//So, the position of the title line is set as a stream pointer. When reading a residue, make the stream pointer go back to this checkpoint.
	std::streampos one_line_before_residue_title;
	one_line_before_residue_title = in_file.tellg();
	std::string line_str;
	getline(in_file, line_str);
	getline(in_file, line_str);
	std::stringstream temp_stream;
	std::string resname,intx,kform;
	getline(in_file, line_str);
	temp_stream << line_str;
	temp_stream >> resname >> intx >> kform;
	if (std::find(query_residue_names.begin(), query_residue_names.end(), resname) != query_residue_names.end() )
	{
		in_file.seekg(one_line_before_residue_title);  //go back to one line before residue title, so the rest of the codes works correctly
	    this->addResidue(std::make_unique<PrepResidue>(in_file));
		one_line_before_residue_title = in_file.tellg();
	}
	resname.clear();
	intx.clear();
	kform.clear();
	while (getline(in_file, line_str))
	{
		if (line_str.find("STOP") != std::string::npos)
		{
			break;
		}
		if (line_str.find("DONE") != std::string::npos)
		{
			one_line_before_residue_title = in_file.tellg();
			getline(in_file, line_str);
			getline(in_file, line_str);
			getline(in_file, line_str);
			std::stringstream line_stream;
			line_stream << line_str;
			line_stream >> resname >> intx >> kform;
		}
		if (std::find(query_residue_names.begin(), query_residue_names.end(), resname) != query_residue_names.end() )
		{
			in_file.seekg(one_line_before_residue_title);  //go back to one line before residue title, so the rest of the codes works correctly
		    this->addResidue(std::make_unique<PrepResidue>(in_file));
			one_line_before_residue_title = in_file.tellg();
		}
		resname.clear();
		intx.clear();
		kform.clear();
	}
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
	stream << std::endl
			<< std::endl;
	for(auto &residue: this->getResidues())
	{
		residue->Write(stream);
	}
	stream << "STOP";
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void PrepFile::Print(std::ostream& out)
{
	for(auto &residue : this->getResidues() )
	{
		out << "**********************************************************************************" << std::endl;
		residue->Print(out);
	}
}
