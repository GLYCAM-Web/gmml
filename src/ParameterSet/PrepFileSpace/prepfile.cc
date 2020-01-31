#include <fstream>
#include <iostream>
#include <iomanip>

#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"

using PrepFileSpace::PrepFile;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
PrepFile::PrepFile(const std::string& prep_file)
{
    path_ = prep_file;
    std::ifstream in_file(prep_file.c_str());
    if(in_file.is_open())
    {
        ReadAllResidues(in_file);
        in_file.close();            /// Close the prep files
    }
    else
    {
        throw PrepFileProcessingException(__LINE__, "Prep file not found");
    }
}

PrepFile::PrepFile(const std::string& prep_file, std::vector<std::string>& query_residue_names)
{
    path_ = prep_file;
    std::ifstream in_file(prep_file.c_str());
    if(in_file.is_open())
    {
        ReadOnlyQueryResidues(in_file, query_residue_names);
        in_file.close();            /// Close the prep files
    }
    else
    {
        throw PrepFileProcessingException(__LINE__, "Prep file not found");
    }
}

PrepFile::PrepFile()
{
    path_ = "";
    residues_ = ResidueMap();
}

PrepFile::~PrepFile()
{
    residues_.clear();
    residues_ = ResidueMap();
}
//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
PrepFile::ResidueMap& PrepFile::GetResidues()
{
    return residues_;
}

std::vector<std::string> PrepFile::GetAllResidueNames()
{
    std::vector<std::string> residue_names;
    for(PrepFile::ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++){
        std::string residue_name = (*it).first;
        residue_names.push_back(residue_name);
    }
    return residue_names;
}

gmml::ResidueNameMap PrepFile::GetAllResidueNamesMap()
{
    gmml::ResidueNameMap residue_names = gmml::ResidueNameMap();
    for(PrepFile::ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++){
        std::string residue_name = (*it).first;
        residue_names[residue_name] = residue_name;
    }
    return residue_names;
}

std::vector<std::string> PrepFile::GetAllAtomNamesOfResidue(std::string residue_name)
{
    std::vector<std::string> atom_names_of_residue;
    ResidueMap residue_map = GetResidues();
    PrepFileSpace::PrepFileResidue* prep_file_residue = residue_map[residue_name];
    std::vector<PrepFileAtom*> atoms = prep_file_residue->GetAtoms();
    for(std::vector<PrepFileAtom*>::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        PrepFileAtom* atom = (*it);
        atom_names_of_residue.push_back(atom->GetName());
    }
    return atom_names_of_residue;
}

//////////////////////////////////////////////////////////
//                         MUTATORS                     //
//////////////////////////////////////////////////////////
void PrepFile::SetPath(std::string path)
{
    path_ = path;
}

void PrepFile::SetResidues(ResidueMap residues)
{
    residues_.clear();
    for(ResidueMap::iterator it = residues.begin(); it != residues.end(); it++)
    {
        PrepFileSpace::PrepFileResidue* residue = (*it).second;
        std::string residue_name = (*it).first;
        residues_[residue_name] = residue;
    }
}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
void PrepFile::ReadAllResidues(std::ifstream &in_file)
{
    std::string header1, header2;
    getline(in_file, header1);
    getline(in_file, header2);


    PrepFileSpace::PrepFileResidue *residue = ProcessResidueSection(in_file);
    while (residue != NULL)
    {
        residues_[residue->name_] = residue;
        residue = ProcessResidueSection(in_file);
    }
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

    PrepFileSpace::PrepFileResidue *residue;
    if (std::find(query_residue_names.begin(), query_residue_names.end(), resname) != query_residue_names.end() ) 
    {
	in_file.seekg(one_line_before_residue_title);  //go back to one line before residue title, so the rest of the codes works correctly 
        residue = ProcessResidueSection(in_file);
        residues_[residue->name_] = residue;
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
            residue = ProcessResidueSection(in_file);
            residues_[residue->name_] = residue;
	    one_line_before_residue_title = in_file.tellg();
	}
        resname.clear();
	intx.clear();
	kform.clear();
    }
}

PrepFileSpace::PrepFileResidue* PrepFile::ProcessResidueSection(std::ifstream &in_file)
{
    PrepFileSpace::PrepFileResidue* residue = new PrepFileSpace::PrepFileResidue();
    residue = residue->LoadFromStream(in_file);
    return residue;
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
        throw PrepFileProcessingException(__LINE__,"File could not be created");
    }
    try
    {
        this->BuildPrepFile(out_file);
    }
    catch(...)
    {
        out_file.close();
    }
}
void PrepFile::BuildPrepFile(std::ofstream &stream)
{
    stream << std::endl
           << std::endl;
    for(ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        PrepFileSpace::PrepFileResidue* residue = (*it).second;
        stream << residue->GetTitle() << std::endl << std::endl
               << std::left << std::setw(4) << residue->GetName() << " " << std::right << std::setw(3) << residue->GetStringFormatOfCoordinateType() << " "
               << std::setw(1) << residue->GetOutputFormat() << std::endl
               << residue->GetStringFormatOfGeometryType() << " " << residue->GetStringFormatOfDummyAtomOmission() << " "
               << residue->GetDummyAtomType() << " " << residue->GetStringFormatOfDummyAtomPosition() << std::endl
               << std::right << std::setw(8) << std::fixed << std::setprecision(3) << residue->GetCharge() << std::endl;
        PrepFileSpace::PrepFileAtomVector atoms = residue->GetAtoms();
        for(PrepFileSpace::PrepFileAtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
        {
            PrepFileAtom* atom = (*it1);
            stream << std::right << std::setw(2) << atom->GetIndex() << " " << std::left << std::setw(4) << atom->GetName() << " " << std::left << std::setw(3) << atom->GetType() << " "
                   << std::setw(1) << atom->GetStringFormatOfTopologicalType() << " " << std::right << std::setw(2) << atom->GetBondIndex() << " " << std::right << std::setw(2) << atom->GetAngleIndex() << " "
                   << std::right << std::setw(2) << atom->GetDihedralIndex() << " ";
            if(atom->GetBondLength() != gmml::dNotSet)
                stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << atom->GetBondLength() << " ";
            else
                stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << 0.0 << " ";
            if(atom->GetAngle() != gmml::dNotSet)
                   stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << atom->GetAngle() << " ";
            else
                stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << 0.0 << " ";
            if(atom->GetDihedral() != gmml::dNotSet)
                stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << atom->GetDihedral();
            else
                stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << 0.0;
            stream << "    " << std::right << std::setw(8) << std::fixed << std::setprecision(4) << atom->GetCharge() << std::endl;
        }
        stream << std::endl;
        PrepFileSpace::PrepFileResidue::DihedralVector dihedrals = residue->GetImproperDihedrals();
        if (dihedrals.size() != 0)
        {
            stream << "IMPROPER" << std::endl;
            for(PrepFileSpace::PrepFileResidue::DihedralVector::iterator it2 = dihedrals.begin(); it2 != dihedrals.end(); it2++)
            {
                PrepFileSpace::PrepFileResidue::Dihedral dihedral = (*it2);
                stream << std::left << std::setw(4) << dihedral.at(0) << std::left << std::setw(4) << dihedral.at(1) << std::left << std::setw(4) << dihedral.at(2) << std::left << std::setw(4) << dihedral.at(3) << std::endl;
            }
            stream << std::endl;
        }
        PrepFileSpace::PrepFileResidue::Loop loops = residue->GetLoops();
        if(loops.size() != 0)
        {
            stream << "LOOP" << std::endl;
            for(PrepFileSpace::PrepFileResidue::Loop::iterator it3 = loops.begin(); it3 != loops.end(); it3++)
            {
                int first_atom_index = (*it3).first;
                int second_atom_index = (*it3).second;
                stream << residue->GetAtomNameByIndex(first_atom_index) << " " << residue->GetAtomNameByIndex(second_atom_index) << std::endl;
            }
        }
        stream << std::endl
               << "DONE" << std::endl;
    }
    stream << "STOP";
}



//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void PrepFile::Print(std::ostream& out)
{
    for(ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        out << "**********************************************************************************" << std::endl;
        it->second->Print(out);
    }
}
