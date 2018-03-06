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
    std::ifstream in_file;
    if(std::ifstream(prep_file.c_str()))
        in_file.open(prep_file.c_str());
    else
    {
        throw PrepFileProcessingException(__LINE__, "Prep file not found");
    }
    Read(in_file);
    in_file.close();            /// Close the prep files
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
void PrepFile::Read(std::ifstream &in_file)
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
        PrepFileSpace::PrepFileResidue::PrepFileAtomVector atoms = residue->GetAtoms();
        for(PrepFileSpace::PrepFileResidue::PrepFileAtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
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
