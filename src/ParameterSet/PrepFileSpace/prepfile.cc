#include <fstream>
#include <iostream>
#include <iomanip>

#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"

using namespace std;
using namespace gmml;
using namespace PrepFileSpace;

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
    in_file.close();            /// Close the parameter files
}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
PrepFile::ResidueMap& PrepFile::GetResidues()
{
    return residues_;
}

vector<string> PrepFile::GetAllResidueNames()
{    
    vector<string> residue_names;
    for(PrepFile::ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++){
        string residue_name = (*it).first;
        residue_names.push_back(residue_name);
    }
    return residue_names;
}

ResidueNameMap PrepFile::GetAllResidueNamesMap()
{
    ResidueNameMap residue_names = ResidueNameMap();
    for(PrepFile::ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++){
        string residue_name = (*it).first;
        residue_names[residue_name] = residue_name;
    }
    return residue_names;
}

vector<string> PrepFile::GetAllAtomNamesOfResidue(string residue_name)
{
    vector<string> atom_names_of_residue;
    ResidueMap residue_map = GetResidues();
    PrepFileResidue* prep_file_residue = residue_map[residue_name];
    vector<PrepFileAtom*> atoms = prep_file_residue->GetAtoms();
    for(vector<PrepFileAtom*>::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        PrepFileAtom* atom = (*it);
        atom_names_of_residue.push_back(atom->GetName());
    }
    return atom_names_of_residue;
}

AtomNameMap PrepFile::GetAllAtomNamesOfResidueMap(string residue_name)
{
    AtomNameMap atom_names_of_residue = AtomNameMap();
    ResidueMap residue_map = GetResidues();
    PrepFileResidue* prep_file_residue = residue_map[residue_name];
    vector<PrepFileAtom*> atoms = prep_file_residue->GetAtoms();
    for(vector<PrepFileAtom*>::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        PrepFileAtom* atom = (*it);
        atom_names_of_residue[atom->GetName()] = atom->GetName();
    }
    return atom_names_of_residue;
}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
void PrepFile::Read(ifstream &in_file)
{
    string header1, header2;
    getline(in_file, header1);
    getline(in_file, header2);

    PrepFileResidue *residue = ProcessResidueSection(in_file);
    while (residue != NULL)
    {
        residues_[residue->name_] = residue;
        residue = ProcessResidueSection(in_file);
    }
}

PrepFileResidue* PrepFile::ProcessResidueSection(ifstream &in_file)
{
    PrepFileResidue* residue = new PrepFileResidue();
    residue = residue->LoadFromStream(in_file);
    return residue;
}
void PrepFile::Write(const string& prep_file)
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
void PrepFile::BuildPrepFile(ofstream &stream)
{
    stream << endl
           << endl;
    for(ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        PrepFileResidue* residue = (*it).second;
        stream << residue->GetTitle() << endl << endl
               << left << setw(4) << residue->GetName() << " " << right << setw(3) << residue->GetStringFormatOfCoordinateType() << " " << setw(1) << residue->GetOutputFormat() << endl
               << residue->GetStringFormatOfGeometryType() << " " << residue->GetStringFormatOfDummyAtomOmission() << " " << residue->GetDummyAtomType() << " " << residue->GetStringFormatOfDummyAtomPosition() << endl
               << right << setw(8) << fixed << setprecision(3) << residue->GetCharge() << endl;
        PrepFileResidue::PrepFileAtomVector atoms = residue->GetAtoms();
        for(PrepFileResidue::PrepFileAtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
        {
            PrepFileAtom* atom = (*it1);
            stream << right << setw(2) << atom->GetIndex() << " " << left << setw(4) << atom->GetName() << " " << left << setw(3) << atom->GetType() << " "
                   << setw(1) << atom->GetStringFormatOfTopologicalType() << " " << right << setw(2) << atom->GetBondIndex() << " " << right << setw(2) << atom->GetAngleIndex() << " "
                   << right << setw(2) << atom->GetDihedralIndex() << " " << right << setw(8) << fixed << setprecision(3) << atom->GetBondLength() << " "
                   << right << setw(8) << fixed << setprecision(3) << atom->GetAngle() << " " << right << setw(8) << fixed << setprecision(3) << atom->GetDihedral()
                   << "    " << right << setw(8) << fixed << setprecision(4) << atom->GetCharge() << endl;
        }
        stream << endl;
        PrepFileResidue::DihedralVector dihedrals = residue->GetImproperDihedrals();
        if (dihedrals.size() != 0)
        {
            stream << "IMPROPER" << endl;
            for(PrepFileResidue::DihedralVector::iterator it2 = dihedrals.begin(); it2 != dihedrals.end(); it2++)
            {
                PrepFileResidue::Dihedral dihedral = (*it2);
                stream << left << setw(4) << dihedral.at(0) << left << setw(4) << dihedral.at(1) << left << setw(4) << dihedral.at(2) << left << setw(4) << dihedral.at(3) << endl
                       << endl;

            }
        }
        PrepFileResidue::Loop loops = residue->GetLoops();
        if(loops.size() != 0)
        {
            stream << "LOOP" << endl;
            for(PrepFileResidue::Loop::iterator it3 = loops.begin(); it3 != loops.end(); it3++)
            {
                int first_atom_index = (*it3).first;
                int second_atom_index = (*it3).second;
                stream << residue->GetAtomNameByIndex(first_atom_index) << " " << residue->GetAtomNameByIndex(second_atom_index) << endl;
            }
        }
        stream << endl
               << "DONE" << endl;
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
        out << "**********************************************************************************" << endl;
        it->second->Print(out);
    }
}

