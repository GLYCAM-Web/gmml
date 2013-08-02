#include <fstream>
#include <iostream>
#include <iomanip>

#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"

using namespace std;
using namespace gmml;
using namespace PrepFileSpace;

//////////////////////////////////// CONSTRUCTOR ///////////////////////////////////////
PrepFile::PrepFile(const std::string& prep_file)
{
    path_ = prep_file;
    std::ifstream in_file;
    try
    {
        in_file.open(prep_file.c_str());
    }
    catch(...)
    {
        throw PrepFileProcessingException(__LINE__,"File not found");
    }
    Read(in_file);
    in_file.close();            // Close the parameter files
}

///////////////////////////////////////// ACCESSOR ////////////////////////////////////////
PrepFile::ResidueMap& PrepFile::GetResidues()
{
    return residues_;
}

///////////////////////////////////////// FUNCTIONS ///////////////////////////////////////
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

////////////////////////////////// DISPLAY FUNCTION ////////////////////////////////////
void PrepFile::Print(std::ostream& out)
{
    for(ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        out << "**********************************************************************************" << endl;
        it->second->Print(out);
    }
}

