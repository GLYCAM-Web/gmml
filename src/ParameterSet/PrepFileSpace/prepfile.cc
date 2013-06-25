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
    try
    {
        std::ifstream in_file(prep_file.c_str());
        Read(in_file);
        in_file.close();            // Close the parameter files
    }
    catch(...)
    {
        throw PrepFileProcessingException(__LINE__,"File not found");
    }
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

    PrepFileResidue *residue = NULL;
    while ((residue = ProcessResidueSection(in_file)) != NULL)
        residues_[residue->name_] = residue;
}

PrepFileResidue* PrepFile::ProcessResidueSection(ifstream &in_file)
{

}

