#include <iostream>
#include <fstream>
#include "../../includes/utils.hpp"
#include "../../includes/FileSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../includes/FileSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../includes/FileSet/PdbqtFileSpace/pdbqtfileprocessingexception.hpp"

using namespace PdbqtFileSpace;
using namespace gmml;
using namespace std;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtFile::PdbqtFile()
{
    path_ = "GMM-Generated";
    models_ = NULL;
}

PdbqtFile::PdbqtFile(const string &pdbqt_file)
{
    path_ = pdbqt_file;
    models_ = NULL;

    ifstream in_file;
    if(ifstream(pdbqt_file.c_str()))
    {
        cout << "Opening PDBQT file ..." << endl;
        in_file.open(pdbqt_file.c_str());
    }
    else
    {
        throw PdbqtFileProcessingException(__LINE__, "PDB file not found");
    }
    if(!Read(in_file))
    {
        throw PdbqtFileProcessingException(__LINE__, "Reading PDB file exception");
    }
    in_file.close();            /// Close the pdbqt files
}

PdbqtFile* PdbqtFile::LoadPdbqtFile()
{
    PdbqtFile* pdbqt = new PdbqtFile();
    return pdbqt;
}

PdbqtFile* PdbqtFile::LoadPdbqtFile(const std::string &pdbqt_file)
{
    PdbqtFile* pdbqt = new PdbqtFile(pdbqt_file);
    return pdbqt;
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbqtFile::GetPath()
{
    return path_;
}

PdbqtModelCard* PdbqtFile::GetModels()
{
    return models_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtFile::SetPath(string pdbqt_path)
{
    path_ = pdbqt_path;
}

void PdbqtFile::SetModels(PdbqtModelCard *models)
{
    models_ = new PdbqtModelCard();
    models_ = models;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////
bool PdbqtFile::Read(ifstream &in_file)
{
    if(!this->ParseCards(in_file))
        return false;
}

bool PdbqtFile::ParseCards(ifstream &in_stream)
{
    string line;

    /// Unable to read file
    if (!getline(in_stream, line))
    {
        cout << "Wrong input file format" << endl;
        throw PdbqtFileProcessingException("Error reading file");
    }
    string record_name = Split(line, " ").at(0);
    if(record_name.compare("MODEL") == 0)
    {
        if(!ParseModelCard(in_stream, line))
            return false;
    }
    else
    {
        cout << "Wrong input file format" << endl;
        return false;
    }
    return true;
}

bool PdbqtFile::ParseModelCard(ifstream &stream, string &line)
{
    stringstream stream_block;
    stream_block << line << endl;
    if(!getline(stream, line))
    {
        cout << "Model card corruption" << endl;
        cout << "Wrong input file format" << endl;
        return false;
    }
    string record_name = Split(line, " ").at(0);
    record_name = Trim(record_name);

    while(record_name.compare("MODEL") == 0 || record_name.compare("COMPND") == 0 || record_name.compare("REMARK") == 0
          || record_name.compare("ROOT") == 0 || record_name.compare("ATOM") == 0 || record_name.compare("ENDROOT") == 0
          || record_name.compare("BRANCH") == 0 || record_name.compare("ENDBRANCH") == 0 || record_name.compare("HETATM") == 0
          || record_name.compare("TORSDOF") == 0 || record_name.compare("ENDMDL") == 0)
    {
        stream_block << line << endl;
        if(getline(stream, line))
        {
            record_name = line.substr(0,6);
            record_name = Trim(record_name);
        }
        else
        {
            cout << "Model card corruption" << endl;
            cout << "Wrong input file format" << endl;
            return false;
        }
    }

    // Model card
//    models_ = new PdbqtModelCard(stream_block);
    return true;
}
