#include <iostream>
#include <fstream>
#include "../../../includes/utils.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtfileprocessingexception.hpp"

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
PdbqtFile::PdbqtResidueAtomsMap PdbqtFile::GetAllAtomsInOrder(vector<string> &key_order)
{
    PdbqtFile::PdbqtResidueAtomsMap residue_atom_map;
    map<string, bool> inserted_residues;
    PdbqtModelCard::PdbqtModelMap models = models_->GetModels();
    PdbqtModel* model = (*models.begin()).second;
    PdbqtModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbqtAtomCard* atom_card = residue_set->GetAtoms();
    PdbqtAtomCard::PdbqtAtomMap atoms = atom_card->GetAtoms();
    for(PdbqtAtomCard::PdbqtAtomMap::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
    {
        PdbqtAtom* atom = (*it2).second;
        string residue_name = atom->GetAtomResidueName();
        char chain_id = atom->GetAtomChainId();
        int sequence_number = atom->GetAtomResidueSequenceNumber();
        char insertion_code = atom->GetAtomInsertionCode();
        char alternate_location = atom->GetAtomAlternateLocation();
        stringstream ss;
        ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
        string key = ss.str();
        if(!inserted_residues[key])
        {
            residue_atom_map[key] = new vector<PdbqtAtom*>();
            inserted_residues[key] = true;
            key_order.push_back(key);
        }
        residue_atom_map[key]->push_back(atom);

    }
    return residue_atom_map;
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
            record_name = Split(line, " ").at(0);
            record_name = Trim(record_name);
        }
        else
        {
            if(record_name.compare("ENDMDL") == 0)
                break;
            else
            {
                cout << "Model card corruption" << endl;
                cout << "Wrong input file format" << endl;
                return false;
            }
        }
    }
    // Model card
    models_ = new PdbqtModelCard(stream_block);
    return true;
}

void PdbqtFile::Print(ostream &out)
{
    if(models_ != NULL)
        models_->Print(out);
}
