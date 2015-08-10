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
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtcompoundcard.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtremarkcard.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtrootcard.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqtbranchcard.hpp"
#include "../../../includes/FileSet/PdbqtFileSpace/pdbqttorsionaldofcard.hpp"

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
        gmml::log(__LINE__, __FILE__,  gmml::INF,"Opening PDBQT file ...");
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
    gmml::log(__LINE__, __FILE__,  gmml::INF,"End of file");
    cout << "End of file" << endl;
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
        gmml::log(__LINE__, __FILE__,  gmml::ERR,"Wrong input file format");
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
        gmml::log(__LINE__, __FILE__,  gmml::ERR,"Wrong input file format");
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
        gmml::log(__LINE__, __FILE__,  gmml::ERR,"Model card corruption");
        cout << "Model card corruption" << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR,"Wrong input file format");
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
                gmml::log(__LINE__, __FILE__,  gmml::ERR,"Model card corruption");
                cout << "Model card corruption" << endl;
                gmml::log(__LINE__, __FILE__,  gmml::ERR,"Wrong input file format");
                cout << "Wrong input file format" << endl;
                return false;
            }
        }
    }
    // Model card
    models_ = new PdbqtModelCard(stream_block);
    return true;
}

void PdbqtFile::Write(const std::string& pdbqt_file)
{
    std::ofstream out_file;
    try
    {
        out_file.open(pdbqt_file.c_str());
    }
    catch(...)
    {
        throw PdbqtFileProcessingException(__LINE__,"File could not be created");
    }
    try
    {
        if(this->models_ != NULL)
        {
            this->ResolveModelCard(out_file);
        }
    }
    catch(...)
    {
        out_file.close();            /// Close the pdbqt file
    }
}

void PdbqtFile::ResolveModelCard(std::ofstream& stream)
{
    PdbqtModelCard::PdbqtModelMap models = models_->GetModels();
    for(PdbqtModelCard::PdbqtModelMap::iterator it = models.begin(); it != models.end(); it++)
    {
        PdbqtModel* model = (*it).second;
        stream << models_->GetRecordName();
        if(model->GetModelSerialNumber() != iNotSet)
            stream << " " << model->GetModelSerialNumber() << endl;
        else
            stream << " " << endl;

        PdbqtCompoundCard* compound_card = model->GetModelCompoundCard();
        if(compound_card != NULL)
        {
            if(compound_card->GetValue().compare("") != 0)
                stream << compound_card->GetRecordName() << " " << compound_card->GetValue() << endl;
        }

        PdbqtModel::RemarkCardVector remarks = model->GetRemarks();
        for(PdbqtModel::RemarkCardVector::iterator it1 = remarks.begin(); it1 != remarks.end(); it1++)
        {
            PdbqtRemarkCard* remark_card = (*it1);
            stream << remark_card->GetRecordName() << " " << remark_card->GetValue() << endl;
        }

        PdbqtModelResidueSet* model_residue_set = model->GetModelResidueSet();
        PdbqtRootCard* root_card = model_residue_set->GetRoots();
        PdbqtAtomCard* atom_card = root_card->GetRootAtoms();
        if(atom_card != NULL)
        {
            stream << root_card->GetRecordName() << endl;
            ResolveRootCard(stream, atom_card);
            stream << "ENDROOT" << endl;

            PdbqtModelResidueSet::BranchCardVector branches = model_residue_set->GetBranches();
            ResolveBranchCards(stream, branches);
        }
        PdbqtModel::TorsionalDoFCardVector tors_dof_vector = model->GetTorsionalDoFCards();
        for(PdbqtModel::TorsionalDoFCardVector::iterator it1 = tors_dof_vector.begin(); it1 != tors_dof_vector.end(); it1++)
        {
            PdbqtTorsionalDoFCard* tor_dof = (*it1);
            stream << tor_dof->GetRecordName() << " " << tor_dof->GetNumberofTorsionalDoF() << endl;
        }
        stream << "ENDMDL" << endl;

    }
}
void PdbqtFile::ResolveBranchCards(ofstream& stream, PdbqtModelResidueSet::BranchCardVector branches)
{
    for(PdbqtModelResidueSet::BranchCardVector::iterator it = branches.begin(); it != branches.end(); it++)
    {
        PdbqtBranchCard* branch_card = (*it);
        stream << branch_card->GetRecordName() << " ";
        if(branch_card->GetSolidAtomSerialNumber() != iNotSet)
            stream << right << setw(3) << branch_card->GetSolidAtomSerialNumber() << " ";
        else
            stream << right << setw(4) << " ";
        if(branch_card->GetRotatbleAtomSerialNumber() != iNotSet)
            stream << right << setw(3) << branch_card->GetRotatbleAtomSerialNumber() << endl;
        else
            stream << right << setw(3) << " " << endl;

        PdbqtAtomCard* atom_card = branch_card->GetRotatableAtomSet();
        if(atom_card != NULL)
            ResolveRootCard(stream, atom_card);

        PdbqtModelResidueSet::BranchCardVector sub_branches = branch_card->GetSubBranches();
        if(sub_branches.size() != 0)
            ResolveBranchCards(stream, sub_branches);

        stream << "ENDBRANCH" << " " << right << setw(3) << branch_card->GetSolidAtomSerialNumber() << " " << right << setw(3) << branch_card->GetRotatbleAtomSerialNumber() << endl;
    }
}

void PdbqtFile::ResolveRootCard(ofstream& stream, PdbqtAtomCard* atom_card)
{
    PdbqtAtomCard::PdbqtAtomMap atom_map = atom_card->GetAtoms();
    for(PdbqtAtomCard::PdbqtAtomMap::iterator it = atom_map.begin(); it != atom_map.end(); it++)
    {
        PdbqtAtom* atom = (*it).second;
        stream << left << setw(6) << atom->GetType();
        if(atom->GetAtomSerialNumber() != iNotSet)
            stream << right << setw(5) << atom->GetAtomSerialNumber();
        else
            stream << right << setw(5) << " ";
        stream << left << setw(1) << " "
               << left << setw(4) << atom->GetAtomName();
        if(atom->GetAtomAlternateLocation() == BLANK_SPACE)
            stream << left << setw(1) << ' ';
        else
            stream << left << setw(1) << atom->GetAtomAlternateLocation();
        stream << right << setw(3) << atom->GetAtomResidueName()
               << left << setw(1) << " ";
        if(atom->GetAtomChainId() == BLANK_SPACE)
            stream << left << setw(1) << ' ';
        else
            stream << left << setw(1) << atom->GetAtomChainId();
        if(atom->GetAtomResidueSequenceNumber() != iNotSet)
            stream << right << setw(4) << atom->GetAtomResidueSequenceNumber();
        else
            stream << right << setw(4) << " ";
        if(atom->GetAtomInsertionCode() == BLANK_SPACE)
            stream << left << setw(1) <<  ' ';
        else
            stream << left << setw(1) << atom->GetAtomInsertionCode();

        stream << left << setw(3) << " ";
        if(atom->GetAtomOrthogonalCoordinate().CompareTo(Geometry::Coordinate(dNotSet, dNotSet, dNotSet)) == false)
        {
            stream << right << setw(8) << fixed << setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetX()
                   << right << setw(8) << fixed << setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetY()
                   << right << setw(8) << fixed << setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetZ();
        }
        else
        {
            stream << right << setw(8) << " "
                   << right << setw(8) << " "
                   << right << setw(8) << " ";
        }
        if(atom->GetAtomOccupancy() != dNotSet)
            stream << right << setw(6) << fixed << setprecision(2) << atom->GetAtomOccupancy();
        else
            stream << right << setw(6) << " ";
        if(atom->GetAtomTempretureFactor() != dNotSet)
            stream << right << setw(6) << fixed << setprecision(2) << atom->GetAtomTempretureFactor();
        else
            stream << right << setw(6) << " ";
        stream << right << setw(4) << " ";
        if(atom->GetAtomCharge() != dNotSet)
            stream << right << setw(6) << fixed << setprecision(3) << atom->GetAtomCharge();
        else
            stream << right << setw(6) << " ";
        stream <<  " " << left << setw(1) <<  atom->GetAtomType() << endl;
//        stream << setw(1) << " " << endl;
    }
}

void PdbqtFile::Print(ostream &out)
{
    if(models_ != NULL)
        models_->Print(out);
}
