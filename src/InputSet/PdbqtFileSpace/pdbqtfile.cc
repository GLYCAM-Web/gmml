#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "../../../includes/utils.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtfileprocessingexception.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtcompoundcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtremarkcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtrootcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtbranchcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqttorsionaldofcard.hpp"

using PdbqtFileSpace::PdbqtFile;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtFile::PdbqtFile()
{
    path_ = "GMML-Generated";
    models_ = NULL;
}

PdbqtFile::PdbqtFile(const std::string &pdbqt_file)
{
    path_ = pdbqt_file;
    models_ = NULL;

    std::ifstream in_file;
    if(std::ifstream(pdbqt_file.c_str()))
    {
        gmml::log(__LINE__, __FILE__,  gmml::INF,"Opening PDBQT file ...");
//        std::cout << "Opening PDBQT file ..." << std::endl;
        in_file.open(pdbqt_file.c_str());
    }
    else
    {
        throw PdbqtFileProcessingException(__LINE__, "PDBQT file not found");
    }
    if(!Read(in_file))
    {
        throw PdbqtFileProcessingException(__LINE__, "Reading PDBQT file exception");
    }
    gmml::log(__LINE__, __FILE__,  gmml::INF,"End of file");
//    std::cout << "End of file" << std::endl;
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
std::string PdbqtFile::GetPath()
{
    return path_;
}

PdbqtFileSpace::PdbqtModelCard* PdbqtFile::GetModels()
{
    return models_;
}
PdbqtFile::PdbqtResidueAtomsMap PdbqtFile::GetAllAtomsInOrder(std::vector<std::string> &key_order)
{
    PdbqtFile::PdbqtResidueAtomsMap residue_atom_map;
    std::map<std::string, bool> inserted_residues;
    PdbqtFileSpace::PdbqtModelCard::PdbqtModelMap models = models_->GetModels();
    PdbqtModel* model = (*models.begin()).second;
    PdbqtModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbqtAtomCard* atom_card = residue_set->GetAtoms();
    PdbqtAtomCard::PdbqtAtomMap atoms = atom_card->GetAtoms();
    for(PdbqtAtomCard::PdbqtAtomMap::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
    {
        PdbqtAtom* atom = (*it2).second;
        std::string residue_name = atom->GetAtomResidueName();
        char chain_id = atom->GetAtomChainId();
        int sequence_number = atom->GetAtomResidueSequenceNumber();
        char insertion_code = atom->GetAtomInsertionCode();
        char alternate_location = atom->GetAtomAlternateLocation();
        std::stringstream ss;
        ss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
        std::string key = ss.str();
        if(!inserted_residues[key])
        {
            residue_atom_map[key] = new std::vector<PdbqtAtom*>();
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
void PdbqtFile::SetPath(std::string pdbqt_path)
{
    path_ = pdbqt_path;
}

void PdbqtFile::SetModels(PdbqtFileSpace::PdbqtModelCard *models)
{
    models_ = new PdbqtFileSpace::PdbqtModelCard();
    models_ = models;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////
bool PdbqtFile::Read(std::ifstream &in_file)
{
    if(!this->ParseCards(in_file))
        return false;
    return true;
}

bool PdbqtFile::ParseCards(std::ifstream &in_stream)
{

    std::string line,record_name;

    while (getline(in_stream, line)){
        record_name = gmml::Split(line, " ").at(0);
        if(record_name.compare("MODEL") == 0 || record_name.compare("ROOT") == 0 || record_name.compare("BRANCH") == 0 
		|| record_name.compare("ATOM") == 0 || record_name.compare("HETATM") == 0)
        {
	    int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
	    in_stream.seekg(offset, in_stream.cur);//Go back one line
            if(!ParseModelCard(in_stream, line))
                return false;
        }
    }
    

    return true;
}

bool PdbqtFile::ParseModelCard(std::ifstream &stream, std::string &line)
{
    //TODO:catch downstream errors and return false. 
    while (getline(stream, line)){

        if (line.find("MODEL") != std::string::npos || line.find("COMPND") != std::string::npos || line.find("REMARK") != std::string::npos 
	    || line.find("ROOT") != std::string::npos || line.find("ATOM") != std::string::npos || line.find("ENDROOT") != std::string::npos 
	    || line.find("BRANCH") != std::string::npos || line.find("ENDBRANCH") != std::string::npos || line.find("HETATM") != std::string::npos 
	    || line.find("TORSDOF") != std::string::npos || line.find("ENDMDL") != std::string::npos){

	    int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
	    stream.seekg(offset, stream.cur);//Go back one line
	    
    	    // Model card
    	    models_ = new PdbqtFileSpace::PdbqtModelCard(stream);
	    
	}
    }
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
    PdbqtFileSpace::PdbqtModelCard::PdbqtModelMap models = models_->GetModels();
    for(PdbqtFileSpace::PdbqtModelCard::PdbqtModelMap::iterator it = models.begin(); it != models.end(); it++)
    {
        PdbqtModel* model = (*it).second;
        stream << models_->GetRecordName();
        if(model->GetModelSerialNumber() != gmml::iNotSet)
            stream << " " << model->GetModelSerialNumber() << std::endl;
        else
            stream << " " << std::endl;

        PdbqtCompoundCard* compound_card = model->GetModelCompoundCard();
        if(compound_card != NULL)
        {
            if(compound_card->GetValue().compare("") != 0)
                stream << compound_card->GetRecordName() << " " << compound_card->GetValue() << std::endl;
        }

        PdbqtModel::RemarkCardVector remarks = model->GetRemarks();
        for(PdbqtModel::RemarkCardVector::iterator it1 = remarks.begin(); it1 != remarks.end(); it1++)
        {
            PdbqtRemarkCard* remark_card = (*it1);
            stream << remark_card->GetRecordName() << " " << remark_card->GetValue() << std::endl;
        }

        PdbqtModelResidueSet* model_residue_set = model->GetModelResidueSet();
        PdbqtRootCard* root_card = model_residue_set->GetRoots();
        PdbqtAtomCard* atom_card = root_card->GetRootAtoms();
        if(atom_card != NULL)
        {
            stream << root_card->GetRecordName() << std::endl;
            ResolveRootCard(stream, atom_card);
            stream << "ENDROOT" << std::endl;

            PdbqtModelResidueSet::BranchCardVector branches = model_residue_set->GetBranches();
            ResolveBranchCards(stream, branches);
        }
        PdbqtModel::TorsionalDoFCardVector tors_dof_vector = model->GetTorsionalDoFCards();
        for(PdbqtModel::TorsionalDoFCardVector::iterator it1 = tors_dof_vector.begin(); it1 != tors_dof_vector.end(); it1++)
        {
            PdbqtTorsionalDoFCard* tor_dof = (*it1);
            stream << tor_dof->GetRecordName() << " " << tor_dof->GetNumberofTorsionalDoF() << std::endl;
        }
        stream << "ENDMDL" << std::endl;

    }
}
void PdbqtFile::ResolveBranchCards(std::ofstream& stream, PdbqtModelResidueSet::BranchCardVector branches)
{
    for(PdbqtModelResidueSet::BranchCardVector::iterator it = branches.begin(); it != branches.end(); it++)
    {
        PdbqtBranchCard* branch_card = (*it);
        stream << branch_card->GetRecordName() << " ";
        if(branch_card->GetSolidAtomSerialNumber() != gmml::iNotSet)
            stream << std::right << std::setw(3) << branch_card->GetSolidAtomSerialNumber() << " ";
        else
            stream << std::right << std::setw(4) << " ";
        if(branch_card->GetRotatbleAtomSerialNumber() != gmml::iNotSet)
            stream << std::right << std::setw(3) << branch_card->GetRotatbleAtomSerialNumber() << std::endl;
        else
            stream << std::right << std::setw(3) << " " << std::endl;

        PdbqtAtomCard* atom_card = branch_card->GetRotatableAtomSet();
        if(atom_card != NULL)
            ResolveRootCard(stream, atom_card);

        PdbqtModelResidueSet::BranchCardVector sub_branches = branch_card->GetSubBranches();
        if(sub_branches.size() != 0)
            ResolveBranchCards(stream, sub_branches);

        stream << "ENDBRANCH" << " " << std::right << std::setw(3) << branch_card->GetSolidAtomSerialNumber() << " " << std::right << std::setw(3) << branch_card->GetRotatbleAtomSerialNumber() << std::endl;
    }
}

void PdbqtFile::ResolveRootCard(std::ofstream& stream, PdbqtAtomCard* atom_card)
{
    PdbqtAtomCard::PdbqtAtomMap atom_map = atom_card->GetAtoms();
    for(PdbqtAtomCard::PdbqtAtomMap::iterator it = atom_map.begin(); it != atom_map.end(); it++)
    {
        PdbqtAtom* atom = (*it).second;
        stream << std::left << std::setw(6) << atom->GetType();
        if(atom->GetAtomSerialNumber() != gmml::iNotSet)
            stream << std::right << std::setw(5) << atom->GetAtomSerialNumber();
        else
            stream << std::right << std::setw(5) << " ";
        stream << std::left << std::setw(1) << " "
               << std::left << std::setw(4) << atom->GetAtomName();
        if(atom->GetAtomAlternateLocation() == gmml::BLANK_SPACE)
            stream << std::left << std::setw(1) << ' ';
        else
            stream << std::left << std::setw(1) << atom->GetAtomAlternateLocation();
        stream << std::right << std::setw(3) << atom->GetAtomResidueName()
               << std::left << std::setw(1) << " ";
        if(atom->GetAtomChainId() == gmml::BLANK_SPACE)
            stream << std::left << std::setw(1) << ' ';
        else
            stream << std::left << std::setw(1) << atom->GetAtomChainId();
        if(atom->GetAtomResidueSequenceNumber() != gmml::iNotSet)
            stream << std::right << std::setw(4) << atom->GetAtomResidueSequenceNumber();
        else
            stream << std::right << std::setw(4) << " ";
        if(atom->GetAtomInsertionCode() == gmml::BLANK_SPACE)
            stream << std::left << std::setw(1) <<  ' ';
        else
            stream << std::left << std::setw(1) << atom->GetAtomInsertionCode();

        stream << std::left << std::setw(3) << " ";
        if(atom->GetAtomOrthogonalCoordinate().CompareTo(GeometryTopology::Coordinate(gmml::dNotSet, gmml::dNotSet, gmml::dNotSet)) == false)
        {
            stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetX()
                   << std::right << std::setw(8) << std::fixed << std::setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetY()
                   << std::right << std::setw(8) << std::fixed << std::setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetZ();
        }
        else
        {
            stream << std::right << std::setw(8) << " "
                   << std::right << std::setw(8) << " "
                   << std::right << std::setw(8) << " ";
        }
        if(atom->GetAtomOccupancy() != gmml::dNotSet)
            stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << atom->GetAtomOccupancy();
        else
            stream << std::right << std::setw(6) << " ";
        if(atom->GetAtomTempretureFactor() != gmml::dNotSet)
            stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << atom->GetAtomTempretureFactor();
        else
            stream << std::right << std::setw(6) << " ";
        stream << std::right << std::setw(4) << " ";
        if(atom->GetAtomCharge() != gmml::dNotSet)
            stream << std::right << std::setw(6) << std::fixed << std::setprecision(3) << atom->GetAtomCharge();
        else
            stream << std::right << std::setw(6) << " ";
        stream <<  " " << std::left << std::setw(1) <<  atom->GetAtomType() << std::endl;
//        stream << std::setw(1) << " " << std::endl;
    }
}

void PdbqtFile::Print(std::ostream &out)
{
    if(models_ != NULL)
        models_->Print(out);
}
