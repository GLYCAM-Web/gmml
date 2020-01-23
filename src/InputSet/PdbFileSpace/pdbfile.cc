// Author: Alireza Khatamian

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <exception>
#include <cctype>
#include <string>
#include "../../../includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheadercard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbobsoletesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbobsoletecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbtitlesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsplitsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbcaveatsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbcompoundsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbcompoundspecification.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsourcesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsourcecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbkeywordssection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbexperimentaldatasection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbnummodelcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodeltypesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbauthorsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbrevisiondatasection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbrevisiondatacard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsupersededentriessection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsupersededentriescard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbjournalsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbremarksection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbdatabasereferencesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbdatabasereference.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsequenceadvancedsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsequenceadvancedcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbresiduesequencesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbresiduesequencecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbresiduemodificationsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbresiduemodificationcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogensection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogencard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogennamesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogennamecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogensynonymsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogensynonymcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbformulasection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbformulacard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbhelixsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbhelixcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbhelixresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsheetsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsheetcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsheetstrand.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsheetstrandresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbdisulfidebondsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbdisulfideresiduebond.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbdisulfideresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinksection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbcispeptidesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbcispeptidecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsitesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsitecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsiteresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbcrystallographiccard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdboriginxncard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdboriginxnsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbscalencard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbscalensection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmatrixnsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmatrixncard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbconnectsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmastercard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbresidue.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/GeometryTopology/coordinate.hpp"
#include "../../../includes/MolecularModeling/assembly.hpp"
#include "../../../includes/Glycan/ontologyvocabulary.hpp"

using PdbFileSpace::PdbFile;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbFile::PdbFile()
{
    path_ = "GMML-Generated";
    header_ = NULL;
    obsolete_ = NULL;
    title_ = NULL;
    split_ = NULL;
    caveat_ = NULL;
    compound_ = NULL;
    source_ = NULL;
    keywords_ = NULL;
    experimental_data_ = NULL;
    number_of_models_ = NULL;
    model_type_ = NULL;
    author_ = NULL;
    revision_data_ = NULL;
    superseded_entries_ = NULL;
    journal_ = NULL;
    remark_cards_ = NULL;
    database_reference_ = NULL;
    sequence_advanced_ = NULL;
    residues_sequence_ = NULL;
    residue_modification_cards_ = NULL;
    heterogen_cards_ = NULL;
    heterogen_name_cards_ = NULL;
    heterogen_synonym_cards_ = NULL;
    formulas_ = NULL;
    helix_cards_ = NULL;
    sheet_cards_ = NULL;
    disulfide_bonds_ = NULL;
    link_cards_ = NULL;
    cis_peptide_ = NULL;
    site_cards_ = NULL;
    crystallography_ = NULL;
    origins_ = NULL;
    scales_ = NULL;
    matrices_ = NULL;
    models_ = NULL;
    connectivities_ = NULL;
    serial_number_mapping_ = PdbFile::PdbSerialNumberMapping();
    sequence_number_mapping_ = PdbFile::PdbSequenceNumberMapping();
    master_ = NULL;
}

PdbFile::PdbFile(const std::string &pdb_file)
{
    path_ = pdb_file;
    header_ = NULL;
    obsolete_ = NULL;
    title_ = NULL;
    split_ = NULL;
    caveat_ = NULL;
    compound_ = NULL;
    source_ = NULL;
    keywords_ = NULL;
    experimental_data_ = NULL;
    number_of_models_ = NULL;
    model_type_ = NULL;
    author_ = NULL;
    revision_data_ = NULL;
    superseded_entries_ = NULL;
    journal_ = NULL;
    remark_cards_ = NULL;
    database_reference_ = NULL;
    sequence_advanced_ = NULL;
    residues_sequence_ = NULL;
    residue_modification_cards_ = NULL;
    heterogen_cards_ = NULL;
    heterogen_name_cards_ = NULL;
    heterogen_synonym_cards_ = NULL;
    formulas_ = NULL;
    helix_cards_ = NULL;
    sheet_cards_ = NULL;
    disulfide_bonds_ = NULL;
    link_cards_ = NULL;
    cis_peptide_ = NULL;
    site_cards_ = NULL;
    crystallography_ = NULL;
    origins_ = NULL;
    scales_ = NULL;
    matrices_ = NULL;
    models_ = NULL;
    connectivities_ = NULL;
    serial_number_mapping_ = PdbFile::PdbSerialNumberMapping();
    sequence_number_mapping_ = PdbFile::PdbSequenceNumberMapping();
    master_ = NULL;
    std::ifstream in_file;
    if(std::ifstream(pdb_file.c_str()))
    {
        gmml::log(__LINE__, __FILE__,  gmml::INF, "Opening PDB file ...");
        // std::cout << "Opening PDB file ..." << std::endl;
        in_file.open(pdb_file.c_str());
    }
    else
    {
        throw PdbFileProcessingException(__LINE__, "PDB file not found");
    }

    std::string line = "";
    std::string temp = "";
    std::stringstream ss;
    while(!in_file.eof())
    {
        if(!getline(in_file, line))
            break;
        else
        {
            temp = line.substr(0,6);
            temp = gmml::Trim(temp);
            if(temp.find("END") != std::string::npos || temp.compare("END") == 0)
                break;
            else if(!line.empty())
                ss << line << std::endl;
        }
    }
    in_file.close();
    //I think this erases PDb files sometimes, as they only have "END"
    //Dave 2/1/19
    // if(temp.find("END") == std::string::npos || temp.compare("END") != 0)
    // {
    //     std::ofstream out_file;
    //     out_file.open(pdb_file.c_str());
    //     out_file << ss.str() << "END";
    //     out_file.close();
    // }
    // else
    // {
    //     std::ofstream out_file;
    //     out_file.open(pdb_file.c_str());
    //     out_file << ss.str() << temp;
    //     out_file.close();
    // }
    in_file.open(pdb_file.c_str());
    if(!Read(in_file))
    {
        throw PdbFileProcessingException(__LINE__, "Reading PDB file exception");
    }
    if(header_ != NULL)
    {
      std::string PDBname = header_->GetIdentifierCode();
      gmml::log(__LINE__, __FILE__,  gmml::INF, PDBname);
    }

    in_file.close();            /// Close the pdb files
}

PdbFile::PdbFile(std::stringstream& atomStream)
{
  path_ = "";
  header_ = NULL;
  obsolete_ = NULL;
  title_ = NULL;
  split_ = NULL;
  caveat_ = NULL;
  compound_ = NULL;
  source_ = NULL;
  keywords_ = NULL;
  experimental_data_ = NULL;
  number_of_models_ = NULL;
  model_type_ = NULL;
  author_ = NULL;
  revision_data_ = NULL;
  superseded_entries_ = NULL;
  journal_ = NULL;
  remark_cards_ = NULL;
  database_reference_ = NULL;
  sequence_advanced_ = NULL;
  residues_sequence_ = NULL;
  residue_modification_cards_ = NULL;
  heterogen_cards_ = NULL;
  heterogen_name_cards_ = NULL;
  heterogen_synonym_cards_ = NULL;
  formulas_ = NULL;
  helix_cards_ = NULL;
  sheet_cards_ = NULL;
  disulfide_bonds_ = NULL;
  link_cards_ = NULL;
  cis_peptide_ = NULL;
  site_cards_ = NULL;
  crystallography_ = NULL;
  origins_ = NULL;
  scales_ = NULL;
  matrices_ = NULL;
  models_ = NULL;
  connectivities_ = NULL;
  serial_number_mapping_ = PdbFile::PdbSerialNumberMapping();
  sequence_number_mapping_ = PdbFile::PdbSequenceNumberMapping();
  master_ = NULL;

  if(!Read(atomStream))
  {
    throw PdbFileProcessingException(__LINE__, "Reading atom stringstream failed");
  }

}

PdbFile* PdbFile::LoadPdbFile()
{
    PdbFile* pdb = new PdbFile();
    return pdb;
}

PdbFile* PdbFile::LoadPdbFile(const std::string &pdb_file)
{
    PdbFile* pdb = new PdbFile(pdb_file);
    return pdb;
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbFile::GetPath()
{
    return path_;
}

PdbFileSpace::PdbHeaderCard* PdbFile::GetHeader()
{
    return header_;
}

PdbFileSpace::PdbObsoleteSection* PdbFile::GetObsoleteCards()
{
    return obsolete_;
}

PdbFileSpace::PdbTitleSection* PdbFile::GetTitle()
{
    return title_;
}

PdbFileSpace::PdbSplitSection* PdbFile::GetSplit()
{
    return split_;
}

PdbFileSpace::PdbCaveatSection* PdbFile::GetCaveat()
{
    return caveat_;
}

PdbFileSpace::PdbCompoundSection* PdbFile::GetCompound()
{
    return compound_;
}

PdbFileSpace::PdbSourceSection* PdbFile::GetSourceCards()
{
    return source_;
}

PdbFileSpace::PdbKeywordsSection* PdbFile::GetKeywords()
{
    return keywords_;
}

PdbFileSpace::PdbExperimentalDataSection* PdbFile::GetExperimentalData()
{
    return experimental_data_;
}

PdbFileSpace::PdbNumModelCard* PdbFile::GetNumberOfModels()
{
    return number_of_models_;
}

PdbFileSpace::PdbModelTypeSection* PdbFile::GetModelType()
{
    return model_type_;
}

PdbFileSpace::PdbAuthorSection* PdbFile::GetAuthor()
{
    return author_;
}

PdbFileSpace::PdbRevisionDataSection* PdbFile::GetRevisionDataCards()
{
    return revision_data_;
}

PdbFileSpace::PdbSupersededEntriesSection* PdbFile::GetSupersededEntriesCards()
{
    return superseded_entries_;
}

PdbFileSpace::PdbJournalSection* PdbFile::GetJournal()
{
    return journal_;
}

PdbFileSpace::PdbRemarkSection* PdbFile::GetRemarks()
{
    return remark_cards_;
}

PdbFileSpace::PdbDatabaseReferenceSection* PdbFile::GetDatabaseReferences()
{
    return database_reference_;
}

PdbFileSpace::PdbSequenceAdvancedSection* PdbFile::GetSequenceAdvanced()
{
    return sequence_advanced_;
}

PdbFileSpace::PdbResidueSequenceSection* PdbFile::GetResiduesSequence()
{
    return residues_sequence_;
}

PdbFileSpace::PdbResidueModificationSection* PdbFile::GetResidueModification()
{
    return residue_modification_cards_;
}

PdbFileSpace::PdbHeterogenSection* PdbFile::GetHeterogenCards()
{
    return heterogen_cards_;
}

PdbFileSpace::PdbHeterogenNameSection* PdbFile::GetHeterogenNameCards()
{
    return heterogen_name_cards_;
}

PdbFileSpace::PdbHeterogenSynonymSection* PdbFile::GetHeterogenSynonymCards()
{
    return heterogen_synonym_cards_;
}

PdbFileSpace::PdbFormulaSection* PdbFile::GetFormulaCards()
{
    return formulas_;
}

PdbFileSpace::PdbHelixSection* PdbFile::GetHelixCards()
{
    return helix_cards_;
}

PdbFileSpace::PdbSheetSection* PdbFile::GetSheets()
{
    return sheet_cards_;
}

PdbFileSpace::PdbDisulfideBondSection* PdbFile::GetDisulfideBonds()
{
    return disulfide_bonds_;
}

PdbFileSpace::PdbLinkSection* PdbFile::GetResidueLinkCards()
{
    return link_cards_;
}

PdbFileSpace::PdbCISPeptideSection* PdbFile::GetCISPeptide()
{
    return cis_peptide_;
}

PdbFileSpace::PdbSiteSection* PdbFile::GetSites()
{
    return site_cards_;
}

PdbFileSpace::PdbCrystallographicCard* PdbFile::GetCrystallography()
{
    return crystallography_;
}

PdbFileSpace::PdbOriginXnSection* PdbFile::GetOrigins()
{
    return origins_;
}

PdbFileSpace::PdbScaleNSection* PdbFile::GetScales()
{
    return scales_;
}

PdbFileSpace::PdbMatrixNSection* PdbFile::GetMatrices()
{
    return matrices_;
}

PdbFileSpace::PdbModelSection* PdbFile::GetModels()
{
    return models_;
}

PdbFileSpace::PdbConnectSection* PdbFile::GetConnectivities()
{
    return connectivities_;
}

PdbFile::PdbSerialNumberMapping PdbFile::GetSerialNumberMapping()
{
    return serial_number_mapping_;
}

PdbFile::PdbSerialNumberMapping PdbFile::GetSequenceNumberMapping()
{
    return sequence_number_mapping_;
}

PdbFile::PdbPairVectorAtomNamePositionFlag PdbFile::GetAllResidueNames()
{
    PdbPairVectorAtomNamePositionFlag residue_names;
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    PdbModelCard* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
    for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomSection* atom_card = (*it1);
        PdbAtomSection::PdbAtomCardOrderVector atoms = atom_card->GetOrderedAtomCards();
        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            PdbFileSpace::PdbAtomCard* atom = (*it2);
            unsigned int dist = distance(atoms.begin(), it2);
            std::string atom_residue_name = atom->GetAtomResidueName();
            if(dist == 0)
            {
                std::pair<std::string, std::string> pair_residue_position = std::make_pair(atom_residue_name, "S");
                if(find(residue_names.begin(), residue_names.end(), pair_residue_position) == residue_names.end())
                    residue_names.push_back(pair_residue_position);
            }
            else if(dist == atoms.size() - 1)
            {
                std::pair<std::string, std::string> pair_residue_position = std::make_pair(atom_residue_name, "E");
                if(find(residue_names.begin(), residue_names.end(), pair_residue_position) == residue_names.end())
                    residue_names.push_back(pair_residue_position);
            }
            std::pair<std::string, std::string> pair_residue_position = std::make_pair(atom_residue_name, " ");
            if(find(residue_names.begin(), residue_names.end(), pair_residue_position) == residue_names.end())
                residue_names.push_back(pair_residue_position);
        }
    }
    PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
    {
        PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
        PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
        for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
        {
            PdbFileSpace::PdbAtomCard* atom = (*it2);
            std::string atom_residue_name = atom->GetAtomResidueName();
            std::pair<std::string, std::string> pair_residue_position = std::make_pair(atom_residue_name, " ");
            if(find(residue_names.begin(), residue_names.end(), pair_residue_position) == residue_names.end())
            {
                residue_names.push_back(pair_residue_position);
            }
        }
    }
    return residue_names;
}

PdbFile::PdbPairVectorAtomNamePositionFlag PdbFile::GetAllResidueNamesFromAtomSection()
{
    PdbPairVectorAtomNamePositionFlag residue_names;
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    PdbModelCard* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
    for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomSection* atom_card = (*it1);
        PdbAtomSection::PdbAtomCardOrderVector atoms = atom_card->GetOrderedAtomCards();
        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            PdbFileSpace::PdbAtomCard* atom = (*it2);
            unsigned int dist = distance(atoms.begin(), it2);
            std::string atom_residue_name = atom->GetAtomResidueName();
            if(dist == 0)
            {
                std::pair<std::string, std::string> pair_residue_position = std::make_pair(atom_residue_name, "S");
                if(find(residue_names.begin(), residue_names.end(), pair_residue_position) == residue_names.end())
                    residue_names.push_back(pair_residue_position);
            }
            else if(dist == atoms.size() - 1)
            {
                std::pair<std::string, std::string> pair_residue_position = std::make_pair(atom_residue_name, "E");
                if(find(residue_names.begin(), residue_names.end(), pair_residue_position) == residue_names.end())
                    residue_names.push_back(pair_residue_position);
            }
            else
            {
                std::pair<std::string, std::string> pair_residue_position = std::make_pair(atom_residue_name, " ");
                if(find(residue_names.begin(), residue_names.end(), pair_residue_position) == residue_names.end())
                    residue_names.push_back(pair_residue_position);
            }
        }
    }
    return residue_names;

}

PdbFile::PdbResidueVector PdbFile::GetAllResidues()
{
    PdbFile::PdbResidueVector residues;
    std::map<std::string, bool> inserted_residues;
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    PdbModelCard* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
    for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomSection* atom_card = (*it1);
        PdbAtomSection::PdbAtomCardOrderVector atoms = atom_card->GetOrderedAtomCards();
        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            PdbFileSpace::PdbAtomCard* atom = (*it2);
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
                PdbResidue* res = new PdbResidue(residue_name, chain_id, sequence_number, insertion_code, alternate_location);
                residues.push_back(res);
                inserted_residues[key] = true;
            }
        }
    }
    PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
    {
        PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
        PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
        for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
        {
            PdbFileSpace::PdbAtomCard* atom = (*it2);
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
                PdbResidue* res = new PdbResidue(residue_name, chain_id, sequence_number, insertion_code, alternate_location);
                residues.push_back(res);
                inserted_residues[key] = true;
            }
        }
    }
    return residues;
}

PdbFile::PdbResidueVector PdbFile::GetAllResiduesFromAtomSection()
{
    PdbFile::PdbResidueVector residues = PdbFile::PdbResidueVector();
    std::map<std::string, bool> inserted_residues = std::map<std::string, bool>();
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    PdbModelCard* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
    for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomSection* atom_card = (*it1);
        PdbAtomSection::PdbAtomCardOrderVector atoms = atom_card->GetOrderedAtomCards();
        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            PdbFileSpace::PdbAtomCard* atom = (*it2);
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
                PdbResidue* res = new PdbResidue(residue_name, chain_id, sequence_number, insertion_code, alternate_location);
                residues.push_back(res);
                inserted_residues[key] = true;
            }
        }
    }
    return residues;
}

PdbFile::PdbAtomCardVector PdbFile::GetAllAtomsOfResidue(PdbResidue *residue)
{
    std::string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    std::stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    std::string target_key = ss.str();

    PdbAtomCardVector atoms_of_residue;
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    PdbModelCard* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
    for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomSection* atom_card = (*it1);
        PdbAtomSection::PdbAtomCardOrderVector atoms = atom_card->GetOrderedAtomCards();
        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            PdbFileSpace::PdbAtomCard* atom = (*it2);
            std::string residue_name = atom->GetAtomResidueName();
            char chain_id = atom->GetAtomChainId();
            int sequence_number = atom->GetAtomResidueSequenceNumber();
            char insertion_code = atom->GetAtomInsertionCode();
            char alternate_location = atom->GetAtomAlternateLocation();
            std::stringstream sss;
            sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
            std::string key = sss.str();
            if(target_key.compare(key) == 0)
            {
                atoms_of_residue.push_back(atom);
            }
        }
    }
    PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
    {
        PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
        PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
        for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
        {
            PdbFileSpace::PdbAtomCard* atom = (*it2);
            std::string residue_name = atom->GetAtomResidueName();
            char chain_id = atom->GetAtomChainId();
            int sequence_number = atom->GetAtomResidueSequenceNumber();
            char insertion_code = atom->GetAtomInsertionCode();
            char alternate_location = atom->GetAtomAlternateLocation();
            std::stringstream sss;
            sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
            std::string key = sss.str();
            if(target_key.compare(key) == 0)
            {
                atoms_of_residue.push_back(atom);
            }
        }
    }
    return atoms_of_residue;
}

PdbFile::PdbResidueAtomsMap PdbFile::GetAllAtomsOfResidues()
{
    PdbFile::PdbResidueAtomsMap residue_atom_map;
    std::map<std::string, bool> inserted_residues;
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    PdbModelCard* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
    for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomSection* atom_card = (*it1);
        PdbAtomSection::PdbAtomCardOrderVector atoms = atom_card->GetOrderedAtomCards();
        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            PdbFileSpace::PdbAtomCard* atom = (*it2);
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
                residue_atom_map[key] = new std::vector<PdbFileSpace::PdbAtomCard*>();
                inserted_residues[key] = true;
            }
            residue_atom_map[key]->push_back(atom);

        }
    }
    PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
    {
        PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
        PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
        for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
        {
            PdbFileSpace::PdbAtomCard* atom = (*it2);
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
                residue_atom_map[key] = new std::vector<PdbFileSpace::PdbAtomCard*>();
                inserted_residues[key] = true;
            }
            residue_atom_map[key]->push_back(atom);
        }
    }
    return residue_atom_map;
}

PdbFile::PdbResidueAtomsMap PdbFile::GetAllAtomsInOrder(std::vector<std::string>& key_order)
{
    PdbFile::PdbResidueAtomsMap residue_atom_map;
    std::map<std::string, bool> inserted_residues;
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    PdbModelCard* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
    for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
    {
        PdbAtomSection* atom_card = (*it1);
        PdbAtomSection::PdbAtomCardOrderVector atoms = atom_card->GetOrderedAtomCards();
        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
        {
            PdbFileSpace::PdbAtomCard* atom = (*it2);
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
                residue_atom_map[key] = new std::vector<PdbFileSpace::PdbAtomCard*>();
                inserted_residues[key] = true;
                key_order.push_back(key);
            }
            residue_atom_map[key]->push_back(atom);

        }
    }
    PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
    {
        PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
        PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
        for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
        {
            PdbFileSpace::PdbAtomCard* atom = (*it2);
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
                residue_atom_map[key] = new std::vector<PdbFileSpace::PdbAtomCard*>();
                inserted_residues[key] = true;
                key_order.push_back(key);
            }
            residue_atom_map[key]->push_back(atom);
        }
    }
    return residue_atom_map;
}

PdbFileSpace::PdbAtomCard* PdbFile::GetAtomOfResidueByName(PdbResidue *residue, std::string atom_name, PdbFile::PdbResidueAtomsMap residue_atom_map)
{
    std::string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    std::stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    std::string target_key = ss.str();
    PdbAtomCardVector* atoms = residue_atom_map[target_key];

    for(PdbAtomCardVector::iterator it = atoms->begin(); it != atoms->end(); it++)
    {
        PdbFileSpace::PdbAtomCard* atom = (*it);
        if(atom->GetAtomName().compare(atom_name) == 0)
            return atom;
    }
    return NULL;
}

PdbFileSpace::PdbAtomCard* PdbFile::GetAtomOfResidueByName(PdbResidue *residue, std::string atom_name)
{
    PdbAtomCardVector atoms = GetAllAtomsOfResidue(residue);

    for(PdbAtomCardVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        PdbFileSpace::PdbAtomCard* atom = (*it);
        if(atom->GetAtomName().compare(atom_name) == 0)
            return atom;
    }
    return NULL;
}

PdbFileSpace::PdbAtomCard* PdbFile::GetAtomOfResidueByAtomKey(std::string atom_key)
{
    std::vector<std::string> key_tokens = gmml::Split(atom_key, "_");
    int serial_number = gmml::ConvertString<int>(key_tokens.at(1));
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    PdbModelCard* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
    for(PdbModelResidueSet::AtomCardVector::iterator it = atom_cards.begin(); it != atom_cards.end(); it++)
    {
        PdbAtomSection* atom_card = (*it);
        PdbAtomSection::PdbAtomMap atom_map = atom_card->GetAtomCards();
        if(atom_map[serial_number] != NULL)
            return atom_map[serial_number];
    }
    PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it = heterogen_atom_cards.begin(); it != heterogen_atom_cards.end(); it++)
    {
        PdbHeterogenAtomSection* heterogen_atom_card = (*it);
        PdbHeterogenAtomSection::PdbHeterogenAtomCardMap heterogen_atom_map = heterogen_atom_card->GetHeterogenAtomCards();
        if(heterogen_atom_map[serial_number] != NULL)
            return heterogen_atom_map[serial_number];
    }
    return NULL;
}

PdbFileSpace::PdbAtomCard* PdbFile::GetAtomBySerialNumber(int serial_number)
{
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    PdbModelCard* model = (*models.begin()).second;
    PdbModelResidueSet* residue_set = model->GetModelResidueSet();
    PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
    for(PdbModelResidueSet::AtomCardVector::iterator it = atom_cards.begin(); it != atom_cards.end(); it++)
    {
        PdbAtomSection* atom_card = (*it);
        PdbAtomSection::PdbAtomMap atom_map = atom_card->GetAtomCards();
        if(atom_map[serial_number] != NULL)
            return atom_map[serial_number];
    }
    PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it = heterogen_atom_cards.begin(); it != heterogen_atom_cards.end(); it++)
    {
        PdbHeterogenAtomSection* heterogen_atom_card = (*it);
        PdbHeterogenAtomSection::PdbHeterogenAtomCardMap heterogen_atom_map = heterogen_atom_card->GetHeterogenAtomCards();
        if(heterogen_atom_map[serial_number] != NULL)
            return heterogen_atom_map[serial_number];
    }
    return NULL;
}

std::vector<std::string> PdbFile::GetAllAtomNamesOfResidue(PdbResidue *residue, PdbFile::PdbResidueAtomsMap residue_atom_map)
{
    std::string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    std::stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    std::string target_key = ss.str();
    PdbAtomCardVector* atoms = residue_atom_map[target_key];

    std::vector<std::string> atom_names;
    for(PdbAtomCardVector::iterator it = atoms->begin(); it != atoms->end(); it++)
    {
        PdbFileSpace::PdbAtomCard* atom = (*it);
        atom_names.push_back(atom->GetAtomName());
    }
    return atom_names;
}

std::vector<std::string> PdbFile::GetAllAtomNamesOfResidue(PdbResidue *residue)
{
    PdbAtomCardVector atoms = GetAllAtomsOfResidue(residue);

    std::vector<std::string> atom_names;
    for(PdbAtomCardVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        PdbFileSpace::PdbAtomCard* atom = (*it);
        atom_names.push_back(atom->GetAtomName());
    }
    return atom_names;
}

PdbFile::PdbPairVectorTerCardPositions PdbFile::GetAllTerCardPositions(std::vector<std::string> glycam_residue_names)
{
    std::vector<std::pair<char, int> > ter_card_positions = std::vector<std::pair<char, int> >();
    // After residues that has no tails or has more than or equal two tails
    PdbFile::PdbResidueVector residues = this->GetAllResiduesFromAtomSection();
    for(PdbFile::PdbResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        unsigned int dist = distance(residues.begin(), it);
        if(dist != residues.size() - 1)
        {
            PdbResidue* residue = (*it);
            std::string residue_name = residue->GetResidueName();
            if(find(glycam_residue_names.begin(), glycam_residue_names.end(), residue_name) != glycam_residue_names.end())
            {
                // No tail || has more than or equal two tails
                if(residue_name[0] == '0' || isalpha(residue_name[0]))
                {
                    ter_card_positions.push_back(std::make_pair(residue->GetResidueChainId(), residue->GetResidueSequenceNumber() + 1));
                }
            }
        }
    }
    return ter_card_positions;
}

PdbFileSpace::PdbMasterCard* PdbFile::GetMasterCard()
{
  return master_;
}
//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbFile::SetPath(std::string pdb_path)
{
    path_ = pdb_path;
}
void PdbFile::SetHeader(PdbFileSpace::PdbHeaderCard *header)
{
    header_ = header;

}
void PdbFile::SetObsolete(PdbFileSpace::PdbObsoleteSection *obsolete)
{
    obsolete_ = obsolete;
}
void PdbFile::SetTitle(PdbFileSpace::PdbTitleSection *title)
{
    title_ = title;
}
void PdbFile::SetSplit(PdbFileSpace::PdbSplitSection *split)
{
    split_ = split;
}
void PdbFile::SetCaveat(PdbFileSpace::PdbCaveatSection *caveat)
{
    caveat_ = caveat;
}
void PdbFile::SetCompound(PdbFileSpace::PdbCompoundSection *compound)
{
    compound_ = compound;
}
void PdbFile::SetSourceCards(PdbFileSpace::PdbSourceSection *source)
{
    source_ = source;
}
void PdbFile::SetKeywords(PdbFileSpace::PdbKeywordsSection *keywords)
{
    keywords_ = keywords;
}
void PdbFile::SetExperimentalData(PdbFileSpace::PdbExperimentalDataSection *experimental_data)
{
    experimental_data_ = experimental_data;
}
void PdbFile::SetNumberOfModels(PdbFileSpace::PdbNumModelCard *number_of_models)
{
    number_of_models_ = number_of_models;
}
void PdbFile::SetModelType(PdbFileSpace::PdbModelTypeSection *model_type)
{
    model_type_ = model_type;
}
void PdbFile::SetAuthor(PdbFileSpace::PdbAuthorSection *author)
{
    author_ = author;
}
void PdbFile::SetRevisionDataCards(PdbFileSpace::PdbRevisionDataSection *revision_data)
{
    revision_data_ = revision_data;
}
void PdbFile::SetSupersededEntriesCards(PdbFileSpace::PdbSupersededEntriesSection *superseded_entries)
{
    superseded_entries_ = superseded_entries;
}
void PdbFile::SetJournal(PdbFileSpace::PdbJournalSection *journal)
{
    journal_ = journal;
}
void PdbFile::SetRemarks(PdbFileSpace::PdbRemarkSection *remark_cards)
{
    remark_cards_ = remark_cards;
}
void PdbFile::SetDatabaseReferences(PdbFileSpace::PdbDatabaseReferenceSection *database_reference)
{
    database_reference_ = database_reference;
}
void PdbFile::SetSequenceAdvanced(PdbFileSpace::PdbSequenceAdvancedSection *sequence_advanced)
{
    sequence_advanced_ = sequence_advanced;
}
void PdbFile::SetResiduesSequence(PdbFileSpace::PdbResidueSequenceSection *residues_sequence)
{
    residues_sequence_ = residues_sequence;
}
void PdbFile::SetResidueModification(PdbFileSpace::PdbResidueModificationSection *residue_modification_cards)
{
    residue_modification_cards_ = residue_modification_cards;
}
void PdbFile::SetHeterogens(PdbFileSpace::PdbHeterogenSection *heterogen_cards)
{
    heterogen_cards_ = heterogen_cards;
}
void PdbFile::SetHeterogensName(PdbFileSpace::PdbHeterogenNameSection *heterogens_name)
{
    heterogen_name_cards_ = heterogens_name;
}
void PdbFile::SetHeterogenSynonyms(PdbFileSpace::PdbHeterogenSynonymSection *heterogen_synonym_cards)
{
    heterogen_synonym_cards_ = heterogen_synonym_cards;
}
void PdbFile::SetFormulas(PdbFileSpace::PdbFormulaSection *formulas)
{
    formulas_ = formulas;
}
void PdbFile::SetHelixes(PdbFileSpace::PdbHelixSection *helixes)
{
    helix_cards_ = helixes;
}
void PdbFile::SetSheets(PdbFileSpace::PdbSheetSection *sheet_cards)
{
    sheet_cards_ = sheet_cards;
}
void PdbFile::SetDisulfideBonds(PdbFileSpace::PdbDisulfideBondSection *disulfide_bonds)
{
    disulfide_bonds_ = disulfide_bonds;
}
void PdbFile::SetLinks(PdbFileSpace::PdbLinkSection *links)
{
    link_cards_ = links;
}
void PdbFile::SetCISPeptide(PdbFileSpace::PdbCISPeptideSection *cis_peptide)
{
    cis_peptide_ = cis_peptide;
}
void PdbFile::SetSites(PdbFileSpace::PdbSiteSection *site_cards)
{
    site_cards_ = site_cards;
}
void PdbFile::SetCrystallography(PdbFileSpace::PdbCrystallographicCard *crystallography)
{
    crystallography_ = crystallography;
}
void PdbFile::SetOrigins(PdbFileSpace::PdbOriginXnSection *origins)
{
    origins_ = origins;
}
void PdbFile::SetScales(PdbFileSpace::PdbScaleNSection *scales)
{
    scales_ = scales;
}
void PdbFile::SetMatrices(PdbFileSpace::PdbMatrixNSection *matrices)
{
    matrices_ = matrices;
}
void PdbFile::SetModels(PdbFileSpace::PdbModelSection *models)
{
    models_ = models;
}
void PdbFile::SetConnectivities(PdbFileSpace::PdbConnectSection *connectivities)
{
    connectivities_ = connectivities;
}
void PdbFile::SetSerialNumberMapping(PdbSerialNumberMapping serial_number_mapping)
{
    serial_number_mapping_.clear();
    for(PdbSerialNumberMapping::iterator it = serial_number_mapping.begin(); it != serial_number_mapping.end(); it++)
    {
        int source_serial = (*it).first;
        int mapped_serial = (*it).second;
        serial_number_mapping_[source_serial] = mapped_serial;
    }
}
void PdbFile::SetSequenceNumberMapping(PdbSequenceNumberMapping sequence_number_mapping)
{
    sequence_number_mapping_.clear();
    for(PdbSequenceNumberMapping::iterator it = sequence_number_mapping.begin(); it != sequence_number_mapping.end(); it++)
    {
        int mapped_sequence = (*it).first;
        int source_sequence = (*it).second;
        sequence_number_mapping_[mapped_sequence] = source_sequence;
    }
}

void PdbFile::DeleteResidue(PdbResidue *residue)
{
    std::string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    std::stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    std::string target_key = ss.str();

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin();it != models.end(); it++)
    {
        PdbModelCard* model = (*it).second;
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards = PdbModelResidueSet::AtomCardVector();
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms = PdbAtomSection::PdbAtomMap();
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back((atom));
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }

            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards = PdbModelResidueSet::HeterogenAtomCardVector();
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        (*it).second = model;
    }
    models_->SetModels(models);
    this->UpdateConnectCard();
}

void PdbFile::DeleteResidues(PdbResidueVector target_residues)
{
    std::map<std::string, PdbResidue* > pdb_residue_map = std::map<std::string, PdbResidue* >();
    for(PdbResidueVector::iterator it = target_residues.begin(); it != target_residues.end(); it++)
    {
        PdbResidue* residue = (*it);
        std::string target_residue_name = residue->GetResidueName();
        char target_residue_chain_id = residue->GetResidueChainId();
        int target_residue_sequence_number = residue->GetResidueSequenceNumber();
        char target_residue_insertion_code = residue->GetResidueInsertionCode();
        char target_residue_alternate_location = residue->GetResidueAlternateLocation();
        std::stringstream ss;
        ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
        std::string target_key = ss.str();
        pdb_residue_map[target_key] = residue;
    }

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin();it != models.end(); it++)
    {
        PdbModelCard* model = (*it).second;
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards = PdbModelResidueSet::AtomCardVector();
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms = PdbAtomSection::PdbAtomMap();
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_residue_map.find(key) == pdb_residue_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards = PdbModelResidueSet::HeterogenAtomCardVector();
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_residue_map.find(key) == pdb_residue_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        (*it).second = model;
    }
    models_->SetModels(models);
    this->UpdateConnectCard();
}

void PdbFile::DeleteResidueWithTheGivenModelNumber(PdbResidue *residue, int model_number)
{
    std::string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    std::stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    std::string target_key = ss.str();

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards = PdbModelResidueSet::AtomCardVector();
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms = PdbAtomSection::PdbAtomMap();
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards = PdbModelResidueSet::HeterogenAtomCardVector();
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        models[model_number] = model;

        models_->SetModels(models);
        this->UpdateConnectCard();
    }
}

void PdbFile::DeleteResiduesWithTheGivenModelNumber(PdbResidueVector target_residues, int model_number)
{
    std::map<std::string, PdbResidue* > pdb_residue_map = std::map<std::string, PdbResidue* >();
    for(PdbResidueVector::iterator it = target_residues.begin(); it != target_residues.end(); it++)
    {
        PdbResidue* residue = (*it);
        std::string target_residue_name = residue->GetResidueName();
        char target_residue_chain_id = residue->GetResidueChainId();
        int target_residue_sequence_number = residue->GetResidueSequenceNumber();
        char target_residue_insertion_code = residue->GetResidueInsertionCode();
        char target_residue_alternate_location = residue->GetResidueAlternateLocation();
        std::stringstream ss;
        ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
        std::string target_key = ss.str();
        pdb_residue_map[target_key] = residue;
    }

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards = PdbModelResidueSet::AtomCardVector();
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms = PdbAtomSection::PdbAtomMap();
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_residue_map.find(key) == pdb_residue_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards = PdbModelResidueSet::HeterogenAtomCardVector();
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_residue_map.find(key) == pdb_residue_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                    {
                        serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        models[model_number] = model;

        models_->SetModels(models);
        this->UpdateConnectCard();
    }
}

void PdbFile::DeleteAtom(PdbFileSpace::PdbAtomCard* target_atom)
{
    std::string target_residue_name = target_atom->GetAtomResidueName();
    char target_residue_chain_id = target_atom->GetAtomChainId();
    int target_residue_sequence_number = target_atom->GetAtomResidueSequenceNumber();
    char target_residue_insertion_code = target_atom->GetAtomInsertionCode();
    char target_residue_alternate_location = target_atom->GetAtomAlternateLocation();
    std::stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    std::string target_key = ss.str();

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin();it != models.end(); it++)
    {
        PdbModelCard* model = (*it).second;
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;

        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_atoms[serial_number] = atom;
                        updated_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_heterogen_atoms[serial_number] = atom;
                        updated_heterogen_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        (*it).second = model;
    }
    models_->SetModels(models);
    this->UpdateConnectCard();
}

void PdbFile::DeleteAtoms(PdbAtomCardVector target_atoms)
{
    std::map<std::string, PdbFileSpace::PdbAtomCard* > pdb_atom_map = std::map<std::string, PdbFileSpace::PdbAtomCard* >();
    for(PdbAtomCardVector::iterator it = target_atoms.begin(); it != target_atoms.end(); it++)
    {
        PdbFileSpace::PdbAtomCard* target_atom = (*it);
        std::string target_residue_name = target_atom->GetAtomResidueName();
        char target_residue_chain_id = target_atom->GetAtomChainId();
        int target_residue_sequence_number = target_atom->GetAtomResidueSequenceNumber();
        char target_residue_insertion_code = target_atom->GetAtomInsertionCode();
        char target_residue_alternate_location = target_atom->GetAtomAlternateLocation();
        std::stringstream ss;
        ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
        std::string target_key = ss.str();
        pdb_atom_map[target_key] = target_atom;
    }

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin();it != models.end(); it++)
    {
        PdbModelCard* model = (*it).second;
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;

        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_atom_map.find(key) == pdb_atom_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    PdbFileSpace::PdbAtomCard* target_atom = pdb_atom_map[key];
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_atoms[serial_number] = atom;
                        updated_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_atom_map.find(key) == pdb_atom_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    PdbFileSpace::PdbAtomCard* target_atom = pdb_atom_map[key];
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_heterogen_atoms[serial_number] = atom;
                        updated_heterogen_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        (*it).second = model;
    }
    models_->SetModels(models);
    this->UpdateConnectCard();
}

void PdbFile::DeleteAtomWithTheGivenModelNumber(PdbFileSpace::PdbAtomCard* target_atom, int model_number)
{
    std::string target_residue_name = target_atom->GetAtomResidueName();
    char target_residue_chain_id = target_atom->GetAtomChainId();
    int target_residue_sequence_number = target_atom->GetAtomResidueSequenceNumber();
    char target_residue_insertion_code = target_atom->GetAtomInsertionCode();
    char target_residue_alternate_location = target_atom->GetAtomAlternateLocation();
    std::stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    std::string target_key = ss.str();

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;

        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_atoms[serial_number] = atom;
                        updated_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) != 0)
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_heterogen_atoms[serial_number] = atom;
                        updated_heterogen_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        models[model_number] = model;
        models_->SetModels(models);
        this->UpdateConnectCard();
    }
}

void PdbFile::DeleteAtomsWithTheGivenModelNumber(PdbAtomCardVector target_atoms, int model_number)
{
    std::map<std::string, PdbFileSpace::PdbAtomCard* > pdb_atom_map = std::map<std::string, PdbFileSpace::PdbAtomCard* >();
    for(PdbAtomCardVector::iterator it = target_atoms.begin(); it != target_atoms.end(); it++)
    {
        PdbFileSpace::PdbAtomCard* target_atom = (*it);
        std::string target_residue_name = target_atom->GetAtomResidueName();
        char target_residue_chain_id = target_atom->GetAtomChainId();
        int target_residue_sequence_number = target_atom->GetAtomResidueSequenceNumber();
        char target_residue_insertion_code = target_atom->GetAtomInsertionCode();
        char target_residue_alternate_location = target_atom->GetAtomAlternateLocation();
        std::stringstream ss;
        ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
        std::string target_key = ss.str();
        pdb_atom_map[target_key] = target_atom;
    }

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;

        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_atom_map.find(key) == pdb_atom_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_atoms[serial_number] = atom;
                    updated_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    PdbFileSpace::PdbAtomCard* target_atom = pdb_atom_map[key];
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_atoms[serial_number] = atom;
                        updated_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
            serial_number++;
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(pdb_atom_map.find(key) == pdb_atom_map.end())
                {
                    serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                    atom->SetAtomSerialNumber(serial_number);
                    updated_heterogen_atoms[serial_number] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                    serial_number++;
                }
                else
                {
                    PdbFileSpace::PdbAtomCard* target_atom = pdb_atom_map[key];
                    std::string atom_name = atom->GetAtomName();
                    if(atom_name.compare(target_atom->GetAtomName()) != 0)
                    {
                        serial_number_mapping_[atom->GetAtomSerialNumber()] = serial_number;
                        atom->SetAtomSerialNumber(serial_number);
                        updated_heterogen_atoms[serial_number] = atom;
                        updated_heterogen_atoms_vector.push_back(atom);
                        serial_number++;
                    }
                    else
                    {
                        if(serial_number_mapping_.find(atom->GetAtomSerialNumber()) != serial_number_mapping_.end())
                        {
                            serial_number_mapping_.erase(atom->GetAtomSerialNumber());
                        }
                    }
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        models[model_number] = model;
        models_->SetModels(models);
        this->UpdateConnectCard();
    }
}

void PdbFile::UpdateResidueName(PdbResidue *residue, std::string updated_residue_name)
{
    std::string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    std::stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    std::string target_key = ss.str();

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin();it != models.end(); it++)
    {
        PdbModelCard* model = (*it).second;
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) == 0)
                {
                    atom->SetAtomResidueName(updated_residue_name);
                    updated_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_atoms_vector.push_back(atom);
                }
                else
                {
                    updated_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_atoms_vector.push_back(atom);
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) == 0)
                {
                    atom->SetAtomResidueName(updated_residue_name);
                    updated_heterogen_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                }
                else
                {
                    updated_heterogen_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        (*it).second = model;
    }
    models_->SetModels(models);
}

void PdbFile::UpdateResidueNameWithTheGivenModelNumber(PdbResidue *residue, std::string updated_residue_name, int model_number)
{
    std::string target_residue_name = residue->GetResidueName();
    char target_residue_chain_id = residue->GetResidueChainId();
    int target_residue_sequence_number = residue->GetResidueSequenceNumber();
    char target_residue_insertion_code = residue->GetResidueInsertionCode();
    char target_residue_alternate_location = residue->GetResidueAlternateLocation();
    std::stringstream ss;
    ss << target_residue_name << "_" << target_residue_chain_id << "_" << target_residue_sequence_number << "_" << target_residue_insertion_code << "_" << target_residue_alternate_location;
    std::string target_key = ss.str();

    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) == 0)
                {
                    atom->SetAtomResidueName(updated_residue_name);
                    updated_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_atoms_vector.push_back(atom);
                }
                else
                {
                    updated_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_atoms_vector.push_back(atom);
                }
            }
            atom_card->SetAtomCards(updated_atoms);
            atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(atom_card);
        }
        residue_set->SetAtomCards(updated_atom_cards);
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                std::string residue_name = atom->GetAtomResidueName();
                char chain_id = atom->GetAtomChainId();
                int heterogen_atom_sequence_number = atom->GetAtomResidueSequenceNumber();
                char insertion_code = atom->GetAtomInsertionCode();
                char alternate_location = atom->GetAtomAlternateLocation();
                std::stringstream sss;
                sss << residue_name << "_" << chain_id << "_" << heterogen_atom_sequence_number << "_" << insertion_code << "_" << alternate_location;
                std::string key = sss.str();
                if(target_key.compare(key) == 0)
                {
                    atom->SetAtomResidueName(updated_residue_name);
                    updated_heterogen_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                }
                else
                {
                    updated_heterogen_atoms[atom->GetAtomSerialNumber()] = atom;
                    updated_heterogen_atoms_vector.push_back(atom);
                }
            }
            heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(heterogen_atom_card);
        }
        residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);

        model->SetModelResidueSet(residue_set);
        models[model_number] = model;

        models_->SetModels(models);
    }
}

void PdbFile::InsertResidueBefore(PdbAtomSection* residue)
{
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    PdbFileSpace::PdbModelSection::PdbModelCardMap updated_models;
    for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin(); it != models.end(); it++)
    {
        PdbModelCard* model = (*it).second;
        PdbModelCard* updated_model = new PdbModelCard();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        int sequence_number = 1;
        int offset = 0;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection* updated_atom_card = new PdbAtomSection();
            updated_atom_card->SetRecordName(atom_card->GetRecordName());
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();
            bool located = true;

            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_of_residue = residue->GetOrderedAtomCards();
            PdbFileSpace::PdbAtomCard* first_atom_in_residue = (*(ordered_atoms_of_residue.begin()));
            char residue_chain_id = first_atom_in_residue->GetAtomChainId();
            int residue_sequence_number = first_atom_in_residue->GetAtomResidueSequenceNumber();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);

                if(residue_chain_id == atom->GetAtomChainId() && residue_sequence_number == atom->GetAtomResidueSequenceNumber())
                {
                    if(it2 != ordered_atoms.begin())
                    {
                        PdbFileSpace::PdbAtomCard* previous_atom = (*(--it2));
                        it2++;
                        if(previous_atom->GetAtomResidueSequenceNumber() != atom->GetAtomResidueSequenceNumber())
                        {
                            int diff = atom->GetAtomResidueSequenceNumber() - previous_atom->GetAtomResidueSequenceNumber();
                            if(diff == 1)
                                sequence_number++;
                            else
                                sequence_number += diff;
                        }
                    }
                    else
                        sequence_number = atom->GetAtomResidueSequenceNumber() + offset;
                    if(located)
                    {
                        // TODO: update coordinates with respect to it2
                        GeometryTopology::CoordinateVector coordinate_set = GeometryTopology::CoordinateVector();
                        GeometryTopology::Coordinate* base_coordinate = new GeometryTopology::Coordinate((*it2)->GetAtomOrthogonalCoordinate());
                        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                            coordinate_set.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
                        base_coordinate->TranslateAll(coordinate_set, gmml::BOND_LENGTH, -1);
                        base_coordinate->RotateAngularAll(coordinate_set, 180.0, -1);
                        base_coordinate->RotateTorsionalAll(coordinate_set, 180.0, -1);
                        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                        {
                            PdbFileSpace::PdbAtomCard* atom_of_residue = (*it3);
                            int index = distance(ordered_atoms_of_residue.begin(), it3);
                            PdbFileSpace::PdbAtomCard* new_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom_of_residue->GetAtomName(),atom_of_residue->GetAtomAlternateLocation(),
                                                            atom_of_residue->GetAtomResidueName(),atom_of_residue->GetAtomChainId(), sequence_number,
                                                            atom_of_residue->GetAtomInsertionCode(), coordinate_set.at(index),
                                                            atom_of_residue->GetAtomOccupancy(), atom_of_residue->GetAtomTempretureFactor(),
                                                            atom_of_residue->GetAtomElementSymbol(), atom_of_residue->GetAtomCharge(), atom->GetAlternateAtomCards());
                            updated_atoms[serial_number] = new_atom;
                            updated_atoms_vector.push_back(new_atom);
                            serial_number++;
                        }
                        sequence_number_mapping_[sequence_number] = gmml::iNotSet;
                        sequence_number++;
                        offset++;
                        located = false;
                    }
                    if(!located)
                    {
                        PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                            atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                            atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(),
                                                            atom->GetAlternateAtomCards());
                        updated_atoms[serial_number] = updated_atom;
                        updated_atoms_vector.push_back(updated_atom);
                        serial_number++;
                        if(sequence_number != atom->GetAtomResidueSequenceNumber())
                            sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                    }
                }
                else
                {
                    if(it2 != ordered_atoms.begin())
                    {

                        PdbFileSpace::PdbAtomCard* previous_atom = (*(--it2));
                        it2++;
                        if(previous_atom->GetAtomResidueSequenceNumber() != atom->GetAtomResidueSequenceNumber())
                        {
                            int diff = atom->GetAtomResidueSequenceNumber() - previous_atom->GetAtomResidueSequenceNumber();
                            if(diff == 1)
                                sequence_number++;
                            else
                                sequence_number += diff;
                        }
                    }
                    else
                        sequence_number = atom->GetAtomResidueSequenceNumber() + offset;
                    PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(), atom->GetAlternateAtomCards());
                    updated_atoms[serial_number] = updated_atom;
                    updated_atoms_vector.push_back(updated_atom);
                    serial_number++;
                    if(sequence_number != atom->GetAtomResidueSequenceNumber())
                        sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                }
            }
            updated_atom_card->SetAtomCards(updated_atoms);
            updated_atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(updated_atom_card);
            serial_number++;
        }
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection* updated_heterogen_atom_card = new PdbHeterogenAtomSection();
            updated_heterogen_atom_card->SetRecordName(heterogen_atom_card->GetRecordName());
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* heterogen_atom = (*it2);
                PdbFileSpace::PdbAtomCard* updated_heterogen_atom = new PdbFileSpace::PdbAtomCard(serial_number, heterogen_atom->GetAtomName(),heterogen_atom->GetAtomAlternateLocation(),
                                                              heterogen_atom->GetAtomResidueName(), heterogen_atom->GetAtomChainId(), heterogen_atom->GetAtomResidueSequenceNumber(),
                                                              heterogen_atom->GetAtomInsertionCode(), heterogen_atom->GetAtomOrthogonalCoordinate(),
                                                              heterogen_atom->GetAtomOccupancy(), heterogen_atom->GetAtomTempretureFactor(),
                                                              heterogen_atom->GetAtomElementSymbol(), heterogen_atom->GetAtomCharge(), heterogen_atom->GetAlternateAtomCards());
                updated_heterogen_atoms[serial_number] = updated_heterogen_atom;
                updated_heterogen_atoms_vector.push_back(updated_heterogen_atom);
                serial_number++;
            }
            updated_heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            updated_heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(updated_heterogen_atom_card);
        }
        updated_residue_set->SetAtomCards(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;
    }
    models_->SetModels(updated_models);
}

void PdbFile::InsertResidueBeforeWithTheGivenModelNumber(PdbAtomSection* residue, int model_number)
{
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard* model = models[model_number];
        PdbFileSpace::PdbModelSection::PdbModelCardMap updated_models;
        PdbModelCard* updated_model = new PdbModelCard();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        int sequence_number = 1;
        int offset = 0;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection* updated_atom_card = new PdbAtomSection();
            updated_atom_card->SetRecordName(atom_card->GetRecordName());
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();
            bool located = true;

            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_of_residue = residue->GetOrderedAtomCards();
            PdbFileSpace::PdbAtomCard* first_atom_in_residue = (*(ordered_atoms_of_residue.begin()));
            char residue_chain_id = first_atom_in_residue->GetAtomChainId();
            int residue_sequence_number = first_atom_in_residue->GetAtomResidueSequenceNumber();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);

                if(residue_chain_id == atom->GetAtomChainId() && residue_sequence_number == atom->GetAtomResidueSequenceNumber())
                {
                    if(it2 != ordered_atoms.begin())
                    {
                        PdbFileSpace::PdbAtomCard* previous_atom = (*(--it2));
                        it2++;
                        if(previous_atom->GetAtomResidueSequenceNumber() != atom->GetAtomResidueSequenceNumber())
                        {
                            int diff = atom->GetAtomResidueSequenceNumber() - previous_atom->GetAtomResidueSequenceNumber();
                            if(diff == 1)
                                sequence_number++;
                            else
                                sequence_number += diff;
                        }
                    }
                    else
                        sequence_number = atom->GetAtomResidueSequenceNumber() + offset;
                    if(located)
                    {
                        // TODO: update coordinates with respect to it2
                        GeometryTopology::CoordinateVector coordinate_set = GeometryTopology::CoordinateVector();
                        GeometryTopology::Coordinate* base_coordinate = new GeometryTopology::Coordinate((*it2)->GetAtomOrthogonalCoordinate());
                        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                            coordinate_set.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
                        base_coordinate->TranslateAll(coordinate_set, gmml::BOND_LENGTH, -1);
                        base_coordinate->RotateAngularAll(coordinate_set, 180.0, -1);
                        base_coordinate->RotateTorsionalAll(coordinate_set, 180.0, -1);
                        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                        {
                            PdbFileSpace::PdbAtomCard* atom_of_residue = (*it3);
                            int index = distance(ordered_atoms_of_residue.begin(), it3);
                            PdbFileSpace::PdbAtomCard* new_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom_of_residue->GetAtomName(),atom_of_residue->GetAtomAlternateLocation(),
                                                            atom_of_residue->GetAtomResidueName(),atom_of_residue->GetAtomChainId(), sequence_number,
                                                            atom_of_residue->GetAtomInsertionCode(), coordinate_set.at(index),
                                                            atom_of_residue->GetAtomOccupancy(), atom_of_residue->GetAtomTempretureFactor(),
                                                            atom_of_residue->GetAtomElementSymbol(), atom_of_residue->GetAtomCharge(), atom->GetAlternateAtomCards());
                            updated_atoms[serial_number] = new_atom;
                            updated_atoms_vector.push_back(new_atom);
                            serial_number++;
                        }
                        sequence_number_mapping_[sequence_number] = gmml::iNotSet;
                        sequence_number++;
                        offset++;
                        located = false;
                    }
                    if(!located)
                    {
                        PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                            atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                            atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(),
                                                            atom->GetAlternateAtomCards());
                        updated_atoms[serial_number] = updated_atom;
                        updated_atoms_vector.push_back(updated_atom);
                        serial_number++;
                        if(sequence_number != atom->GetAtomResidueSequenceNumber())
                            sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                    }
                }
                else
                {
                    if(it2 != ordered_atoms.begin())
                    {

                        PdbFileSpace::PdbAtomCard* previous_atom = (*(--it2));
                        it2++;
                        if(previous_atom->GetAtomResidueSequenceNumber() != atom->GetAtomResidueSequenceNumber())
                        {
                            int diff = atom->GetAtomResidueSequenceNumber() - previous_atom->GetAtomResidueSequenceNumber();
                            if(diff == 1)
                                sequence_number++;
                            else
                                sequence_number += diff;
                        }
                    }
                    else
                        sequence_number = atom->GetAtomResidueSequenceNumber() + offset;
                    PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(), atom->GetAlternateAtomCards());
                    updated_atoms[serial_number] = updated_atom;
                    updated_atoms_vector.push_back(updated_atom);
                    serial_number++;
                    if(sequence_number != atom->GetAtomResidueSequenceNumber())
                        sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                }
            }
            updated_atom_card->SetAtomCards(updated_atoms);
            updated_atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(updated_atom_card);
            serial_number++;
        }
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection* updated_heterogen_atom_card = new PdbHeterogenAtomSection();
            updated_heterogen_atom_card->SetRecordName(heterogen_atom_card->GetRecordName());
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector  updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* heterogen_atom = (*it2);
                PdbFileSpace::PdbAtomCard* updated_heterogen_atom = new PdbFileSpace::PdbAtomCard(serial_number, heterogen_atom->GetAtomName(),heterogen_atom->GetAtomAlternateLocation(),
                                                              heterogen_atom->GetAtomResidueName(), heterogen_atom->GetAtomChainId(), heterogen_atom->GetAtomResidueSequenceNumber(),
                                                              heterogen_atom->GetAtomInsertionCode(), heterogen_atom->GetAtomOrthogonalCoordinate(),
                                                              heterogen_atom->GetAtomOccupancy(), heterogen_atom->GetAtomTempretureFactor(),
                                                              heterogen_atom->GetAtomElementSymbol(), heterogen_atom->GetAtomCharge(),
                                                              heterogen_atom->GetAlternateAtomCards());
                updated_heterogen_atoms[serial_number] = updated_heterogen_atom;
                updated_heterogen_atoms_vector.push_back(updated_heterogen_atom);
                serial_number++;
            }
            updated_heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            updated_heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(updated_heterogen_atom_card);
        }
        updated_residue_set->SetAtomCards(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;

        models_->SetModels(updated_models);
    }
}

void PdbFile::InsertResidueAfter(PdbAtomSection* residue)
{
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    PdbFileSpace::PdbModelSection::PdbModelCardMap updated_models;
    for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin(); it != models.end(); it++)
    {
        PdbModelCard* model = (*it).second;
        PdbModelCard* updated_model = new PdbModelCard();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        int sequence_number = 1;
        int offset = 0;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection* updated_atom_card = new PdbAtomSection();
            updated_atom_card->SetRecordName(atom_card->GetRecordName());
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();
            bool located = false;

            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_of_residue = residue->GetOrderedAtomCards();
            PdbFileSpace::PdbAtomCard* first_atom_in_residue = (*(ordered_atoms_of_residue.begin()));
            char residue_chain_id = first_atom_in_residue->GetAtomChainId();
            int residue_sequence_number = first_atom_in_residue->GetAtomResidueSequenceNumber();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);

                if(residue_chain_id == atom->GetAtomChainId() && residue_sequence_number == atom->GetAtomResidueSequenceNumber())
                {
                    if(it2 != ordered_atoms.begin())
                    {
                        PdbFileSpace::PdbAtomCard* previous_atom = (*(--it2));
                        it2++;
                        if(previous_atom->GetAtomResidueSequenceNumber() != atom->GetAtomResidueSequenceNumber())
                        {
                            int diff = atom->GetAtomResidueSequenceNumber() - previous_atom->GetAtomResidueSequenceNumber();
                            if(diff == 1)
                                sequence_number++;
                            else
                                sequence_number += diff;
                        }
                    }
                    else
                    {
                        sequence_number = atom->GetAtomResidueSequenceNumber() + offset;
                    }
                    PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(), atom->GetAlternateAtomCards());
                    updated_atoms[serial_number] = updated_atom;
                    updated_atoms_vector.push_back(updated_atom);
                    serial_number++;
                    if(sequence_number != atom->GetAtomResidueSequenceNumber())
                        sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                    located = true;
                }
                else
                {
                    if(located)
                    {
                        // TODO: update coordinates with respect to it2
                        GeometryTopology::CoordinateVector coordinate_set = GeometryTopology::CoordinateVector();
                        GeometryTopology::Coordinate* base_coordinate = new GeometryTopology::Coordinate((*it2)->GetAtomOrthogonalCoordinate());
                        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                            coordinate_set.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
                        base_coordinate->TranslateAll(coordinate_set, gmml::BOND_LENGTH, 1);
                        base_coordinate->RotateAngularAll(coordinate_set, 180.0, 1);
                        base_coordinate->RotateTorsionalAll(coordinate_set, 180.0, 1);
                        sequence_number++;
                        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                        {
                            PdbFileSpace::PdbAtomCard* atom_of_residue = (*it3);
                            int index = distance(ordered_atoms_of_residue.begin(), it3);
                            PdbFileSpace::PdbAtomCard* new_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom_of_residue->GetAtomName(),atom_of_residue->GetAtomAlternateLocation(),
                                                            atom_of_residue->GetAtomResidueName(),atom_of_residue->GetAtomChainId(), sequence_number,
                                                            atom_of_residue->GetAtomInsertionCode(), coordinate_set.at(index),
                                                            atom_of_residue->GetAtomOccupancy(), atom_of_residue->GetAtomTempretureFactor(),
                                                            atom_of_residue->GetAtomElementSymbol(), atom_of_residue->GetAtomCharge(), atom->GetAlternateAtomCards());
                            updated_atoms[serial_number] = new_atom;
                            updated_atoms_vector.push_back(new_atom);
                            serial_number++;
                        }
                        sequence_number_mapping_[sequence_number] = gmml::iNotSet;
                        sequence_number++;
                        offset++;
                        located = false;
                    }
                    if(!located)
                    {
                        if(it2 != ordered_atoms.begin())
                        {
                            PdbFileSpace::PdbAtomCard* previous_atom = (*(--it2));
                            it2++;
                            if(previous_atom->GetAtomResidueSequenceNumber() != atom->GetAtomResidueSequenceNumber())
                            {
                                int diff = atom->GetAtomResidueSequenceNumber() - previous_atom->GetAtomResidueSequenceNumber();
                                if(diff == 1)
                                    sequence_number++;
                                else
                                    sequence_number += diff;
                            }
                        }
                        else
                        {
                            sequence_number = atom->GetAtomResidueSequenceNumber() + offset;
                        }
                        PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                            atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                            atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(),
                                                            atom->GetAlternateAtomCards());
                        updated_atoms[serial_number] = updated_atom;
                        updated_atoms_vector.push_back(updated_atom);
                        serial_number++;
                        if(sequence_number != atom->GetAtomResidueSequenceNumber())
                            sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                    }
                }
                if(it2 == --ordered_atoms.end() && located)
                {
                    // TODO: update coordinates with respect to it2
                    GeometryTopology::CoordinateVector coordinate_set = GeometryTopology::CoordinateVector();
                    GeometryTopology::Coordinate* base_coordinate = new GeometryTopology::Coordinate((*it2)->GetAtomOrthogonalCoordinate());
                    for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                        coordinate_set.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
                    base_coordinate->TranslateAll(coordinate_set, gmml::BOND_LENGTH, 1);
                    base_coordinate->RotateAngularAll(coordinate_set, 180.0, 1);
                    base_coordinate->RotateTorsionalAll(coordinate_set, 180.0, 1);
                    sequence_number++;
                    for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                    {
                        PdbFileSpace::PdbAtomCard* atom_of_residue = (*it3);
                        int index = distance(ordered_atoms_of_residue.begin(), it3);
                        PdbFileSpace::PdbAtomCard* new_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom_of_residue->GetAtomName(),atom_of_residue->GetAtomAlternateLocation(),
                                                        atom_of_residue->GetAtomResidueName(),atom_of_residue->GetAtomChainId(), sequence_number,
                                                        atom_of_residue->GetAtomInsertionCode(), coordinate_set.at(index),
                                                        atom_of_residue->GetAtomOccupancy(), atom_of_residue->GetAtomTempretureFactor(),
                                                        atom_of_residue->GetAtomElementSymbol(), atom_of_residue->GetAtomCharge(), atom_of_residue->GetAlternateAtomCards());
                        updated_atoms[serial_number] = new_atom;
                        updated_atoms_vector.push_back(new_atom);
                        serial_number++;
                    }
                    sequence_number_mapping_[sequence_number] = gmml::iNotSet;
                    located = false;
                    sequence_number++;
                    offset++;
                }
            }
            updated_atom_card->SetAtomCards(updated_atoms);
            updated_atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(updated_atom_card);
            serial_number++;
        }
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection* updated_heterogen_atom_card = new PdbHeterogenAtomSection();
            updated_heterogen_atom_card->SetRecordName(heterogen_atom_card->GetRecordName());
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* heterogen_atom = (*it2);
                PdbFileSpace::PdbAtomCard* updated_heterogen_atom = new PdbFileSpace::PdbAtomCard(serial_number, heterogen_atom->GetAtomName(),heterogen_atom->GetAtomAlternateLocation(),
                                                              heterogen_atom->GetAtomResidueName(), heterogen_atom->GetAtomChainId(), heterogen_atom->GetAtomResidueSequenceNumber(),
                                                              heterogen_atom->GetAtomInsertionCode(), heterogen_atom->GetAtomOrthogonalCoordinate(),
                                                              heterogen_atom->GetAtomOccupancy(), heterogen_atom->GetAtomTempretureFactor(),
                                                              heterogen_atom->GetAtomElementSymbol(), heterogen_atom->GetAtomCharge(), heterogen_atom->GetAlternateAtomCards());
                updated_heterogen_atoms[serial_number] = updated_heterogen_atom;
                updated_heterogen_atoms_vector.push_back(updated_heterogen_atom);
                serial_number++;
            }
            updated_heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            updated_heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(updated_heterogen_atom_card);
        }
        updated_residue_set->SetAtomCards(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;
    }
    models_->SetModels(updated_models);
}

void PdbFile::InsertResidueAfterWithTheGivenModelNumber(PdbAtomSection* residue, int model_number)
{
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard* model = models[model_number];
        PdbFileSpace::PdbModelSection::PdbModelCardMap updated_models;
        PdbModelCard* updated_model = new PdbModelCard();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        int sequence_number = 1;
        int offset = 0;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection* updated_atom_card = new PdbAtomSection();
            updated_atom_card->SetRecordName(atom_card->GetRecordName());
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap updated_atoms;
            PdbAtomSection::PdbAtomCardOrderVector updated_atoms_vector = PdbAtomSection::PdbAtomCardOrderVector();
            bool located = false;

            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_of_residue = residue->GetOrderedAtomCards();
            PdbFileSpace::PdbAtomCard* first_atom_in_residue = (*(ordered_atoms_of_residue.begin()));
            char residue_chain_id = first_atom_in_residue->GetAtomChainId();
            int residue_sequence_number = first_atom_in_residue->GetAtomResidueSequenceNumber();

            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);

                if(residue_chain_id == atom->GetAtomChainId() && residue_sequence_number == atom->GetAtomResidueSequenceNumber())
                {
                    if(it2 != ordered_atoms.begin())
                    {
                        PdbFileSpace::PdbAtomCard* previous_atom = (*(--it2));
                        it2++;
                        if(previous_atom->GetAtomResidueSequenceNumber() != atom->GetAtomResidueSequenceNumber())
                        {
                            int diff = atom->GetAtomResidueSequenceNumber() - previous_atom->GetAtomResidueSequenceNumber();
                            if(diff == 1)
                                sequence_number++;
                            else
                                sequence_number += diff;
                        }
                    }
                    else
                    {
                        sequence_number = atom->GetAtomResidueSequenceNumber() + offset;
                    }
                    PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(), atom->GetAlternateAtomCards());
                    updated_atoms[serial_number] = updated_atom;
                    updated_atoms_vector.push_back(updated_atom);
                    serial_number++;
                    if(sequence_number != atom->GetAtomResidueSequenceNumber())
                        sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                    located = true;
                }
                else
                {
                    if(located)
                    {
                        // TODO: update coordinates with respect to it2
                        GeometryTopology::CoordinateVector coordinate_set = GeometryTopology::CoordinateVector();
                        GeometryTopology::Coordinate* base_coordinate = new GeometryTopology::Coordinate((*it2)->GetAtomOrthogonalCoordinate());
                        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                            coordinate_set.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
                        base_coordinate->TranslateAll(coordinate_set, gmml::BOND_LENGTH, 1);
                        base_coordinate->RotateAngularAll(coordinate_set, 180.0, 1);
                        base_coordinate->RotateTorsionalAll(coordinate_set, 180.0, 1);
                        sequence_number++;
                        for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                        {
                            PdbFileSpace::PdbAtomCard* atom_of_residue = (*it3);
                            int index = distance(ordered_atoms_of_residue.begin(), it3);
                            PdbFileSpace::PdbAtomCard* new_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom_of_residue->GetAtomName(),atom_of_residue->GetAtomAlternateLocation(),
                                                            atom_of_residue->GetAtomResidueName(),atom_of_residue->GetAtomChainId(), sequence_number,
                                                            atom_of_residue->GetAtomInsertionCode(), coordinate_set.at(index),
                                                            atom_of_residue->GetAtomOccupancy(), atom_of_residue->GetAtomTempretureFactor(),
                                                            atom_of_residue->GetAtomElementSymbol(), atom_of_residue->GetAtomCharge(), atom->GetAlternateAtomCards());
                            updated_atoms[serial_number] = new_atom;
                            updated_atoms_vector.push_back(new_atom);
                            serial_number++;
                        }
                        sequence_number++;
                        offset++;
                        located = false;
                    }
                    if(!located)
                    {
                        if(it2 != ordered_atoms.begin())
                        {
                            PdbFileSpace::PdbAtomCard* previous_atom = (*(--it2));
                            it2++;
                            if(previous_atom->GetAtomResidueSequenceNumber() != atom->GetAtomResidueSequenceNumber())
                            {
                                int diff = atom->GetAtomResidueSequenceNumber() - previous_atom->GetAtomResidueSequenceNumber();
                                if(diff == 1)
                                    sequence_number++;
                                else
                                    sequence_number += diff;
                            }
                        }
                        else
                        {
                            sequence_number = atom->GetAtomResidueSequenceNumber() + offset;
                        }
                        PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                            atom->GetAtomChainId(), sequence_number, atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                            atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(),
                                                            atom->GetAlternateAtomCards());
                        updated_atoms[serial_number] = updated_atom;
                        updated_atoms_vector.push_back(updated_atom);
                        serial_number++;
                        if(sequence_number != atom->GetAtomResidueSequenceNumber())
                            sequence_number_mapping_[sequence_number] = atom->GetAtomResidueSequenceNumber();
                    }
                }
                if(it2 == --ordered_atoms.end() && located)
                {
                    // TODO: update coordinates with respect to it2
                    GeometryTopology::CoordinateVector coordinate_set = GeometryTopology::CoordinateVector();
                    GeometryTopology::Coordinate* base_coordinate = new GeometryTopology::Coordinate((*it2)->GetAtomOrthogonalCoordinate());
                    for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                        coordinate_set.push_back(new GeometryTopology::Coordinate((*it3)->GetAtomOrthogonalCoordinate()));
                    base_coordinate->TranslateAll(coordinate_set, gmml::BOND_LENGTH, 1);
                    base_coordinate->RotateAngularAll(coordinate_set, 180.0, 1);
                    base_coordinate->RotateTorsionalAll(coordinate_set, 180.0, 1);
                    sequence_number++;
                    for(PdbAtomSection::PdbAtomCardOrderVector::iterator it3 = ordered_atoms_of_residue.begin(); it3 != ordered_atoms_of_residue.end(); it3++)
                    {
                        PdbFileSpace::PdbAtomCard* atom_of_residue = (*it3);
                        int index = distance(ordered_atoms_of_residue.begin(), it3);
                        PdbFileSpace::PdbAtomCard* new_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom_of_residue->GetAtomName(),atom_of_residue->GetAtomAlternateLocation(),
                                                        atom_of_residue->GetAtomResidueName(),atom_of_residue->GetAtomChainId(), sequence_number,
                                                        atom_of_residue->GetAtomInsertionCode(), coordinate_set.at(index),
                                                        atom_of_residue->GetAtomOccupancy(), atom_of_residue->GetAtomTempretureFactor(),
                                                        atom_of_residue->GetAtomElementSymbol(), atom_of_residue->GetAtomCharge(), atom->GetAlternateAtomCards());
                        updated_atoms[serial_number] = new_atom;
                        updated_atoms_vector.push_back(new_atom);
                        serial_number++;
                    }
                    sequence_number_mapping_[sequence_number] = gmml::iNotSet;
                    located = false;
                    sequence_number++;
                    offset++;
                }
            }
            updated_atom_card->SetAtomCards(updated_atoms);
            updated_atom_card->SetOrderedAtomCards(updated_atoms_vector);
            updated_atom_cards.push_back(updated_atom_card);
            serial_number++;
        }
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        PdbModelResidueSet::HeterogenAtomCardVector updated_heterogen_atom_cards;
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection* updated_heterogen_atom_card = new PdbHeterogenAtomSection();
            updated_heterogen_atom_card->SetRecordName(heterogen_atom_card->GetRecordName());
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap updated_heterogen_atoms;
            PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector updated_heterogen_atoms_vector = PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* heterogen_atom = (*it2);
                PdbFileSpace::PdbAtomCard* updated_heterogen_atom = new PdbFileSpace::PdbAtomCard(serial_number, heterogen_atom->GetAtomName(),heterogen_atom->GetAtomAlternateLocation(),
                                                              heterogen_atom->GetAtomResidueName(), heterogen_atom->GetAtomChainId(), heterogen_atom->GetAtomResidueSequenceNumber(),
                                                              heterogen_atom->GetAtomInsertionCode(), heterogen_atom->GetAtomOrthogonalCoordinate(),
                                                              heterogen_atom->GetAtomOccupancy(), heterogen_atom->GetAtomTempretureFactor(),
                                                              heterogen_atom->GetAtomElementSymbol(), heterogen_atom->GetAtomCharge(),
                                                              heterogen_atom->GetAlternateAtomCards());
                updated_heterogen_atoms[serial_number] = updated_heterogen_atom;
                updated_heterogen_atoms_vector.push_back(updated_heterogen_atom);
                serial_number++;
            }
            updated_heterogen_atom_card->SetHeterogenAtoms(updated_heterogen_atoms);
            updated_heterogen_atom_card->SetOrderedHeterogenAtoms(updated_heterogen_atoms_vector);
            updated_heterogen_atom_cards.push_back(updated_heterogen_atom_card);
        }
        updated_residue_set->SetAtomCards(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(updated_heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;
        models_->SetModels(updated_models);
    }
}
void PdbFile::SplitAtomCardOfModelCard(char split_point_chain_id, int split_point_sequence_number)
{
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    PdbFileSpace::PdbModelSection::PdbModelCardMap updated_models;
    for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin(); it != models.end(); it++)
    {
        PdbModelCard* model = (*it).second;
        PdbModelCard* updated_model = new PdbModelCard();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection* updated_atom_card_first_part = new PdbAtomSection();
            PdbAtomSection* updated_atom_card_second_part = new PdbAtomSection();
            updated_atom_card_first_part->SetRecordName(atom_card->GetRecordName());
            updated_atom_card_second_part->SetRecordName(atom_card->GetRecordName());
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap atoms_first_part;
            PdbAtomSection::PdbAtomMap atoms_second_part;
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_first_part = PdbAtomSection::PdbAtomCardOrderVector();
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_second_part = PdbAtomSection::PdbAtomCardOrderVector();

            bool located = false;
            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                char residue_chain_id = atom->GetAtomChainId();
                int residue_sequence_number = atom->GetAtomResidueSequenceNumber();

                if(residue_chain_id == split_point_chain_id && residue_sequence_number == split_point_sequence_number)
                {
                    located = true;
                }
                if(located)
                {
                    serial_number++;
                    PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), atom->GetAtomResidueSequenceNumber(), atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(), atom->GetAlternateAtomCards());
                    atoms_second_part[serial_number] = updated_atom;
                    ordered_atoms_second_part.push_back(updated_atom);
                }
                if(!located)
                {
                    PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), atom->GetAtomResidueSequenceNumber(), atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(), atom->GetAlternateAtomCards());
                    atoms_first_part[serial_number] = updated_atom;
                    ordered_atoms_first_part.push_back(updated_atom);
                    serial_number++;
                }
            }
            updated_atom_card_first_part->SetAtomCards(atoms_first_part);
            updated_atom_card_first_part->SetOrderedAtomCards(ordered_atoms_first_part);
            updated_atom_card_second_part->SetAtomCards(atoms_second_part);
            updated_atom_card_second_part->SetOrderedAtomCards(ordered_atoms_second_part);
            if(ordered_atoms_first_part.size() > 0)
            {
                updated_atom_cards.push_back(updated_atom_card_first_part);
                serial_number++;
            }
            if(ordered_atoms_second_part.size() > 0)
            {
                updated_atom_cards.push_back(updated_atom_card_second_part);
                serial_number++;
            }
        }
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();

        updated_residue_set->SetAtomCards(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;
    }
    models_->SetModels(updated_models);
}

void PdbFile::SplitAtomCardOfModelCardWithTheGivenModelNumber(char split_point_chain_id, int split_point_sequence_number, int model_number)
{
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbFileSpace::PdbModelSection::PdbModelCardMap updated_models;
        PdbModelCard* model = models[model_number];
        PdbModelCard* updated_model = new PdbModelCard();
        updated_model->SetModelSerialNumber(model->GetModelSerialNumber());
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet* updated_residue_set = new PdbModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        PdbModelResidueSet::AtomCardVector updated_atom_cards;
        int serial_number = 1;
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            PdbAtomSection* updated_atom_card_first_part = new PdbAtomSection();
            PdbAtomSection* updated_atom_card_second_part = new PdbAtomSection();
            updated_atom_card_first_part->SetRecordName(atom_card->GetRecordName());
            updated_atom_card_second_part->SetRecordName(atom_card->GetRecordName());
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
            PdbAtomSection::PdbAtomMap atoms_first_part;
            PdbAtomSection::PdbAtomMap atoms_second_part;
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_first_part = PdbAtomSection::PdbAtomCardOrderVector();
            PdbAtomSection::PdbAtomCardOrderVector ordered_atoms_second_part = PdbAtomSection::PdbAtomCardOrderVector();

            bool located = false;
            for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2);
                char residue_chain_id = atom->GetAtomChainId();
                int residue_sequence_number = atom->GetAtomResidueSequenceNumber();

                if(residue_chain_id == split_point_chain_id && residue_sequence_number == split_point_sequence_number)
                {
                    located = true;
                }
                if(located)
                {
                    serial_number++;
                    PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), atom->GetAtomResidueSequenceNumber(), atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(), atom->GetAlternateAtomCards());
                    atoms_second_part[serial_number] = updated_atom;
                    ordered_atoms_second_part.push_back(updated_atom);
                }
                if(!located)
                {
                    PdbFileSpace::PdbAtomCard* updated_atom = new PdbFileSpace::PdbAtomCard(serial_number, atom->GetAtomName(),atom->GetAtomAlternateLocation(), atom->GetAtomResidueName(),
                                                        atom->GetAtomChainId(), atom->GetAtomResidueSequenceNumber(), atom->GetAtomInsertionCode(), atom->GetAtomOrthogonalCoordinate(),
                                                        atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(), atom->GetAtomElementSymbol(), atom->GetAtomCharge(), atom->GetAlternateAtomCards());
                    atoms_first_part[serial_number] = updated_atom;
                    ordered_atoms_first_part.push_back(updated_atom);
                    serial_number++;
                }
            }
            updated_atom_card_first_part->SetAtomCards(atoms_first_part);
            updated_atom_card_first_part->SetOrderedAtomCards(ordered_atoms_first_part);
            updated_atom_card_second_part->SetAtomCards(atoms_second_part);
            updated_atom_card_second_part->SetOrderedAtomCards(ordered_atoms_second_part);
            if(ordered_atoms_first_part.size() > 0)
            {
                updated_atom_cards.push_back(updated_atom_card_first_part);
                serial_number++;
            }
            if(ordered_atoms_second_part.size() > 0)
            {
                updated_atom_cards.push_back(updated_atom_card_second_part);
                serial_number++;
            }
        }
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();

        updated_residue_set->SetAtomCards(updated_atom_cards);
        updated_residue_set->SetHeterogenAtoms(heterogen_atom_cards);
        updated_model->SetModelResidueSet(updated_residue_set);
        updated_models[updated_model->GetModelSerialNumber()] = updated_model;
        models_->SetModels(updated_models);
    }
}

void PdbFile::UpdateConnectCard()
{
    if(connectivities_ != NULL)
    {
        PdbFileSpace::PdbConnectSection::BondedAtomsSerialNumbersMap bondedAtomsSerialNumbersMap = connectivities_->GetBondedAtomsSerialNumbers();
        PdbFileSpace::PdbConnectSection::BondedAtomsSerialNumbersMap new_bonded_atoms_serial_numbers_map = PdbFileSpace::PdbConnectSection::BondedAtomsSerialNumbersMap();
        for(PdbFileSpace::PdbConnectSection::BondedAtomsSerialNumbersMap::iterator it = bondedAtomsSerialNumbersMap.begin(); it != bondedAtomsSerialNumbersMap.end(); it++)
        {
            int source_serial_number = (*it).first;
            std::vector<int> bonded_serial_numbers = (*it).second;
            if(serial_number_mapping_.find(source_serial_number) != serial_number_mapping_.end())
            {
                int new_source_serial_number = serial_number_mapping_[source_serial_number];
                new_bonded_atoms_serial_numbers_map[new_source_serial_number] = std::vector<int>();
                for(std::vector<int>::iterator it1 = bonded_serial_numbers.begin(); it1 != bonded_serial_numbers.end(); it1++)
                {
                    int bonded_serial_number = (*it1);
                    if(serial_number_mapping_.find(bonded_serial_number) != serial_number_mapping_.end())
                    {
                        int new_bonded_serial_number = serial_number_mapping_[bonded_serial_number];
                        new_bonded_atoms_serial_numbers_map[new_source_serial_number].push_back(new_bonded_serial_number);
                    }
                }
            }
        }
        connectivities_->SetBondedAtomsSerialNumbers(new_bonded_atoms_serial_numbers_map);
    }
}
void PdbFile::SetMasterCard(PdbFileSpace::PdbMasterCard *master)
{
    master_ = master;
}
//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////
bool PdbFile::Read(std::ifstream &in_file)
{
    if(!this->ParseCards(in_file))
        return false;
	return true;
}

bool PdbFile::Read(std::stringstream &atomstream)
{
  if(!this->ParseAtomStream(atomstream))
    return false;
  return true;
}

bool PdbFile::ParseAtomStream(std::stringstream &atomstream)
{
  //This function is to take an input of just atom cards in a stringstream and create a PdbFile object

  // std::string line;
  // getline(atomstream, line);
  // line = gmml::ExpandLine(line, gmml::iPdbLineLength);
  // std::string record_name = line.substr(0,6);
  // record_name = gmml::Trim(record_name);
  // atomstream.seekg(0, atomstream.beg);
  models_ = new PdbFileSpace::PdbModelSection(atomstream);

  if(models_ == NULL)
  {
    return false;
  }
  return true;
}

bool PdbFile::ParseCards(std::ifstream &in_stream)
{
    std::string line;

    /// Unable to read file
    if (!getline(in_stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format");
//        std::cout << "Wrong input file format" << std::endl;
        throw PdbFileProcessingException("Error reading file");
    }

    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("HEADER") == 0)
    {
        if(!ParseHeaderCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("OBSLTE") == 0)
    {
        if(!ParseObsoleteSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("TITLE") == 0)
    {
        if(!ParseTitleSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("SPLIT") == 0)
    {
        if(!ParseSplitSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("CAVEAT") == 0)
    {
        if(!ParseCaveatSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("COMPND") == 0)
    {
        if(!ParseCompoundSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("SOURCE") == 0)
    {
        if(!ParseSourceSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("KEYWDS") == 0)
    {
        if(!ParseKeywordsSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("EXPDTA") == 0)
    {
        if(!ParseExperimentalDataSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("NUMMDL") == 0)
    {
        if(!ParseNumModelCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("MDLTYP") == 0)
    {
        if(!ParseModelTypeSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("AUTHOR") == 0)
    {
        if(!ParseAuthorSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("REVDAT") == 0)
    {
        if(!ParseRevisionDataSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("SPRSDE") == 0)
    {
        if(!ParseSupersededEntriesSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("JRNL") == 0)
    {
        if(!ParseJournalSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("REMARK") == 0)
    {
        if(!ParseRemarkSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.find("DBREF") != std::string::npos)
    {
        if(!ParseDatabaseReferenceSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("SEQADV") == 0)
    {
        if(!ParseSequenceAdvancedSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("SEQRES") == 0)
    {
        if(!ParseResidueSequenceSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("MODRES") == 0)
    {
        if(!ParseResidueModificationSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("HET") == 0)
    {
        if(!ParseHeterogenSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("HETNAM") == 0)
    {
        if(!ParseHeterogenNameSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("HETSYN") == 0)
    {
        if(!ParseHeterogenSynonymSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("FORMUL") == 0)
    {
        if(!ParseFormulaSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("HELIX") == 0)
    {
        if(!ParseHelixSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("SHEET") == 0)
    {
        if(!ParseSheetSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("SSBOND") == 0)
    {
        if(!ParseDisulfideBondSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("LINK") == 0)
    {
        if(!ParseLinkSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("CISPEP") == 0)
    {
        if(!ParseCISPeptideSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("SITE") == 0)
    {
        if(!ParseSiteSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("CRYST1") == 0)
    {
        if(!ParseCrystallographyCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.find("ORIGX") != std::string::npos)
    {
        if(!ParseOriginCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.find("SCALE") != std::string::npos)
    {
        if(!ParseScaleCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("TVECT") == 0)
    {
      while(record_name.compare("TVECT") == 0)
      {
        getline(in_stream, line);//skip for now. TODO figure out if we care about having this data
        record_name = line.substr(0,6);
        record_name = gmml::Trim(record_name);
      }
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.find("MTRIX") != std::string::npos)
    {
        if(!ParseMatrixSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("MODEL") == 0 || record_name.compare("ATOM") == 0 || record_name.compare("HETATM") == 0)
    {
        if(!ParseModelSection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.compare("CONECT") == 0)
    {
        if(!ParseConnectivitySection(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    if(record_name.compare("MODEL") == 0)
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Multiple connect card between model cards --> Unexpected entry");
//        std::cout << "Multiple connect card between model cards --> Unexpected entry" << std::endl;
        return false;
    }
    if(record_name.compare("MASTER") == 0)
    {
        if(!ParseMasterCard(in_stream, line))
            return false;
    }
    record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    if(record_name.find("END") != std::string::npos || record_name.compare("END") == 0)
    {
        if(!ParseEndCard(in_stream, line))
            return false;
        return true;
    }
    else
    {
        // gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format");
        // std::cout << "Wrong input file format" << std::endl;
        std::stringstream ss;
        ss << record_name << " is an Unknown record name.";
        gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
        std::cerr << ss.str() << std::endl;
        // return false;
    }
    return true;
}

bool PdbFile::ParseHeaderCard(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Header card corupption");
//        std::cout << "Header card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format");
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("HEADER") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Header card corruption");
//            std::cout << "Header card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }

    header_ = new PdbFileSpace::PdbHeaderCard(stream_block);
    return true;
}

bool PdbFile::ParseObsoleteSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Obsolete card corruption");
//        std::cout << "Obsolete card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format");
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("OBSLTE") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Obsolete card corruption");
//            std::cout << "Obsolete card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format");
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }

    obsolete_= new PdbFileSpace::PdbObsoleteSection(stream_block);
    // obsolete_->Print();
    return true;
}

bool PdbFile::ParseTitleSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Title card corruption");
//        std::cout << "Title card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format");
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("TITLE") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Title card corruption");
//            std::cout << "Title card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format");
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }

    title_ = new PdbFileSpace::PdbTitleSection(stream_block);
    return true;
}

bool PdbFile::ParseSplitSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Split card corruption" );
//        std::cout << "Split card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("SPLIT") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Split card corruption" );
//            std::cout << "Split card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    split_ = new PdbFileSpace::PdbSplitSection(stream_block);
    // split_->Print();
    return true;
}

bool PdbFile::ParseCaveatSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Caveat card corruption" );
//        std::cout << "Caveat card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("CAVEAT") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Caveat card corruption" );
//            std::cout << "Caveat card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    caveat_ = new PdbFileSpace::PdbCaveatSection(stream_block);
    // caveat_->Print();
    return true;
}

bool PdbFile::ParseCompoundSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Compound card corruption" );
//        std::cout << "Compound card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("COMPND") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Compound card corruption" );
//            std::cout << "Compound card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }

    compound_ = new PdbFileSpace::PdbCompoundSection(stream_block);
    return true;
}

bool PdbFile::ParseSourceSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Source card corruption" );
//        std::cout << "Source card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("SOURCE") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Source card corruption" );
//            std::cout << "Source card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    source_ = new PdbFileSpace::PdbSourceSection(stream_block);
    // source_->Print();
    return true;
}

bool PdbFile::ParseKeywordsSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Keyword card corruption" );
//        std::cout << "Keyword card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("KEYWDS") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Keyword card corruption" );
//            std::cout << "Keyword card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    keywords_ =  new PdbFileSpace::PdbKeywordsSection(stream_block);
    // keywords_->Print();
    return true;
}

bool PdbFile::ParseExperimentalDataSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Experimental data card corruption" );
//        std::cout << "Experimental data card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("EXPDTA") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Experimental data card corruption" );
//            std::cout << "Experimental data card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    experimental_data_ = new PdbFileSpace::PdbExperimentalDataSection(stream_block);
    // experimental_data_->Print();
    return true;
}

bool PdbFile::ParseNumModelCard(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Number of model card corruption" );
//        std::cout << "Number of model card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("NUMMDL") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Number of model card corruption" );
//            std::cout << "Number of model card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    number_of_models_ = new PdbFileSpace::PdbNumModelCard(stream_block);
    return true;
}

bool PdbFile::ParseModelTypeSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Model type card corruption" );
//        std::cout << "Model type card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("MDLTYP") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Model type card corruption" );
//            std::cout << "Model type card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    model_type_ = new PdbFileSpace::PdbModelTypeSection(stream_block);
    return true;
}

bool PdbFile::ParseAuthorSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Author card corruption" );
//        std::cout << "Author card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("AUTHOR") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Author card corruption" );
//            std::cout << "Author card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    author_ = new PdbFileSpace::PdbAuthorSection(stream_block);
    // author_->Print();
    return true;
}

bool PdbFile::ParseRevisionDataSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Revision data card corruption" );
//        std::cout << "Revision data card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("REVDAT") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Revision data card corruption" );
//            std::cout << "Revision data card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    revision_data_ = new PdbFileSpace::PdbRevisionDataSection(stream_block);
    // revision_data_->Print();
    return true;
}

bool PdbFile::ParseSupersededEntriesSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Superseded entries card corruption" );
//        std::cout << "Superseded entries card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("SPRSDE") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Superseded entries card corruption" );
//            std::cout << "Superseded entries card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    superseded_entries_ = new PdbFileSpace::PdbSupersededEntriesSection(stream_block);
    // superseded_entries_->Print();
    return true;
}

bool PdbFile::ParseJournalSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Journal card corruption" );
//        std::cout << "Journal card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("JRNL") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Journal card corruption" );
//            std::cout << "Journal card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    journal_ = new PdbFileSpace::PdbJournalSection(stream_block);
    // journal_->Print();
    return true;
}

bool PdbFile::ParseRemarkSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Remark card corruption" );
//        std::cout << "Remark card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("REMARK") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Remark card corruption" );
//            std::cout << "Remark card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    remark_cards_ = new PdbFileSpace::PdbRemarkSection(stream_block);
    // std::cout << remark_cards_->GetResolution();
    return true;
}

bool PdbFile::ParseDatabaseReferenceSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "database reference card corruption" );
//        std::cout << "database reference card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("DBREF") == 0 || record_name.compare("DBREF1") == 0 ||record_name.compare("DBREF2") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "database reference card corruption" );
//            std::cout << "database reference card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    database_reference_ = new PdbFileSpace::PdbDatabaseReferenceSection(stream_block);
    // database_reference_->Print();
    return true;
}

bool PdbFile::ParseSequenceAdvancedSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Sequence advanced card corruption" );
//        std::cout << "Sequence advanced card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("SEQADV") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Sequence advanced card corruption" );
//            std::cout << "Sequence advanced card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    sequence_advanced_ = new PdbFileSpace::PdbSequenceAdvancedSection(stream_block);
    // sequence_advanced_->Print();
    return true;
}

bool PdbFile::ParseResidueSequenceSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Sequence residue card corruption" );
//        std::cout << "Sequence residue card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("SEQRES") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Sequence residue card corruption" );
//            std::cout << "Sequence residue card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    residues_sequence_ = new PdbFileSpace::PdbResidueSequenceSection(stream_block);
    return true;
}

bool PdbFile::ParseResidueModificationSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Modification residue card corruption" );
//        std::cout << "Modification residue card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("MODRES") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Modification residue card corruption" );
//            std::cout << "Modification residue card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    residue_modification_cards_ = new PdbFileSpace::PdbResidueModificationSection(stream_block);
    return true;
}

bool PdbFile::ParseHeterogenSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Heterogen card corruption" );
//        std::cout << "Heterogen card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("HET") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Heterogen card corruption" );
//            std::cout << "Heterogen card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    heterogen_cards_ = new PdbFileSpace::PdbHeterogenSection(stream_block);
    return true;
}

bool PdbFile::ParseHeterogenNameSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Heterogen name card corruption" );
//        std::cout << "Heterogen name card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("HETNAM") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Heterogen name card corruption" );
//            std::cout << "Heterogen name card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    heterogen_name_cards_ = new PdbFileSpace::PdbHeterogenNameSection(stream_block);
    return true;
}

bool PdbFile::ParseHeterogenSynonymSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Heterogen synonym card corruption" );
//        std::cout << "Heterogen synonym card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("HETSYN") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Heterogen synonym card corruption" );
//            std::cout << "Heterogen synonym card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    heterogen_synonym_cards_ = new PdbFileSpace::PdbHeterogenSynonymSection(stream_block);
    return true;
}

bool PdbFile::ParseFormulaSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Formula card corruption" );
//        std::cout << "Formula card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("FORMUL") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Formula card corruption" );
//            std::cout << "Formula card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    formulas_ = new PdbFileSpace::PdbFormulaSection(stream_block);
    return true;
}

bool PdbFile::ParseHelixSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Helix card corruption" );
//        std::cout << "Helix card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("HELIX") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Helix card corruption" );
//            std::cout << "Helix card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    helix_cards_ = new PdbFileSpace::PdbHelixSection(stream_block);
    return true;
}

bool PdbFile::ParseSheetSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Sheet card corruption" );
//        std::cout << "Sheet card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("SHEET") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Sheet card corruption" );
//            std::cout << "Sheet card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    sheet_cards_ = new PdbFileSpace::PdbSheetSection(stream_block);
    // sheet_cards_->Print();
    return true;
}

bool PdbFile::ParseDisulfideBondSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Disulfide bond card corruption" );
//        std::cout << "Disulfide bond card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("SSBOND") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Disulfide bond card corruption" );
//            std::cout << "Disulfide bond card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    disulfide_bonds_ = new PdbFileSpace::PdbDisulfideBondSection(stream_block);
    return true;
}

bool PdbFile::ParseLinkSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Link card corruption" );
//        std::cout << "Link card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("LINK") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Link card corruption" );
//            std::cout << "Link card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    link_cards_ = new PdbFileSpace::PdbLinkSection(stream_block);
    return true;
}

bool PdbFile::ParseCISPeptideSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "CIS peptide card corruption" );
//        std::cout << "CIS peptide card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("CISPEP") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "CIS peptide card corruption" );
//            std::cout << "CIS peptide card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    cis_peptide_ = new PdbFileSpace::PdbCISPeptideSection(stream_block);
    // cis_peptide_->Print();
    return true;
}

bool PdbFile::ParseSiteSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Site card corruption" );
//        std::cout << "Site card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("SITE") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Site card corruption" );
//            std::cout << "Site card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    site_cards_ = new PdbFileSpace::PdbSiteSection(stream_block);
    return true;
}

bool PdbFile::ParseCrystallographyCard(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Crystallography card corruption" );
//        std::cout << "Crystallography card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("CRYST1") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Crystallography card corruption" );
//            std::cout << "Crystallography card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    crystallography_ = new PdbFileSpace::PdbCrystallographicCard(stream_block);
    return true;
}

bool PdbFile::ParseOriginCard(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Origin card corruption" );
//        std::cout << "Origin card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,5);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("ORIGX") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,5);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Origin card corruption" );
//            std::cout << "Origin card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    origins_ = new PdbFileSpace::PdbOriginXnSection(stream_block);
    return true;
}

bool PdbFile::ParseScaleCard(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Scale card corruption" );
//        std::cout << "Scale card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,5);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("SCALE") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,5);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Scale card corruption" );
//            std::cout << "Scale card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    scales_ = new PdbFileSpace::PdbScaleNSection(stream_block);
    return true;
}

bool PdbFile::ParseMatrixSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Matrix card corruption" );
//        std::cout << "Matrix card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,5);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("MTRIX") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,5);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Matrix card corruption" );
//            std::cout << "Matrix card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    matrices_ = new PdbFileSpace::PdbMatrixNSection(stream_block);
    return true;
}

bool PdbFile::ParseModelSection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Model card corruption" );
//        std::cout << "Model card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("MODEL") == 0 || record_name.compare("ATOM") == 0 || record_name.compare("ANISOU") == 0
          || record_name.compare("TER") == 0 || record_name.compare("HETATM") == 0 || record_name.compare("ENDMDL") == 0)
        //          || record_name.find("TER") != std::string::npos || record_name.find("ENDMDL") != std::string::npos)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Model card corruption" );
//            std::cout << "Model card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    // Model card
       // gmml::log(__LINE__, __FILE__,  gmml::ERR, stream_block.str());
    models_ = new PdbFileSpace::PdbModelSection(stream_block);
    return true;
}

bool PdbFile::ParseConnectivitySection(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Connectivity card corruption" );
//        std::cout << "Connectivity card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("CONECT") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Connectivity card corruption" );
//            std::cout << "Connectivity card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    connectivities_ = new PdbFileSpace::PdbConnectSection(stream_block);
    return true;
}

bool PdbFile::ParseMasterCard(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Master card corruption" );
//        std::cout << "Master card corruption" << std::endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//        std::cout << "Wrong input file format" << std::endl;
        return false;
    }
    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);

    while(record_name.compare("MASTER") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Master card corruption" );
//            std::cout << "Master card corruption" << std::endl;
            gmml::log(__LINE__, __FILE__,  gmml::ERR, "Wrong input file format" );
//            std::cout << "Wrong input file format" << std::endl;
            return false;
        }
    }
    master_ = new PdbFileSpace::PdbMasterCard(stream_block);
    return true;
}

bool PdbFile::ParseEndCard(std::ifstream& stream, std::string& line)
{
    std::stringstream stream_block;
    stream_block << line << std::endl;
    if(!getline(stream, line))
    {
        gmml::log(__LINE__, __FILE__,  gmml::INF, "End of file" );
        // std::cout << "End of file" << std::endl;
        return true;
    }

    line = gmml::ExpandLine(line, gmml::iPdbLineLength);
    std::string record_name = line.substr(0,6);
    record_name = gmml::Trim(record_name);
    while(record_name.find("END") != std::string::npos || record_name.compare("END") == 0)
    {
        stream_block << line << std::endl;
        if(getline(stream, line))
        {
            line = gmml::ExpandLine(line, gmml::iPdbLineLength);
            record_name = line.substr(0,6);
            record_name = gmml::Trim(record_name);
        }
        else
        {
            gmml::log(__LINE__, __FILE__,  gmml::INF, "End of file" );
            return true;
        }
    }
    // end_ = new PdbEndCard(stream_block);
    return true;
}
void PdbFile::WriteToStringstream(std::ostringstream& pdbstream)
{
  this->ResolveCards(pdbstream);
}
void PdbFile::Write(const std::string& pdb_file)
{
    std::ofstream out_file;
    try
    {
        out_file.open(pdb_file.c_str());
    }
    catch(...)
    {
        throw PdbFileProcessingException(__LINE__,"File could not be created");
    }
    try
    {
        this->ResolveCards(out_file);
    }
    catch(...)
    {
        out_file.close();            /// Close the parameter files
    }
}

void PdbFile::WriteWithTheGivenModelNumber(const std::string& pdb_file, int model_number)
{
    std::ofstream out_file;
    try
    {
        out_file.open(pdb_file.c_str());
    }
    catch(...)
    {
        throw PdbFileProcessingException(__LINE__,"File could not be created");
    }
    try
    {
        this->ResolveCardsWithTheGivenModelNumber(out_file, model_number);
    }
    catch(...)
    {
        out_file.close();            /// Close the parameter files
    }
}

void PdbFile::ResolveCards(std::ostream& out_stream)
{
    if(this->header_ != NULL)
    {
        this->ResolveHeaderCard(out_stream);
    }
    if(this->obsolete_ != NULL)
    {
        this->ResolveObsoleteCards(out_stream);
    }
    if(this->title_ != NULL)
    {
        this->ResolveTitleCards(out_stream);
    }
    if(this->split_ != NULL)
    {
        this->ResolveSplitCards(out_stream);
    }if(this->caveat_ != NULL)
    {
        this->ResolveCaveatCards(out_stream);
    }
    if(this->compound_ != NULL)
    {
        this->ResolveCompoundCards(out_stream);
    }
    if(this->source_ != NULL)
    {
        this->ResolveSourceCards(out_stream);
    }
    if(this->keywords_ != NULL)
    {
        this->ResolveKeywordCards(out_stream);
    }
    if(this->experimental_data_ != NULL)
    {
        this->ResolveExperimentalDataCards(out_stream);
    }
    if(this->number_of_models_ != NULL)
    {
        this->ResolveCompoundCards(out_stream);
    }
    if(this->model_type_ != NULL)
    {
        this->ResolveNumModelCard(out_stream);
    }
    if(this->author_ != NULL)
    {
        this->ResolveAuthorCards(out_stream);
    }
    if(this->revision_data_ != NULL)
    {
        this->ResolveRevisionDataCards(out_stream);
    }
    if(this->superseded_entries_ != NULL)
    {
        this->ResolveSupersededEntriesCards(out_stream);
    }
    if(this->journal_ != NULL)
    {
        this->ResolveJournalCards(out_stream);
    }
    if(this->remark_cards_ != NULL)
    {
        this->ResolveRemarkCards(out_stream);
    }
    if(this->database_reference_ != NULL)
    {
        this->ResolveDatabaseReferenceCards(out_stream);
    }
    if(this->sequence_advanced_ != NULL)
    {
        this->ResolveSequenceAdvancedCards(out_stream);
    }
    if(this->residues_sequence_ != NULL)
    {
        this->ResolveSequenceResidueCards(out_stream);
    }
    if(this->residue_modification_cards_ != NULL)
    {
        this->ResolveModificationResidueCards(out_stream);
    }
    if(this->heterogen_cards_ != NULL)
    {
        this->ResolveHeterogenCards(out_stream);
    }
    if(this->heterogen_name_cards_ != NULL)
    {
        this->ResolveHeterogenNameCards(out_stream);
    }
    if(this->heterogen_synonym_cards_ != NULL)
    {
        this->ResolveHeterogenSynonymCards(out_stream);
    }
    if(this->formulas_ != NULL)
    {
        this->ResolveFormulaCards(out_stream);
    }
    if(this->helix_cards_ != NULL)
    {
        this->ResolveHelixCards(out_stream);
    }
    if(this->sheet_cards_ != NULL)
    {
        this->ResolveSheetCards(out_stream);
    }
    if(this->disulfide_bonds_ != NULL)
    {
        this->ResolveDisulfideBondCards(out_stream);
    }
    if(this->link_cards_ != NULL)
    {
        this->ResolveLinkCards(out_stream);
    }
    if(this->cis_peptide_ != NULL)
    {
        this->ResolveCISPeptideCards(out_stream);
    }
    if(this->site_cards_ != NULL)
    {
        this->ResolveSiteCards(out_stream);
    }
    if(this->crystallography_ != NULL)
    {
        this->ResolveCrystallographyCard(out_stream);
    }
    if(this->origins_ != NULL)
    {
        this->ResolveOriginCard(out_stream);
    }
    if(this->scales_ != NULL)
    {
        this->ResolveScaleCard(out_stream);
    }
    if(this->matrices_ != NULL)
    {
        this->ResolveMatrixCards(out_stream);
    }
    if(this->models_ != NULL)
    {
        this->ResolveModelCards(out_stream);
    }
    if(this->connectivities_ != NULL)
    {
        this->ResolveConnectivityCards(out_stream);
    }
    if(this->master_ != NULL)
    {
        this->ResolveMasterCards(out_stream);
    }
    this->ResolveEndCard(out_stream);
}

void PdbFile::ResolveCardsWithTheGivenModelNumber(std::ofstream& out_stream, int model_number)
{
    if(this->header_ != NULL)
    {
        this->ResolveHeaderCard(out_stream);
    }
    if(this->obsolete_ != NULL)
    {
        this->ResolveObsoleteCards(out_stream);
    }
    if(this->title_ != NULL)
    {
        this->ResolveTitleCards(out_stream);
    }
    if(this->split_ != NULL)
    {
        this->ResolveSplitCards(out_stream);
    }if(this->caveat_ != NULL)
    {
        this->ResolveCaveatCards(out_stream);
    }
    if(this->compound_ != NULL)
    {
        this->ResolveCompoundCards(out_stream);
    }
    if(this->source_ != NULL)
    {
        this->ResolveSourceCards(out_stream);
    }
    if(this->keywords_ != NULL)
    {
        this->ResolveKeywordCards(out_stream);
    }
    if(this->experimental_data_ != NULL)
    {
        this->ResolveExperimentalDataCards(out_stream);
    }
    if(this->number_of_models_ != NULL)
    {
        this->ResolveCompoundCards(out_stream);
    }
    if(this->model_type_ != NULL)
    {
        this->ResolveNumModelCard(out_stream);
    }
    if(this->author_ != NULL)
    {
        this->ResolveAuthorCards(out_stream);
    }
    if(this->revision_data_ != NULL)
    {
        this->ResolveRevisionDataCards(out_stream);
    }
    if(this->superseded_entries_ != NULL)
    {
        this->ResolveSupersededEntriesCards(out_stream);
    }
    if(this->journal_ != NULL)
    {
        this->ResolveJournalCards(out_stream);
    }
    if(this->remark_cards_ != NULL)
    {
        this->ResolveRemarkCards(out_stream);
    }
    if(this->database_reference_ != NULL)
    {
        this->ResolveDatabaseReferenceCards(out_stream);
    }
    if(this->sequence_advanced_ != NULL)
    {
        this->ResolveSequenceAdvancedCards(out_stream);
    }
    if(this->residues_sequence_ != NULL)
    {
        this->ResolveSequenceResidueCards(out_stream);
    }
    if(this->residue_modification_cards_ != NULL)
    {
        this->ResolveModificationResidueCards(out_stream);
    }
    if(this->heterogen_cards_ != NULL)
    {
        this->ResolveHeterogenCards(out_stream);
    }
    if(this->heterogen_name_cards_ != NULL)
    {
        this->ResolveHeterogenNameCards(out_stream);
    }
    if(this->heterogen_synonym_cards_ != NULL)
    {
        this->ResolveHeterogenSynonymCards(out_stream);
    }
    if(this->formulas_ != NULL)
    {
        this->ResolveFormulaCards(out_stream);
    }
    if(this->helix_cards_ != NULL)
    {
        this->ResolveHelixCards(out_stream);
    }
    if(this->sheet_cards_ != NULL)
    {
        this->ResolveSheetCards(out_stream);
    }
    if(this->disulfide_bonds_ != NULL)
    {
        this->ResolveDisulfideBondCards(out_stream);
    }
    if(this->link_cards_ != NULL)
    {
        this->ResolveLinkCards(out_stream);
    }
    if(this->cis_peptide_ != NULL)
    {
        this->ResolveCISPeptideCards(out_stream);
    }
    if(this->site_cards_ != NULL)
    {
        this->ResolveSiteCards(out_stream);
    }
    if(this->crystallography_ != NULL)
    {
        this->ResolveCrystallographyCard(out_stream);
    }
    if(this->origins_ != NULL)
    {
        this->ResolveOriginCard(out_stream);
    }
    if(this->scales_ != NULL)
    {
        this->ResolveScaleCard(out_stream);
    }
    if(this->matrices_ != NULL)
    {
        this->ResolveMatrixCards(out_stream);
    }
    if(this->models_ != NULL)
    {
        this->ResolveModelCardWithTheGivenModelNumber(out_stream, model_number);
    }
    if(this->connectivities_ != NULL)
    {
        this->ResolveConnectivityCards(out_stream);
    }
    if(this->master_ != NULL)
    {
        this->ResolveMasterCards(out_stream);
    }
    this->ResolveEndCard(out_stream);
}

void PdbFile::ResolveHeaderCard(std::ostream& stream)
{
    stream << std::left << std::setw(6) << header_->GetRecordName()
           << std::left << std::setw(4) << " "
           << std::left << std::setw(40) << header_->GetClassification()
           << std::left << std::setw(9) << header_->GetDepositionDate()
           << std::left << std::setw(3) << " "
           << std::right << std::setw(4) << header_->GetIdentifierCode()
           << std::left << std::setw(14) << " "
           << std::endl;
}

void PdbFile::ResolveObsoleteCards(std::ostream& stream)
{
    stream << std::left << std::setw(6) << obsolete_->GetRecordName()
          << std::left << std::setw(2) << " "
          << std::left << std::setw(2) << obsolete_->GetContinuation()
          << std::left << std::setw(1) << " "
          << std::left << std::setw(9) << obsolete_->GetReplacementDate()
          << std::left << std::setw(1) << " ";
          std::vector<std::string> identifier_codes = obsolete_->GetIdentifierCodes();
          for(unsigned int i = 0; i < identifier_codes.size(); i++)
          {
            stream << std::left << std::setw(4) << identifier_codes[i]
                  << std::left << std::setw(5) << "      ";
          }
          stream << std::left << std::setw(14) << " "
          << std::endl;
}

void PdbFile::ResolveTitleCards(std::ostream& stream)
{
    const int MAX_TITLE_LENGTH_IN_LINE = 70;
    stream << std::left << std::setw(6) << title_->GetRecordName()
           << std::left << std::setw(2) << " ";
    if((int)title_->GetTitle().length() > MAX_TITLE_LENGTH_IN_LINE)
    {
        stream << std::right << std::setw(2) << " "
               << std::left << std::setw(70) << title_->GetTitle().substr(0,MAX_TITLE_LENGTH_IN_LINE)
               << std::endl;

        int counter = ceil((double)(title_->GetTitle().length()) / MAX_TITLE_LENGTH_IN_LINE);
        for(int i = 2; i <= counter; i++)
        {
            if(i != counter)
            {
                stream << std::left << std::setw(6) << title_->GetRecordName()
                       << std::left << std::setw(2) << " "
                       << std::right << std::setw(2) << i
                       << std::left << std::setw(70) << title_->GetTitle().substr(MAX_TITLE_LENGTH_IN_LINE*(i-1), MAX_TITLE_LENGTH_IN_LINE)
                       << std::endl;
            }
            else
            {
                stream << std::left << std::setw(6) << title_->GetRecordName()
                       << std::left << std::setw(2) << " "
                       << std::right << std::setw(2) << i
                       << std::left << std::setw(70) << title_->GetTitle().substr(MAX_TITLE_LENGTH_IN_LINE*(i-1), title_->GetTitle().length()-MAX_TITLE_LENGTH_IN_LINE*(i-1))
                       << std::endl;
            }
        }
    }
    else
    {
        stream << std::right << std::setw(2) << " "
               << std::left << std::setw(70) << title_->GetTitle()
               << std::endl;
    }
}

void PdbFile::ResolveSplitCards(std::ostream& stream)
{
     const int MAX_SPLIT_ID_IN_LINE = 79;
     stream << std::left << std::setw(6) << split_->GetRecordName()
            << std::left << std::setw(2) << " ";
     if((int)split_->GetSplit().length() > MAX_SPLIT_ID_IN_LINE)
     {
         stream << std::left << std::setw(79) << split_->GetSplit().substr(0,MAX_SPLIT_ID_IN_LINE)
                << std::endl;

         int counter = ceil((double)(split_->GetSplit().length()) / MAX_SPLIT_ID_IN_LINE);
         for(int i = 2; i <= counter; i++)
         {
             if(i != counter)
             {
                 stream << std::left << std::setw(6) << split_->GetRecordName()
                        << std::left << std::setw(2) << " "
                        << std::right << std::setw(2) << i
                        << std::left << std::setw(79) << split_->GetSplit().substr(MAX_SPLIT_ID_IN_LINE*(i-1), MAX_SPLIT_ID_IN_LINE)
                        << std::endl;
             }
             else
             {
                 stream << std::left << std::setw(6) << split_->GetRecordName()
                        << std::left << std::setw(2) << " "
                        << std::right << std::setw(2) << i
                        << std::left << std::setw(79) << split_->GetSplit().substr(MAX_SPLIT_ID_IN_LINE*(i-1), split_->GetSplit().length()-MAX_SPLIT_ID_IN_LINE*(i-1))
                        << std::endl;
             }
         }
     }
     else
     {
         stream << std::right << std::setw(2) << " "
                << std::left << std::setw(79) << split_->GetSplit()
                << std::endl;
     }
}

void PdbFile::ResolveCaveatCards(std::ostream& stream)
{
    const int MAX_CAVEAT_LENGTH_IN_LINE = 70;
    stream << std::left << std::setw(6) << caveat_->GetRecordName()
           << std::left << std::setw(2) << " "
           << std::left << std::setw(4) << header_->GetIdentifierCode();
    if((int)caveat_->GetCaveat().length() > MAX_CAVEAT_LENGTH_IN_LINE)
    {
        stream << std::left << std::setw(70) << caveat_->GetCaveat().substr(0,MAX_CAVEAT_LENGTH_IN_LINE)
               << std::endl;

        int counter = ceil((double)(caveat_->GetCaveat().length()) / MAX_CAVEAT_LENGTH_IN_LINE);
        for(int i = 2; i <= counter; i++)
        {
            if(i != counter)
            {
                stream << std::left << std::setw(6) << caveat_->GetRecordName()
                       << std::left << std::setw(2) << " "
                       << std::right << std::setw(2) << i
                       << std::left << std::setw(70) << caveat_->GetCaveat().substr(MAX_CAVEAT_LENGTH_IN_LINE*(i-1), MAX_CAVEAT_LENGTH_IN_LINE)
                       << std::endl;
            }
            else
            {
                stream << std::left << std::setw(6) << caveat_->GetRecordName()
                       << std::left << std::setw(2) << " "
                       << std::right << std::setw(2) << i
                       << std::left << std::setw(70) << caveat_->GetCaveat().substr(MAX_CAVEAT_LENGTH_IN_LINE*(i-1), caveat_->GetCaveat().length()-MAX_CAVEAT_LENGTH_IN_LINE*(i-1))
                       << std::endl;
            }
        }
    }
    else
    {
        stream << std::right << std::setw(2) << " "
               << std::left << std::setw(70) << caveat_->GetCaveat()
               << std::endl;
    }
}

void PdbFile::ResolveCompoundCards(std::ostream& stream)
{
    const int MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE = 70;
    stream << std::left << std::setw(6) << compound_->GetRecordName()
           << std::left << std::setw(1) << " "
           << std::right << std::setw(3) << " ";

    PdbFileSpace::PdbCompoundSection::PdbCompoundSpecificationMap compound_specification_map = compound_->GetCompoundSpecifications();
    if((*(compound_specification_map.begin())).second->GetMoleculeId() != "")
    {
        std::stringstream ss;
        ss << "MOL_ID: " << (*(compound_specification_map.begin())).second->GetMoleculeId() << ";";
        stream << std::left << std::setw(70) << ss.str() << std::endl;
    }
    else
    {
        std::stringstream ss;
        ss << " UNKNOWN;";
        stream << std::left << std::setw(70) << ss.str() << std::endl;
    }
    bool first = true;
    int counter = 2;
    for(PdbFileSpace::PdbCompoundSection::PdbCompoundSpecificationMap::iterator it = compound_specification_map.begin(); it != compound_specification_map.end(); it++)
    {
        PdbCompoundSpecification* compound_specification = (*it).second;
        if(!first)
        {
            if(compound_specification->GetMoleculeId() != "")
            {
                std::stringstream ss;
                ss << " MOL_ID: " << compound_specification->GetMoleculeId() << ";";
                stream << std::left << std::setw(6) << compound_->GetRecordName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << counter
                       << std::left << std::setw(70) << ss.str() << std::endl;
                counter++;
            }
            else
            {
                std::stringstream ss;
                ss << " UNKNOWN;";
                stream << std::left << std::setw(6) << compound_->GetRecordName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << counter
                       << std::left << std::setw(70) << ss.str() << std::endl;
                counter++;
            }
        }

        /// Molecule name specification
        if(compound_specification->GetMoleculeName() != "")
        {
            std::stringstream molecule_name;
            molecule_name << " MOLECULE: " << compound_specification->GetMoleculeName() << ";";
            int length = molecule_name.str().length();

            if(length <= MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE)
            {
                std::stringstream ss;
                ss << molecule_name.str();
                stream << std::left << std::setw(6) << compound_->GetRecordName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << counter
                       << std::left << std::setw(70) << ss.str() << std::endl;
                counter++;
            }
            else
            {
                int number_of_lines = ceil((double)(length) / MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);

                for(int i = 1; i <= number_of_lines; i++)
                {
                    if(i != number_of_lines)
                    {
                        std::stringstream ss;
                        ss << molecule_name.str().substr((i-1)*(MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE), MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << std::left << std::setw(6) << compound_->GetRecordName()
                               << std::left << std::setw(1) << " "
                               << std::right << std::setw(3) << counter
                               << std::left << std::setw(70) << ss.str()
                               << std::endl;
                        counter++;
                    }
                    else
                    {
                        std::stringstream ss;
                        ss << molecule_name.str().substr((i-1)*(MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE), length - (i-1)*(MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE));
                        stream << std::left << std::setw(6) << compound_->GetRecordName()
                               << std::left << std::setw(1) << " "
                               << std::right << std::setw(3) << counter
                               << std::left << std::setw(70) << ss.str()
                               << std::endl;
                        counter++;
                    }
                }

            }
        }

        /// Molecule chain ids specification
        if(compound_specification->GetChainIds().size() > 0)
        {
            std::vector<std::string> chain_ids = compound_specification->GetChainIds();
            std::stringstream chain_id;
            chain_id << " CHAIN: ";
            for(std::vector<std::string>::iterator it1 = chain_ids.begin(); it1 != chain_ids.end(); it1++)
            {
                if(it1 < chain_ids.end()-1)
                    chain_id << (*it1) << ",";
                else
                    chain_id << (*it1);
            }
            chain_id << ";";
            int length = chain_id.str().length();
            if(length <= MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE)
            {
                std::stringstream ss;
                ss << chain_id.str();
                stream << std::left << std::setw(6) << compound_->GetRecordName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << counter
                       << std::left << std::setw(70) << ss.str() << std::endl;
                counter++;
            }
            else
            {
                int number_of_lines = ceil((double)(length) / MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                for(int i = 1; i <= number_of_lines; i++)
                {
                    if(i != number_of_lines)
                    {
                        std::stringstream ss;
                        ss << chain_id.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << std::left << std::setw(6) << compound_->GetRecordName()
                               << std::left << std::setw(1) << " "
                               << std::right << std::setw(3) << counter
                               << std::left << std::setw(70) << ss.str() << std::endl;
                        counter++;
                    }
                    else
                    {
                        std::stringstream ss;
                        ss << chain_id.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, length - (i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << std::left << std::setw(6) << compound_->GetRecordName()
                               << std::left << std::setw(1) << " "
                               << std::right << std::setw(3) << counter
                               << std::left << std::setw(70) << ss.str() << std::endl;
                        counter++;
                    }
                }
            }
        }

        /// Fragment specification
        if(compound_specification->GetFragment() != "")
        {
            std::stringstream fragment;
            fragment << " FRAGMENT: " << compound_specification->GetFragment() << ";";
            int length = fragment.str().length();
            if(length <= MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE)
            {
                std::stringstream ss;
                ss << ss.str();
                stream << std::left << std::setw(6) << compound_->GetRecordName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << counter
                       << std::left << std::setw(70) << ss.str() << std::endl;
                counter++;
            }
            else
            {
                int number_of_lines = ceil((double)(length) / MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                for(int i = 1; i <= number_of_lines; i++)
                {
                    if(i != number_of_lines)
                    {
                        std::stringstream ss;
                        ss << ss.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << std::left << std::setw(6) << compound_->GetRecordName()
                               << std::left << std::setw(1) << " "
                               << std::right << std::setw(3) << counter
                               << std::left << std::setw(70) << ss.str() << std::endl;
                        counter++;
                    }
                    else
                    {
                        std::stringstream ss;
                        ss << ss.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, length - (i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << std::left << std::setw(6) << compound_->GetRecordName()
                               << std::left << std::setw(1) << " "
                               << std::right << std::setw(3) << counter
                               << std::left << std::setw(70) << ss.str() << std::endl;
                        counter++;
                    }
                }
            }
        }

        /// Molecule synonyms specification
        if(compound_specification->GetMoleculeSynonyms().size() > 0)
        {
            std::vector<std::string> molecule_synonyms = compound_specification->GetMoleculeSynonyms();
            std::stringstream synonyms;
            synonyms << " SYNONYM: ";
            for(std::vector<std::string>::iterator it1 = molecule_synonyms.begin(); it1 != molecule_synonyms.end(); it1++)
            {
                if(it1 < molecule_synonyms.end()-1)
                    synonyms << (*it1) << ",";
                else
                    synonyms << (*it1);
            }
            synonyms << ";";
            int length = synonyms.str().length();
            if(length <= MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE)
            {
                std::stringstream ss;
                ss << synonyms.str();
                stream << std::left << std::setw(6) << compound_->GetRecordName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << counter
                       << std::left << std::setw(70) << ss.str() << std::endl;
                counter++;
            }
            else
            {
                int number_of_lines = ceil((double)(length) / MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                for(int i = 1; i <= number_of_lines; i++)
                {
                    if(i != number_of_lines)
                    {
                        std::stringstream ss;
                        ss << synonyms.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << std::left << std::setw(6) << compound_->GetRecordName()
                               << std::left << std::setw(1) << " "
                               << std::right << std::setw(3) << counter
                               << std::left << std::setw(70) << ss.str() << std::endl;
                        counter++;
                    }
                    else
                    {
                        std::stringstream ss;
                        ss << synonyms.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, length - (i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << std::left << std::setw(6) << compound_->GetRecordName()
                               << std::left << std::setw(1) << " "
                               << std::right << std::setw(3) << counter
                               << std::left << std::setw(70) << ss.str() << std::endl;
                        counter++;
                    }
                }
            }
        }

        /// Enzyme commission numbers specification
        if(compound_specification->GetEnzymeCommissionNumbers().size() > 0)
        {
            std::vector<std::string> enzyme_commission_numbers = compound_specification->GetEnzymeCommissionNumbers();
            std::stringstream commission_numbers;
            commission_numbers << " EC: ";
            for(std::vector<std::string>::iterator it1 = enzyme_commission_numbers.begin(); it1 != enzyme_commission_numbers.end(); it1++)
            {
                if(it1 < enzyme_commission_numbers.end()-1)
                    commission_numbers << (*it1) << ",";
                else
                    commission_numbers << (*it1);
            }
            commission_numbers << ";";
            int length = commission_numbers.str().length();
            if(length <= MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE)
            {
                std::stringstream ss;
                ss << commission_numbers.str();
                stream << std::left << std::setw(6) << compound_->GetRecordName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << counter
                       << std::left << std::setw(70) << ss.str() << std::endl;
                counter++;
            }
            else
            {
                int number_of_lines = ceil((double)(length) / MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                for(int i = 1; i <= number_of_lines; i++)
                {
                    if(i != number_of_lines)
                    {
                        std::stringstream ss;
                        ss << commission_numbers.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << std::left << std::setw(6) << compound_->GetRecordName()
                               << std::left << std::setw(1) << " "
                               << std::right << std::setw(3) << counter
                               << std::left << std::setw(70) << ss.str() << std::endl;
                        counter++;
                    }
                    else
                    {
                        std::stringstream ss;
                        ss << commission_numbers.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, length - (i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << std::left << std::setw(6) << compound_->GetRecordName()
                               << std::left << std::setw(1) << " "
                               << std::right << std::setw(3) << counter
                               << std::left << std::setw(70) << ss.str() << std::endl;
                        counter++;
                    }
                }
            }
        }

        /// Engineered specification
        if(compound_specification->GetIsEngineered())
        {
            std::stringstream ss;
            ss << " ENGINEERED: YES";
            stream << std::left << std::setw(6) << compound_->GetRecordName()
                   << std::left << std::setw(1) << " "
                   << std::right << std::setw(3) << counter
                   << std::left << std::setw(70) << ss.str() << std::endl;
            counter++;
        }

        /// Mutation specification
        if(compound_specification->GetHasMutation())
        {
            std::stringstream ss;
            ss << " MUTATION: YES;";
            stream << std::left << std::setw(6) << compound_->GetRecordName()
                   << std::left << std::setw(1) << " "
                   << std::right << std::setw(3) << counter
                   << std::left << std::setw(70) << ss.str() << std::endl;
            counter++;
        }

        /// Other comments specification
        if(compound_specification->GetComments() != "")
        {
            std::stringstream comments;
            comments << " OTHER_DETAILS: " << compound_specification->GetComments() << ";";
            int length = comments.str().length();
            if(length <= MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE)
            {
                std::stringstream ss;
                ss << comments.str();
                stream << std::left << std::setw(6) << compound_->GetRecordName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << counter
                       << std::left << std::setw(70) << ss.str() << std::endl;
                counter++;
            }
            else
            {
                int number_of_lines = ceil((double)(length) / MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                for(int i = 1; i <= number_of_lines; i++)
                {
                    if(i != number_of_lines)
                    {
                        std::stringstream ss;
                        ss << comments.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << std::left << std::setw(6) << compound_->GetRecordName()
                               << std::left << std::setw(1) << " "
                               << std::right << std::setw(3) << counter
                               << std::left << std::setw(70) << ss.str() << std::endl;
                        counter++;
                    }
                    else
                    {
                        std::stringstream ss;
                        ss << comments.str().substr((i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE, length - (i-1)*MAX_LENGTH_OF_COMPOUND_SPEC_IN_LINE);
                        stream << std::left << std::setw(6) << compound_->GetRecordName()
                               << std::left << std::setw(1) << " "
                               << std::right << std::setw(3) << counter
                               << std::left << std::setw(70) << ss.str() << std::endl;
                        counter++;
                    }
                }
            }
        }

        first = false;

    }

}

void PdbFile::ResolveSourceCards(std::ostream& stream)
{
  int SOURCE_COUNT = 1;
  SourceCardVector source_cards = source_->GetSourceCards();
  for (SourceCardVector::iterator it = source_cards.begin(); it != source_cards.end(); it++)
  {
    if (SOURCE_COUNT ==1)
    {
      stream << std::left << std::setw(6) << (*it)->GetRecordName()
             << std::left << std::setw(4) << " "
             << std::left << std::setw(69) << (*it)->GetToken() +":"+ (*it)->GetValue()
             << std::endl;
      SOURCE_COUNT ++;
    }
    else
    {
      stream << std::left << std::setw(6) << (*it)->GetRecordName()
             << std::right << std::setw(4) << SOURCE_COUNT
             << std::left << std::setw(68) << (*it)->GetToken() +":"+ (*it)->GetValue()
             << std::endl;
      SOURCE_COUNT ++;
    }
  }
}

void PdbFile::ResolveKeywordCards(std::ostream& stream)
{
  const int MAX_KEYWORDS_LENGTH_IN_LINE = 70;
  stream << std::left << std::setw(6) << keywords_->GetRecordName()
         << std::left << std::setw(2) << " ";
  if((int)keywords_->GetKeywords().length() > MAX_KEYWORDS_LENGTH_IN_LINE)
  {
      stream << std::left << std::setw(2) << " "
             << std::left << std::setw(70) << keywords_->GetKeywords().substr(0,MAX_KEYWORDS_LENGTH_IN_LINE)
             << std::endl;

      int counter = ceil((double)(keywords_->GetKeywords().length()) / MAX_KEYWORDS_LENGTH_IN_LINE);
      for(int i = 2; i <= counter; i++)
      {
          if(i != counter)
          {
              stream << std::left << std::setw(6) << keywords_->GetRecordName()
                     << std::left << std::setw(2) << " "
                     << std::right << std::setw(2) << i
                     << std::left << std::setw(70) << keywords_->GetKeywords().substr(MAX_KEYWORDS_LENGTH_IN_LINE*(i-1), MAX_KEYWORDS_LENGTH_IN_LINE)
                     << std::endl;
          }
          else
          {
              stream << std::left << std::setw(6) << keywords_->GetRecordName()
                     << std::left << std::setw(2) << " "
                     << std::right << std::setw(2) << i
                     << std::left << std::setw(70) << keywords_->GetKeywords().substr(MAX_KEYWORDS_LENGTH_IN_LINE*(i-1), keywords_->GetKeywords().length()-MAX_KEYWORDS_LENGTH_IN_LINE*(i-1))
                     << std::endl;
          }
      }
  }
  else
  {
      stream << std::right << std::setw(2) << " "
             << std::left << std::setw(70) << keywords_->GetKeywords()
             << std::endl;
  }
}

void PdbFile::ResolveExperimentalDataCards(std::ostream& stream)
{
  const int MAX_EXPDTA_LENGTH_IN_LINE = 70;
  stream << std::left << std::setw(6) << experimental_data_->GetRecordName()
         << std::left << std::setw(2) << " ";
  if((int)experimental_data_->GetExperimentalData().length() > MAX_EXPDTA_LENGTH_IN_LINE)
  {
      stream << std::left << std::setw(70) << experimental_data_->GetExperimentalData().substr(0,MAX_EXPDTA_LENGTH_IN_LINE)
             << std::endl;

      int counter = ceil((double)(experimental_data_->GetExperimentalData().length()) / MAX_EXPDTA_LENGTH_IN_LINE);
      for(int i = 2; i <= counter; i++)
      {
          if(i != counter)
          {
              stream << std::left << std::setw(6) << experimental_data_->GetRecordName()
                     << std::left << std::setw(2) << " "
                     << std::right << std::setw(2) << i
                     << std::left << std::setw(70) << experimental_data_->GetExperimentalData().substr(MAX_EXPDTA_LENGTH_IN_LINE*(i-1), MAX_EXPDTA_LENGTH_IN_LINE)
                     << std::endl;
          }
          else
          {
              stream << std::left << std::setw(6) << experimental_data_->GetRecordName()
                     << std::left << std::setw(2) << " "
                     << std::right << std::setw(2) << i
                     << std::left << std::setw(70) << experimental_data_->GetExperimentalData().substr(MAX_EXPDTA_LENGTH_IN_LINE*(i-1), experimental_data_->GetExperimentalData().length()-MAX_EXPDTA_LENGTH_IN_LINE*(i-1))
                     << std::endl;
          }
      }
  }
  else
  {
      stream << std::right << std::setw(2) << " "
             << std::left << std::setw(70) << experimental_data_->GetExperimentalData()
             << std::endl;
  }
}

void PdbFile::ResolveNumModelCard(std::ostream& stream)
{
    stream << std::left << std::setw(6) << number_of_models_->GetRecordName()
           << std::left << std::setw(4) << " ";
    if(number_of_models_->GetNumberOfModels() != gmml::iNotSet)
        stream << std::right << std::setw(4) << number_of_models_->GetNumberOfModels();
    else
        stream << std::right << std::setw(4) << " ";
    stream << std::left << std::setw(66) << " "
           << std::endl;
}

void PdbFile::ResolveModelTypeCards(std::ostream& stream)
{
    stream << std::left << std::setw(6) << model_type_->GetRecordName()
           << std::left << std::setw(2) << " ";
    std::stringstream ss;
    for(std::vector<std::string>::iterator it = model_type_->GetComments().begin(); it != model_type_->GetComments().end(); it++)
    {
        if(it != model_type_->GetComments().end() - 1)
        {
            ss << (*it) << "; ";
        }
        else
        {
            ss << (*it);
        }
    }
    if(ss.str().length() > 70)
    {
        stream << std::right << std::setw(2) << " "
               << std::left << std::setw(70) << ss.str().substr(0,70)
               << std::endl;

        int counter = ceil((double)(ss.str().length()) / 70);
        for(int i = 2; i <= counter; i++)
        {
            if(i != counter)
            {
                stream << std::left << std::setw(6) << model_type_->GetRecordName()
                       << std::left << std::setw(2) << " "
                       << std::right << std::setw(2) << i
                       << std::left << std::setw(70) << ss.str().substr(70*(i-1), 70)
                       << std::endl;
            }
            else
            {
                stream << std::left << std::setw(6) << model_type_->GetRecordName()
                       << std::left << std::setw(2) << " "
                       << std::right << std::setw(2) << i
                       << std::left << std::setw(70) << ss.str().substr(70*(i-1), ss.str().length() - (i-1)*70)
                       << std::endl;
            }
        }
    }
    else
    {
        stream << std::right << std::setw(2) << " "
               << std::left << std::setw(70) << ss.str()
               << std::endl;
    }
}

void PdbFile::ResolveAuthorCards(std::ostream& stream)
{
  const int MAX_AUTHOR_LENGTH_IN_LINE = 70;
  stream << std::left << std::setw(6) << author_->GetRecordName()
         << std::left << std::setw(2) << " ";
  if((int)author_->GetAuthor().length() > MAX_AUTHOR_LENGTH_IN_LINE)
  {
      stream << std::left << std::setw(70) << author_->GetAuthor().substr(0,MAX_AUTHOR_LENGTH_IN_LINE)
             << std::endl;

      int counter = ceil((double)(author_->GetAuthor().length()) / MAX_AUTHOR_LENGTH_IN_LINE);
      for(int i = 2; i <= counter; i++)
      {
          if(i != counter)
          {
              stream << std::left << std::setw(6) << author_->GetRecordName()
                     << std::left << std::setw(2) << " "
                     << std::right << std::setw(2) << i
                     << std::left << std::setw(70) << author_->GetAuthor().substr(MAX_AUTHOR_LENGTH_IN_LINE*(i-1), MAX_AUTHOR_LENGTH_IN_LINE)
                     << std::endl;
          }
          else
          {
              stream << std::left << std::setw(6) << author_->GetRecordName()
                     << std::left << std::setw(2) << " "
                     << std::right << std::setw(2) << i
                     << std::left << std::setw(70) << author_->GetAuthor().substr(MAX_AUTHOR_LENGTH_IN_LINE*(i-1), author_->GetAuthor().length()-MAX_AUTHOR_LENGTH_IN_LINE*(i-1))
                     << std::endl;
          }
      }
  }
  else
  {
      stream << std::right << std::setw(2) << " "
             << std::left << std::setw(70) << author_->GetAuthor()
             << std::endl;
  }
}

void PdbFile::ResolveRevisionDataCards(std::ostream& stream)
{
  int REVDAT_COUNT = 1;
  int REVDAT_NUM = 0;
  RevisionDataCardVector revision_data_cards = revision_data_->GetRevisionDataCards();
  for (RevisionDataCardVector::iterator it = revision_data_cards.begin(); it != revision_data_cards.end(); it++)
  {
    if (REVDAT_NUM != (*it)->GetModificationNumber())
    {
      REVDAT_COUNT = 1;
    }
    if (REVDAT_COUNT == 1)
    {
      stream << std::left << std::setw(6) << (*it)->GetRecordName()
             << std::left << std::setw(1) << " "
             << std::right << std::setw(3) << (*it)->GetModificationNumber()
             << std::right << std::setw(3) << " "
             << std::left << std::setw(9) << (*it)->GetModificationDate()
             << std::left << std::setw(1) << " "
             << std::left << std::setw(4) << (*it)->GetModificationID()
             << std::left << std::setw(4) << " "
             << std::left << std::setw(1) << (*it)->GetModificationType()
             << std::left << std::setw(7) << " "
             << std::left << std::setw(27) << (*it)->GetModificationDetails()
             << std::endl;
    }
    else if (REVDAT_NUM == (*it)->GetModificationNumber())
    {
      stream << std::left << std::setw(6) << (*it)->GetRecordName()
             << std::left << std::setw(1) << " "
             << std::right << std::setw(3) << (*it)->GetModificationNumber()
             << std::right << std::setw(2) << REVDAT_COUNT
             << std::left << std::setw(1) << " "
             << std::left << std::setw(9) << (*it)->GetModificationDate()
             << std::left << std::setw(1) << " "
             << std::left << std::setw(4) << (*it)->GetModificationID()
             << std::left << std::setw(4) << " "
             << std::left << std::setw(1) << (*it)->GetModificationType()
             << std::left << std::setw(7) << " "
             << std::left << std::setw(27) << (*it)->GetModificationDetails()
             << std::endl;
      REVDAT_COUNT ++;
    }

    REVDAT_NUM = (*it)->GetModificationNumber();
  }

}

void PdbFile::ResolveSupersededEntriesCards(std::ostream& stream)
{
	stream << "";
}

void PdbFile::ResolveJournalCards(std::ostream& stream)
{
  const int MAX_JOURNAL_LENGTH_IN_LINE = 67;
  stream << std::left << std::setw(6) << journal_->GetRecordName()
         << std::left << std::setw(6) << " ";
  if((int)journal_->GetText().length() > MAX_JOURNAL_LENGTH_IN_LINE)
  {
      stream << std::left << std::setw(67) << journal_->GetText().substr(0,MAX_JOURNAL_LENGTH_IN_LINE)
             << std::endl;

      int counter = ceil((double)(journal_->GetText().length()) / MAX_JOURNAL_LENGTH_IN_LINE);
      for(int i = 2; i <= counter; i++)
      {
          if(i != counter)
          {
              stream << std::left << std::setw(6) << journal_->GetRecordName()
                     << std::left << std::setw(6) << " "
                     << std::left << std::setw(67) << journal_->GetText().substr(MAX_JOURNAL_LENGTH_IN_LINE*(i-1), MAX_JOURNAL_LENGTH_IN_LINE)
                     << std::endl;
          }
          else
          {
              stream << std::left << std::setw(6) << journal_->GetRecordName()
                     << std::left << std::setw(6) << " "
                     << std::left << std::setw(67) << journal_->GetText().substr(MAX_JOURNAL_LENGTH_IN_LINE*(i-1), journal_->GetText().length()-MAX_JOURNAL_LENGTH_IN_LINE*(i-1))
                     << std::endl;
          }
      }
  }
  else
  {
      stream << std::right << std::setw(6) << " "
             << std::left << std::setw(67) << journal_->GetText()
             << std::endl;
  }
}

void PdbFile::ResolveRemarkCards(std::ostream& stream)
{
  stream << std::left << remark_cards_->GetRemarks();
}

void PdbFile::ResolveDatabaseReferenceCards(std::ostream& stream)
{
  DatabaseReferenceVector database_reference_cards = database_reference_->GetDatabaseReferences();
  for (DatabaseReferenceVector::iterator it = database_reference_cards.begin(); it != database_reference_cards.end(); it++)
  {
    if ((*it)->GetRecordName() == "DBREF ")
    {
      stream << std::left << std::setw(6) << (*it)->GetRecordName()
             << std::left << std::setw(1) << " "
             << std::left << std::setw(4) << (*it)->GetIDCode()
             << std::left << std::setw(1) << " "
             << std::left << std::setw(1) << (*it)->GetChainID()
             << std::left << std::setw(1) << " "
             << std::right << std::setw(4) << (*it)->GetSeqBegin()
             << std::right << std::setw(1) << (*it)->GetInsertBegin()
             << std::left << std::setw(1) << " "
             << std::right << std::setw(4) << (*it)->GetSeqEnd()
             << std::right << std::setw(1) << (*it)->GetInsertEnd()
             << std::left << std::setw(1) << " "
             << std::left << std::setw(6) << (*it)->GetDatabase()
             << std::left << std::setw(1) << " "
             << std::left << std::setw(8) << (*it)->GetDatabaseAccession()
             << std::left << std::setw(1) << " "
             << std::left << std::setw(12) << (*it)->GetDatabaseIDCode()
             << std::left << std::setw(1) << " "
             << std::right << std::setw(5) << (*it)->GetDatabaseSeqBegin()
             << std::right << std::setw(1) << (*it)->GetDatabaseInsBegin()
             << std::left << std::setw(1) << " "
             << std::right << std::setw(5) << (*it)->GetDatabaseSeqEnd()
             << std::right << std::setw(1) << (*it)->GetDatabaseInsEnd()
             << std::endl;
    }
    else if (((*it)->GetRecordName() == "DBREF1"))
    {
      stream << std::left << std::setw(6) << (*it)->GetRecordName()
             << std::left << std::setw(1) << " "
             << std::left << std::setw(4) << (*it)->GetIDCode()
             << std::left << std::setw(1) << " "
             << std::left << std::setw(1) << (*it)->GetChainID()
             << std::left << std::setw(1) << " "
             << std::right << std::setw(4) << (*it)->GetSeqBegin()
             << std::right << std::setw(1) << (*it)->GetInsertBegin()
             << std::left << std::setw(1) << " "
             << std::right << std::setw(4) << (*it)->GetSeqEnd()
             << std::right << std::setw(1) << (*it)->GetInsertEnd()
             << std::left << std::setw(1) << " "
             << std::left << std::setw(6) << (*it)->GetDatabase()
             << std::left << std::setw(16) << " "
             << std::left << std::setw(15) << (*it)->GetDatabaseIDCode()
             << std::endl
             << std::left << std::setw(6) << "DBREF2"
             << std::left << std::setw(1) << " "
             << std::left << std::setw(4) << (*it)->GetIDCode()
             << std::left << std::setw(1) << " "
             << std::left << std::setw(1) << (*it)->GetChainID()
             << std::left << std::setw(6) << " "
             << std::left << std::setw(22) << (*it)->GetDatabaseAccession()
             << std::left << std::setw(5) << " "
             << std::right << std::setw(10) << (*it)->GetDatabaseSeqBegin()
             << std::left << std::setw(2) << " "
             << std::right << std::setw(10) << (*it)->GetDatabaseSeqEnd()
             << std::endl;
    }

  }
}

void PdbFile::ResolveSequenceAdvancedCards(std::ostream& stream)
{
  SequenceAdvancedCardVector sequence_advanced_cards = sequence_advanced_->GetSequenceAdvancedCards();
  for (SequenceAdvancedCardVector::iterator it = sequence_advanced_cards.begin(); it != sequence_advanced_cards.end(); it++)
  {
    stream << std::left << std::setw(6) << (*it)->GetRecordName()
           << std::left << std::setw(1) << " "
           << std::right << std::setw(3) << (*it)->GetIdentifierCode()
           << std::left << std::setw(1) << " "
           << std::left << std::setw(3) << (*it)->GetResidueName()
           << std::left << std::setw(1) << " "
           << std::left << std::setw(1) << (*it)->GetChainId()
           << std::left << std::setw(1) << " "
           << std::right << std::setw(4) << (*it)->GetSequenceNumber()
           << std::left << std::setw(1) << (*it)->GetInsertionCode()
           << std::left << std::setw(1) << " "
           << std::left << std::setw(4) << (*it)->GetDatabase()
           << std::left << std::setw(1) << " "
           << std::left << std::setw(9) << (*it)->GetDatabaseAccession()
           << std::left << std::setw(1) << " "
           << std::right << std::setw(3) << (*it)->GetDatabaseResidue()
           << std::left << std::setw(1) << " "
           << std::right << std::setw(5) << (*it)->GetDatabaseSequence()
           << std::left << std::setw(1) << " "
           << std::left << std::setw(30) << (*it)->GetConflict()
           << std::endl;
  }
}

void PdbFile::ResolveSequenceResidueCards(std::ostream& stream)
{
    PdbFileSpace::PdbResidueSequenceSection::ResidueSequenceCardMap residue_sequence_map = residues_sequence_->GetResidueSequenceChain();
    for(PdbFileSpace::PdbResidueSequenceSection::ResidueSequenceCardMap::iterator it = residue_sequence_map.begin(); it != residue_sequence_map.end(); it++)
    {

        PdbResidueSequenceCard* residue_sequence = (*it).second;
        int serial_number = 1;
        const int MAX_RESIDUE_IN_SINGLE_LINE = 13;
        std::vector<std::string> residue_names = residue_sequence->GetResidueNames();
        if(residue_sequence->GetNumberOfResidues() <= MAX_RESIDUE_IN_SINGLE_LINE)
        {
            std::stringstream ss;
            for(std::vector<std::string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
            {
                ss << std::right << std::setw(1) << " "
                   << std::right << std::setw(3) << (*it1);
            }
            for(int i = 0; i < MAX_RESIDUE_IN_SINGLE_LINE - residue_sequence->GetNumberOfResidues(); i++)
            {
                ss << std::right << std::setw(1) << " "
                   << std::right << std::setw(3) << " ";
            }
            stream << std::left << std::setw(6) << residues_sequence_->GetRecordName()
                   << std::left << std::setw(1) << " "
                   << std::right << std::setw(3) << serial_number
                   << std::right << std::setw(1) << " "
                   << std::right << std::setw(1) << residue_sequence->GetChainId()
                   << std::right << std::setw(1) << " ";
            if(residue_sequence->GetNumberOfResidues() != gmml::iNotSet)
                stream << std::right << std::setw(4) << residue_sequence->GetNumberOfResidues();
            else
                stream << std::right << std::setw(4) << " ";
            stream << std::right << std::setw(1) << " "
                   << std::right << std::setw(52) << ss.str()
                   << std::right << std::setw(10) << " "
                   << std::endl;
        }
        else
        {
            int number_of_lines = ceil((double)(residue_sequence->GetNumberOfResidues()) / MAX_RESIDUE_IN_SINGLE_LINE);
            for(int i = 0; i < number_of_lines; i++)
            {
                std::stringstream ss;
                if(i != number_of_lines - 1)
                {
                    for(std::vector<std::string>::iterator it1 = residue_names.begin() + i * MAX_RESIDUE_IN_SINGLE_LINE;
                        it1 != residue_names.begin() + (i+1) * MAX_RESIDUE_IN_SINGLE_LINE; it1++)
                    {
                        ss << std::right << std::setw(1) << " "
                           << std::right << std::setw(3) << (*it1);
                    }
                }
                else
                {
                    for(std::vector<std::string>::iterator it1 = residue_names.begin() + i * MAX_RESIDUE_IN_SINGLE_LINE;
                        it1 != residue_names.end(); it1++)
                    {
                        ss << std::right << std::setw(1) << " "
                           << std::right << std::setw(3) << (*it1);
                    }
                    for(int i = 0; i < MAX_RESIDUE_IN_SINGLE_LINE - residue_sequence->GetNumberOfResidues() % MAX_RESIDUE_IN_SINGLE_LINE; i++)
                    {
                        ss << std::right << std::setw(1) << " "
                           << std::right << std::setw(3) << " ";
                    }
                }
                stream << std::left << std::setw(6) << residues_sequence_->GetRecordName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << serial_number
                       << std::right << std::setw(1) << " "
                       << std::right << std::setw(1) << residue_sequence->GetChainId()
                       << std::right << std::setw(1) << " ";
                if(residue_sequence->GetNumberOfResidues() != gmml::iNotSet)
                    stream << std::right << std::setw(4) << residue_sequence->GetNumberOfResidues();
                else
                    stream << std::right << std::setw(4) << " ";
                stream << std::right << std::setw(1) << " "
                       << std::right << std::setw(52) << ss.str()
                       << std::right << std::setw(10) << " "
                       << std::endl;
                serial_number++;
            }
        }

    }
}

void PdbFile::ResolveModificationResidueCards(std::ostream& stream)
{
    PdbFileSpace::PdbResidueModificationSection::ResidueModificationCardMap residue_modification_cards_map = residue_modification_cards_->GetResidueModificationCards();
    for(PdbFileSpace::PdbResidueModificationSection::ResidueModificationCardMap::iterator it = residue_modification_cards_map.begin(); it != residue_modification_cards_map.end(); it++)
    {
        PdbResidueModificationCard* residue_modification_cards = (*it).second;
        stream << std::left << std::setw(6) << residue_modification_cards_->GetRecordName()
               << std::left << std::setw(1) << " "
               << std::right << std::setw(4) << residue_modification_cards->GetIdCode()
               << std::left << std::setw(1) << " "
               << std::right << std::setw(3) << residue_modification_cards->GetResidueName()
               << std::left << std::setw(1) << " "
               << std::right << std::setw(1) << residue_modification_cards->GetChainId()
               << std::left << std::setw(1) << " ";
        if(residue_modification_cards->GetSequenceNumber() != gmml::iNotSet)
            stream << std::right << std::setw(4) << residue_modification_cards->GetSequenceNumber();
        else
            stream << std::right << std::setw(4) << " ";
        stream << std::right << std::setw(1) << residue_modification_cards->GetInsertionCode()
               << std::left << std::setw(1) << " "
               << std::right << std::setw(3) << residue_modification_cards->GetStandardResidueName()
               << std::left << std::setw(2) << " "
               << std::left << std::setw(41) << residue_modification_cards->GetDscr()
               << std::left << std::setw(10) << " "
               << std::endl;
    }
}

void PdbFile::ResolveHeterogenCards(std::ostream& stream)
{
    PdbFileSpace::PdbHeterogenSection::HeterogenCardMap heterogen_map = heterogen_cards_->GetHeterogenCards();
    for(PdbFileSpace::PdbHeterogenSection::HeterogenCardMap::iterator it = heterogen_map.begin(); it != heterogen_map.end(); it++)
    {

        PdbHeterogenCard* heterogen = (*it).second;
        stream << std::left << std::setw(6) << heterogen_cards_->GetRecordName()
               << std::left << std::setw(1) << " "
               << std::right << std::setw(3) << heterogen->GetHeterogenId()
               << std::left << std::setw(2) << " "
               << std::right << std::setw(1) << heterogen->GetChainId();
        if(heterogen->GetSequenceNumber() != gmml::iNotSet)
            stream << std::right << std::setw(4) << heterogen->GetSequenceNumber();
        else
            stream << std::right << std::setw(4) << " ";
        stream << std::right << std::setw(1) << heterogen->GetInsertionCode()
               << std::left << std::setw(2) << " ";
        if(heterogen->GetNumberOfHeterogenAtoms() != gmml::iNotSet)
            stream << std::right << std::setw(5) << heterogen->GetNumberOfHeterogenAtoms();
        else
            stream << std::right << std::setw(5) << " ";
        stream << std::left << std::setw(5) << " "
               << std::left << std::setw(40) << heterogen->GetDscr()
               << std::left << std::setw(10) << " "
               << std::endl;
    }
}

void PdbFile::ResolveHeterogenNameCards(std::ostream& stream)
{
    PdbFileSpace::PdbHeterogenNameSection::HeterogenNameCardMap heterogen_name_map = heterogen_name_cards_->GetHeterogenNameCards();
    for(PdbFileSpace::PdbHeterogenNameSection::HeterogenNameCardMap::iterator it = heterogen_name_map.begin(); it != heterogen_name_map.end(); it++)
    {
        PdbHeterogenNameCard* heterogen_name = (*it).second;
        const int MAX_NAME_LENGTH_IN_LINE = 55;
        if((int)heterogen_name->GetHeterogenName().length() > MAX_NAME_LENGTH_IN_LINE)
        {
            stream << std::left << std::setw(6) << heterogen_name_cards_->GetRecordName()
                   << std::left << std::setw(2) << " "
                   << std::right << std::setw(2) << " "
                   << std::left << std::setw(1) << " "
                   << std::right << std::setw(3) << heterogen_name->GetHeterogenIdentifier()
                   << std::left << std::setw(1) << " "
                   << std::left << std::setw(55) << heterogen_name->GetHeterogenName().substr(0,MAX_NAME_LENGTH_IN_LINE)
                   << std::left << std::setw(10) << " "
                   << std::endl;
            int counter = ceil((double)(heterogen_name->GetHeterogenName().length()) / MAX_NAME_LENGTH_IN_LINE);
            for(int i = 2; i <= counter; i++)
            {
                if(i != counter)
                {
                    stream << std::left << std::setw(6) << heterogen_name_cards_->GetRecordName()
                           << std::left << std::setw(2) << " "
                           << std::right << std::setw(2) << i
                           << std::left << std::setw(1) << " "
                           << std::right << std::setw(3) << heterogen_name->GetHeterogenIdentifier()
                           << std::left << std::setw(1) << " "
                           << std::left << std::setw(55) << heterogen_name->GetHeterogenName().substr(MAX_NAME_LENGTH_IN_LINE*(i-1),MAX_NAME_LENGTH_IN_LINE)
                           << std::left << std::setw(10) << " "
                           << std::endl;
                }
                else
                {
                    stream << std::left << std::setw(6) << heterogen_name_cards_->GetRecordName()
                           << std::left << std::setw(2) << " "
                           << std::right << std::setw(2) << i
                           << std::left << std::setw(1) << " "
                           << std::right << std::setw(3) << heterogen_name->GetHeterogenIdentifier()
                           << std::left << std::setw(1) << " "
                           << std::left << std::setw(55) << heterogen_name->GetHeterogenName().substr(MAX_NAME_LENGTH_IN_LINE*(i-1),heterogen_name->GetHeterogenName().length()-MAX_NAME_LENGTH_IN_LINE*(i-1))
                           << std::left << std::setw(10) << " "
                           << std::endl;
                }
            }
        }
        else
        {
            stream << std::left << std::setw(6) << heterogen_name_cards_->GetRecordName()
                   << std::left << std::setw(2) << " "
                   << std::right << std::setw(2) << " "
                   << std::left << std::setw(1) << " "
                   << std::right << std::setw(3) << heterogen_name->GetHeterogenIdentifier()
                   << std::left << std::setw(1) << " "
                   << std::left << std::setw(55) << heterogen_name->GetHeterogenName()
                   << std::left << std::setw(10) << " "
                   << std::endl;
        }
    }
}

void PdbFile::ResolveHeterogenSynonymCards(std::ostream& stream)
{
    PdbFileSpace::PdbHeterogenSynonymSection::HeterogenSynonymCardMap heterogen_synonym_map = heterogen_synonym_cards_->GetHeterogensSynonymCards();
    for(PdbFileSpace::PdbHeterogenSynonymSection::HeterogenSynonymCardMap::iterator it = heterogen_synonym_map.begin(); it != heterogen_synonym_map.end(); it++)
    {
        PdbHeterogenSynonymCard* heterogen_synonym_card = (*it).second;
        std::stringstream ss;
        std::vector<std::string> synonyms = heterogen_synonym_card->GetHeterogenSynonymCards();
        for(std::vector<std::string>::iterator it = synonyms.begin(); it != synonyms.end(); it++)
        {
            if(it != synonyms.end() - 1)
            {
                ss << (*it) << "; ";
            }
            else
            {
                ss << (*it) << ";";
            }
        }
        const int MAX_SYNONYM_LENGTH_IN_LINE = 55;
        if((int)ss.str().length() > MAX_SYNONYM_LENGTH_IN_LINE)
        {
            stream << std::left << std::setw(6) << heterogen_synonym_cards_->GetRecordName()
                   << std::left << std::setw(2) << " "
                   << std::right << std::setw(2) << " "
                   << std::left << std::setw(1) << " "
                   << std::right << std::setw(3) << heterogen_synonym_card->GetHeterogenIdentifier()
                   << std::left << std::setw(1) << " "
                   << std::left << std::setw(55) << ss.str().substr(0,MAX_SYNONYM_LENGTH_IN_LINE)
                   << std::left << std::setw(10) << " "
                   << std::endl;
            int counter = ceil((double)(ss.str().length()) / MAX_SYNONYM_LENGTH_IN_LINE);
            for(int i = 2; i <= counter; i++)
            {
                if(i != counter)
                {
                    stream << std::left << std::setw(6) << heterogen_synonym_cards_->GetRecordName()
                           << std::left << std::setw(2) << " "
                           << std::right << std::setw(2) << i
                           << std::left << std::setw(1) << " "
                           << std::right << std::setw(3) << heterogen_synonym_card->GetHeterogenIdentifier()
                           << std::left << std::setw(1) << " "
                           << std::left << std::setw(55) << ss.str().substr(MAX_SYNONYM_LENGTH_IN_LINE*(i-1),MAX_SYNONYM_LENGTH_IN_LINE)
                           << std::left << std::setw(10) << " "
                           << std::endl;
                }
                else
                {
                    stream << std::left << std::setw(6) << heterogen_synonym_cards_->GetRecordName()
                           << std::left << std::setw(2) << " "
                           << std::right << std::setw(2) << i
                           << std::left << std::setw(1) << " "
                           << std::right << std::setw(3) << heterogen_synonym_card->GetHeterogenIdentifier()
                           << std::left << std::setw(1) << " "
                           << std::left << std::setw(55) << ss.str().substr(MAX_SYNONYM_LENGTH_IN_LINE*(i-1),ss.str().length()-MAX_SYNONYM_LENGTH_IN_LINE*(i-1))
                           << std::left << std::setw(10) << " "
                           << std::endl;
                }
            }
        }
        else
        {
            stream << std::left << std::setw(6) << heterogen_synonym_cards_->GetRecordName()
                   << std::left << std::setw(2) << " "
                   << std::right << std::setw(2) << " "
                   << std::left << std::setw(1) << " "
                   << std::right << std::setw(3) << heterogen_synonym_card->GetHeterogenIdentifier()
                   << std::left << std::setw(1) << " "
                   << std::left << std::setw(55) << ss.str()
                   << std::left << std::setw(10) << " "
                   << std::endl;
        }
    }
}

void PdbFile::ResolveFormulaCards(std::ostream& stream)
{
    PdbFileSpace::PdbFormulaSection::FormulaCardMap formula_map = formulas_->GetFormulaCards();
    for(PdbFileSpace::PdbFormulaSection::FormulaCardMap::iterator it = formula_map.begin(); it != formula_map.end(); it++)
    {
        PdbFormulaCard* formula = (*it).second;
        const int MAX_NAME_LENGTH_IN_LINE = 51;
        if((int)formula->GetChemicalFormula().length() > MAX_NAME_LENGTH_IN_LINE)
        {
            stream << std::left << std::setw(6) << formulas_->GetRecordName()
                   << std::left << std::setw(2) << " ";
            if(formula->GetComponentNumber() != gmml::iNotSet)
                stream << std::right << std::setw(2) << formula->GetComponentNumber();
            else
                stream << std::right << std::setw(2) << " ";
            stream << std::left << std::setw(2) << " "
                   << std::right << std::setw(3) << formula->GetHeterogenIdentifier()
                   << std::left << std::setw(1) << " "
                   << std::right << std::setw(2) << " "
                   << std::left << std::setw(1) << " "
                   << std::left << std::setw(51) << formula->GetChemicalFormula().substr(0,MAX_NAME_LENGTH_IN_LINE)
                   << std::left << std::setw(10) << " "
                   << std::endl;
            int counter = ceil((double)(formula->GetChemicalFormula().length()) / MAX_NAME_LENGTH_IN_LINE);
            for(int i = 2; i <= counter; i++)
            {
                if(i != counter)
                {
                    stream << std::left << std::setw(6) << formulas_->GetRecordName()
                           << std::left << std::setw(2) << " ";
                    if(formula->GetComponentNumber() != gmml::iNotSet)
                        stream << std::right << std::setw(2) << formula->GetComponentNumber();
                    else
                        stream << std::right << std::setw(2) << " ";
                    stream << std::left << std::setw(2) << " "
                           << std::right << std::setw(3) << formula->GetHeterogenIdentifier()
                           << std::left << std::setw(1) << " "
                           << std::right << std::setw(2) << i
                           << std::left << std::setw(1) << " "
                           << std::left << std::setw(51) << formula->GetChemicalFormula().substr(MAX_NAME_LENGTH_IN_LINE*(i-1),MAX_NAME_LENGTH_IN_LINE)
                           << std::left << std::setw(10) << " "
                           << std::endl;
                }
                else
                {
                    stream << std::left << std::setw(6) << formulas_->GetRecordName()
                           << std::left << std::setw(2) << " ";
                    if(formula->GetComponentNumber()!= gmml::iNotSet)
                        stream << std::right << std::setw(2) << formula->GetComponentNumber();
                    else
                        stream << std::right << std::setw(2) << " ";
                    stream << std::left << std::setw(2) << " "
                           << std::right << std::setw(3) << formula->GetHeterogenIdentifier()
                           << std::left << std::setw(1) << " "
                           << std::right << std::setw(2) << i
                           << std::left << std::setw(1) << " "
                           << std::left << std::setw(51) << formula->GetChemicalFormula().substr(MAX_NAME_LENGTH_IN_LINE*(i-1),formula->GetChemicalFormula().length()-MAX_NAME_LENGTH_IN_LINE*(i-1))
                           << std::left << std::setw(10) << " "
                           << std::endl;
                }
            }
        }
        else
        {
            stream << std::left << std::setw(6) << formulas_->GetRecordName()
                   << std::left << std::setw(2) << " ";
            if(formula->GetComponentNumber() != gmml::iNotSet)
                stream << std::right << std::setw(2) << formula->GetComponentNumber();
            else
                stream << std::right << std::setw(2) << " ";
            stream << std::left << std::setw(2) << " "
                   << std::right << std::setw(3) << formula->GetHeterogenIdentifier()
                   << std::left << std::setw(1) << " "
                   << std::right << std::setw(2) << " "
                   << std::left << std::setw(1) << " "
                   << std::left << std::setw(51) << formula->GetChemicalFormula().substr(0,MAX_NAME_LENGTH_IN_LINE)
                   << std::left << std::setw(10) << " "
                   << std::endl;
        }
    }
}

void PdbFile::ResolveHelixCards(std::ostream& stream)
{
    PdbFileSpace::PdbHelixSection::HelixCardMap helix_map = helix_cards_->GetHelixCards();
    int counter = helix_map.size();
    int serial_number = 1;
    while(serial_number <= counter)
    {
        for(PdbFileSpace::PdbHelixSection::HelixCardMap::iterator it = helix_map.begin(); it != helix_map.end(); it++)
        {

            PdbHelixCard* helix = (*it).second;
            PdbHelixCard::HelixResidueVector helix_residues = helix->GetHelixResidues();
            if(helix->GetHelixSerialNumber() == serial_number)
            {
                stream << std::left << std::setw(6) << helix_cards_->GetRecordName()
                       << std::left << std::setw(1) << " ";
                if(helix->GetHelixSerialNumber() != gmml::iNotSet)
                    stream << std::right << std::setw(3) << helix->GetHelixSerialNumber();
                else
                    stream << std::right << std::setw(3) << " ";
                stream << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << helix->GetHelixId()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << helix_residues.at(0)->GetResidueName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(1) << helix_residues.at(0)->GetResidueChainId()
                       << std::left << std::setw(1) << " ";
                if(helix_residues.at(0)->GetResidueSequenceNumber() != gmml::iNotSet)
                    stream << std::right << std::setw(4) << helix_residues.at(0)->GetResidueSequenceNumber();
                else
                    stream << std::right << std::setw(4) << " ";
                stream << std::right << std::setw(1) << helix_residues.at(0)->GetResidueInsertionCode()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << helix_residues.at(1)->GetResidueName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(1) << helix_residues.at(1)->GetResidueChainId()
                       << std::left << std::setw(1) << " ";
                if(helix_residues.at(1)->GetResidueSequenceNumber() != gmml::iNotSet)
                    stream << std::right << std::setw(4) << helix_residues.at(1)->GetResidueSequenceNumber();
                else
                    stream << std::right << std::setw(4) << " ";
                stream << std::right << std::setw(1) << helix_residues.at(1)->GetResidueInsertionCode();
                if(helix->GetHelixClass() != UnknownHelix)
                    stream << std::right << std::setw(2) << helix->GetHelixClass();
                else
                    stream << std::right << std::setw(2) << " ";
                stream << std::left << std::setw(30) << helix->GetComment()
                       << std::left << std::setw(1) << " ";
                if(helix->GetHelixLength() != gmml::dNotSet)
                    stream << std::right << std::setw(5) << helix->GetHelixLength();
                else
                    stream << std::right << std::setw(5) << " ";
                stream << std::left << std::setw(4) << " "
                       << std::endl;
                break;
            }
        }
        serial_number++;
    }
}

void PdbFile::ResolveSheetCards(std::ostream& stream)
{
    PdbFileSpace::PdbSheetSection::SheetCardMap sheet_map = sheet_cards_->GetSheets();
    for(PdbFileSpace::PdbSheetSection::SheetCardMap::iterator it = sheet_map.begin(); it != sheet_map.end(); it++)
    {
        PdbSheetCard* sheet = (*it).second;
        PdbSheetCard::SheetStrandVector strands = sheet->GetStrands();
        int serial_number = 1;
        for(PdbSheetCard::SheetStrandVector::iterator it1 = strands.begin(); it1 != strands.end(); it1++)
        {
            PdbSheetStrand* strand = (*it1);
            PdbSheetStrand::SheetStrandResidueVector strand_residues = strand->GetStrandResidues();
            if(strand->GetSense() != 0)
            {
                stream << std::left << std::setw(6) << sheet_cards_->GetRecordName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << serial_number
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << sheet->GetSheetId();
                if(sheet->GetNumberOfStrands() != gmml::iNotSet)
                    stream << std::right << std::setw(2) << sheet->GetNumberOfStrands();
                else
                    stream << std::right << std::setw(2) << " ";
                stream << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << strand_residues.at(0)->GetResidueName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(1) << strand_residues.at(0)->GetResidueChainId();
                if(strand_residues.at(0)->GetResidueSequenceNumber() != gmml::iNotSet)
                    stream << std::right << std::setw(4) << strand_residues.at(0)->GetResidueSequenceNumber();
                else
                    stream << std::right << std::setw(4) << " ";
                stream << std::right << std::setw(1) << strand_residues.at(0)->GetResidueInsertionCode()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << strand_residues.at(1)->GetResidueName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(1) << strand_residues.at(1)->GetResidueChainId();
                if(strand_residues.at(1)->GetResidueSequenceNumber() != gmml::iNotSet)
                    stream << std::right << std::setw(4) << strand_residues.at(1)->GetResidueSequenceNumber();
                else
                    stream << std::right << std::setw(4) << " ";
                stream << std::right << std::setw(1) << strand_residues.at(1)->GetResidueInsertionCode();
                if(strand->GetSense() != UnknownStrand)
                    stream << std::right << std::setw(2) << strand->GetSense();
                else
                    stream << std::right << std::setw(2) << " ";
                stream << std::left << std::setw(2) << " "
                       << std::left << std::setw(3) << strand->GetCurrentAtom()
                       << std::right << std::setw(3) << strand_residues.at(2)->GetResidueName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(1) << strand_residues.at(2)->GetResidueChainId();
                if(strand_residues.at(2)->GetResidueSequenceNumber() != gmml::iNotSet)
                    stream << std::right << std::setw(4) << strand_residues.at(2)->GetResidueSequenceNumber();
                else
                    stream << std::right << std::setw(4) << " ";
                stream << std::right << std::setw(1) << strand_residues.at(2)->GetResidueInsertionCode()
                       << std::left << std::setw(1) << " "
                       << std::left << std::setw(4) << strand->GetPreviousAtom()
                       << std::right << std::setw(3) << strand_residues.at(3)->GetResidueName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(1) << strand_residues.at(3)->GetResidueChainId();
                if(strand_residues.at(3)->GetResidueSequenceNumber() != gmml::iNotSet)
                    stream << std::right << std::setw(4) << strand_residues.at(3)->GetResidueSequenceNumber();
                else
                    stream << std::right << std::setw(4) << " ";
                stream << std::right << std::setw(1) << strand_residues.at(3)->GetResidueInsertionCode()
                       << std::left << std::setw(10) << " "
                       << std::endl;
            }
            else
            {
                stream << std::left << std::setw(6) << sheet_cards_->GetRecordName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << serial_number
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << sheet->GetSheetId();
                if(sheet->GetNumberOfStrands() != gmml::iNotSet)
                    stream << std::right << std::setw(2) << sheet->GetNumberOfStrands();
                else
                    stream << std::right << std::setw(2) << " ";
                stream << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << strand_residues.at(0)->GetResidueName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(1) << strand_residues.at(0)->GetResidueChainId();
                if(strand_residues.at(0)->GetResidueSequenceNumber() != gmml::iNotSet)
                    stream << std::right << std::setw(4) << strand_residues.at(0)->GetResidueSequenceNumber();
                else
                    stream << std::right << std::setw(4) << " ";
                stream << std::right << std::setw(1) << strand_residues.at(0)->GetResidueInsertionCode()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(3) << strand_residues.at(1)->GetResidueName()
                       << std::left << std::setw(1) << " "
                       << std::right << std::setw(1) << strand_residues.at(1)->GetResidueChainId();
                if(strand_residues.at(1)->GetResidueSequenceNumber() != gmml::iNotSet)
                    stream << std::right << std::setw(4) << strand_residues.at(1)->GetResidueSequenceNumber();
                else
                    stream << std::right << std::setw(4) << " ";
                stream << std::right << std::setw(1) << strand_residues.at(1)->GetResidueInsertionCode();
                if(strand->GetSense() != UnknownStrand)
                    stream << std::right << std::setw(2) << strand->GetSense();
                else
                    stream << std::right << std::setw(2) << " ";
                stream << std::left << std::setw(40) << " "
                       << std::endl;
            }
            serial_number++;
        }

    }
}

void PdbFile::ResolveDisulfideBondCards(std::ostream& stream)
{
    PdbFileSpace::PdbDisulfideBondSection::DisulfideResidueBondMap disulfide_bond_map = disulfide_bonds_->GetDisulfideResidueBonds();
    for(PdbFileSpace::PdbDisulfideBondSection::DisulfideResidueBondMap::iterator it = disulfide_bond_map.begin(); it != disulfide_bond_map.end(); it++)
    {
        PdbDisulfideResidueBond* disulfide_bonds = (*it).second;
        PdbDisulfideResidueBond::DisulfideResidueVector disulfide_bonds_residues = disulfide_bonds->GetResidues();
        stream << std::left << std::setw(6) << disulfide_bonds_->GetRecordName()
               << std::left << std::setw(1) << " ";
        if(disulfide_bonds->GetSerialNumber() != gmml::iNotSet)
            stream << std::right << std::setw(3) << disulfide_bonds->GetSerialNumber();
        else
            stream << std::right << std::setw(3) << " ";
        stream << std::left << std::setw(1) << " "
               << std::right << std::setw(3) << disulfide_bonds_residues.at(0)->GetResidueName()
               << std::left << std::setw(1) << " "
               << std::right << std::setw(1) << disulfide_bonds_residues.at(0)->GetResidueChainId()
               << std::left << std::setw(1) << " ";
        if(disulfide_bonds_residues.at(0)->GetResidueSequenceNumber() != gmml::iNotSet)
            stream << std::right << std::setw(4) << disulfide_bonds_residues.at(0)->GetResidueSequenceNumber();
        else
            stream << std::right << std::setw(4) << " ";
        stream << std::right << std::setw(1) << disulfide_bonds_residues.at(0)->GetResidueInsertionCode()
               << std::left << std::setw(3) << " "
               << std::right << std::setw(3) << disulfide_bonds_residues.at(1)->GetResidueName()
               << std::left << std::setw(1) << " "
               << std::right << std::setw(1) << disulfide_bonds_residues.at(1)->GetResidueChainId()
               << std::left << std::setw(1) << " ";
        if(disulfide_bonds_residues.at(1)->GetResidueSequenceNumber() != gmml::iNotSet)
            stream << std::right << std::setw(4) << disulfide_bonds_residues.at(1)->GetResidueSequenceNumber();
        else
            stream << std::right << std::setw(4) << " ";
        stream << std::right << std::setw(1) << disulfide_bonds_residues.at(1)->GetResidueInsertionCode()
               << std::left << std::setw(23) << " ";
        if(disulfide_bonds_residues.at(0)->GetSymmetryOperator() != gmml::iNotSet)
            stream << std::right << std::setw(6) << disulfide_bonds_residues.at(0)->GetSymmetryOperator();
        else
            stream << std::right << std::setw(6) << " ";
        stream << std::left << std::setw(1) << " ";
        if(disulfide_bonds_residues.at(1)->GetSymmetryOperator() != gmml::iNotSet)
            stream << std::right << std::setw(6) << disulfide_bonds_residues.at(1)->GetSymmetryOperator();
        else
            stream << std::right << std::setw(6) << " ";
        stream << std::left << std::setw(1) << " ";
        if(disulfide_bonds->GetBondLength() != gmml::dNotSet)
            stream << std::right << std::setw(5) << std::fixed << std::setprecision(2) << disulfide_bonds->GetBondLength();
        else
            stream << std::right << std::setw(5) << " ";
        stream << std::left << std::setw(2) << " "
               << std::endl;
    }
}

void PdbFile::ResolveLinkCards(std::ostream& stream)
{
    PdbFileSpace::PdbLinkSection::LinkCardVector links = link_cards_->GetResidueLinkCards();
    for(PdbFileSpace::PdbLinkSection::LinkCardVector::iterator it = links.begin(); it != links.end(); it++)
    {
        PdbLinkCard* link = (*it);
        PdbLinkCard::LinkResidueVector link_residues = link->GetResidues();
        stream << std::left << std::setw(6) << link_cards_->GetRecordName()
               << std::left << std::setw(6) << " "
               << std::left << std::setw(4) << link_residues.at(0)->GetAtomName();
        if(link_residues.at(0)->GetAlternateLocationIndicator() != gmml::BLANK_SPACE)
            stream << std::right << std::setw(1) << link_residues.at(0)->GetAlternateLocationIndicator();
        else
            stream << std::right << std::setw(1) << " ";
        stream << std::right << std::setw(3) << link_residues.at(0)->GetResidueName()
               << std::left << std::setw(1) << " ";
        if(link_residues.at(0)->GetResidueChainId() != gmml::BLANK_SPACE)
            stream << std::right << std::setw(1) << link_residues.at(0)->GetResidueChainId();
        else
            stream << std::right << std::setw(1) << " ";
        if(link_residues.at(0)->GetResidueSequenceNumber() != gmml::iNotSet)
            stream << std::right << std::setw(4) << link_residues.at(0)->GetResidueSequenceNumber();
        else
            stream << std::right << std::setw(4) << " ";
        if(link_residues.at(0)->GetResidueInsertionCode() != gmml::BLANK_SPACE)
            stream << std::right << std::setw(1) << link_residues.at(0)->GetResidueInsertionCode();
        else
            stream << std::right << std::setw(1) << " ";

        stream << std::left << std::setw(15) << " "
               << std::left << std::setw(4) << link_residues.at(1)->GetAtomName();
        if(link_residues.at(1)->GetAlternateLocationIndicator() != gmml::BLANK_SPACE)
            stream << std::right << std::setw(1) << link_residues.at(1)->GetAlternateLocationIndicator();
        else
            stream << std::right << std::setw(1) << " ";
        stream << std::right << std::setw(3) << link_residues.at(1)->GetResidueName()
               << std::left << std::setw(1) << " ";
        if(link_residues.at(1)->GetResidueChainId() != gmml::BLANK_SPACE)
            stream << std::right << std::setw(1) << link_residues.at(1)->GetResidueChainId();
        else
            stream << std::right << std::setw(1) << " ";
        if(link_residues.at(1)->GetResidueSequenceNumber() != gmml::iNotSet)
            stream << std::right << std::setw(4) << link_residues.at(1)->GetResidueSequenceNumber();
        else
            stream << std::right << std::setw(4) << " ";
        if(link_residues.at(1)->GetResidueInsertionCode() != gmml::BLANK_SPACE)
            stream << std::right << std::setw(1) << link_residues.at(1)->GetResidueInsertionCode();
        else
            stream << std::right << std::setw(1) << " ";
        stream << std::left << std::setw(2) << " ";
        if(link_residues.at(0)->GetSymmetryOperator() != gmml::iNotSet)
            stream << std::right << std::setw(6) << link_residues.at(0)->GetSymmetryOperator();
        else
            stream << std::right << std::setw(6) << " ";
        stream << std::left << std::setw(1) << " ";
        if(link_residues.at(1)->GetSymmetryOperator() != gmml::iNotSet)
            stream << std::right << std::setw(6) << link_residues.at(1)->GetSymmetryOperator();
        else
            stream << std::right << std::setw(6) << " ";
        stream << std::left << std::setw(1) << " ";
        if(link->GetLinkLength() != gmml::dNotSet)
            stream << std::right << std::setw(5) << std::fixed << std::setprecision(2) << link->GetLinkLength();
        else
            stream << std::right << std::setw(5) << " ";
        stream << std::left << std::setw(2) << " "
               << std::endl;
    }
}

void PdbFile::ResolveCISPeptideCards(std::ostream& stream)
{
  CISPeptideCardVector cis_peptide_cards = cis_peptide_->GetCISPeptideCards();
  for (CISPeptideCardVector::iterator it = cis_peptide_cards.begin(); it != cis_peptide_cards.end(); it++)
  {
    stream << std::left << std::setw(6) << (*it)->GetRecordName()
           << std::left << std::setw(1) << " "
           << std::right << std::setw(3) << (*it)->GetSerialNumber()
           << std::left << std::setw(1) << " "
           << std::left << std::setw(3) << (*it)->GetPeptide1ResidueName()
           << std::left << std::setw(1) << " "
           << std::left << std::setw(1) << (*it)->GetPeptide1ChainId()
           << std::left << std::setw(1) << " "
           << std::right << std::setw(4) << (*it)->GetPeptide1SequenceNumber()
           << std::left << std::setw(1) << (*it)->GetPeptide1InsertionCode()
           << std::left << std::setw(3) << " "
           << std::left << std::setw(3) << (*it)->GetPeptide2ResidueName()
           << std::left << std::setw(1) << " "
           << std::left << std::setw(1) << (*it)->GetPeptide2ChainId()
           << std::left << std::setw(1) << " "
           << std::right << std::setw(4) << (*it)->GetPeptide2SequenceNumber()
           << std::left << std::setw(1) << (*it)->GetPeptide2InsertionCode()
           << std::left << std::setw(7) << " "
           << std::right << std::setw(3) << (*it)->GetModelNumber()
           << std::left << std::setw(7) << " "
           << std::right << std::setw(6) << (*it)->GetMeasure()
           << std::endl;
  }
}

void PdbFile::ResolveSiteCards(std::ostream& stream)
{
    PdbFileSpace::PdbSiteSection::PdbSiteCardMap site_map = site_cards_->GetResidueSiteCards();
    for(PdbFileSpace::PdbSiteSection::PdbSiteCardMap::iterator it = site_map.begin(); it != site_map.end(); it++)
    {
        PdbSiteCard* site = (*it).second;
        PdbSiteCard::SiteResidueVector site_residues = site->GetResidues();
        const int MAX_RESIDUE_IN_LINE = 4;
        const int RESIDUE_LENGHT_IN_LINE = 11;
        int number_of_residues = site_residues.size();

        int sequence_number = 1;
        if(number_of_residues > MAX_RESIDUE_IN_LINE)
        {
            int number_of_lines = ceil((double)(number_of_residues) / MAX_RESIDUE_IN_LINE);
            while(sequence_number <= number_of_lines)
            {
                if(sequence_number != number_of_lines)
                {
                    std::stringstream ss;
                    for(PdbSiteCard::SiteResidueVector::iterator it1 = site_residues.begin()+(sequence_number - 1)*MAX_RESIDUE_IN_LINE;
                        it1 != site_residues.begin()+(sequence_number)*MAX_RESIDUE_IN_LINE; it1++)
                    {
                        PdbSiteResidue* residue = (*it1);
                        ss << std::left << std::setw(1) << " "
                           << std::right << std::setw(3) << residue->GetResidueName()
                           << std::left << std::setw(1) << " "
                           << std::right << std::setw(1) << residue->GetResidueChainId();
                        if(residue->GetResidueSequenceNumber() != gmml::iNotSet)
                            ss << std::right << std::setw(4) << residue->GetResidueSequenceNumber();
                        else
                            ss << std::right << std::setw(4) << " ";
                        ss << std::right << std::setw(1) << residue->GetResidueInsertionCode();
                    }
                    ss << std::left << std::setw(19) << " ";
                    stream << std::left << std::setw(6) << site_cards_->GetRecordName()
                           << std::left << std::setw(1) << " "
                           << std::right << std::setw(3) << sequence_number
                           << std::left << std::setw(1) << " "
                           << std::right << std::setw(3) << site->GetSiteId()
                           << std::left << std::setw(1) << " ";
                    if(site->GetNumberOfResidues() != gmml::iNotSet)
                        stream << std::right << std::setw(2) << site->GetNumberOfResidues();
                    else
                        stream << std::right << std::setw(2) << " ";
                    stream << std::left << std::setw(63) << ss.str()
                           << std::endl;
                }
                else
                {
                    std::stringstream ss;
                    for(PdbSiteCard::SiteResidueVector::iterator it1 = site_residues.begin()+(sequence_number - 1)*MAX_RESIDUE_IN_LINE;
                        it1 != site_residues.end(); it1++)
                    {
                        PdbSiteResidue* residue = (*it1);
                        ss << std::left << std::setw(1) << " "
                           << std::right << std::setw(3) << residue->GetResidueName()
                           << std::left << std::setw(1) << " "
                           << std::right << std::setw(1) << residue->GetResidueChainId();
                        if(residue->GetResidueSequenceNumber() != gmml::iNotSet)
                            ss << std::right << std::setw(4) << residue->GetResidueSequenceNumber();
                        else
                            ss << std::right << std::setw(4) << " ";
                        ss << std::right << std::setw(1) << residue->GetResidueInsertionCode();
                    }
                    if((sequence_number*MAX_RESIDUE_IN_LINE-number_of_residues)*RESIDUE_LENGHT_IN_LINE != 0)
                        ss << std::left << std::setw((sequence_number*MAX_RESIDUE_IN_LINE-number_of_residues)*RESIDUE_LENGHT_IN_LINE) << " ";
                    ss << std::left << std::setw(19) << " ";
                    stream << std::left << std::setw(6) << site_cards_->GetRecordName()
                           << std::left << std::setw(1) << " "
                           << std::right << std::setw(3) << sequence_number
                           << std::left << std::setw(1) << " "
                           << std::right << std::setw(3) << site->GetSiteId()
                           << std::left << std::setw(1) << " ";
                    if(site->GetNumberOfResidues() != gmml::iNotSet)
                        stream << std::right << std::setw(2) << site->GetNumberOfResidues();
                    else
                        stream << std::right << std::setw(2) << " ";
                    stream << std::left << std::setw(63) << ss.str()
                           << std::endl;
                }
                sequence_number++;
            }

        }
        else
        {
            std::stringstream ss;
            for(PdbSiteCard::SiteResidueVector::iterator it1 = site_residues.begin(); it1 != site_residues.end(); it1++)
            {
                PdbSiteResidue* residue = (*it1);
                ss << std::left << std::setw(1) << " "
                   << std::right << std::setw(3) << residue->GetResidueName()
                   << std::left << std::setw(1) << " "
                   << std::right << std::setw(1) << residue->GetResidueChainId();
                if(residue->GetResidueSequenceNumber() != gmml::iNotSet)
                    ss << std::right << std::setw(4) << residue->GetResidueSequenceNumber();
                else
                    ss << std::right << std::setw(4) << " ";
                ss << std::right << std::setw(1) << residue->GetResidueInsertionCode();
            }
            if((MAX_RESIDUE_IN_LINE-number_of_residues)*RESIDUE_LENGHT_IN_LINE != 0)
                ss << std::left << std::setw((MAX_RESIDUE_IN_LINE-number_of_residues)*RESIDUE_LENGHT_IN_LINE) << " ";
            ss << std::left << std::setw(19) << " ";

            stream << std::left << std::setw(6) << site_cards_->GetRecordName()
                   << std::left << std::setw(1) << " "
                   << std::right << std::setw(3) << sequence_number
                   << std::left << std::setw(1) << " "
                   << std::right << std::setw(3) << site->GetSiteId()
                   << std::left << std::setw(1) << " ";
            if(site->GetNumberOfResidues() != gmml::iNotSet)
                stream << std::right << std::setw(2) << site->GetNumberOfResidues();
            else
                stream << std::right << std::setw(2) << " ";
            stream << std::left << std::setw(63) << ss.str()
                   << std::endl;
            sequence_number++;
        }
    }
}

void PdbFile::ResolveCrystallographyCard(std::ostream& stream)
{
    stream << std::left << std::setw(6) << crystallography_->GetRecordName();
    if(crystallography_->GetA() != gmml::dNotSet)
        stream << std::right << std::setw(9) << std::fixed << std::setprecision(3) << crystallography_->GetA();
    else
        stream << std::right << std::setw(9) << " ";
    if(crystallography_->GetB() != gmml::dNotSet)
        stream << std::right << std::setw(9) << std::fixed << std::setprecision(3) << crystallography_->GetB();
    else
        stream << std::right << std::setw(9) << " ";
    if(crystallography_->GetC() != gmml::dNotSet)
        stream << std::right << std::setw(9) << std::fixed << std::setprecision(3) << crystallography_->GetC();
    else
        stream << std::right << std::setw(9) << " ";
    if(crystallography_->GetAlpha() != gmml::dNotSet)
        stream << std::right << std::setw(7) << std::fixed << std::setprecision(2) << crystallography_->GetAlpha();
    else
        stream << std::right << std::setw(7) << " ";
    if(crystallography_->GetBeta() != gmml::dNotSet)
        stream << std::right << std::setw(7) << std::fixed << std::setprecision(2) << crystallography_->GetBeta();
    else
        stream << std::right << std::setw(7) << " ";
    if(crystallography_->GetGamma() != gmml::dNotSet)
        stream << std::right << std::setw(7) << std::fixed << std::setprecision(2) << crystallography_->GetGamma();
    else
        stream << std::right << std::setw(7) << " ";
    stream << std::left << std::setw(1) << " "
           << std::left << std::setw(11) << crystallography_->GetSpaceGroup();
    if(crystallography_->GetZValue() != gmml::iNotSet)
        stream << std::right << std::setw(4) << crystallography_->GetZValue();
    else
        stream << std::right << std::setw(4) << " ";
    stream << std::left << std::setw(10) << " "
           << std::endl;
}

void PdbFile::ResolveOriginCard(std::ostream& stream)
{
    PdbFileSpace::PdbOriginXnSection::OriginXnCardVector origins = origins_->GetOriginXN();
    for(PdbFileSpace::PdbOriginXnSection::OriginXnCardVector::iterator it = origins.begin(); it != origins.end(); it++)
    {
        PdbOriginXnCard* origin = (*it);
        std::stringstream ss;
        ss << origin->GetRecordName() << origin->GetN();
        stream << std::left << std::setw(6) << ss.str()
               << std::left << std::setw(4) << " ";
        if(origin->GetOrigin().CompareTo(GeometryTopology::Coordinate(gmml::dNotSet, gmml::dNotSet, gmml::dNotSet)) == false)
        {
            stream << std::right << std::setw(10) << std::fixed << std::setprecision(6) << origin->GetOrigin().GetX()
                   << std::right << std::setw(10) << std::fixed << std::setprecision(6) << origin->GetOrigin().GetY()
                   << std::right << std::setw(10) << std::fixed << std::setprecision(6) << origin->GetOrigin().GetZ();
        }
        else
        {
            stream << std::right << std::setw(10) << " "
                   << std::right << std::setw(10) << " "
                   << std::right << std::setw(10) << " ";
        }

        stream << std::left << std::setw(5) << " ";
        if(origin->GetT() != gmml::dNotSet)
            stream << std::right << std::setw(10) << std::fixed << std::setprecision(5) << origin->GetT();
        else
            stream << std::right << std::setw(10) << " ";
        stream << std::left << std::setw(25) << " "
               << std::endl;
    }
}

void PdbFile::ResolveScaleCard(std::ostream& stream)
{
    PdbFileSpace::PdbScaleNSection::ScaleNCardVector scales = scales_->GetScaleNCard();
    for(PdbFileSpace::PdbScaleNSection::ScaleNCardVector::iterator it = scales.begin(); it != scales.end(); it++)
    {
        PdbScaleNCard* scale = (*it);
        std::stringstream ss;
        ss << scale->GetRecordName() << scale->GetN();
        stream << std::left << std::setw(6) << ss.str()
               << std::left << std::setw(4) << " ";
        if(scale->GetScaleVector().CompareTo(GeometryTopology::Coordinate(gmml::dNotSet, gmml::dNotSet, gmml::dNotSet)) == false)
        {
            stream << std::right << std::setw(10) << std::fixed << std::setprecision(6) << scale->GetScaleVector().GetX()
                   << std::right << std::setw(10) << std::fixed << std::setprecision(6) << scale->GetScaleVector().GetY()
                   << std::right << std::setw(10) << std::fixed << std::setprecision(6) << scale->GetScaleVector().GetZ();
        }
        else
        {
            stream << std::right << std::setw(10) << " "
                   << std::right << std::setw(10) << " "
                   << std::right << std::setw(10) << " ";
        }
        stream << std::left << std::setw(5) << " ";
        if(scale->GetU() != gmml::dNotSet)
            stream << std::right << std::setw(10) << std::fixed << std::setprecision(5) << scale->GetU();
        else
            stream << std::right << std::setw(10) << " ";
        stream << std::left << std::setw(25) << " "
               << std::endl;
    }
}

void PdbFile::ResolveMatrixCards(std::ostream& stream)
{
    PdbFileSpace::PdbMatrixNSection::MatrixNVectorVector matrices = matrices_->GetMatrixN();
    int number_of_matrix_entries = matrices.at(0).size();
    for(int i = 0; i < number_of_matrix_entries; i++)
    {
        for(unsigned int j = 0; j < 3; j++)
        {
            PdbFileSpace::PdbMatrixNSection::MatrixNVector matrix_vector = matrices.at(j);
            PdbMatrixNCard* matrix = matrix_vector.at(i);
            std::stringstream ss;
            ss << matrix->GetRecordName() << matrix->GetN();
            stream << std::left << std::setw(6) << ss.str()
                   << std::left << std::setw(1) << " ";
            if(matrix->GetSerialNumber() != gmml::iNotSet)
                stream << std::right << std::setw(3) << matrix->GetSerialNumber();
            else
                stream << std::right << std::setw(3) << " ";
            if(matrix->GetTransformationVector().CompareTo(GeometryTopology::Coordinate(gmml::dNotSet, gmml::dNotSet, gmml::dNotSet)) == false)
            {
                stream << std::right << std::setw(10) << std::fixed << std::setprecision(6) << matrix->GetTransformationVector().GetX()
                       << std::right << std::setw(10) << std::fixed << std::setprecision(6) << matrix->GetTransformationVector().GetY()
                       << std::right << std::setw(10) << std::fixed << std::setprecision(6) << matrix->GetTransformationVector().GetZ();
            }
            else
            {
                stream << std::right << std::setw(10) << " "
                       << std::right << std::setw(10) << " "
                       << std::right << std::setw(10) << " ";
            }
            stream << std::left << std::setw(5) << " ";
            if(matrix->GetV() != gmml::dNotSet)
                stream << std::right << std::setw(10) << std::fixed << std::setprecision(5) << matrix->GetV();
            else
                stream << std::right << std::setw(10) << " ";
            stream << std::left << std::setw(4) << " ";
            if(matrix->GetIGiven() != gmml::iNotSet)
                stream << std::right << std::setw(1) << matrix->GetIGiven();
            else
                stream << std::right << std::setw(1) << " ";
            stream << std::left << std::setw(20) << " "
                   << std::endl;
        }
    }

}

void PdbFile::ResolveModelCards(std::ostream& stream)
{
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    int number_of_models = models.size();
    if(number_of_models == 1)
    {
        for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin(); it != models.end(); it++)
        {
            PdbModelCard* model = (*it).second;
            PdbModelResidueSet* residue_set = model->GetModelResidueSet();
            PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
            for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
            {
                PdbAtomSection* atom_card = (*it1);
                int serial_number = 0;
                std::string residue_name = "";
                char chain_id = ' ';
                int residue_sequence_number = 0;
                char insertion_code = ' ';
                PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
                int atoms_size = ordered_atoms.size();
                for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
                {
                    PdbFileSpace::PdbAtomCard* atom = (*it2);
                    stream << std::left << std::setw(6) << atom_card->GetRecordName();
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
                    stream << std::left << std::setw(10) << " "
                           << std::right << std::setw(2) << atom->GetAtomElementSymbol()
                           << std::left << std::setw(2) << atom->GetAtomCharge()
                           << std::endl;
                    serial_number = atom->GetAtomSerialNumber();
                    residue_name = atom->GetAtomResidueName();
                    chain_id = atom->GetAtomChainId();
                    residue_sequence_number = atom->GetAtomResidueSequenceNumber();
                }
                if(atoms_size != 0)
                {
                    stream << std::left << std::setw(6) << "TER";
                    if(serial_number != gmml::iNotSet)
                        stream << std::right << std::setw(5) << (serial_number+1);
                    else
                        stream << std::right << std::setw(5) << " ";
                    stream << std::left << std::setw(6) << " "
                           << std::right << std::setw(3) << residue_name
                           << std::left << std::setw(1) << " ";
                    if(chain_id == gmml::BLANK_SPACE)
                        stream << std::left << std::setw(1) << " ";
                    else
                        stream << std::left << std::setw(1) << chain_id;
                    if(residue_sequence_number != gmml::iNotSet)
                        stream << std::right << std::setw(4) << residue_sequence_number;
                    else
                        stream << std::right << std::setw(4) << " ";
                    stream << std::left << std::setw(1) << insertion_code
                           << std::left << std::setw(53) << " "
                           << std::endl;
                }
            }
            PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
            for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
            {
                PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
                PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
                for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
                {
                    PdbFileSpace::PdbAtomCard* heterogen_atom = (*it2);
                    stream << std::left << std::setw(6) << heterogen_atom_card->GetRecordName();
                    if(heterogen_atom->GetAtomSerialNumber() != gmml::iNotSet)
                        stream << std::right << std::setw(5) << heterogen_atom->GetAtomSerialNumber();
                    else
                        stream << std::right << std::setw(5) << " ";
                    stream << std::left << std::setw(1) << " "
                           << std::left << std::setw(4) << heterogen_atom->GetAtomName();
                    //stream << std::left << std::setw(2) << " "
                           //<< std::left << std::setw(3) << heterogen_atom->GetAtomName(); //Temporarily commenting out the correct line to write out an incorrect file, which ADT needs.
                    if(heterogen_atom->GetAtomAlternateLocation() == gmml::BLANK_SPACE)
                        stream << std::left << std::setw(1) << ' ';
                    else
                        stream << std::left << std::setw(1) << heterogen_atom->GetAtomAlternateLocation();
                    stream << std::right << std::setw(3) << heterogen_atom->GetAtomResidueName()
                           << std::left << std::setw(1) << " ";
                    if(heterogen_atom->GetAtomChainId() == gmml::BLANK_SPACE)
                        stream << std::left << std::setw(1) << ' ';
                    else
                        stream << std::left << std::setw(1) << heterogen_atom->GetAtomChainId();
                    if(heterogen_atom->GetAtomResidueSequenceNumber() != gmml::iNotSet)
                        stream << std::right << std::setw(4) << heterogen_atom->GetAtomResidueSequenceNumber();
                    else
                        stream << std::right << std::setw(4) << " ";
                    if(heterogen_atom->GetAtomInsertionCode() == gmml::BLANK_SPACE)
                        stream << std::left << std::setw(1) << ' ';
                    else
                        stream << std::left << std::setw(1) << heterogen_atom->GetAtomInsertionCode();
                    stream << std::left << std::setw(3) << " ";
                    if(heterogen_atom->GetAtomOrthogonalCoordinate().CompareTo(GeometryTopology::Coordinate(gmml::dNotSet, gmml::dNotSet, gmml::dNotSet)) == false)
                    {
                        stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetX()
                               << std::right << std::setw(8) << std::fixed << std::setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetY()
                               << std::right << std::setw(8) << std::fixed << std::setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetZ();
                    }
                    else
                    {
                        stream << std::right << std::setw(8) << " "
                               << std::right << std::setw(8) << " "
                               << std::right << std::setw(8) << " ";
                    }
                    if(heterogen_atom->GetAtomOccupancy() != gmml::dNotSet)
                        stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << heterogen_atom->GetAtomOccupancy();
                    else
                        stream << std::right << std::setw(6) << " ";
                    if(heterogen_atom->GetAtomTempretureFactor() != gmml::dNotSet)
                        stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << heterogen_atom->GetAtomTempretureFactor();
                    else
                        stream << std::right << std::setw(6) << " ";
                    stream << std::left << std::setw(10) << " "
                           << std::right << std::setw(2) << heterogen_atom->GetAtomElementSymbol()
                           << std::left << std::setw(2) << heterogen_atom->GetAtomCharge()
                           << std::endl;
                }
            }
        }
    }
    else
    {
        for(PdbFileSpace::PdbModelSection::PdbModelCardMap::iterator it = models.begin(); it != models.end(); it++)
        {
            PdbModelCard* model = (*it).second;
            stream << std::left << std::setw(6) << models_->GetRecordName()
                   << std::left << std::setw(4) << " ";
            if(model->GetModelSerialNumber() != gmml::iNotSet)
                stream << std::right << std::setw(4) << model->GetModelSerialNumber();
            else
                stream << std::right << std::setw(4) << " ";
            stream << std::left << std::setw(66) << " "
                   << std::endl;
            PdbModelResidueSet* residue_set = model->GetModelResidueSet();
            PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
            for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
            {
                PdbAtomSection* atom_card = (*it1);
                int serial_number = 0;
                std::string residue_name = "";
                char chain_id = ' ';
                int residue_sequence_number = 0;
                char insertion_code = ' ';
                PdbAtomSection::PdbAtomCardOrderVector ordered_atoms = atom_card->GetOrderedAtomCards();
                int atoms_size = ordered_atoms.size();
                for(PdbAtomSection::PdbAtomCardOrderVector::iterator it2 = ordered_atoms.begin(); it2 != ordered_atoms.end(); it2++)
                {
                    PdbFileSpace::PdbAtomCard* atom = (*it2);
                    stream << std::left << std::setw(6) << atom_card->GetRecordName();
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
                        stream << std::left << std::setw(1) << ' ';
                    else
                        stream << std::left << std::setw(1) << atom->GetAtomInsertionCode();
                    stream << std::left << std::setw(3) << " ";
                    if(atom->GetAtomOrthogonalCoordinate().CompareTo(GeometryTopology::Coordinate(gmml::dNotSet, gmml::dNotSet, gmml::dNotSet)) == false)
                        stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetX()
                               << std::right << std::setw(8) << std::fixed << std::setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetY()
                               << std::right << std::setw(8) << std::fixed << std::setprecision(3) << atom->GetAtomOrthogonalCoordinate().GetZ();
                    else
                        stream << std::right << std::setw(8) << " "
                               << std::right << std::setw(8) << " "
                               << std::right << std::setw(8) << " ";
                    if(atom->GetAtomOccupancy() != gmml::dNotSet)
                        stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << atom->GetAtomOccupancy();
                    else
                        stream << std::right << std::setw(6) << " ";
                    if(atom->GetAtomTempretureFactor() != gmml::dNotSet)
                        stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << atom->GetAtomTempretureFactor();
                    else
                        stream << std::right << std::setw(6) << " ";
                    stream << std::left << std::setw(10) << " "
                           << std::right << std::setw(2) << atom->GetAtomElementSymbol()
                           << std::left << std::setw(2) << atom->GetAtomCharge()
                           << std::endl;
                    serial_number = atom->GetAtomSerialNumber();
                    residue_name = atom->GetAtomResidueName();
                    chain_id = atom->GetAtomChainId();
                    residue_sequence_number = atom->GetAtomResidueSequenceNumber();
                }
                if(atoms_size != 0)
                {
                    stream << std::left << std::setw(6) << "TER";
                    if(serial_number != gmml::iNotSet)
                        stream << std::right << std::setw(5) << (serial_number+1);
                    else
                        stream << std::right << std::setw(5) << " ";
                    stream << std::left << std::setw(6) << " "
                           << std::right << std::setw(3) << residue_name
                           << std::left << std::setw(1) << " ";
                    if(chain_id == gmml::BLANK_SPACE)
                        stream << std::left << std::setw(1) << " ";
                    else
                        stream << std::left << std::setw(1) << chain_id;
                    if(residue_sequence_number != gmml::iNotSet)
                        stream << std::right << std::setw(4) << residue_sequence_number;
                    else
                        stream << std::right << std::setw(4) << " ";
                    stream << std::left << std::setw(1) << insertion_code
                           << std::left << std::setw(53) << " "
                           << std::endl;
                }
            }
            PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
            for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
            {
                PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
                PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector ordered_heterogen_atoms = heterogen_atom_card->GetOrderedHeterogenAtomCards();
                for(PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector::iterator it2 = ordered_heterogen_atoms.begin(); it2 != ordered_heterogen_atoms.end(); it2++)
                {
                    PdbFileSpace::PdbAtomCard* heterogen_atom = (*it2);
                    stream << std::left << std::setw(6) << heterogen_atom_card->GetRecordName();
                    if(heterogen_atom->GetAtomSerialNumber() != gmml::iNotSet)
                        stream << std::right << std::setw(5) << heterogen_atom->GetAtomSerialNumber();
                    else
                        stream << std::right << std::setw(5) << " ";
                    stream << std::left << std::setw(1) << " "
                           << std::left << std::setw(4) << heterogen_atom->GetAtomName();
                    if(heterogen_atom->GetAtomAlternateLocation() == gmml::BLANK_SPACE)
                        stream << std::left << std::setw(1) << ' ';
                    else
                        stream << std::left << std::setw(1) << heterogen_atom->GetAtomAlternateLocation();
                    stream << std::right << std::setw(3) << heterogen_atom->GetAtomResidueName()
                           << std::left << std::setw(1) << " ";
                    if(heterogen_atom->GetAtomChainId() == gmml::BLANK_SPACE)
                        stream << std::left << std::setw(1) << ' ';
                    else
                        stream << std::left << std::setw(1) << heterogen_atom->GetAtomChainId();
                    if(heterogen_atom->GetAtomResidueSequenceNumber() != gmml::iNotSet)
                        stream << std::right << std::setw(4) << heterogen_atom->GetAtomResidueSequenceNumber();
                    else
                        stream << std::right << std::setw(4) << " ";
                    if(heterogen_atom->GetAtomInsertionCode() == gmml::BLANK_SPACE)
                        stream << std::left << std::setw(1) << ' ';
                    else
                        stream << std::left << std::setw(1) << heterogen_atom->GetAtomInsertionCode();
                    stream << std::left << std::setw(3) << " ";
                    if(heterogen_atom->GetAtomOrthogonalCoordinate().CompareTo(GeometryTopology::Coordinate(gmml::dNotSet, gmml::dNotSet, gmml::dNotSet)) == false)
                        stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetX()
                               << std::right << std::setw(8) << std::fixed << std::setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetY()
                               << std::right << std::setw(8) << std::fixed << std::setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetZ();
                    else
                        stream << std::right << std::setw(8) << " "
                               << std::right << std::setw(8) << " "
                               << std::right << std::setw(8) << " ";
                    if(heterogen_atom->GetAtomOccupancy() != gmml::dNotSet)
                        stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << heterogen_atom->GetAtomOccupancy();
                    else
                        stream << std::right << std::setw(6) << " ";
                    if(heterogen_atom->GetAtomTempretureFactor() != gmml::dNotSet)
                        stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << heterogen_atom->GetAtomTempretureFactor();
                    else
                        stream << std::right << std::setw(6) << " ";
                    stream << std::left << std::setw(10) << " "
                           << std::right << std::setw(2) << heterogen_atom->GetAtomElementSymbol()
                           << std::left << std::setw(2) << heterogen_atom->GetAtomCharge()
                           << std::endl;
                }
            }
            stream << std::left << std::setw(6) << "ENDMDL"
                   << std::left << std::setw(74) << " "
                   << std::endl;
        }
    }
}

void PdbFile::ResolveModelCardWithTheGivenModelNumber(std::ostream& stream, int model_number)
{
    PdbFileSpace::PdbModelSection::PdbModelCardMap models = models_->GetModels();
    if(models.size() != 0)
    {
        PdbModelCard* model = models[model_number];
        PdbModelResidueSet* residue_set = model->GetModelResidueSet();
        PdbModelResidueSet::AtomCardVector atom_cards = residue_set->GetAtomCards();
        for(PdbModelResidueSet::AtomCardVector::iterator it1 = atom_cards.begin(); it1 != atom_cards.end(); it1++)
        {
            PdbAtomSection* atom_card = (*it1);
            int serial_number = 0;
            std::string residue_name = "";
            char chain_id = ' ';
            int residue_sequence_number = 0;
            char insertion_code = ' ';
            PdbAtomSection::PdbAtomMap atoms = atom_card->GetAtomCards();
            int atoms_size = atoms.size();
            for(PdbAtomSection::PdbAtomMap::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* atom = (*it2).second;
                stream << std::left << std::setw(6) << atom_card->GetRecordName();
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
                    stream << std::left << std::setw(1) << ' ';
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
                stream << std::left << std::setw(10) << " "
                       << std::right << std::setw(2) << atom->GetAtomElementSymbol()
                       << std::left << std::setw(2) << atom->GetAtomCharge()
                       << std::endl;
                serial_number = atom->GetAtomSerialNumber();
                residue_name = atom->GetAtomResidueName();
                chain_id = atom->GetAtomChainId();
                residue_sequence_number = atom->GetAtomResidueSequenceNumber();
                std::vector<PdbFileSpace::PdbAtomCard*> alternate_atom_cards = atom->GetAlternateAtomCards();
                if(alternate_atom_cards.size() != 0)
                {
                  for(std::vector<PdbFileSpace::PdbAtomCard*>::iterator thisAltAtom = alternate_atom_cards.begin(); thisAltAtom != alternate_atom_cards.end(); thisAltAtom++)
                  {
                    atom = *thisAltAtom;
                    stream << std::left << std::setw(6) << atom_card->GetRecordName();
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
                        stream << std::left << std::setw(1) << ' ';
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
                    stream << std::left << std::setw(10) << " "
                           << std::right << std::setw(2) << atom->GetAtomElementSymbol()
                           << std::left << std::setw(2) << atom->GetAtomCharge()
                           << std::endl;
                    serial_number = atom->GetAtomSerialNumber();
                    residue_name = atom->GetAtomResidueName();
                    chain_id = atom->GetAtomChainId();
                  }
                }
            }
            if(atoms_size != 0)
            {
                stream << std::left << std::setw(6) << "TER";
                if(serial_number != gmml::iNotSet)
                    stream << std::right << std::setw(5) << (serial_number+1);
                else
                    stream << std::right << std::setw(5) << " ";
                stream << std::left << std::setw(6) << " "
                       << std::right << std::setw(3) << residue_name
                       << std::left << std::setw(1) << " ";
                if(chain_id == gmml::BLANK_SPACE)
                    stream << std::left << std::setw(1) << " ";
                else
                    stream << std::left << std::setw(1) << chain_id;
                if(residue_sequence_number != gmml::iNotSet)
                    stream << std::right << std::setw(4) << residue_sequence_number;
                else
                    stream << std::right << std::setw(4) << " ";
                stream << std::left << std::setw(1) << insertion_code
                       << std::left << std::setw(53) << " "
                       << std::endl;
            }
        }
        PdbModelResidueSet::HeterogenAtomCardVector heterogen_atom_cards = residue_set->GetHeterogenAtomCards();
        for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it1 = heterogen_atom_cards.begin(); it1 != heterogen_atom_cards.end(); it1++)
        {
            PdbHeterogenAtomSection* heterogen_atom_card = (*it1);
            PdbHeterogenAtomSection::PdbHeterogenAtomCardMap heterogen_atoms = heterogen_atom_card->GetHeterogenAtomCards();
            for(PdbHeterogenAtomSection::PdbHeterogenAtomCardMap::iterator it2 = heterogen_atoms.begin(); it2 != heterogen_atoms.end(); it2++)
            {
                PdbFileSpace::PdbAtomCard* heterogen_atom = (*it2).second;
                stream << std::left << std::setw(6) << heterogen_atom_card->GetRecordName();
                if(heterogen_atom->GetAtomSerialNumber() != gmml::iNotSet)
                    stream << std::right << std::setw(5) << heterogen_atom->GetAtomSerialNumber();
                else
                    stream << std::right << std::setw(5) << " ";
                stream << std::left << std::setw(1) << " "
                       << std::left << std::setw(4) << heterogen_atom->GetAtomName();
                if(heterogen_atom->GetAtomAlternateLocation() == gmml::BLANK_SPACE)
                    stream << std::left << std::setw(1) << ' ';
                else
                    stream << std::left << std::setw(1) << heterogen_atom->GetAtomAlternateLocation();
                stream << std::right << std::setw(3) << heterogen_atom->GetAtomResidueName()
                       << std::left << std::setw(1) << " ";
                if(heterogen_atom->GetAtomChainId() == gmml::BLANK_SPACE)
                    stream << std::left << std::setw(1) << ' ';
                else
                    stream << std::left << std::setw(1) << heterogen_atom->GetAtomChainId();
                if(heterogen_atom->GetAtomResidueSequenceNumber() != gmml::iNotSet)
                    stream << std::right << std::setw(4) << heterogen_atom->GetAtomResidueSequenceNumber();
                else
                    stream << std::right << std::setw(4) << " ";
                if(heterogen_atom->GetAtomInsertionCode() == gmml::BLANK_SPACE)
                    stream << std::left << std::setw(1) << ' ';
                else
                    stream << std::left << std::setw(1) << heterogen_atom->GetAtomInsertionCode();
                stream << std::left << std::setw(3) << " ";
                if(heterogen_atom->GetAtomOrthogonalCoordinate().CompareTo(GeometryTopology::Coordinate(gmml::dNotSet, gmml::dNotSet, gmml::dNotSet)) == false)
                {
                    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetX()
                           << std::right << std::setw(8) << std::fixed << std::setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetY()
                           << std::right << std::setw(8) << std::fixed << std::setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetZ();
                }
                else
                {
                    stream << std::right << std::setw(8) << " "
                           << std::right << std::setw(8) << " "
                           << std::right << std::setw(8) << " ";
                }
                if(heterogen_atom->GetAtomOccupancy() != gmml::dNotSet)
                    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << heterogen_atom->GetAtomOccupancy();
                else
                    stream << std::right << std::setw(6) << " ";
                if(heterogen_atom->GetAtomTempretureFactor() != gmml::dNotSet)
                    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << heterogen_atom->GetAtomTempretureFactor();
                else
                    stream << std::right << std::setw(6) << " ";
                stream << std::left << std::setw(10) << " "
                       << std::right << std::setw(2) << heterogen_atom->GetAtomElementSymbol()
                       << std::left << std::setw(2) << heterogen_atom->GetAtomCharge()
                       << std::endl;

                std::vector<PdbFileSpace::PdbAtomCard*> alternate_atom_cards = heterogen_atom->GetAlternateAtomCards();
                if(alternate_atom_cards.size() != 0)
                {
                  for(std::vector<PdbFileSpace::PdbAtomCard*>::iterator thisAltAtom = alternate_atom_cards.begin(); thisAltAtom != alternate_atom_cards.end(); thisAltAtom++)
                  {
                    heterogen_atom = (*thisAltAtom);
                    stream << std::left << std::setw(6) << heterogen_atom_card->GetRecordName();
                    if(heterogen_atom->GetAtomSerialNumber() != gmml::iNotSet)
                        stream << std::right << std::setw(5) << heterogen_atom->GetAtomSerialNumber();
                    else
                        stream << std::right << std::setw(5) << " ";
                    stream << std::left << std::setw(1) << " "
                           << std::left << std::setw(4) << heterogen_atom->GetAtomName();
                    if(heterogen_atom->GetAtomAlternateLocation() == gmml::BLANK_SPACE)
                        stream << std::left << std::setw(1) << ' ';
                    else
                        stream << std::left << std::setw(1) << heterogen_atom->GetAtomAlternateLocation();
                    stream << std::right << std::setw(3) << heterogen_atom->GetAtomResidueName()
                           << std::left << std::setw(1) << " ";
                    if(heterogen_atom->GetAtomChainId() == gmml::BLANK_SPACE)
                        stream << std::left << std::setw(1) << ' ';
                    else
                        stream << std::left << std::setw(1) << heterogen_atom->GetAtomChainId();
                    if(heterogen_atom->GetAtomResidueSequenceNumber() != gmml::iNotSet)
                        stream << std::right << std::setw(4) << heterogen_atom->GetAtomResidueSequenceNumber();
                    else
                        stream << std::right << std::setw(4) << " ";
                    if(heterogen_atom->GetAtomInsertionCode() == gmml::BLANK_SPACE)
                        stream << std::left << std::setw(1) << ' ';
                    else
                        stream << std::left << std::setw(1) << heterogen_atom->GetAtomInsertionCode();
                    stream << std::left << std::setw(3) << " ";
                    if(heterogen_atom->GetAtomOrthogonalCoordinate().CompareTo(GeometryTopology::Coordinate(gmml::dNotSet, gmml::dNotSet, gmml::dNotSet)) == false)
                    {
                        stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetX()
                               << std::right << std::setw(8) << std::fixed << std::setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetY()
                               << std::right << std::setw(8) << std::fixed << std::setprecision(3) << heterogen_atom->GetAtomOrthogonalCoordinate().GetZ();
                    }
                    else
                    {
                        stream << std::right << std::setw(8) << " "
                               << std::right << std::setw(8) << " "
                               << std::right << std::setw(8) << " ";
                    }
                    if(heterogen_atom->GetAtomOccupancy() != gmml::dNotSet)
                        stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << heterogen_atom->GetAtomOccupancy();
                    else
                        stream << std::right << std::setw(6) << " ";
                    if(heterogen_atom->GetAtomTempretureFactor() != gmml::dNotSet)
                        stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << heterogen_atom->GetAtomTempretureFactor();
                    else
                        stream << std::right << std::setw(6) << " ";
                    stream << std::left << std::setw(10) << " "
                           << std::right << std::setw(2) << heterogen_atom->GetAtomElementSymbol()
                           << std::left << std::setw(2) << heterogen_atom->GetAtomCharge()
                           << std::endl;
                  }
                }
            }
        }
    }
}

void PdbFile::ResolveConnectivityCards(std::ostream& stream)
{
    PdbFileSpace::PdbConnectSection::BondedAtomsSerialNumbersMap bonded_atoms = connectivities_->GetBondedAtomsSerialNumbers();
    for(PdbFileSpace::PdbConnectSection::BondedAtomsSerialNumbersMap::iterator it = bonded_atoms.begin(); it != bonded_atoms.end(); it++)
    {
        std::vector<int> bonded_atoms_serial_number = (*it).second;
        int source_atom_serial_number = (*it).first;
        int number_of_bonded_atoms = bonded_atoms_serial_number.size();
        const int MAX_SERIAL_NUMBER_IN_LINE = 4;
        const int SERIAL_NUMBER_LENGTH = 5;
        if(number_of_bonded_atoms <= MAX_SERIAL_NUMBER_IN_LINE)
        {
            stream << std::left << std::setw(6) << connectivities_->GetRecordName();
            if(source_atom_serial_number != gmml::iNotSet)
                stream << std::right << std::setw(5) << source_atom_serial_number;
            else
                stream << std::right << std::setw(5) << " ";
            for(std::vector<int>::iterator it1 = bonded_atoms_serial_number.begin(); it1 != bonded_atoms_serial_number.end(); it1++)
            {
                int serial_number = (*it1);
                if(serial_number != gmml::iNotSet)
                    stream << std::right << std::setw(5) << serial_number;
                else
                    stream << std::right << std::setw(5) << " ";
            }
            if((MAX_SERIAL_NUMBER_IN_LINE-number_of_bonded_atoms)*SERIAL_NUMBER_LENGTH != 0)
                stream << std::left << std::setw((MAX_SERIAL_NUMBER_IN_LINE-number_of_bonded_atoms)*SERIAL_NUMBER_LENGTH) << " "
                       << std::left << std::setw(49) << " "
                       << std::endl;
            else
                stream << std::left << std::setw(49) << " "
                       << std::endl;
        }
        else
        {
            int number_of_lines = ceil((double)(number_of_bonded_atoms) / MAX_SERIAL_NUMBER_IN_LINE);
            for(int i = 1; i <= number_of_lines; i++)
            {
                stream << std::left << std::setw(6) << connectivities_->GetRecordName();
                if(source_atom_serial_number != gmml::iNotSet)
                    stream << std::right << std::setw(5) << source_atom_serial_number;
                else
                    stream << std::right << std::setw(5) << " ";
                if(i != number_of_lines)
                {
                    for(std::vector<int>::iterator it1 = bonded_atoms_serial_number.begin() + (i-1) * MAX_SERIAL_NUMBER_IN_LINE;
                        it1 != bonded_atoms_serial_number.begin() + i * MAX_SERIAL_NUMBER_IN_LINE; it1++)
                    {
                        int serial_number = (*it1);
                        if(serial_number != gmml::iNotSet)
                            stream << std::right << std::setw(5) << serial_number;
                        else
                            stream << std::right << std::setw(5) << " ";
                    }
                    stream << std::left << std::setw(49) << " "
                           << std::endl;
                }
                else
                {
                    for(std::vector<int>::iterator it1 = bonded_atoms_serial_number.begin() + (i-1) * MAX_SERIAL_NUMBER_IN_LINE;
                        it1 != bonded_atoms_serial_number.end(); it1++)
                    {
                        int serial_number = (*it1);
                        if(serial_number != gmml::iNotSet)
                            stream << std::right << std::setw(5) << serial_number;
                        else
                            stream << std::right << std::setw(5) << " ";
                    }
                    if((MAX_SERIAL_NUMBER_IN_LINE-(number_of_bonded_atoms-(i-1)*MAX_SERIAL_NUMBER_IN_LINE))*SERIAL_NUMBER_LENGTH != 0)
                        stream << std::left << std::setw((MAX_SERIAL_NUMBER_IN_LINE-number_of_bonded_atoms)*SERIAL_NUMBER_LENGTH) << " "
                               << std::left << std::setw(49) << " "
                               << std::endl;
                    else
                        stream << std::left << std::setw(49) << " "
                               << std::endl;
                }
            }
        }
    }
}

void PdbFile::ResolveMasterCards(std::ostream& stream)
{
  stream << std::left << std::setw(6) << master_->GetRecordName()
         << std::right << std::setw(4) << " "
         << std::right << std::setw(5) << master_->GetNumRemark()
         << std::right << std::setw(5) << 0
         << std::right << std::setw(5) << master_->GetNumHet()
         << std::right << std::setw(5) << master_->GetNumHelix()
         << std::right << std::setw(5) << master_->GetNumSheet()
         << std::right << std::setw(5) << 0
         << std::right << std::setw(5) << master_->GetNumSite()
         << std::right << std::setw(5) << master_->GetNumXForm()
         << std::right << std::setw(5) << master_->GetNumCoord()
         << std::right << std::setw(5) << master_->GetNumTer()
         << std::right << std::setw(5) << master_->GetNumConnect()
         << std::right << std::setw(5) << master_->GetNumSeq()
         << std::endl;
}

void PdbFile::ResolveEndCard(std::ostream& stream)
{
    stream << std::left << std::setw(6) << "END" << std::left << std::setw(74) << " " << std::endl;
}

void PdbFile::PrintOntology(std::stringstream& ont_stream)
{
  //Match formatting of Ontology
  std::stringstream uri;
  uri << Ontology::ONT_PREFIX;
  if(header_ != NULL)
  {
    uri << header_->GetIdentifierCode();
  }
  else
  {
    std::string file = gmml::Split(path_.substr(path_.find_last_of('/') + 1), ".").at(0);
    // std::transform(file.begin(), file.end(),file.begin(), ::toupper);
    uri << file;
  }
  std::string uriStr = uri.str();
  std::transform(uriStr.begin(), uriStr.end(), uriStr.begin(), ::tolower);

  gmml::AddLiteral( uriStr, Ontology::TYPE, Ontology::PDB, ont_stream );

  //Return PDB_ID
  if(header_ != NULL)
    gmml::AddLiteral( uriStr, Ontology::id, this->header_->GetIdentifierCode(), ont_stream );

  //Return Protein Acession Number
  //TODO add check to make sure that it is the Uniprot database reference
  if(database_reference_ != NULL)
    gmml::AddLiteral( uriStr, Ontology::hasProteinID, this->database_reference_->GetUniprotIDs(), ont_stream );

  //Return Title
  if(title_ != NULL)
    gmml::AddLiteral( uriStr, Ontology::hasTitle, this->title_->GetTitle(), ont_stream );

  //Return Authors
  if(author_ != NULL)
    gmml::AddLiteral( uriStr, Ontology::hasAuthors, this->author_->GetAuthor(), ont_stream );


  if(journal_ != NULL)
  {
      //Return JOURNAL
      gmml::AddLiteral( uriStr, Ontology::hasJournal, this->journal_->GetReference(), ont_stream );

      //Return DOI
      gmml::AddLiteral( uriStr, Ontology::hasDOI, this->journal_->GetDOI(), ont_stream );

      //Return PMID
      gmml::AddLiteral( uriStr, Ontology::hasPMID, this->journal_->GetPMID(), ont_stream );
  }



  if(remark_cards_ != NULL)
  {
    //Return Resolution
    gmml::AddDecimal( uriStr, Ontology::hasResolution, this->remark_cards_->GetResolution(), ont_stream );

    //Return B Factor
    gmml::AddDecimal( uriStr, Ontology::hasBFactor, this->remark_cards_->GetBFactor(), ont_stream );
  }
}
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

void PdbFile::Print(std::ostream &out)
{
    if(header_ != NULL)
    {
        out << "******************************* HEADER *******************************" << std::endl;
        header_->Print(out);
    }
    if(title_ != NULL)
    {
        out << "******************************** TITLE *******************************" << std::endl;
        title_->Print(out);
    }
    if(split_ != NULL)
    {
        out << "******************************** SPLIT *******************************" << std::endl;
        split_->Print(out);
    }
    if(caveat_ != NULL)
    {
        out << "******************************** CAVEAT *******************************" << std::endl;
        caveat_->Print(out);
    }
    if(compound_ != NULL)
    {
        out << "****************************** COMPOUND ******************************" << std::endl;
        compound_->Print(out);
    }
    if(source_ != NULL)
    {
        out << "******************************** SOURCE *******************************" << std::endl;
        source_->Print(out);
    }
    if(keywords_ != NULL)
    {
        out << "******************************** KEYWORDS *******************************" << std::endl;
        keywords_->Print(out);
    }
    if(experimental_data_ != NULL)
    {
        out << "*************************** EXPERIMENTAL DATA ****************************" << std::endl;
        experimental_data_->Print(out);
    }
    if(number_of_models_ != NULL)
    {
        out << "************************** NUMBER OF MODELS **************************" << std::endl;
        number_of_models_->Print(out);
    }
    if(model_type_ != NULL)
    {
        out << "***************************** MODEL TYPE *****************************" << std::endl;
        model_type_->Print(out);
    }
    if(author_ != NULL)
    {
        out << "******************************** AUTHOR ******************************" << std::endl;
        author_->Print(out);
    }
    if(revision_data_ != NULL)
    {
        out << "**************************** REVISION DATA ***************************" << std::endl;
        revision_data_->Print(out);
    }
    if(superseded_entries_ != NULL)
    {
        out << "************************** SUPERSEDED ENTRIES *************************" << std::endl;
        superseded_entries_->Print(out);
    }
    if(journal_ != NULL)
    {
        out << "******************************** JOURNAL *******************************" << std::endl;
        journal_->Print(out);
    }
    if(remark_cards_ != NULL)
    {
        out << "****************************** REMARKS *******************************" << std::endl;
        remark_cards_->Print(out);
    }
    if(database_reference_ != NULL)
    {
        out << "************************** DATABASE REFERENCE *************************" << std::endl;
        database_reference_->Print(out);
    }
    if(sequence_advanced_ != NULL)
    {
        out << "*************************** SEQUENCE ADVANCED **************************" << std::endl;
        sequence_advanced_->Print(out);
    }
    if(residues_sequence_ != NULL)
    {
        out << "************************** RESIDUE SEQUENCE **************************" << std::endl;
        residues_sequence_->Print(out);
    }
    if(residue_modification_cards_ != NULL)
    {
        out << "************************ RESIDUE MODIFICATION ************************" << std::endl;
        residue_modification_cards_->Print(out);
    }
    if(heterogen_cards_ != NULL)
    {
        out << "***************************** HETEROGEN ******************************" << std::endl;
        heterogen_cards_->Print(out);
    }
    if(heterogen_name_cards_ != NULL)
    {
        out << "*************************** HETEROGEN NAME ***************************" << std::endl;
        heterogen_name_cards_->Print(out);
    }
    if(heterogen_synonym_cards_ != NULL)
    {
        out << "************************** HETEROGEN SYNONYM *************************" << std::endl;
        heterogen_synonym_cards_->Print(out);
    }
    if(formulas_ != NULL)
    {
        out << "******************************* FORMULA ******************************" << std::endl;
        formulas_->Print(out);
    }
    if(helix_cards_ != NULL)
    {
        out << "******************************** HELIX *******************************" << std::endl;
        helix_cards_->Print(out);
    }
    if(sheet_cards_ != NULL)
    {
        out << "******************************** SHEET *******************************" << std::endl;
        sheet_cards_->Print(out);
    }
    if(disulfide_bonds_ != NULL)
    {
        out << "*************************** DISULFIDE BOND ***************************" << std::endl;
        disulfide_bonds_->Print(out);
    }
    if(link_cards_ != NULL)
    {
        out << "******************************** LINK ********************************" << std::endl;
        link_cards_->Print(out);
    }
    if(cis_peptide_ != NULL)
    {
        out << "***************************** CIS PEPTIDE ****************************" << std::endl;
        cis_peptide_->Print(out);
    }
    if(site_cards_ != NULL)
    {
        out << "******************************** SITE ********************************" << std::endl;
        site_cards_->Print(out);
    }
    if(crystallography_ != NULL)
    {
        out << "************************** CRYSTALLOGRAPHIC **************************" << std::endl;
        crystallography_->Print(out);
    }
    if(origins_ != NULL)
    {
        out << "******************************* ORIGIN *******************************" << std::endl;
        origins_->Print(out);
    }
    if(scales_ != NULL)
    {
        out << "******************************** SCALE *******************************" << std::endl;
        scales_->Print(out);
    }
    if(matrices_ != NULL)
    {
        out << "******************************* MATRIX *******************************" << std::endl;
        matrices_->Print(out);
    }
    if(models_ != NULL)
    {
        out << "******************************* MODEL ********************************" << std::endl;
        models_->Print(out);
    }
    if(connectivities_ != NULL)
    {
        out << "******************************* CONNECT ******************************" << std::endl;
        connectivities_->Print(out);
    }
    if(master_ != NULL)
    {
        out << "******************************** MASTER *******************************" << std::endl;
        master_->Print(out);
    }
}
