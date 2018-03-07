#include "../../../includes/InputSet/PdbFileSpace/pdbsequenceadvancedcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbSequenceAdvancedCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSequenceAdvancedCard::PdbSequenceAdvancedCard() {}
PdbSequenceAdvancedCard::PdbSequenceAdvancedCard(std::string &line)
{
    if (!gmml::Trim(line).empty())
    {
      record_name_ = line.substr(0, 6);
      id_code_ = line.substr(7, 4);
      residue_name_ = line.substr(12, 3);
      chain_id_ = line.substr(16, 1);
      seq_num_ = gmml::ConvertString<int>(line.substr(18, 4));
      i_code_ = line.substr(22, 1);
      database_ = line.substr(24, 4);
      db_accession_ = line.substr(29, 9);
      db_res_ = line.substr(39, 3);
      db_seq_ = gmml::ConvertString<int>(line.substr(43, 5));
      conflict_ = line.substr(49, 30);
    }
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

std::string PdbSequenceAdvancedCard::GetRecordName()
{
    return record_name_;
}

std::string PdbSequenceAdvancedCard::GetIdentifierCode()
{
    return id_code_;
}

std::string PdbSequenceAdvancedCard::GetResidueName()
{
    return residue_name_;
}
std::string PdbSequenceAdvancedCard::GetChainId()
{
    return chain_id_;
}

int PdbSequenceAdvancedCard::GetSequenceNumber()
{
    return seq_num_;
}
std::string PdbSequenceAdvancedCard::GetInsertionCode()
{
    return i_code_;
}

std::string PdbSequenceAdvancedCard::GetDatabase()
{
    return database_;
}
std::string PdbSequenceAdvancedCard::GetDatabaseAccession()
{
    return db_accession_;
}

std::string PdbSequenceAdvancedCard::GetDatabaseResidue()
{
    return db_res_;
}
int PdbSequenceAdvancedCard::GetDatabaseSequence()
{
    return db_seq_;
}

std::string PdbSequenceAdvancedCard::GetConflict()
{
    return conflict_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbSequenceAdvancedCard::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbSequenceAdvancedCard::SetIdentifierCode(const std::string id_code)
{
     id_code_ = id_code;
}

void PdbSequenceAdvancedCard::SetResidueName(const std::string residue_name)
{
     residue_name_ = residue_name;
}
void PdbSequenceAdvancedCard::SetChainId(const std::string chain_id)
{
     chain_id_ = chain_id;
}

void PdbSequenceAdvancedCard::SetSequenceNumber(int seq_num)
{
     seq_num_ = seq_num;
}
void PdbSequenceAdvancedCard::SetInsertionCode(const std::string i_code)
{
     i_code_ = i_code;
}

void PdbSequenceAdvancedCard::SetDatabase(const std::string database)
{
     database_ = database;
}
void PdbSequenceAdvancedCard::SetDatabaseAccession(const std::string db_accession)
{
     db_accession_ = db_accession;
}

void PdbSequenceAdvancedCard::SetDatabaseResidue(const std::string db_res)
{
     db_res_ = db_res;
}
void PdbSequenceAdvancedCard::SetDatabaseSequence(int db_seq)
{
     db_seq_ = db_seq;
}

void PdbSequenceAdvancedCard::SetConflict(const std::string conflict)
{
     conflict_ = conflict;
}


//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbSequenceAdvancedCard::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_;
    out << "ID Code: " << id_code_;
    out << "Residue Name: " << residue_name_;
    out << "Chain ID: " << chain_id_;
    out << "Sequence Number: " << seq_num_;
    out << "Insertion Code: " << i_code_;
    out << "Database: " << database_;
    out << "Database Accession: " << db_accession_;
    out << "Database Residue Name: " << db_res_;
    out << "Database Sequence Number: " << db_seq_;
    out << "Conflict Comment: " << conflict_;
    out << std::endl;
}
