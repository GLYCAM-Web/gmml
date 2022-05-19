#include "includes/InputSet/PdbFile/SectionClasses/databaseReferenceRecord.hpp"
#include "includes/utils.hpp"
#include "includes/common.hpp"

using pdb::DatabaseReference;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
DatabaseReference::DatabaseReference(std::string &line)
{
    std::string temp = line;
    if (!gmml::Trim(temp).empty())
    {
        record_name_ = line.substr(0,6);

        id_code_ = line.substr(7,4);

        chain_id_ = line.substr(12,1);

        seq_begin_ = gmml::ConvertString<int>(line.substr(14,4));

        insert_begin_ = line.substr(18,1);

        seq_end_ = gmml::ConvertString<int>(line.substr(20,4));

        insert_end_ = line.substr(24,1);

        database_ = line.substr(26,6);

        if (record_name_ == "DBREF ")
        {
            db_accession_ = line.substr(33,8);

            db_id_code_ = line.substr(42,12);

            db_seq_begin_ = gmml::ConvertString<int>(line.substr(55,5));

            db_ins_beg_ = line.substr(60,1);

            db_seq_end_ = gmml::ConvertString<int>(line.substr(62,5));

            db_ins_end_ = line.substr(67,1);
        }
        else if (record_name_ == "DBREF1")
        {
            std::size_t position = line.find("DBREF2");
            if (position!=std::string::npos)
            {
                std::string dbref2 = line.substr(position, 68);

                db_accession_ = dbref2.substr(18,22);

                db_id_code_ = line.substr(47,15);

                db_seq_begin_ = gmml::ConvertString<int>(dbref2.substr(45,10));

                db_ins_beg_ = ' ';

                db_seq_end_ = gmml::ConvertString<int>(line.substr(57,10));

                db_ins_end_ = ' ';
            }
        }
    }
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

std::string DatabaseReference::GetRecordName() const
{
    return record_name_;
}

std::string DatabaseReference::GetIDCode() const
{
    return id_code_;
}

std::string DatabaseReference::GetChainID() const
{
    return chain_id_;
}

int DatabaseReference::GetSeqBegin() const
{
    return seq_begin_;
}

std::string DatabaseReference::GetInsertBegin() const
{
    return insert_begin_;
}

int DatabaseReference::GetSeqEnd() const
{
    return seq_end_;
}

std::string DatabaseReference::GetInsertEnd() const
{
    return insert_end_;
}

std::string DatabaseReference::GetDatabase() const
{
    return database_;
}

std::string DatabaseReference::GetDatabaseAccession() const
{
    return db_accession_;
}

std::string DatabaseReference::GetDatabaseIDCode() const
{
    return db_id_code_;
}

int DatabaseReference::GetDatabaseSeqBegin() const
{
    return db_seq_begin_;
}

std::string DatabaseReference::GetDatabaseInsBegin() const
{
    return db_ins_beg_;
}

int DatabaseReference::GetDatabaseSeqEnd() const
{
    return db_seq_end_;
}

std::string DatabaseReference::GetDatabaseInsEnd() const
{
    return db_ins_end_;
}

std::string DatabaseReference::GetUniprotID() const
{
    if(this->GetDatabase() == "UNP   ")
    {
        return this->GetDatabaseAccession() + " ";
    }
    return "";
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void DatabaseReference::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void DatabaseReference::SetIDCode(const std::string id_code)
{
    id_code_ = id_code;
}

void DatabaseReference::SetChainId(const std::string chain_id)
{
    chain_id_ = chain_id;
}

void DatabaseReference::SetSeqBegin(const int seq_begin)
{
    seq_begin_ = seq_begin;
}

void DatabaseReference::SetInsertBegin(const std::string insert_begin)
{
    insert_begin_ = insert_begin;
}

void DatabaseReference::SetSeqEnd(const int seq_end)
{
    seq_end_ = seq_end;
}

void DatabaseReference::SetInsertEnd(const std::string insert_end)
{
    insert_end_ = insert_end;
}

void DatabaseReference::SetDatabase(const std::string database)
{
    database_ = database;
}

void DatabaseReference::SetDatabaseAccession(const std::string db_accession)
{
    db_accession_ = db_accession;
}

void DatabaseReference::SetDatabaseIDCode(const std::string db_id_code)
{
    db_id_code_ = db_id_code;
}

void DatabaseReference::SetDatabaseSeqBegin(const int db_seq_begin)
{
    db_seq_begin_ = db_seq_begin;
}

void DatabaseReference::SetDatabaseInsBegin(const std::string db_ins_beg)
{
    db_ins_beg_ = db_ins_beg;
}

void DatabaseReference::SetDatabaseSeqEnd(const int db_seq_end)
{
    db_seq_end_ = db_seq_end;
}

void DatabaseReference::SetDatabaseInsEnd(const std::string db_ins_end)
{
    db_ins_end_ = db_ins_end;
}
//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void DatabaseReference::Print(std::ostream &out) const
{
    out << "Record Name: " << this->GetRecordName();
    out << "ID Code: " << this->GetIDCode();
    out << "Chain ID: " << this->GetChainID();
    out << "Sequence Begin: " << this->GetSeqBegin();
    out << "Insert Begin: " << this->GetInsertBegin();
    out << "Sequence End: " << this->GetSeqEnd();
    out << "Insert End: " << this->GetInsertEnd();
    out << "Database: " << this->GetDatabase();
    out << "DB accession: " << this->GetDatabaseAccession();
    out << "DB ID Code: " << this->GetDatabaseIDCode();
    out << "DB Sequence Begin: " << this->GetDatabaseSeqBegin();
    out << "DB Insert Begin: " << this->GetDatabaseInsBegin();
    out << "DB Sequence End: " << this->GetDatabaseSeqEnd();
    out << "DB Insert End: " << this->GetDatabaseInsEnd();
    out << std::endl;
}

void DatabaseReference::Write(std::ostream& stream) const
{
    if (this->GetRecordName() == "DBREF ")
    {
        stream << std::left << std::setw(6) << this->GetRecordName()
               << std::left << std::setw(1) << " "
               << std::left << std::setw(4) << this->GetIDCode()
               << std::left << std::setw(1) << " "
               << std::left << std::setw(1) << this->GetChainID()
               << std::left << std::setw(1) << " "
               << std::right << std::setw(4) << this->GetSeqBegin()
               << std::right << std::setw(1) << this->GetInsertBegin()
               << std::left << std::setw(1) << " "
               << std::right << std::setw(4) << this->GetSeqEnd()
               << std::right << std::setw(1) << this->GetInsertEnd()
               << std::left << std::setw(1) << " "
               << std::left << std::setw(6) << this->GetDatabase()
               << std::left << std::setw(1) << " "
               << std::left << std::setw(8) << this->GetDatabaseAccession()
               << std::left << std::setw(1) << " "
               << std::left << std::setw(12) << this->GetDatabaseIDCode()
               << std::left << std::setw(1) << " "
               << std::right << std::setw(5) << this->GetDatabaseSeqBegin()
               << std::right << std::setw(1) << this->GetDatabaseInsBegin()
               << std::left << std::setw(1) << " "
               << std::right << std::setw(5) << this->GetDatabaseSeqEnd()
               << std::right << std::setw(1) << this->GetDatabaseInsEnd()
               << std::endl;
    }
    else if (this->GetRecordName() == "DBREF1")
    {
        stream << std::left << std::setw(6) << this->GetRecordName()
               << std::left << std::setw(1) << " "
               << std::left << std::setw(4) << this->GetIDCode()
               << std::left << std::setw(1) << " "
               << std::left << std::setw(1) << this->GetChainID()
               << std::left << std::setw(1) << " "
               << std::right << std::setw(4) << this->GetSeqBegin()
               << std::right << std::setw(1) << this->GetInsertBegin()
               << std::left << std::setw(1) << " "
               << std::right << std::setw(4) << this->GetSeqEnd()
               << std::right << std::setw(1) << this->GetInsertEnd()
               << std::left << std::setw(1) << " "
               << std::left << std::setw(6) << this->GetDatabase()
               << std::left << std::setw(16) << " "
               << std::left << std::setw(15) << this->GetDatabaseIDCode()
               << std::endl
               << std::left << std::setw(6) << "DBREF2"
               << std::left << std::setw(1) << " "
               << std::left << std::setw(4) << this->GetIDCode()
               << std::left << std::setw(1) << " "
               << std::left << std::setw(1) << this->GetChainID()
               << std::left << std::setw(6) << " "
               << std::left << std::setw(22) << this->GetDatabaseAccession()
               << std::left << std::setw(5) << " "
               << std::right << std::setw(10) << this->GetDatabaseSeqBegin()
               << std::left << std::setw(2) << " "
               << std::right << std::setw(10) << this->GetDatabaseSeqEnd()
               << std::endl;
    }
}
