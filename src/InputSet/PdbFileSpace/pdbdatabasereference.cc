#include "../../../includes/InputSet/PdbFileSpace/pdbdatabasereference.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace PdbFileSpace;


//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbDatabaseReference::PdbDatabaseReference() {}
PdbDatabaseReference::PdbDatabaseReference(string &line)
{
    std::string temp = line;
    if (!Trim(temp).empty())
    {
      record_name_ = line.substr(0,6);

      id_code_ = line.substr(7,4);

      chain_id_ = line.substr(12,1);

      seq_begin_ = atoi(line.substr(14,4).c_str());

      insert_begin_ = line.substr(18,1);

      seq_end_ = atoi(line.substr(20,4).c_str());

      insert_end_ = line.substr(24,1);

      database_ = line.substr(26,6);

      if (record_name_ == "DBREF ")
      {
        db_accession_ = line.substr(33,8);

        db_id_code_ = line.substr(42,12);

        db_seq_begin_ = atoi(line.substr(55,5).c_str());

        db_ins_beg_ = line.substr(60,1);

        db_seq_end_ = atoi(line.substr(62,5).c_str());

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

             db_seq_begin_ = atoi(dbref2.substr(45,10).c_str());

             db_ins_beg_ = ' ';

             db_seq_end_ = atoi(line.substr(57,10).c_str());

             db_ins_end_ = ' ';
           }
      }
    }
}


//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

string PdbDatabaseReference::GetRecordName()
{
    return record_name_;
}

string PdbDatabaseReference::GetIDCode()
{
    return id_code_;
}

string PdbDatabaseReference::GetChainID()
{
    return chain_id_;
}

int PdbDatabaseReference::GetSeqBegin()
{
    return seq_begin_;
}

string PdbDatabaseReference::GetInsertBegin()
{
    return insert_begin_;
}

int PdbDatabaseReference::GetSeqEnd()
{
    return seq_end_;
}

string PdbDatabaseReference::GetInsertEnd()
{
    return insert_end_;
}

string PdbDatabaseReference::GetDatabase()
{
    return database_;
}

string PdbDatabaseReference::GetDatabaseAccession()
{
    return db_accession_;
}

string PdbDatabaseReference::GetDatabaseIDCode()
{
    return db_id_code_;
}

int PdbDatabaseReference::GetDatabaseSeqBegin()
{
    return db_seq_begin_;
}

string PdbDatabaseReference::GetDatabaseInsBegin()
{
    return db_ins_beg_;
}

int PdbDatabaseReference::GetDatabaseSeqEnd()
{
    return db_seq_end_;
}

string PdbDatabaseReference::GetDatabaseInsEnd()
{
    return db_ins_end_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbDatabaseReference::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

void PdbDatabaseReference::SetIDCode(const string id_code)
{
    id_code_ = id_code;
}

void PdbDatabaseReference::SetChainId(const string chain_id)
{
    chain_id_ = chain_id;
}

void PdbDatabaseReference::SetSeqBegin(const int seq_begin)
{
    seq_begin_ = seq_begin;
}

void PdbDatabaseReference::SetInsertBegin(const string insert_begin)
{
    insert_begin_ = insert_begin;
}

void PdbDatabaseReference::SetSeqEnd(const int seq_end)
{
    seq_end_ = seq_end;
}

void PdbDatabaseReference::SetInsertEnd(const string insert_end)
{
    insert_end_ = insert_end;
}

void PdbDatabaseReference::SetDatabase(const string database)
{
    database_ = database;
}

void PdbDatabaseReference::SetDatabaseAccession(const string db_accession)
{
    db_accession_ = db_accession;
}

void PdbDatabaseReference::SetDatabaseIDCode(const string db_id_code)
{
    db_id_code_ = db_id_code;
}

void PdbDatabaseReference::SetDatabaseSeqBegin(const int db_seq_begin)
{
    db_seq_begin_ = db_seq_begin;
}

void PdbDatabaseReference::SetDatabaseInsBegin(const string db_ins_beg)
{
    db_ins_beg_ = db_ins_beg;
}

void PdbDatabaseReference::SetDatabaseSeqEnd(const int db_seq_end)
{
    db_seq_end_ = db_seq_end;
}

void PdbDatabaseReference::SetDatabaseInsEnd(const string db_ins_end)
{
    db_ins_end_ = db_ins_end;
}


//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void PdbDatabaseReference::Print(ostream &out)
{
    out << "Record Name: " << record_name_;
    out << "ID Code: " << id_code_;
    out << "Chain ID: " << chain_id_;
    out << "Sequence Begin: " << seq_begin_;
    out << "Insert Begin: " << insert_begin_;
    out << "Sequence End: " << seq_end_;
    out << "Insert End: " << insert_end_;
    out << "Database: " << database_;
    out << "DB accession: " << db_accession_;
    out << "DB ID Code: " << db_id_code_;
    out << "DB Sequence Begin: " << db_seq_begin_;
    out << "DB Insert Begin: " << db_ins_beg_;
    out << "DB Sequence End: " << db_seq_end_;
    out << "DB Insert End: " << db_ins_end_;
    out << endl;
}
