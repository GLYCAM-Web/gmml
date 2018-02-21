#include "../../../includes/InputSet/PdbFileSpace/pdbdatabasereferencesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbdatabasereference.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbDatabaseReferenceSection::PdbDatabaseReferenceSection() {}
PdbDatabaseReferenceSection::PdbDatabaseReferenceSection(stringstream &stream_block)
{
    string line;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        string record_name = line.substr(0, 6);
        if(record_name == "DBREF1")
        { //This gets this line and the next line, as DBREF1&2 are two line references
          PdbDatabaseReference* database_reference = new PdbDatabaseReference(line);
          AddDatabaseReferences(database_reference);
          getline(stream_block, line);
          temp = line;
          PdbDatabaseReference* database_reference_2 = new PdbDatabaseReference(line);
          AddDatabaseReferences(database_reference_2);
        }
        else if (record_name == "DBREF ")
        { //DBREF is a one line reference
          PdbDatabaseReference* database_reference = new PdbDatabaseReference(line);
          AddDatabaseReferences(database_reference);
        }
        getline(stream_block, line);
        temp = line;
    }
    }

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbDatabaseReferenceSection::DatabaseReferenceVector PdbDatabaseReferenceSection::GetDatabaseReferences(){
    return database_reference_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbDatabaseReferenceSection::SetDatabaseReferences(DatabaseReferenceVector database_reference){
    database_reference_.clear();
    for(DatabaseReferenceVector::iterator it = database_reference.begin(); it != database_reference.end(); it++)
    {
        database_reference.push_back(*it);
    }
}

void PdbDatabaseReferenceSection::AddDatabaseReferences(PdbDatabaseReference *database_reference)
{
    database_reference_.push_back(database_reference);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbDatabaseReferenceSection::Print(ostream &out)
{
    for(DatabaseReferenceVector::iterator it = database_reference_.begin(); it != database_reference_.end(); it++)
            (*it)->Print(out);
    out << endl;
}
