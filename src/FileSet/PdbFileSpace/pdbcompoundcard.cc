#include <iostream>

#include "../../../includes/FileSet/PdbFileSpace/pdbcompoundcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbcompoundspecification.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbCompoundCard::PdbCompoundCard() : record_name_("COMPND"){}
PdbCompoundCard::PdbCompoundCard(const string& record_name) : record_name_(record_name) {}

PdbCompoundCard::PdbCompoundCard(stringstream& stream_block)
{
    string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    line = Trim(line);
    while (!Trim(line).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            is_record_name_set=true;
        }
        stringstream ss, specification_block;
        ss << line.substr(10,70);
        specification_block << ss.str() << endl;

        getline(stream_block, line);

        while (line.find("MOL_ID") == string::npos && !Trim(line).empty()){
            stringstream sss;
            sss << line.substr(10,70);

            specification_block << sss.str() << endl;
            getline(stream_block, line);
        }
        PdbCompoundSpecification* compound_specification = new PdbCompoundSpecification(specification_block);
        compound_specifications_[compound_specification->GetMoleculeId()] = compound_specification;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbCompoundCard::GetRecordName()
{
    return record_name_;
}

PdbCompoundCard::PdbCompoundSpecificationMap PdbCompoundCard::GetCompoundSpecifications()
{
    return compound_specifications_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbCompoundCard::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
