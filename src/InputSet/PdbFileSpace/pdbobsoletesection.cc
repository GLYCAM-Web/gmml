
#include "../../../includes/InputSet/PdbFileSpace/pdbobsoletecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbobsoletesection.hpp"
#include "../../../includes/utils.hpp"


using namespace std;
using namespace gmml;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbObsoleteSection::PdbObsoleteSection() : record_name_("OBSLTE") {}

PdbObsoleteSection::PdbObsoleteSection(const std::string& record_name, const std::string& continuation,
                  const std::string& replacement_date, const std::string& identifier_codes) :
                  record_name_(record_name), continuation_(continuation), replacement_date_(replacement_date),
                  identifier_codes_(identifier_codes) {}

PdbObsoleteSection::PdbObsoleteSection(stringstream &stream_block)
{
    string line;
    bool is_record_name_set = false;
    bool is_continuation_set = false;
    getline(stream_block, line);
    string temp = line;
    int i = 0;
    while (!Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            Trim(record_name_);
            is_record_name_set=true;
        }
        if(!is_continuation_set){
            continuation_ = line.substr(8,2);
            Trim(continuation_);
            is_continuation_set=true;
        }
        stringstream obsolete_block;
        obsolete_block << line << endl;
        getline(stream_block, line);
        temp = line;
        PdbObsoleteCard* obsolete_card = new PdbObsoleteCard(obsolete_block);
        obsolete_cards_[i] = obsolete_card;
        i++;
    }
    replacement_date_ = obsolete_cards_[0]->GetReplacementDate();
    identifier_codes_ = "";
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbObsoleteSection::GetRecordName()
{
    return record_name_;
}

string PdbObsoleteSection::GetContinuation()
{
    return continuation_;
}

string PdbObsoleteSection::GetReplacementDate()
{
    return replacement_date_;
}

string PdbObsoleteSection::GetIdentifierCodes()
{
    return identifier_codes_;
}

PdbObsoleteSection::ObsoleteCardVector PdbObsoleteSection::GetObsoleteCards()
{
    return obsolete_cards_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbObsoleteSection::SetRecordName(const string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbObsoleteSection::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "=========== Obsolete =============" << endl;
    for(PdbObsoleteSection::ObsoleteCardVector::iterator it = obsolete_cards_.begin(); it != obsolete_cards_.end(); it++)
    {
        out << "Obsolete Card: ";
        (*it)->Print(out);
        out << endl;
    }
}
