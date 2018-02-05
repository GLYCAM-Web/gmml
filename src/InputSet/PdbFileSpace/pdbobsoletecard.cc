
#include "../../../includes/InputSet/PdbFileSpace/pdbobsoletecard.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbObsoleteCard::PdbObsoleteCard() : replacement_date_(""), replacement_entries_() {}
PdbObsoleteCard::PdbObsoleteCard(const string &replacement_date, const vector<std::string> &replacement_entries)
    : replacement_date_(replacement_date), replacement_entries_(replacement_entries) {}

PdbObsoleteCard::PdbObsoleteCard(stringstream& obsolete_block)
{
    string line;
    bool is_replacement_date_set = false;
    stringstream ss;
    getline(obsolete_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(!is_replacement_date_set){
            replacement_date_ = line.substr(12,8);
            Trim(replacement_date_);
            is_replacement_date_set = true;
        }
        ss << line.substr(21,54);

        getline(obsolete_block, line);
        temp = line;
    }
    string all_replacement_entries_ = ss.str();
    std::size_t pos = 0, found;
    while ((found = all_replacement_entries_.find_first_of('     ', pos)) != std::string::npos)
    {
      replacement_entries_.push_back(all_replacement_entries_.substr(pos, found - pos));
      pos = found + 1;
    }
    replacement_entries_.push_back(all_replacement_entries_.substr(pos));
}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbObsoleteCard::GetReplacementDate()
{
    return replacement_date_;
}

vector<std::string> PdbObsoleteCard::GetReplacementEntries()
{
    return replacement_entries_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbObsoleteCard::SetReplacementDate(const string replacement_date)
{
    replacement_date_ = replacement_date;
}

void PdbObsoleteCard::SetReplacementEntries(const vector<std::string> replacement_entries)
{
    replacement_entries_ = replacement_entries;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void  PdbObsoleteCard::Print(ostream &out)
{
    out << "Replacement Date: " << replacement_date_
        << " ";
    for (int i = 0; i < replacement_entries_.size(); i++)
     out << "Replacement Entries: " << replacement_entries_[i] << ", ";
}
