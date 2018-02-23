
#include "../../../includes/InputSet/PdbFileSpace/pdbobsoletecard.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbObsoleteCard::PdbObsoleteCard() : replacement_date_(""), identifier_codes_() {}
PdbObsoleteCard::PdbObsoleteCard(const string &replacement_date,
                                const std::vector<std::string> &identifier_codes) :
                                replacement_date_(replacement_date),
                                identifier_codes_(identifier_codes) {}

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
    string all_identifier_codes_ = ss.str();
    std::size_t pos = 0, found;
    while ((found = all_identifier_codes_.find_first_of('      ', pos)) != std::string::npos)
    {
      identifier_codes_.push_back(all_identifier_codes_.substr(pos, found - pos));
      pos = found + 1;
    }
    identifier_codes_.push_back(all_identifier_codes_.substr(pos));
}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbObsoleteCard::GetReplacementDate()
{
    return replacement_date_;
}

vector<std::string> PdbObsoleteCard::GetIdentifierCodes()
{
    return identifier_codes_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbObsoleteCard::SetReplacementDate(const string replacement_date)
{
    replacement_date_ = replacement_date;
}

void PdbObsoleteCard::SetIdentifierCodes(const vector<std::string> identifier_codes)
{
    identifier_codes_ = identifier_codes;
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
    for (int i = 0; i < identifier_codes_.size(); i++)
    {
      out << "Identifier Codes: " << identifier_codes_[i] << ", ";
    }
}
