
#include "../../../includes/InputSet/PdbFileSpace/pdbobsoletecard.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbObsoleteCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbObsoleteCard::PdbObsoleteCard() : replacement_date_(""), identifier_codes_() {}
PdbObsoleteCard::PdbObsoleteCard(const std::string &replacement_date,
                                const std::vector<std::string> &identifier_codes) :
                                replacement_date_(replacement_date),
                                identifier_codes_(identifier_codes) {}

PdbObsoleteCard::PdbObsoleteCard(std::stringstream& obsolete_block)
{
    std::string line;
    bool is_replacement_date_set = false;
    std::stringstream ss;
    getline(obsolete_block, line);
    std::string temp = line;
    std::cout << line << "\n";
    while (!gmml::Trim(temp).empty())
    {
      if(line.find("OBSLTE") != std::string::npos)
      {
        if(!is_replacement_date_set){
            replacement_date_ = line.substr(12,8);
            gmml::Trim(replacement_date_);
            is_replacement_date_set = true;
        }
        ss << line.substr(21,54);

        getline(obsolete_block, line);
        temp = line;
      }
    }
    std::string all_identifier_codes_ = ss.str();
    std::size_t pos = 0, found;
    while ((found = all_identifier_codes_.find_first_of("      ", pos)) != std::string::npos)
    {
      identifier_codes_.push_back(all_identifier_codes_.substr(pos, found - pos));
      pos = found + 1;
    }
    identifier_codes_.push_back(all_identifier_codes_.substr(pos));
}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbObsoleteCard::GetReplacementDate()
{
    return replacement_date_;
}

std::vector<std::string> PdbObsoleteCard::GetIdentifierCodes()
{
    return identifier_codes_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbObsoleteCard::SetReplacementDate(const std::string replacement_date)
{
    replacement_date_ = replacement_date;
}

void PdbObsoleteCard::SetIdentifierCodes(const std::vector<std::string> identifier_codes)
{
    identifier_codes_ = identifier_codes;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void  PdbObsoleteCard::Print(std::ostream &out)
{
    out << "Replacement Date: " << replacement_date_
        << " ";
    for (unsigned int i = 0; i < identifier_codes_.size(); i++)
    {
      out << "Identifier Codes: " << identifier_codes_[i] << ", ";
    }
}
