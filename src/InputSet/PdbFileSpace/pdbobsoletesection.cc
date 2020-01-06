
#include "../../../includes/InputSet/PdbFileSpace/pdbobsoletecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbobsoletesection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbObsoleteSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbObsoleteSection::PdbObsoleteSection() : record_name_("OBSLTE"),
                                          continuation_(""),
                                          replacement_date_(""),
                                          identifier_codes_()  {}

PdbObsoleteSection::PdbObsoleteSection(const std::string& record_name,
                                        const std::string& continuation,
                                        const std::string& replacement_date,
                                        const std::vector<std::string>& identifier_codes) :
                                        record_name_(record_name),
                                        continuation_(continuation),
                                        replacement_date_(replacement_date),
                                        identifier_codes_(identifier_codes) {}

PdbObsoleteSection::PdbObsoleteSection(std::stringstream &stream_block)
{
    std::string line;
    bool is_record_name_set = false;
    bool is_continuation_set = false;
    getline(stream_block, line);
    std::string temp = line;
    int position = 0;
    while (!gmml::Trim(temp).empty())
    {
      if(line.find("OBSLTE") != std::string::npos)
      {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            gmml::Trim(record_name_);
            is_record_name_set=true;
        }
        if(!is_continuation_set){
            continuation_ = line.substr(8,2);
            gmml::Trim(continuation_);
            is_continuation_set=true;
        }
        std::stringstream obsolete_block;
        obsolete_block << line << "\n";
        std::cout << position;
        PdbObsoleteCard* obsolete_card = new PdbObsoleteCard(obsolete_block);
        std::cout << obsolete_card->GetReplacementDate() << "\n";
        obsolete_cards_.push_back(obsolete_card);
        replacement_date_ = obsolete_cards_[position]->GetReplacementDate();
        std::vector<std::string> IDcodes = obsolete_cards_[position]->GetIdentifierCodes();
        identifier_codes_.insert(identifier_codes_.end(), IDcodes.begin(), IDcodes.end());
        position++;
        getline(stream_block, line);
        temp = line;
      }
    }

}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbObsoleteSection::GetRecordName()
{
    return record_name_;
}

std::string PdbObsoleteSection::GetContinuation()
{
    return continuation_;
}

std::string PdbObsoleteSection::GetReplacementDate()
{
    return replacement_date_;
}

std::vector<std::string> PdbObsoleteSection::GetIdentifierCodes()
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
void PdbObsoleteSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbObsoleteSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "=========== Obsolete =============" << std::endl;
    for(PdbObsoleteSection::ObsoleteCardVector::iterator it = obsolete_cards_.begin(); it != obsolete_cards_.end(); it++)
    {
        out << "Obsolete Card: ";
        (*it)->Print(out);
        out << std::endl;
    }
}
