// Created by: Dave Montgomery

#include "../../../includes/InputSet/PdbFileSpace/pdbremarksection.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbRemarkSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbRemarkSection::PdbRemarkSection() : record_name_("REMARK"){}

PdbRemarkSection::PdbRemarkSection(const std::string& record_name, const std::string& remark_cards)
{
    record_name_ = record_name;
    remark_cards_ = remark_cards;
}

PdbRemarkSection::PdbRemarkSection(std::stringstream &stream_block)
{
    std::string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(line.find("REMARK") != std::string::npos)
        {
            if(!is_record_name_set)
            {
                record_name_ = line.substr(0,6);
                gmml::Trim(record_name_);
                is_record_name_set=true;
            }
            std::stringstream remark_cardstream;
            while(line.find("REMARK") != std::string::npos)
            {
                remark_cardstream << line << std::endl;
                if (line.find("2 RESOLUTION.")!= std::string::npos)
                {
                    std::string tmp_resolution = line.substr(23,7);
                    gmml::Trim(tmp_resolution);
                    // if( tmp_resolution != "NULL" ) {
                    //     this->SetResolution( std::stof( tmp_resolution ) );
                    // }
                    try {
                        this->SetResolution( std::stof( tmp_resolution ) );
                    } catch (const std::invalid_argument& error) {
                        gmml::log(__LINE__, __FILE__, gmml::ERR, "RESOLUTION is not a valid float value. Value:\t" + tmp_resolution);
                    }
                }
                if (line.find("MEAN B VALUE")!= std::string::npos)
                {
                    int start = line.find(":") + 1;
                    std::string tmp_b_factor = line.substr(start,80-start);
                    gmml::Trim( tmp_b_factor );
                    if(tmp_b_factor == "NULL") {
                        this->SetBFactor(0);
                    }
                    try {
                        this->SetBFactor( std::stof( tmp_b_factor ) );
                    } catch(const std::invalid_argument& error) {
                        gmml::log(__LINE__, __FILE__, gmml::ERR, "MEAN B VALUE is not a valid float value. Value:\t" + tmp_b_factor);
                        this->SetBFactor(0);
                    }
                }
                getline(stream_block,line);
                temp = line;
            }
            remark_cards_ = remark_cardstream.str();
        }

    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

std::string PdbRemarkSection::GetRecordName(){
    return record_name_;
}

std::string PdbRemarkSection::GetRemarks(){
    return remark_cards_;
}

float PdbRemarkSection::GetResolution(){
    return resolution_;
}

float PdbRemarkSection::GetBFactor(){
    return b_factor_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbRemarkSection::SetRecordName(const std::string record_name){
    record_name_ = record_name;
}

void PdbRemarkSection::SetRemarks(const std::string remark_cards){
    remark_cards_ = remark_cards;
}

void PdbRemarkSection::SetResolution(const float resolution) {
    this->resolution_ = resolution;
}

void PdbRemarkSection::SetBFactor(const float b_factor) {
    this->b_factor_ = b_factor;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbRemarkSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           remark_cards_ << std::endl;

    // for(PdbRemarkSection::PdbRemarkMap::iterator it = remark_cards_.begin(); it != remark_cards_.end(); it++)
    // {
    //     out << "Remark Serial Number: ";
    //     if((it)->first != iNotSet)
    //         out << (it)->first << std::endl;
    //     else
    //         out << " " << std::endl;
    //     (it)->second->Print();
    //     out << std::endl;
    // }
}
