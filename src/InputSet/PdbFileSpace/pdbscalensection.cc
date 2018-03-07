#include "../../../includes/InputSet/PdbFileSpace/pdbscalensection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbscalencard.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbScaleNSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbScaleNSection::PdbScaleNSection() {}
PdbScaleNSection::PdbScaleNSection(std::stringstream &stream_block)
{
    std::string line;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        std::stringstream ss;
        ss << line << std::endl;
        PdbScaleNCard* scale = new PdbScaleNCard(ss);
        AddScaleNCard(scale);
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbScaleNSection::ScaleNCardVector PdbScaleNSection::GetScaleNCard(){
    return scale_n_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbScaleNSection::SetScaleNCard(ScaleNCardVector scale_n){
    scale_n_.clear();
    for(ScaleNCardVector::iterator it = scale_n.begin(); it != scale_n.end(); it++)
    {
        scale_n_.push_back(*it);
    }
}

void PdbScaleNSection::AddScaleNCard(PdbScaleNCard *scale)
{
    scale_n_.push_back(scale);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbScaleNSection::Print(std::ostream &out)
{
    for(PdbScaleNSection::ScaleNCardVector::iterator it = scale_n_.begin(); it != scale_n_.end(); it++)
    {
        (*it)->Print(out);
    }
    out << std::endl;
}
