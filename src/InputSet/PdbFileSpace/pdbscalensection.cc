#include "../../../includes/InputSet/PdbFileSpace/pdbscalensection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbscalencard.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbScaleNSection::PdbScaleNSection() {}
PdbScaleNSection::PdbScaleNSection(stringstream &stream_block)
{
    string line;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        stringstream ss;
        ss << line << endl;
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
void PdbScaleNSection::Print(ostream &out)
{
    for(PdbScaleNSection::ScaleNCardVector::iterator it = scale_n_.begin(); it != scale_n_.end(); it++)
    {
        (*it)->Print(out);
    }
    out << endl;
}
