#include "../../../includes/InputSet/PdbFileSpace/pdbscalen.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbscalencard.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbScaleNCard::PdbScaleNCard() {}
PdbScaleNCard::PdbScaleNCard(stringstream &stream_block)
{
    string line;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        stringstream ss;
        ss << line << endl;
        PdbScaleN* scale = new PdbScaleN(ss);
        AddScaleN(scale);
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbScaleNCard::ScaleNVector PdbScaleNCard::GetScaleN(){
    return scale_n_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbScaleNCard::SetScaleN(ScaleNVector scale_n){
    scale_n_.clear();
    for(ScaleNVector::iterator it = scale_n.begin(); it != scale_n.end(); it++)
    {
        scale_n_.push_back(*it);
    }
}

void PdbScaleNCard::AddScaleN(PdbScaleN *scale)
{
    scale_n_.push_back(scale);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbScaleNCard::Print(ostream &out)
{
    for(PdbScaleNCard::ScaleNVector::iterator it = scale_n_.begin(); it != scale_n_.end(); it++)
    {
        (*it)->Print(out);
    }
    out << endl;
}
