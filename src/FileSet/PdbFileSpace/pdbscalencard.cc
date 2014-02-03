#include "../../../includes/FileSet/PdbFileSpace/pdbscalen.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbscalencard.hpp"
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
void PdbScaleNCard::SetScaleN(const ScaleNVector scale_n){\
    scale_n_ = scale_n;
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
