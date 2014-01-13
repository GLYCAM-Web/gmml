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
    line = Trim(line);
    while (!Trim(line).empty())
    {
        stringstream ss(line);
        PdbScaleN* scale = new PdbScaleN(ss);
        AddScaleN(scale);
        getline(stream_block, line);
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




