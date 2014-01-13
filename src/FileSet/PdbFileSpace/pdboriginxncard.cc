#include "../../../includes/FileSet/PdbFileSpace/pdboriginxn.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdboriginxncard.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbOriginXnCard::PdbOriginXnCard() {}

PdbOriginXnCard::PdbOriginXnCard(stringstream &stream_block)
{
    string line;
    getline(stream_block, line);
    line = Trim(line);
    while (!Trim(line).empty())
    {
        stringstream ss(line);
        PdbOriginXn* origin = new PdbOriginXn(ss);
        AddOriginXN(origin);
        getline(stream_block, line);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbOriginXnCard::OriginXnVector PdbOriginXnCard::GetOriginXN(){
    return origin_x_n_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbOriginXnCard::SetOriginXN(const OriginXnVector origin_x_n){
    origin_x_n_ = origin_x_n;
}

void PdbOriginXnCard::AddOriginXN(PdbOriginXn *origin)
{
    origin_x_n_.push_back(origin);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////




