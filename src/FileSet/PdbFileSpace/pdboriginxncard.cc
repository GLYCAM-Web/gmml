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
    string temp = line;
    while (!Trim(temp).empty())
    {
        stringstream ss;
        ss << line << endl;
        PdbOriginXn* origin = new PdbOriginXn(ss);
        AddOriginXN(origin);
        getline(stream_block, line);
        temp = line;
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

void PdbOriginXnCard::SetOriginXN(OriginXnVector origin_x_n){
    origin_x_n_.clear();
    for(OriginXnVector::iterator it = origin_x_n.begin(); it != origin_x_n.end(); it++)
    {
        origin_x_n_.push_back(*it);
    }
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
void PdbOriginXnCard::Print(ostream &out)
{
    for(PdbOriginXnCard::OriginXnVector::iterator it = origin_x_n_.begin(); it != origin_x_n_.end(); it++)
    {
        (*it)->Print(out);
    }
    out << endl;
}
