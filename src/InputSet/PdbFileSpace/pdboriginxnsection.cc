#include "../../../includes/InputSet/PdbFileSpace/pdboriginxnsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdboriginxncard.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbOriginXnSection::PdbOriginXnSection() {}

PdbOriginXnSection::PdbOriginXnSection(stringstream &stream_block)
{
    string line;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        stringstream ss;
        ss << line << endl;
        PdbOriginXnCard* origin = new PdbOriginXnCard(ss);
        AddOriginXN(origin);
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbOriginXnSection::OriginXnCardVector PdbOriginXnSection::GetOriginXN(){
    return origin_x_n_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbOriginXnSection::SetOriginXN(OriginXnCardVector origin_x_n){
    origin_x_n_.clear();
    for(OriginXnCardVector::iterator it = origin_x_n.begin(); it != origin_x_n.end(); it++)
    {
        origin_x_n_.push_back(*it);
    }
}

void PdbOriginXnSection::AddOriginXN(PdbOriginXnCard *origin)
{
    origin_x_n_.push_back(origin);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbOriginXnSection::Print(ostream &out)
{
    for(PdbOriginXnSection::OriginXnCardVector::iterator it = origin_x_n_.begin(); it != origin_x_n_.end(); it++)
    {
        (*it)->Print(out);
    }
    out << endl;
}
