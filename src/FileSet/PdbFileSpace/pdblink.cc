#include "../../../includes/FileSet/PdbFileSpace/pdblink.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdblinkresidue.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbLink::PdbLink() {}

PdbLink::PdbLink(string &line)
{
    PdbLinkResidue* residue_1 = new  PdbLinkResidue(line.substr(12, 4), ConvertString<char>(line.substr(16,1)), line.substr(17, 3),
                                                    ConvertString<char>(line.substr(21,1)), ConvertString<int>(line.substr(22,4)),
                                                    ConvertString<char>(line.substr(26,1)), ConvertString<int>(line.substr(59,6)));
    PdbLinkResidue* residue_2 = new  PdbLinkResidue(line.substr(42, 4), ConvertString<char>(line.substr(46,1)), line.substr(47, 3),
                                                    ConvertString<char>(line.substr(51,1)), ConvertString<int>(line.substr(52,4)),
                                                    ConvertString<char>(line.substr(56,1)), ConvertString<int>(line.substr(66,6)));
    residues_.push_back(residue_1);
    residues_.push_back(residue_2);
    link_length_ = ConvertString<double>(line.substr(73,5));

}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbLink::LinkResidueVector PdbLink::GetResidues(){
    return residues_;
}

double PdbLink::GetLinkLength(){
    return link_length_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbLink::SetResidues(const LinkResidueVector residues){
    residues_ = residues;
}

void PdbLink::SetLinkLength(double link_length){
    link_length_ = link_length;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////





