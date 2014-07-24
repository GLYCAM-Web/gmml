#include "../../../includes/FileSet/PdbFileSpace/pdblink.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdblinkresidue.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbLink::PdbLink() {}

PdbLink::PdbLink(string &line)
{

    char temp1, temp2, temp3;
    int temp4, temp5;
    string temp6, temp7;
    temp6 = line.substr(12, 4);
    temp6 = Trim(temp6);
    temp7 = line.substr(17, 3);
    temp7 = Trim(temp7);
    if(line.substr(16,1) == " ")
        temp1 = ' ';
    else
        temp1 = ConvertString<char>(line.substr(16,1));
    if(line.substr(21,1) == " ")
        temp2 = ' ';
    else
        temp2 = ConvertString<char>(line.substr(21,1));
    if(line.substr(26,1) == " ")
        temp3 = ' ';
    else
        temp3 = ConvertString<char>(line.substr(26,1));
    if(line.substr(22,4) == "    ")
        temp4 = iNotSet;
    else
        temp4 = ConvertString<int>(line.substr(22,4));
    if(line.substr(59, 6) == "      ")
        temp5 = iNotSet;
    else
        temp5 = ConvertString<int>(line.substr(59,6));
    PdbLinkResidue* residue_1 = new  PdbLinkResidue(temp6, temp1, temp7, temp2, temp4, temp3, temp5);

    temp6 = line.substr(42,4);
    temp6 = Trim(temp6);
    temp7 = line.substr(47, 3);
    temp7 = Trim(temp7);
    if(line.substr(46,1) == " ")
        temp1 = ' ';
    else
        temp1 = ConvertString<char>(line.substr(46,1));
    if(line.substr(51,1) == " ")
        temp2 = ' ';
    else
        temp2 = ConvertString<char>(line.substr(51,1));
    if(line.substr(56,1) == " ")
        temp3 = ' ';
    else
        temp3 = ConvertString<char>(line.substr(56,1));
    if(line.substr(52, 4) == "    ")
        temp4 = iNotSet;
    else
        temp4 = ConvertString<int>(line.substr(52,4));
    if(line.substr(66, 6) == "      ")
        temp5 = iNotSet;
    else
        temp5 = ConvertString<int>(line.substr(66,6));
    PdbLinkResidue* residue_2 = new  PdbLinkResidue(temp6, temp1, temp7, temp2, temp4, temp3, temp5);
    residues_.push_back(residue_1);
    residues_.push_back(residue_2);
    if(line.substr(73, 5) == "     ")
        link_length_ = dNotSet;
    else
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

void PdbLink::SetResidues(LinkResidueVector residues){
    residues_.clear();
    for(LinkResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        residues_.push_back(*it);
    }
}

void PdbLink::AddResidue(PdbLinkResidue *residue)
{
    residues_.push_back(residue);
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
void PdbLink::Print(ostream &out)
{
    out << "------------- Residues ---------------" << endl;
    for(PdbLink::LinkResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        (*it)->Print(out);
        out << endl;
    }
    out << "Linke Lenght: " << link_length_ << endl << endl;
}
