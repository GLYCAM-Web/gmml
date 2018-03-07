#include "../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbLinkCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbLinkCard::PdbLinkCard() {}

PdbLinkCard::PdbLinkCard(std::string &line)
{

    char temp1, temp2, temp3;
    int temp4, temp5;
    std::string temp6, temp7;
    temp6 = line.substr(12, 4);
    temp6 = gmml::Trim(temp6);
    temp7 = line.substr(17, 3);
    temp7 = gmml::Trim(temp7);
    if(line.substr(16,1) == " ")
        temp1 = ' ';
    else
        temp1 = gmml::ConvertString<char>(line.substr(16,1));
    if(line.substr(21,1) == " ")
        temp2 = ' ';
    else
        temp2 = gmml::ConvertString<char>(line.substr(21,1));
    if(line.substr(26,1) == " ")
        temp3 = ' ';
    else
        temp3 = gmml::ConvertString<char>(line.substr(26,1));
    if(line.substr(22,4) == "    ")
        temp4 = gmml::iNotSet;
    else
        temp4 = gmml::ConvertString<int>(line.substr(22,4));
    if(line.substr(59, 6) == "      ")
        temp5 = gmml::iNotSet;
    else
        temp5 = gmml::ConvertString<int>(line.substr(59,6));
    PdbLinkCardResidue* residue_1 = new  PdbLinkCardResidue(temp6, temp1, temp7, temp2, temp4, temp3, temp5);

    temp6 = line.substr(42,4);
    temp6 = gmml::Trim(temp6);
    temp7 = line.substr(47, 3);
    temp7 = gmml::Trim(temp7);
    if(line.substr(46,1) == " ")
        temp1 = ' ';
    else
        temp1 = gmml::ConvertString<char>(line.substr(46,1));
    if(line.substr(51,1) == " ")
        temp2 = ' ';
    else
        temp2 = gmml::ConvertString<char>(line.substr(51,1));
    if(line.substr(56,1) == " ")
        temp3 = ' ';
    else
        temp3 = gmml::ConvertString<char>(line.substr(56,1));
    if(line.substr(52, 4) == "    ")
        temp4 = gmml::iNotSet;
    else
        temp4 = gmml::ConvertString<int>(line.substr(52,4));
    if(line.substr(66, 6) == "      ")
        temp5 = gmml::iNotSet;
    else
        temp5 = gmml::ConvertString<int>(line.substr(66,6));
    PdbLinkCardResidue* residue_2 = new  PdbLinkCardResidue(temp6, temp1, temp7, temp2, temp4, temp3, temp5);
    residues_.push_back(residue_1);
    residues_.push_back(residue_2);
    if(line.substr(73, 5) == "     ")
        link_length_ = gmml::dNotSet;
    else
        link_length_ = gmml::ConvertString<double>(line.substr(73,5));
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbLinkCard::LinkResidueVector PdbLinkCard::GetResidues(){
    return residues_;
}

double PdbLinkCard::GetLinkLength(){
    return link_length_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbLinkCard::SetResidues(LinkResidueVector residues){
    residues_.clear();
    for(LinkResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        residues_.push_back(*it);
    }
}

void PdbLinkCard::AddResidue(PdbLinkCardResidue *residue)
{
    residues_.push_back(residue);
}

void PdbLinkCard::SetLinkLength(double link_length){
    link_length_ = link_length;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbLinkCard::Print(std::ostream &out)
{
    out << "------------- Residues ---------------" << std::endl;
    for(PdbLinkCard::LinkResidueVector::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        (*it)->Print(out);
        out << std::endl;
    }
    out << "Linke Lenght: " << link_length_ << std::endl << std::endl;
}
