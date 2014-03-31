// Author: Alireza Khatamian

#include "../../../includes/FileSet/PdbFileSpace/pdbhelix.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbhelixresidue.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHelix::PdbHelix() : helix_id_(""), helix_serial_number_(dNotSet), helix_class_(POLYPROLINE), comment_(""), helix_length_(dNotSet) {}

PdbHelix::PdbHelix(const string &helix_id, int helix_serial_number, HelixResidueVector helix_residues,
                   PdbHelixClass helix_class, const string &comment, double helix_length)
    : helix_id_(helix_id), helix_serial_number_(helix_serial_number), helix_class_(helix_class), comment_(comment), helix_length_(helix_length)
{
    helix_residues_.clear();
    for(PdbHelix::HelixResidueVector::const_iterator it = helix_residues.begin(); it != helix_residues.end(); it++)
    {
        helix_residues_.push_back(*it);
    }
}

PdbHelix::PdbHelix(stringstream& stream_block)
{
    string line;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        helix_id_ = line.substr(11, 3);
        Trim(helix_id_);
        if(line.substr(7, 3) == "   ")
            helix_serial_number_ = iNotSet;
        else
            helix_serial_number_ = ConvertString<int>(line.substr(7,3));

        string temp0;
        char temp1, temp2;
        int temp3;
        temp0 = line.substr(15, 3);
        temp0 = Trim(temp0);
        if(line.substr(19,1) == " ")
            temp1 = ' ';
        else
            temp1 = ConvertString<char>(line.substr(19, 1));
        if(line.substr(25,1) == " ")
            temp2 = ' ';
        else
            temp2 = ConvertString<char>(line.substr(25, 1));
        if(line.substr(21, 4) == "    ")
            temp3 = iNotSet;
        else
            temp3 = ConvertString<int>(line.substr(21, 4));
        PdbHelixResidue* initial_residue = new PdbHelixResidue(temp0, temp1, temp3, temp2);

        temp0 = line.substr(27,3);
        if(line.substr(31,1) == " ")
            temp1 = ' ';
        else
            temp1 = ConvertString<char>(line.substr(31, 1));
        if(line.substr(37,1) == " ")
            temp2 = ' ';
        else
            temp2 = ConvertString<char>(line.substr(37, 1));
        if(line.substr(33,4) == "    ")
            temp3 = iNotSet;
        else
            temp3 = ConvertString<int>(line.substr(33, 4));
        PdbHelixResidue* terminal_residue = new PdbHelixResidue(temp0, temp1, temp3, temp2);

        helix_residues_.push_back(initial_residue);
        helix_residues_.push_back(terminal_residue);


        int helix_class;
        if(line.substr(38, 2) == "  ")
            helix_class = iNotSet;
        else
            helix_class = ConvertString<int>(line.substr(38,2));
        switch(helix_class)
        {
            case 1:
                helix_class_ = RIGHT_HANDED_ALPHA;
                break;
            case 2:
                helix_class_ = RIGHT_HANDED_OMEGA;
                break;
            case 3:
                helix_class_ = RIGHT_HANDED_PI;
                break;
            case 4:
                helix_class_ = RIGHT_HANDED_GAMMA;
                break;
            case 5:
                helix_class_ = RIGHT_HANDED_310;
                break;
            case 6:
                helix_class_ = LEFT_HANDED_ALPHA;
                break;
            case 7:
                helix_class_ = LEFT_HANDED_OMEGA_;
                break;
            case 8:
                helix_class_ = LEFT_HANDED_GAMMA_;
                break;
            case 9 :
                helix_class_ = RIBBON_27;
                break;
            case 10:
                helix_class_ = POLYPROLINE;
                break;
            case iNotSet:
                helix_class_ = UnknownHelix;
        }
        comment_ = line.substr(40, 30);
        Trim(comment_);
        if(line.substr(71, 5) == "     ")
            helix_length_ = dNotSet;
        else
            helix_length_ = ConvertString<double>(line.substr(71,5));

        getline(stream_block, line);
        temp = line;
    }
}


//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbHelix::GetHelixId()
{
    return helix_id_;
}

int PdbHelix::GetHelixSerialNumber()
{
    return helix_serial_number_;
}

PdbHelix::HelixResidueVector PdbHelix::GetHelixResidues()
{
    return helix_residues_;
}

PdbHelixClass PdbHelix::GetHelixClass()
{
    return helix_class_;
}

string PdbHelix::GetComment()
{
    return comment_;
}

double PdbHelix::GetHelixLength()
{
    return helix_length_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHelix::SetHelixId(const string helix_id)
{
    helix_id_ = helix_id;
}

void PdbHelix::SetHelixSerialNumber(int helix_serial_number)
{
    helix_serial_number_ = helix_serial_number;
}

void PdbHelix::SetHelixResidues(const HelixResidueVector helix_residues)
{
    helix_residues_.clear();
    for(PdbHelix::HelixResidueVector::const_iterator it = helix_residues.begin(); it != helix_residues.end(); it++)
    {
        helix_residues_.push_back(*it);
    }
}

void PdbHelix::AddHelixResidue(PdbHelixResidue *helix_residue)
{
    helix_residues_.push_back(helix_residue);
}

void PdbHelix::SetHelixClass(PdbHelixClass helix_class)
{
    helix_class_ = helix_class;
}

void PdbHelix::SetComment(const string &comment)
{
    comment_ = comment;
}

void PdbHelix::SetHelixLength(double helix_length)
{
    helix_length_ = helix_length;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHelix::Print(ostream &out)
{
    out << "Helix ID: " << helix_id_
        << ", Helix Serial Number: ";
    if(helix_serial_number_ != iNotSet)
        out << helix_serial_number_;
    else
        out << " ";
    out << endl
        << "=============== Helix Residues ==============" << endl;
    for(PdbHelix::HelixResidueVector::iterator it = helix_residues_.begin(); it != helix_residues_.end(); it++)
    {
        (*it)->Print(out);
    }
    out << "Helix Class: ";
    if(helix_class_ != UnknownHelix)
        out << helix_class_;
    else
        out << " ";
    out << ", Comments: " << comment_
        << "Helix Length: ";
    if(helix_length_ != dNotSet)
        out << helix_length_ ;
    else
        out << " ";
    out << endl << endl;
}
