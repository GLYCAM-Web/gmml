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
PdbHelix::PdbHelix() : helix_id_(""), helix_serial_number_(kNotSet), helix_class_(POLYPROLINE), comment_(""), helix_length_(kNotSet) {}

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
    line = Trim(line);
    while (!Trim(line).empty())
    {
        helix_id_ = line.substr(11, 3);
        helix_serial_number_ = ConvertString<int>(line.substr(7,3));

        PdbHelixResidue* initial_residue = new PdbHelixResidue(line.substr(15, 3), ConvertString<char>(line.substr(19, 1)), ConvertString<int>(line.substr(21, 4)), ConvertString<char>(line.substr(25, 1)));
        PdbHelixResidue* terminal_residue = new PdbHelixResidue(line.substr(27, 3), ConvertString<char>(line.substr(31, 1)), ConvertString<int>(line.substr(33, 4)), ConvertString<char>(line.substr(37, 1)));

        helix_residues_.push_back(initial_residue);
        helix_residues_.push_back(terminal_residue);

        int helix_class = ConvertString<int>(line.substr(38,2));
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
        }
        comment_ = line.substr(40, 30);
        helix_length_ = ConvertString<double>(line.substr(71,5));

        getline(stream_block, line);
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

