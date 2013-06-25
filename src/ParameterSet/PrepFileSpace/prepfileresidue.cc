#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"

using namespace std;
using namespace gmml;
using namespace PrepFileSpace;

/////////////////////////////// CONSTRUCTOR ///////////////////////////////
PrepFileResidue::PrepFileResidue() {}

///////////////////////////// FUNCTION ////////////////////////////////////
PrepFileResidue* PrepFileResidue::LoadFromStream(std::ifstream& in_file)
{
    string line, name, dummy_atom_type;
    istringstream ss;
    CoordinateType coordinate_type;
    OutputFormat output_format;
    GeometryType geometry_type;
    DummyAtomOmission dummy_atom_omission;
    DummyAtomPosition dummy_atom_position;
    PrepFileResidue *residue = new PrepFileResidue();

    getline(in_file, line);             // Read the first line of a residue section
    if (trim(line) == "STOP")           // End of file
        return NULL;

    residue->title_ = line;            // Set title of the residue
    getline(in_file, line);             // Blank line, skip
    getline(in_file, line);             // Read the next line

    ss.str(line);                       // Create an stream from the read line

    // Residue name extraction from the stream of the read line
    name = ExtractResidueName(ss);
    residue->name_ = name;

    // Coordinate type extraction from the stream of the read line
    coordinate_type = ExtractResidueCoordinateType(ss);
    residue->coordinate_type_ = coordinate_type;

    // Output format extraction from the stream of the read line
    output_format = ExtractResidueOutputFormat(ss);
    residue->output_format_ = output_format;

    ss.clear();

    getline(in_file, line);             // Read the next line

    ss.str(line);                       // Create an stream from the read line

    // Geometry type extraction from the stream of the read line
    geometry_type = ExtractResidueGeometryType(ss);
    residue->geometry_type_ = geometry_type;

    // Dummy atom omission extraction from the stream of the read line
    dummy_atom_omission = ExtractResidueDummyAtomOmission(ss);
    residue->dummy_atom_omission_ = dummy_atom_omission;

    // Dummy atom type extraction from the stream of the read line
    ss >> dummy_atom_type;
    residue->dummy_atom_type_ = dummy_atom_type;

    // Dummy atom position extraction from the stream of the read line
    dummy_atom_position = ExtractResidueDummyAtomPosition(ss);
    residue->dummy_atom_position_ = dummy_atom_position;

    getline(in_file, line);             // Read the next line

    // Residue charge extraction from the read line
    residue->charge_ = convert_string<double>(line);

    // Process atoms of the residue
    while (getline(in_file, line) && !trim(line).empty())
    {
        PrepFileAtom *atom = new PrepFileAtom(line);
        residue->atoms_.push_back(atom);
    }

    bool done = false;
    while (!done)
    {
        while (getline(in_file, line) && trim(line).empty())
        {
        }
        switch (ExtractSectionType(line))
        {
            case kSectionLoop:
                //residue->impl_->set_loops(in);
                break;
            case kSectionImproper:
                //residue->impl_->set_improper_dihedrals(in);
                break;
            case kSectionDone:
                done = true;
                break;
            case kSectionOther:
                cout << "Unrecognized section in prep file";
                break;
        }
    }

    return residue;
}

std::string PrepFileResidue::ExtractResidueName(std::istream& ss)
{
    string name;
    ss >> name;
    return name;
}

CoordinateType PrepFileResidue::ExtractResidueCoordinateType(std::istream &ss)
{
    string s;
    ss >> s;
    if (s == "XYZ")
        return kXYZ;
    else
        return kINT;
}

OutputFormat PrepFileResidue::ExtractResidueOutputFormat(std::istream& ss)
{
    int val;
    ss >> val;
    if (val == 1)
        return kBinary;
    else
        return kFormatted;
}

GeometryType PrepFileResidue::ExtractResidueGeometryType(std::istream& ss)
{
    string s;
    ss >> s;
    if (s == "CHANGE")
        return kGeometryChange;
    else
        return kGeometryCorrect;
}

DummyAtomOmission PrepFileResidue::ExtractResidueDummyAtomOmission(istream &ss)
{
    string s;
    ss >> s;
    if (s == "NOMIT")
        return kNomit;
    else
        return kOmit;
}

DummyAtomPosition PrepFileResidue::ExtractResidueDummyAtomPosition(istream &ss)
{
    string s;
    ss >> s;
    if (s == "ALL")
        return kPositionAll;
    else
        return kPositionBeg;
}

SectionType PrepFileResidue::ExtractSectionType(string &line)
{
    if (line == "LOOP")
        return kSectionLoop;
    else if (line == "IMPROPER")
        return kSectionImproper;
    else if (line == "DONE")
        return kSectionDone;
    return kSectionOther;
}

////////////////////////// DISPLAY FUNCTION ///////////////////////////////
void PrepFileResidue::Print(std::ostream& out)
{

}
