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

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
PrepFileResidue::PrepFileResidue() {}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
/// Return the index of an atom in a specific residue by a given name
int PrepFileResidue::GetAtomIndexByName(const std::string& name)
{
    for (unsigned int i = 0; i < atoms_.size(); i++)
        if (name == atoms_[i]->name_)
            return i;
    return -1;
}

///Delaram

std::string PrepFileResidue::GetTitle(){
    return title_;
}

std::string PrepFileResidue::GetName(){
    return name_;
}

CoordinateType PrepFileResidue::GetCoordinateType(){
    return coordinate_type_;
}

OutputFormat PrepFileResidue::GetOutputFormat(){
    return output_format_;
}

GeometryType PrepFileResidue::GetGeometryType(){
    return geometry_type_;
}

DummyAtomOmission PrepFileResidue::GetDummyAtomOmission(){
    return dummy_atom_omission_;
}

std::string PrepFileResidue::GetDummyAtomType(){
    return dummy_atom_type_;
}

DummyAtomPosition PrepFileResidue::GetDummyAtomPosition(){
    return dummy_atom_position_;
}

double PrepFileResidue::GetCharge(){
    return charge_;
}

std::vector<PrepFileAtom*> PrepFileResidue::GetAtoms(){
    return atoms_;
}

std::vector<PrepFileResidue::Dihedral> PrepFileResidue::GetImproperDihedrals(){
    return improper_dihedrals_;
}

PrepFileResidue::Loop PrepFileResidue::GetLoops(){
    return loops_;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////

void PrepFileResidue::SetTitle(const std::string title){
    title_ = title;
}

void PrepFileResidue::SetName(const std::string name){
    name_ = name;
}

void PrepFileResidue::SetCoordinateType(CoordinateType coordinate_type){
    coordinate_type_ = coordinate_type;
}

void PrepFileResidue::SetOutputFormat(OutputFormat output_format){
    output_format_ = output_format;
}

void PrepFileResidue::SetBondIndex(GeometryType geometry_type){
    geometry_type_ = geometry_type;
}

void PrepFileResidue::GetDummyAtomOmission(DummyAtomOmission dummy_atom_omission){
    dummy_atom_omission_ = dummy_atom_omission;
}

void PrepFileResidue::SetDummyAtomType(const std::string dummy_atom_type){
    dummy_atom_type_ = dummy_atom_type;
}

void PrepFileResidue::SetDummyAtomPosition(DummyAtomPosition dummy_atom_position){
    dummy_atom_position_ = dummy_atom_position;
}

void PrepFileResidue::SetCharge(double charge){
    charge_ = charge;
}

void PrepFileResidue::SetAtoms(std::vector<PrepFileAtom*> atoms){
    atoms_ = atoms;
}

void PrepFileResidue::AddAtom(PrepFileAtom* atom){
    atoms_.push_back(atom);
}

void PrepFileResidue::SetImproperDihedrals(std::vector<Dihedral> improper_dihedrals){
    improper_dihedrals_ = improper_dihedrals;
}

void PrepFileResidue::AddImproperDihedral(Dihedral improper_dihedral){
    improper_dihedrals_.push_back(improper_dihedral);
}

void PrepFileResidue::SetLoops(Loop loops){
    loops_ = loops;
}

///Delaram

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
/// Create a new residue from a given stream
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

    getline(in_file, line);             /// Read the first line of a residue section
    if (Trim(line) == "STOP")           /// End of file
        return NULL;

    residue->title_ = line;             /// Set title of the residue
    getline(in_file, line);             /// Blank line, skip
    getline(in_file, line);             /// Read the next line

    ss.str(line);                       /// Create an stream from the read line

    /// Residue name extraction from the stream of the read line
    name = ExtractResidueName(ss);
    residue->name_ = name;

    /// Coordinate type extraction from the stream of the read line
    coordinate_type = ExtractResidueCoordinateType(ss);
    residue->coordinate_type_ = coordinate_type;

    /// Output format extraction from the stream of the read line
    output_format = ExtractResidueOutputFormat(ss);
    residue->output_format_ = output_format;

    ss.clear();

    getline(in_file, line);             /// Read the next line

    ss.str(line);                       /// Create an stream from the read line

    /// Geometry type extraction from the stream of the read line
    geometry_type = ExtractResidueGeometryType(ss);
    residue->geometry_type_ = geometry_type;

    /// Dummy atom omission extraction from the stream of the read line
    dummy_atom_omission = ExtractResidueDummyAtomOmission(ss);
    residue->dummy_atom_omission_ = dummy_atom_omission;

    /// Dummy atom type extraction from the stream of the read line
    ss >> dummy_atom_type;
    residue->dummy_atom_type_ = dummy_atom_type;

    /// Dummy atom position extraction from the stream of the read line
    dummy_atom_position = ExtractResidueDummyAtomPosition(ss);
    residue->dummy_atom_position_ = dummy_atom_position;

    getline(in_file, line);             /// Read the next line

    /// Residue charge extraction from the read line
    residue->charge_ = ConvertString<double>(line);

    /// Process atoms of the residue
    while (getline(in_file, line) && !Trim(line).empty())
    {
        PrepFileAtom *atom = new PrepFileAtom(line);
        residue->atoms_.push_back(atom);
    }

    /// Process the extra sections: IMPROPER, LOOP, DONE
    bool done = false;
    while (!done)
    {
        /// Skip blank lines until to reach to a known section title
        getline(in_file, line);
        while (Trim(line).empty())
        {
            getline(in_file, line);
        }
        /// Does a corresponding action based on the section title
        switch (ExtractSectionType(line))
        {
            case kSectionLoop:
                residue->loops_ = residue->ExtractLoops(in_file);
                break;
            case kSectionImproper:
                residue->improper_dihedrals_ = residue->ExtractImproperDihedral(in_file);
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

/// Return residue name from a stream line which is the first column of the 3rd line in each residue section
std::string PrepFileResidue::ExtractResidueName(std::istream& ss)
{
    string name;
    ss >> name;
    return name;
}

/// Return coordinate type from a stream line which is the 2nd column of the 3rd line in each residue section
CoordinateType PrepFileResidue::ExtractResidueCoordinateType(std::istream &ss)
{
    string s;
    ss >> s;
    if (s == "XYZ")
        return kXYZ;
    else
        return kINT;
}

/// Return output format from a stream line which is the 3rd column of the 3rd line in each residue section
OutputFormat PrepFileResidue::ExtractResidueOutputFormat(std::istream& ss)
{
    int val;
    ss >> val;
    if (val == 1)
        return kBinary;
    else
        return kFormatted;
}

/// Return geometry type from a stream line which is the first column of the 4th line in each residue section
GeometryType PrepFileResidue::ExtractResidueGeometryType(std::istream& ss)
{
    string s;
    ss >> s;
    if (s == "CHANGE")
        return kGeometryChange;
    else
        return kGeometryCorrect;
}

/// Return dummy atom omission from a stream line which is the 2nd column of the 4th line in each residue section
DummyAtomOmission PrepFileResidue::ExtractResidueDummyAtomOmission(istream &ss)
{
    string s;
    ss >> s;
    if (s == "NOMIT")
        return kNomit;
    else
        return kOmit;
}

/// Return dummy atom position from a stream line which is the 4th column of the 4th line in each residue section
DummyAtomPosition PrepFileResidue::ExtractResidueDummyAtomPosition(istream &ss)
{
    string s;
    ss >> s;
    if (s == "ALL")
        return kPositionAll;
    else
        return kPositionBeg;
}

/// Return a corresponding title from a stream line which may appear in each residue section
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

/// Parse the loop section of each residue section and return a loop map
PrepFileResidue::Loop PrepFileResidue::ExtractLoops(ifstream &in_file)
{
    Loop loops;
    string line;
    std::stringstream ss;

    getline(in_file, line);
    while (!Trim(line).empty())                         /// Read file until blank line which determines the end of the section
    {
        ss.clear();
        ss.str(line);                                   /// Create a stream from the read line
        string atom_names[2];
        ss >> atom_names[0] >> atom_names[1];           /// Extract atom names from the stream

        int from = GetAtomIndexByName(atom_names[0]);   /// Extract index of the first atom in the loop
        int to = GetAtomIndexByName(atom_names[1]);     /// Extract index of the second atom in the loop

        if (from == -1 || to == -1)
        {
            /// throw error here, unknown atom names
        }
        loops[from] = to;                               /// Add a new entry into the loop map of the residue
        getline(in_file, line);
    }
    return loops;
}

/// Parse the improper dihedral section of each residue section and return a vector of improper dihedrals
vector<PrepFileResidue::Dihedral> PrepFileResidue::ExtractImproperDihedral(ifstream &in_file)
{
    string line;
    std::stringstream ss;
    vector<Dihedral> dihedrals;
    getline(in_file, line);
    while (!Trim(line).empty())                         /// Read file until blank line which determines the end of the section
    {
        string atom_names[4];
        Dihedral dihedral;
        ss.clear();
        ss.str(line);
        /// Extract improper atom types involving in a dihedral from each line
        ss >> atom_names[0]
           >> atom_names[1]
           >> atom_names[2]
           >> atom_names[3];
        /// Push all atoms into a vector of atom types
        for(int i = 0; i < 4; i++)
        {
            dihedral.push_back(atom_names[i]);
        }
        dihedrals.push_back(dihedral);                  /// Create a new dihedral into the vector of dihedrals
        getline(in_file, line);                         /// Read the next line
    }
    return dihedrals;
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void PrepFileResidue::Print(std::ostream& out)
{
    out << "Title: " << title_ << endl;
    out << setw(10) << "ResName"
        << setw(10) << "CrdType"
        << setw(10) << "Output"
        << setw(10) << "GeoType"
        << setw(15) << "DummyOmission"
        << setw(12) << "DummyType"
        << setw(12) << "DummyPos"
        << setw(10) << "Charge" << endl;
    out << setw(10) << name_;

    if(coordinate_type_ == kINT)
        out << setw(10) << "INT";
    else if(coordinate_type_ == kXYZ)
        out << setw(10) << "XYZ";
    else
        out << setw(10) << "--";

    if(output_format_ == kBinary)
        out << setw(10) << "Binary";
    else if(output_format_ == kFormatted)
        out << setw(10) << "NBinary";
    else
        out << setw(10) << "--";

    if(geometry_type_ == kGeometryCorrect)
        out << setw(10) << "Correct";
    else if(geometry_type_ == kGeometryChange)
        out << setw(10) << "Change";
    else
        out << setw(10) << "--";

    if(dummy_atom_omission_ == kOmit)
        out << setw(15) << "YES";
    else if(dummy_atom_omission_ == kNomit)
        out << setw(15) << "NO";
    else
        out << setw(15) << "--";

    out << setw(12) << dummy_atom_type_;

    if(dummy_atom_position_ == kPositionAll)
        out << setw(12) << "ALL";
    else if (dummy_atom_position_ == kPositionBeg)
        out << setw(12) << "BEG";
    else
        out << setw(12) << "--";

    out << setw(10) << charge_ << endl << endl;

    out << setw(3) << "#"
        << setw(6) << "Name"
        << setw(6) << "Type"
        << setw(3) << "TT"
        << setw(4) << "B#"
        << setw(4) << "A#"
        << setw(4) << "D#"
        << setw(10) << "Bond"
        << setw(10) << "Angle"
        << setw(10) << "Dihedral"
        << setw(10) << "Charge"
        << endl;

    for(vector<PrepFileAtom*>::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        (*it)->Print(out);
    }

    out << endl << "Improper dihedrals" << endl;
    for(vector<Dihedral>::iterator it = improper_dihedrals_.begin(); it != improper_dihedrals_.end(); it++)
    {
        for(Dihedral::iterator it1 = it->begin(); it1 != it->end(); it1++)
        {
            out << setw(6) << (*it1);
        }
        out << endl;
    }

    out << endl << "Loops" << endl;
    for(Loop::iterator it = loops_.begin(); it != loops_.end(); it++)
    {
        out << setw(6) << it->first << setw(6) << it->second << endl;
    }

    out << endl;
}
