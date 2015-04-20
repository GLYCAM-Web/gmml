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
PrepFileResidue::PrepFileResidue() : title_(""), name_(""), coordinate_type_(kINT), output_format_(kFormatted), geometry_type_(kGeometryCorrect),
    dummy_atom_omission_(kOmit), dummy_atom_type_("DU"), dummy_atom_position_(kPositionBeg), charge_(0.0)
{
    atoms_ = PrepFileAtomVector();
    improper_dihedrals_ = DihedralVector();
    loops_ = Loop();
}

PrepFileResidue::~PrepFileResidue()
{
    atoms_.clear();
    atoms_ = PrepFileAtomVector();
    improper_dihedrals_.clear();
    improper_dihedrals_ = DihedralVector();
    loops_.clear();
    loops_ = Loop();
}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
/// Return the index of an atom in a specific residue by a given name
int PrepFileResidue::GetAtomIndexByName(const std::string& name)
{
    for (unsigned int i = 0; i < atoms_.size(); i++)
        if (name.compare(atoms_[i]->GetName()) == 0)
            return atoms_[i]->GetIndex();
    return -1;
}

string PrepFileResidue::GetTitle(){
    return title_;
}

string PrepFileResidue::GetName(){
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

string PrepFileResidue::GetDummyAtomType(){
    return dummy_atom_type_;
}

DummyAtomPosition PrepFileResidue::GetDummyAtomPosition(){
    return dummy_atom_position_;
}

double PrepFileResidue::GetCharge(){
    return charge_;
}

PrepFileResidue::PrepFileAtomVector PrepFileResidue::GetAtoms(){
    return atoms_;
}

PrepFileResidue::DihedralVector PrepFileResidue::GetImproperDihedrals(){
    return improper_dihedrals_;
}

PrepFileResidue::Loop PrepFileResidue::GetLoops(){
    return loops_;
}
PrepFileAtom* PrepFileResidue::GetPrepAtomByName(string atom_name)
{
    for(PrepFileAtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        PrepFileAtom* prep_file_atom = (*it);
        string prep_file_atom_name = prep_file_atom->GetName();
        if(prep_file_atom_name.compare(atom_name) == 0)
            return prep_file_atom;
    }
    return NULL;
}
PrepFileResidue::BondedAtomIndexMap PrepFileResidue::GetBondingsOfResidue()
{
    BondedAtomIndexMap bonded_atoms_map = BondedAtomIndexMap();
//    vector<PrepFileAtom*> stack = vector<PrepFileAtom*>();
//    vector<int> number_of_bonds = vector<int>();

    PrepFileAtomVector parents = this->GetAtomsParentVector();
    for(Loop::iterator it = loops_.begin(); it != loops_.end(); it++)
    {
        int from = (*it).first;
        int to = (*it).second;

        bonded_atoms_map[from].push_back(to);
        bonded_atoms_map[to].push_back(from);
    }

    for(PrepFileAtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        PrepFileAtom* atom = *it;
        int index = distance(atoms_.begin(), it);
        int atom_index = atom->GetIndex();
        int parent_index = parents.at(index)->GetIndex();
        if(atom_index != parent_index)
        {
            bonded_atoms_map[atom_index].push_back(parent_index);
            bonded_atoms_map[parent_index].push_back(atom_index);
        }
    }

  /*
    for(PrepFileAtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        PrepFileAtom* atom = (*it);
        if(stack.empty())
        {
            if(atom->GetTopologicalType() == kTopTypeM || atom->GetTopologicalType() == kTopTypeS
                    || atom->GetTopologicalType() == kTopTypeB || atom->GetTopologicalType() == kTopType3)
            {
                stack.push_back(atom);
                if(atom->GetTopologicalType() == kTopTypeM)
                    number_of_bonds.push_back(4);
                if(atom->GetTopologicalType() == kTopTypeS)
                    number_of_bonds.push_back(2);
                if(atom->GetTopologicalType() == kTopTypeB)
                    number_of_bonds.push_back(3);
                if(atom->GetTopologicalType() == kTopType3)
                    number_of_bonds.push_back(4);
            }
        }
        if(!stack.empty())
        {
            if(atom->GetTopologicalType() == kTopTypeE)
            {
                PrepFileAtom* top_of_stack_atom = stack.at(stack.size()-1);
                bonded_atoms_map[top_of_stack_atom->GetIndex()].push_back(atom->GetIndex());
                bonded_atoms_map[atom->GetIndex()].push_back(top_of_stack_atom->GetIndex());
                number_of_bonds.at(number_of_bonds.size()-1)--;
                if(number_of_bonds.at(number_of_bonds.size()-1) == 0)
                {
                    number_of_bonds.pop_back();
                    stack.pop_back();
                }
            }
            if(atom->GetTopologicalType() == kTopTypeM)
            {
                PrepFileAtom* top_of_stack_atom = stack.at(stack.size()-1);
                if(top_of_stack_atom->GetTopologicalType() == kTopTypeM)
                {
                    if(top_of_stack_atom->GetType().compare(this->GetDummyAtomType()) == 0)
                    {
                        bonded_atoms_map[top_of_stack_atom->GetIndex()].push_back(atom->GetIndex());
                        bonded_atoms_map[atom->GetIndex()].push_back(top_of_stack_atom->GetIndex());
                        stack.pop_back();
                        number_of_bonds.pop_back();
                        stack.push_back(atom);
                        number_of_bonds.push_back(4);
                    }
                    else
                    {
                        bonded_atoms_map[top_of_stack_atom->GetIndex()].push_back(atom->GetIndex());
                        bonded_atoms_map[atom->GetIndex()].push_back(top_of_stack_atom->GetIndex());
                        number_of_bonds.at(number_of_bonds.size()-1)--;
                        if(number_of_bonds.at(number_of_bonds.size()-1) == 0)
                        {
                            number_of_bonds.pop_back();
                            stack.pop_back();
                        }
                        stack.push_back(atom);
                        number_of_bonds.push_back(3);
                    }
                }
                if(top_of_stack_atom->GetTopologicalType() == kTopTypeS || top_of_stack_atom->GetTopologicalType() == kTopTypeB
                        || top_of_stack_atom->GetTopologicalType() == kTopType3)
                {
                    bonded_atoms_map[top_of_stack_atom->GetIndex()].push_back(atom->GetIndex());
                    bonded_atoms_map[atom->GetIndex()].push_back(top_of_stack_atom->GetIndex());
                    number_of_bonds.at(number_of_bonds.size()-1)--;
                    if(number_of_bonds.at(number_of_bonds.size()-1) == 0)
                    {
                        number_of_bonds.pop_back();
                        stack.pop_back();
                    }
                    stack.push_back(atom);
                    number_of_bonds.push_back(3);
                }
            }
            if(atom->GetTopologicalType() == kTopTypeS)
            {
                PrepFileAtom* top_of_stack_atom = stack.at(stack.size()-1);
                bonded_atoms_map[top_of_stack_atom->GetIndex()].push_back(atom->GetIndex());
                bonded_atoms_map[atom->GetIndex()].push_back(top_of_stack_atom->GetIndex());
                number_of_bonds.at(number_of_bonds.size()-1)--;
                if(number_of_bonds.at(number_of_bonds.size()-1) == 0)
                {
                    number_of_bonds.pop_back();
                    stack.pop_back();
                }
                stack.push_back(atom);
                number_of_bonds.push_back(1);
            }
            if(atom->GetTopologicalType() == kTopTypeB)
            {
                PrepFileAtom* top_of_stack_atom = stack.at(stack.size()-1);
                bonded_atoms_map[top_of_stack_atom->GetIndex()].push_back(atom->GetIndex());
                bonded_atoms_map[atom->GetIndex()].push_back(top_of_stack_atom->GetIndex());
                number_of_bonds.at(number_of_bonds.size()-1)--;
                if(number_of_bonds.at(number_of_bonds.size()-1) == 0)
                {
                    number_of_bonds.pop_back();
                    stack.pop_back();
                }
                stack.push_back(atom);
                number_of_bonds.push_back(2);
            }
            if(atom->GetTopologicalType() == kTopType3)
            {
                PrepFileAtom* top_of_stack_atom = stack.at(stack.size()-1);
                bonded_atoms_map[top_of_stack_atom->GetIndex()].push_back(atom->GetIndex());
                bonded_atoms_map[atom->GetIndex()].push_back(top_of_stack_atom->GetIndex());
                number_of_bonds.at(number_of_bonds.size()-1)--;
                if(number_of_bonds.at(number_of_bonds.size()-1) == 0)
                {
                    number_of_bonds.pop_back();
                    stack.pop_back();
                }
                stack.push_back(atom);
                number_of_bonds.push_back(3);
            }
        }
    }
*/

    return bonded_atoms_map;
}
string PrepFileResidue::GetAtomNameByIndex(int atom_index)
{
    for(PrepFileAtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        PrepFileAtom* prep_file_atom = (*it);
        if(prep_file_atom->GetIndex() == atom_index)
            return prep_file_atom->GetName();
    }
    return NULL;
}
string PrepFileResidue::GetStringFormatOfCoordinateType(CoordinateType coordinate_type)
{
    switch(coordinate_type)
    {
        case kINT:
            return "INT";
        case kXYZ:
            return "XYZ";
        default:
            return "";
    }
}
string PrepFileResidue::GetStringFormatOfCoordinateType()
{
    switch(coordinate_type_)
    {
        case kINT:
            return "INT";
        case kXYZ:
            return "XYZ";
        default:
            return "";
    }
}
string PrepFileResidue::GetStringFormatOfOutputFormat(OutputFormat output_format)
{
    switch(output_format)
    {
        case kFormatted:
            return "FRM";
        case kBinary:
            return "BIN";
        default:
            return "";
    }
}
string PrepFileResidue::GetStringFormatOfOutputFormat()
{
    switch(output_format_)
    {
        case kFormatted:
            return "FRM";
        case kBinary:
            return "BIN";
        default:
            return "";
    }
}
string PrepFileResidue::GetStringFormatOfGeometryType(GeometryType geometry_type)
{
    switch(geometry_type)
    {
        case kGeometryCorrect:
            return "CORRECT";
        case kGeometryChange:
            return "CHANGE";
        default:
            return "";
    }
}
string PrepFileResidue::GetStringFormatOfGeometryType()
{
    switch(geometry_type_)
    {
        case kGeometryCorrect:
            return "CORRECT";
        case kGeometryChange:
            return "CHANGE";
        default:
            return "";
    }
}
string PrepFileResidue::GetStringFormatOfDummyAtomPosition(DummyAtomPosition dummy_atom_position)
{
    switch(dummy_atom_position)
    {
        case kPositionAll:
            return "ALL";
        case kPositionBeg:
            return "BEG";
        default:
            return "";
    }
}
string PrepFileResidue::GetStringFormatOfDummyAtomPosition()
{
    switch(dummy_atom_position_)
    {
        case kPositionAll:
            return "ALL";
        case kPositionBeg:
            return "BEG";
        default:
            return "";
    }
}
string PrepFileResidue::GetStringFormatOfDummyAtomOmission(DummyAtomOmission dummy_atom_omission)
{
    switch(dummy_atom_omission)
    {
        case kOmit:
            return "OMIT";
        case kNomit:
            return "NOMIT";
        default:
            return "";
    }
}
string PrepFileResidue::GetStringFormatOfDummyAtomOmission()
{
    switch(dummy_atom_omission_)
    {
        case kOmit:
            return "OMIT";
        case kNomit:
            return "NOMIT";
        default:
            return "";
    }
}
string PrepFileResidue::GetStringFormatOfSectionType(SectionType section_type)
{
    switch(section_type)
    {
        case kSectionLoop:
            return "SectionLoop";
        case kSectionImproper:
            return "SectionImproper";
        case kSectionDone:
            return "kSectionDone";
        case kSectionOther:
            return "SectionOther";
        default:
            return "";
    }
}
CoordinateType PrepFileResidue::GetCoordinateTypeFromString(string coordinate_type)
{
    if(coordinate_type.compare("INT") == 0)
        return kINT;
    if(coordinate_type.compare("XYZ") == 0)
        return kXYZ;
    else
        return kINT;
}
OutputFormat PrepFileResidue::GetOutputFormatFromString(string output_format)
{
    if(output_format.compare("Formatted") == 0)
        return kFormatted;
    if(output_format.compare("Binary") == 0)
        return kBinary;
    else
        return kBinary;
}
GeometryType PrepFileResidue::GetGeometryTypeFromString(string geometry_type)
{
    if(geometry_type.compare("GeometryCorrect") == 0)
        return kGeometryCorrect;
    if(geometry_type.compare("GeometryChange") == 0)
        return kGeometryChange;
    else
        return kGeometryCorrect;
}
DummyAtomPosition PrepFileResidue::GetDummyAtomPositionFromString(string dummy_atom_position)
{
    if(dummy_atom_position.compare("PositionAll") == 0)
        return kPositionAll;
    if(dummy_atom_position.compare("PositionBeg") == 0)
        return kPositionBeg;
    else
        return kPositionBeg;
}
DummyAtomOmission PrepFileResidue::GetDummyAtomOmissionFromString(string dummy_atom_omission)
{
    if(dummy_atom_omission.compare("Omit") == 0)
        return kOmit;
    if(dummy_atom_omission.compare("Nomit") == 0)
        return kNomit;
    else
        return kOmit;
}
SectionType PrepFileResidue::GetSectionTypeFromString(string section_type)
{
    if(section_type.compare("SectionLoop") == 0)
        return kSectionLoop;
    if(section_type.compare("SectionImproper") == 0)
        return kSectionImproper;
    if(section_type.compare("SectionDone") == 0)
        return kSectionDone;
    if(section_type.compare("SectionOther") == 0)
        return kSectionOther;
    else
        return kSectionOther;
}
PrepFileAtom* PrepFileResidue::GetPrepAtomByAtomName(string atom_name)
{
    for(PrepFileResidue::PrepFileAtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        PrepFileAtom* atom = (*it);
        if(atom->GetName().compare(atom_name) == 0)
            return atom;
    }
    return NULL;
}

PrepFileResidue::PrepFileAtomVector PrepFileResidue::GetAtomsParentVector()
{
    PrepFileAtomVector parents = PrepFileAtomVector();
    PrepFileAtomVector stack = PrepFileAtomVector();
    vector<int> neighbors = vector<int>();
    for(PrepFileAtomVector::iterator it = this->atoms_.begin(); it != this->atoms_.end(); it++)
    {
        parents.push_back(*it);
        neighbors.push_back(0);
    }
    for(Loop::iterator it = loops_.begin(); it != loops_.end(); it++)
    {
        int from = (*it).first;
        int to = (*it).second;

        neighbors.at(from - 1) = 1;
        neighbors.at(to - 1) = 1;
    }
    for(PrepFileAtomVector::iterator it = this->atoms_.begin(); it != this->atoms_.end(); it++)
    {
        PrepFileAtom* atom = *it;
        int index = distance(atoms_.begin(), it);
        if(stack.empty())
        {
            switch(atom->GetTopologicalType())
            {
                case kTopTypeM:
                case kTopType4:
                case kTopType3:
                    if(neighbors.at(index) == 0)
                        stack.push_back(atom);
                    stack.push_back(atom);
                    stack.push_back(atom);
                    stack.push_back(atom);
                    break;
                case kTopTypeB:
                    if(neighbors.at(index) == 0)
                        stack.push_back(atom);
                    stack.push_back(atom);
                    stack.push_back(atom);
                    break;
                case kTopTypeS:
                    if(neighbors.at(index) == 0)
                        stack.push_back(atom);
                    stack.push_back(atom);
                    break;
                case kTopTypeE:
                    if(neighbors.at(index) == 0)
                        stack.push_back(atom);
                    break;
            }
        }
        else
        {
            switch(atom->GetTopologicalType())
            {
                case kTopTypeM:
                case kTopType4:
                case kTopType3:
                    parents.at(index) = stack.at(stack.size() - 1);
                    stack.pop_back();
                    if(neighbors.at(index) == 0)
                        stack.push_back(atom);
                    stack.push_back(atom);
                    stack.push_back(atom);
                    break;
                case kTopTypeB:
                    parents.at(index) = stack.at(stack.size() - 1);
                    stack.pop_back();
                    stack.push_back(atom);
                    stack.push_back(atom);
                    break;
                case kTopTypeS:
                    parents.at(index) = stack.at(stack.size() - 1);
                    stack.pop_back();
                    stack.push_back(atom);
                    break;
                case kTopTypeE:
                    parents.at(index) = stack.at(stack.size() - 1);
                    stack.pop_back();
                    break;
            }
        }
    }
    return parents;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////

void PrepFileResidue::SetTitle(const string title){
    title_ = title;
}

void PrepFileResidue::SetName(const string name){
    name_ = name;
}

void PrepFileResidue::SetCoordinateType(CoordinateType coordinate_type){
    coordinate_type_ = coordinate_type;
}

void PrepFileResidue::SetOutputFormat(OutputFormat output_format){
    output_format_ = output_format;
}

void PrepFileResidue::SetGeometryType(GeometryType geometry_type){
    geometry_type_ = geometry_type;
}

void PrepFileResidue::SetDummyAtomOmission(DummyAtomOmission dummy_atom_omission){
    dummy_atom_omission_ = dummy_atom_omission;
}

void PrepFileResidue::SetDummyAtomType(const string dummy_atom_type){
    dummy_atom_type_ = dummy_atom_type;
}

void PrepFileResidue::SetDummyAtomPosition(DummyAtomPosition dummy_atom_position){
    dummy_atom_position_ = dummy_atom_position;
}

void PrepFileResidue::SetCharge(double charge){
    charge_ = charge;
}

void PrepFileResidue::SetAtoms(PrepFileAtomVector atoms){
    atoms_.clear();
    for(PrepFileAtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        atoms_.push_back(*it);
    }
}

void PrepFileResidue::AddAtom(PrepFileAtom* atom){
    atoms_.push_back(atom);
}

void PrepFileResidue::SetImproperDihedrals(DihedralVector improper_dihedrals){
    improper_dihedrals_.clear();
    for(DihedralVector::iterator it = improper_dihedrals.begin(); it != improper_dihedrals.end(); it++)
    {
        improper_dihedrals_.push_back(*it);
    }
}

void PrepFileResidue::AddImproperDihedral(Dihedral improper_dihedral){
    improper_dihedrals_.push_back(improper_dihedral);
}

void PrepFileResidue::SetLoops(Loop loops){
    loops_ = loops;
}


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
    if (Trim(line).find("STOP") != string::npos)           /// End of file
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
//                     FUNCTIONS                        //
//////////////////////////////////////////////////////////
 double PrepFileResidue::CalculatePrepResidueCharge()
 {
    PrepFileAtomVector atoms = GetAtoms();
    double residue_charge = 0.0;
    for(PrepFileAtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        PrepFileAtom* atom = (*it);
        residue_charge += atom->GetCharge();
    }
    return residue_charge;
 }

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void PrepFileResidue::Print(std::ostream& out)
{
    BondedAtomIndexMap bonded_atoms_map = this->GetBondingsOfResidue();
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
        << setw(10) << "Bonded"
        << endl;

    for(vector<PrepFileAtom*>::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        (*it)->Print(out);
        vector<int> bonded_atoms = bonded_atoms_map[(*it)->GetIndex()];
        out << "\t";
        for(unsigned int i = 0; i < bonded_atoms.size(); i++)
        {
            if(i != bonded_atoms.size() - 1)
                out << this->GetAtomNameByIndex(bonded_atoms.at(i)) << ", ";
            else
                out << this->GetAtomNameByIndex(bonded_atoms.at(i));
        }
        out << endl;

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
        out << setw(6) << this->GetAtomNameByIndex(it->first) << setw(6) << this->GetAtomNameByIndex(it->second) << endl;
    }

    out << endl;
}
