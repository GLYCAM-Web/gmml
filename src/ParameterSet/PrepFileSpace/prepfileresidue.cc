#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/MolecularModeling/residue.hpp"

using PrepFileSpace::PrepFileResidue;
using MolecularModeling::Residue;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
PrepFileResidue::PrepFileResidue() : title_(""), name_(""), coordinate_type_(PrepFileSpace::kINT), output_format_(PrepFileSpace::kFormatted), geometry_type_(PrepFileSpace::kGeometryCorrect),
    dummy_atom_omission_(PrepFileSpace::kOmit), dummy_atom_type_("DU"), dummy_atom_position_(PrepFileSpace::kPositionBeg), charge_(0.0)
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

std::string PrepFileResidue::GetTitle(){
    return title_;
}

std::string PrepFileResidue::GetName(){
    return name_;
}

PrepFileSpace::CoordinateType PrepFileResidue::GetCoordinateType(){
    return coordinate_type_;
}

PrepFileSpace::OutputFormat PrepFileResidue::GetOutputFormat(){
    return output_format_;
}

PrepFileSpace::GeometryType PrepFileResidue::GetGeometryType(){
    return geometry_type_;
}

PrepFileSpace::DummyAtomOmission PrepFileResidue::GetDummyAtomOmission(){
    return dummy_atom_omission_;
}

std::string PrepFileResidue::GetDummyAtomType(){
    return dummy_atom_type_;
}

PrepFileSpace::DummyAtomPosition PrepFileResidue::GetDummyAtomPosition(){
    return dummy_atom_position_;
}

double PrepFileResidue::GetCharge(){
    return charge_;
}

PrepFileSpace::PrepFileAtomVector PrepFileResidue::GetAtoms(){
    return atoms_;
}

PrepFileResidue::DihedralVector PrepFileResidue::GetImproperDihedrals(){
    return improper_dihedrals_;
}

PrepFileResidue::Loop PrepFileResidue::GetLoops(){
    return loops_;
}
PrepFileSpace::PrepFileAtom* PrepFileResidue::GetPrepAtomByName(std::string atom_name)
{
    for(PrepFileSpace::PrepFileAtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        PrepFileSpace::PrepFileAtom* prep_file_atom = (*it);
        std::string prep_file_atom_name = prep_file_atom->GetName();
        if(prep_file_atom_name.compare(atom_name) == 0)
            return prep_file_atom;
    }
    return NULL;
}
PrepFileResidue::BondedAtomIndexMap PrepFileResidue::GetBondingsOfResidue()
{
    BondedAtomIndexMap bonded_atoms_map = BondedAtomIndexMap();
//    std::vector<PrepFileSpace::PrepFileAtom*> stack = std::vector<PrepFileSpace::PrepFileAtom*>();
//    std::vector<int> number_of_bonds = std::vector<int>();

    PrepFileAtomVector parents = this->GetAtomsParentVector();
    for(Loop::iterator it = loops_.begin(); it != loops_.end(); it++)
    {
        int from = (*it).first;
        int to = (*it).second;

        bonded_atoms_map[from].push_back(to);
        bonded_atoms_map[to].push_back(from);
    }

    for(PrepFileSpace::PrepFileAtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        PrepFileSpace::PrepFileAtom* atom = *it;
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
    for(PrepFileSpace::PrepFileAtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        PrepFileSpace::PrepFileAtom* atom = (*it);
        if(stack.empty())
        {
            if(atom->GetTopologicalType() == gmml::kTopTypeM || atom->GetTopologicalType() == gmml::kTopTypeS
                    || atom->GetTopologicalType() == gmml::kTopTypeB || atom->GetTopologicalType() == gmml::kTopType3)
            {
                stack.push_back(atom);
                if(atom->GetTopologicalType() == gmml::kTopTypeM)
                    number_of_bonds.push_back(4);
                if(atom->GetTopologicalType() == gmml::kTopTypeS)
                    number_of_bonds.push_back(2);
                if(atom->GetTopologicalType() == gmml::kTopTypeB)
                    number_of_bonds.push_back(3);
                if(atom->GetTopologicalType() == gmml::kTopType3)
                    number_of_bonds.push_back(4);
            }
        }
        if(!stack.empty())
        {
            if(atom->GetTopologicalType() == gmml::kTopTypeE)
            {
                PrepFileSpace::PrepFileAtom* top_of_stack_atom = stack.at(stack.size()-1);
                bonded_atoms_map[top_of_stack_atom->GetIndex()].push_back(atom->GetIndex());
                bonded_atoms_map[atom->GetIndex()].push_back(top_of_stack_atom->GetIndex());
                number_of_bonds.at(number_of_bonds.size()-1)--;
                if(number_of_bonds.at(number_of_bonds.size()-1) == 0)
                {
                    number_of_bonds.pop_back();
                    stack.pop_back();
                }
            }
            if(atom->GetTopologicalType() == gmml::kTopTypeM)
            {
                PrepFileSpace::PrepFileAtom* top_of_stack_atom = stack.at(stack.size()-1);
                if(top_of_stack_atom->GetTopologicalType() == gmml::kTopTypeM)
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
                if(top_of_stack_atom->GetTopologicalType() == gmml::kTopTypeS || top_of_stack_atom->GetTopologicalType() == gmml::kTopTypeB
                        || top_of_stack_atom->GetTopologicalType() == gmml::kTopType3)
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
            if(atom->GetTopologicalType() == gmml::kTopTypeS)
            {
                PrepFileSpace::PrepFileAtom* top_of_stack_atom = stack.at(stack.size()-1);
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
            if(atom->GetTopologicalType() == gmml::kTopTypeB)
            {
                PrepFileSpace::PrepFileAtom* top_of_stack_atom = stack.at(stack.size()-1);
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
            if(atom->GetTopologicalType() == gmml::kTopType3)
            {
                PrepFileSpace::PrepFileAtom* top_of_stack_atom = stack.at(stack.size()-1);
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
std::string PrepFileResidue::GetAtomNameByIndex(int atom_index)
{
    for(PrepFileSpace::PrepFileAtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        PrepFileSpace::PrepFileAtom* prep_file_atom = (*it);
        if(prep_file_atom->GetIndex() == atom_index)
            return prep_file_atom->GetName();
    }
    return NULL;
}
std::string PrepFileResidue::GetStringFormatOfCoordinateType(PrepFileSpace::CoordinateType coordinate_type)
{
    switch(coordinate_type)
    {
        case PrepFileSpace::kINT:
            return "INT";
        case PrepFileSpace::kXYZ:
            return "XYZ";
        default:
            return "";
    }
}
std::string PrepFileResidue::GetStringFormatOfCoordinateType()
{
    switch(coordinate_type_)
    {
        case PrepFileSpace::kINT:
            return "INT";
        case PrepFileSpace::kXYZ:
            return "XYZ";
        default:
            return "";
    }
}
std::string PrepFileResidue::GetStringFormatOfOutputFormat(PrepFileSpace::OutputFormat output_format)
{
    switch(output_format)
    {
        case PrepFileSpace::kFormatted:
            return "FRM";
        case PrepFileSpace::kBinary:
            return "BIN";
        default:
            return "";
    }
}
std::string PrepFileResidue::GetStringFormatOfOutputFormat()
{
    switch(output_format_)
    {
        case PrepFileSpace::kFormatted:
            return "FRM";
        case PrepFileSpace::kBinary:
            return "BIN";
        default:
            return "";
    }
}
std::string PrepFileResidue::GetStringFormatOfGeometryType(PrepFileSpace::GeometryType geometry_type)
{
    switch(geometry_type)
    {
        case PrepFileSpace::kGeometryCorrect:
            return "CORRECT";
        case PrepFileSpace::kGeometryChange:
            return "CHANGE";
        default:
            return "";
    }
}
std::string PrepFileResidue::GetStringFormatOfGeometryType()
{
    switch(geometry_type_)
    {
        case PrepFileSpace::kGeometryCorrect:
            return "CORRECT";
        case PrepFileSpace::kGeometryChange:
            return "CHANGE";
        default:
            return "";
    }
}
std::string PrepFileResidue::GetStringFormatOfDummyAtomPosition(DummyAtomPosition dummy_atom_position)
{
    switch(dummy_atom_position)
    {
        case PrepFileSpace::kPositionAll:
            return "ALL";
        case PrepFileSpace::kPositionBeg:
            return "BEG";
        default:
            return "";
    }
}
std::string PrepFileResidue::GetStringFormatOfDummyAtomPosition()
{
    switch(dummy_atom_position_)
    {
        case PrepFileSpace::kPositionAll:
            return "ALL";
        case PrepFileSpace::kPositionBeg:
            return "BEG";
        default:
            return "";
    }
}
std::string PrepFileResidue::GetStringFormatOfDummyAtomOmission(PrepFileSpace::DummyAtomOmission dummy_atom_omission)
{
    switch(dummy_atom_omission)
    {
        case PrepFileSpace::kOmit:
            return "OMIT";
        case PrepFileSpace::kNomit:
            return "NOMIT";
        default:
            return "";
    }
}
std::string PrepFileResidue::GetStringFormatOfDummyAtomOmission()
{
    switch(dummy_atom_omission_)
    {
        case PrepFileSpace::kOmit:
            return "OMIT";
        case PrepFileSpace::kNomit:
            return "NOMIT";
        default:
            return "";
    }
}
std::string PrepFileResidue::GetStringFormatOfSectionType(SectionType section_type)
{
    switch(section_type)
    {
        case PrepFileSpace::kSectionLoop:
            return "SectionLoop";
        case PrepFileSpace::kSectionImproper:
            return "SectionImproper";
        case PrepFileSpace::kSectionDone:
            return "PrepFileSpace::kSectionDone";
        case PrepFileSpace::kSectionOther:
            return "SectionOther";
        default:
            return "";
    }
}
PrepFileSpace::CoordinateType PrepFileResidue::GetCoordinateTypeFromString(std::string coordinate_type)
{
    if(coordinate_type.compare("INT") == 0)
        return PrepFileSpace::kINT;
    if(coordinate_type.compare("XYZ") == 0)
        return PrepFileSpace::kXYZ;
    else
        return PrepFileSpace::kINT;
}
PrepFileSpace::OutputFormat PrepFileResidue::GetOutputFormatFromString(std::string output_format)
{
    if(output_format.compare("Formatted") == 0)
        return PrepFileSpace::kFormatted;
    if(output_format.compare("Binary") == 0)
        return PrepFileSpace::kBinary;
    else
        return PrepFileSpace::kBinary;
}
PrepFileSpace::GeometryType PrepFileResidue::GetGeometryTypeFromString(std::string geometry_type)
{
    if(geometry_type.compare("GeometryCorrect") == 0)
        return PrepFileSpace::kGeometryCorrect;
    if(geometry_type.compare("GeometryChange") == 0)
        return PrepFileSpace::kGeometryChange;
    else
        return PrepFileSpace::kGeometryCorrect;
}
PrepFileSpace::DummyAtomPosition PrepFileResidue::GetDummyAtomPositionFromString(std::string dummy_atom_position)
{
    if(dummy_atom_position.compare("PositionAll") == 0)
        return PrepFileSpace::kPositionAll;
    if(dummy_atom_position.compare("PositionBeg") == 0)
        return PrepFileSpace::kPositionBeg;
    else
        return PrepFileSpace::kPositionBeg;
}
PrepFileSpace::DummyAtomOmission PrepFileResidue::GetDummyAtomOmissionFromString(std::string dummy_atom_omission)
{
    if(dummy_atom_omission.compare("Omit") == 0)
        return PrepFileSpace::kOmit;
    if(dummy_atom_omission.compare("Nomit") == 0)
        return PrepFileSpace::kNomit;
    else
        return PrepFileSpace::kOmit;
}
PrepFileSpace::SectionType PrepFileResidue::GetSectionTypeFromString(std::string section_type)
{
    if(section_type.compare("SectionLoop") == 0)
        return PrepFileSpace::kSectionLoop;
    if(section_type.compare("SectionImproper") == 0)
        return PrepFileSpace::kSectionImproper;
    if(section_type.compare("SectionDone") == 0)
        return PrepFileSpace::kSectionDone;
    if(section_type.compare("SectionOther") == 0)
        return PrepFileSpace::kSectionOther;
    else
        return PrepFileSpace::kSectionOther;
}
PrepFileSpace::PrepFileAtom* PrepFileResidue::GetPrepAtomByAtomName(std::string atom_name)
{
    for(PrepFileSpace::PrepFileAtomVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        PrepFileSpace::PrepFileAtom* atom = (*it);
        if(atom->GetName().compare(atom_name) == 0)
            return atom;
    }
    return NULL;
}

PrepFileSpace::PrepFileAtomVector PrepFileResidue::GetAtomsParentVector()
{
    PrepFileAtomVector parents = PrepFileAtomVector();
    PrepFileAtomVector stack = PrepFileAtomVector();
    std::vector<int> neighbors = std::vector<int>();
    for(PrepFileSpace::PrepFileAtomVector::iterator it = this->atoms_.begin(); it != this->atoms_.end(); it++)
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
    for(PrepFileSpace::PrepFileAtomVector::iterator it = this->atoms_.begin(); it != this->atoms_.end(); it++)
    {
        PrepFileSpace::PrepFileAtom* atom = *it;
        int index = distance(atoms_.begin(), it);
        if(stack.empty())
        {
            switch(atom->GetTopologicalType())
            {
                case gmml::kTopTypeM:
                case gmml::kTopType4:
                case gmml::kTopType3:
                    if(neighbors.at(index) == 0)
                        stack.push_back(atom);
                    stack.push_back(atom);
                    stack.push_back(atom);
                    stack.push_back(atom);
                    break;
                case gmml::kTopTypeB:
                    if(neighbors.at(index) == 0)
                        stack.push_back(atom);
                    stack.push_back(atom);
                    stack.push_back(atom);
                    break;
                case gmml::kTopTypeS:
                    if(neighbors.at(index) == 0)
                        stack.push_back(atom);
                    stack.push_back(atom);
                    break;
                case gmml::kTopTypeE:
                    if(neighbors.at(index) == 0)
                        stack.push_back(atom);
                    break;
            }
        }
        else
        {
            switch(atom->GetTopologicalType())
            {
                case gmml::kTopTypeM:
                case gmml::kTopType4:
                case gmml::kTopType3:
                    parents.at(index) = stack.at(stack.size() - 1);
                    stack.pop_back();
                    if(neighbors.at(index) == 0)
                        stack.push_back(atom);
                    stack.push_back(atom);
                    stack.push_back(atom);
                    break;
                case gmml::kTopTypeB:
                    parents.at(index) = stack.at(stack.size() - 1);
                    stack.pop_back();
                    stack.push_back(atom);
                    stack.push_back(atom);
                    break;
                case gmml::kTopTypeS:
                    parents.at(index) = stack.at(stack.size() - 1);
                    stack.pop_back();
                    stack.push_back(atom);
                    break;
                case gmml::kTopTypeE:
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

void PrepFileResidue::SetTitle(const std::string title){
    title_ = title;
}

void PrepFileResidue::SetName(const std::string name){
    name_ = name;
}

void PrepFileResidue::SetCoordinateType(PrepFileSpace::CoordinateType coordinate_type){
    coordinate_type_ = coordinate_type;
}

void PrepFileResidue::SetOutputFormat(PrepFileSpace::OutputFormat output_format){
    output_format_ = output_format;
}

void PrepFileResidue::SetGeometryType(PrepFileSpace::GeometryType geometry_type){
    geometry_type_ = geometry_type;
}

void PrepFileResidue::SetDummyAtomOmission(PrepFileSpace::DummyAtomOmission dummy_atom_omission){
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

void PrepFileResidue::SetAtoms(PrepFileAtomVector atoms){
    atoms_.clear();
    for(PrepFileSpace::PrepFileAtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        atoms_.push_back(*it);
    }
}

void PrepFileResidue::AddAtom(PrepFileSpace::PrepFileAtom* atom){
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
    std::string line, name, dummy_atom_type;
    std::istringstream ss;
    PrepFileSpace::CoordinateType coordinate_type;
    PrepFileSpace::OutputFormat output_format;
    PrepFileSpace::GeometryType geometry_type;
    PrepFileSpace::DummyAtomOmission dummy_atom_omission;
    DummyAtomPosition dummy_atom_position;
    PrepFileResidue *residue = new PrepFileResidue();

    getline(in_file, line);             /// Read the first line of a residue section
    if (gmml::Trim(line).find("STOP") != std::string::npos)           /// End of file
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
    residue->charge_ = gmml::ConvertString<double>(line);

    /// Process atoms of the residue
    while (getline(in_file, line) && !gmml::Trim(line).empty())
    {
        PrepFileSpace::PrepFileAtom *atom = new PrepFileSpace::PrepFileAtom(line);
        residue->atoms_.push_back(atom);
    }

    /// Process the extra sections: IMPROPER, LOOP, DONE
    bool done = false;
    while (!done)
    {
        /// Skip blank lines until to reach to a known section title
        getline(in_file, line);
        while (gmml::Trim(line).empty())
        {
            getline(in_file, line);
        }
        /// Does a corresponding action based on the section title
        switch (ExtractSectionType(line))
        {
            case PrepFileSpace::kSectionLoop:
                residue->loops_ = residue->ExtractLoops(in_file);
                break;
            case PrepFileSpace::kSectionImproper:
                residue->improper_dihedrals_ = residue->ExtractImproperDihedral(in_file);
                break;
            case PrepFileSpace::kSectionDone:
                done = true;
                break;
            case PrepFileSpace::kSectionOther:
//                std::cout << "Unrecognized section in prep file";
                gmml::log(__LINE__, __FILE__,  gmml::WAR, "Unrecognized section in prep file" );
                break;
        }
    }
    return residue;
}

/// Return residue name from a stream line which is the first column of the 3rd line in each residue section
std::string PrepFileResidue::ExtractResidueName(std::istream& ss)
{
    std::string name;
    ss >> name;
    return name;
}

/// Return coordinate type from a stream line which is the 2nd column of the 3rd line in each residue section
PrepFileSpace::CoordinateType PrepFileResidue::ExtractResidueCoordinateType(std::istream &ss)
{
    std::string s;
    ss >> s;
    if (s == "XYZ")
        return PrepFileSpace::kXYZ;
    else
        return PrepFileSpace::kINT;
}

/// Return output format from a stream line which is the 3rd column of the 3rd line in each residue section
PrepFileSpace::OutputFormat PrepFileResidue::ExtractResidueOutputFormat(std::istream& ss)
{
    int val;
    ss >> val;
    if (val == 1)
        return PrepFileSpace::kBinary;
    else
        return PrepFileSpace::kFormatted;
}

/// Return geometry type from a stream line which is the first column of the 4th line in each residue section
PrepFileSpace::GeometryType PrepFileResidue::ExtractResidueGeometryType(std::istream& ss)
{
    std::string s;
    ss >> s;
    if (s == "CHANGE")
        return PrepFileSpace::kGeometryChange;
    else
        return PrepFileSpace::kGeometryCorrect;
}

/// Return dummy atom omission from a stream line which is the 2nd column of the 4th line in each residue section
PrepFileSpace::DummyAtomOmission PrepFileResidue::ExtractResidueDummyAtomOmission(std::istream &ss)
{
    std::string s;
    ss >> s;
    if (s == "NOMIT")
        return PrepFileSpace::kNomit;
    else
        return PrepFileSpace::kOmit;
}

/// Return dummy atom position from a stream line which is the 4th column of the 4th line in each residue section
PrepFileSpace::DummyAtomPosition PrepFileResidue::ExtractResidueDummyAtomPosition(std::istream &ss)
{
    std::string s;
    ss >> s;
    if (s == "ALL")
        return PrepFileSpace::kPositionAll;
    else
        return PrepFileSpace::kPositionBeg;
}

/// Return a corresponding title from a stream line which may appear in each residue section
PrepFileSpace::SectionType PrepFileResidue::ExtractSectionType(std::string &line)
{
    if (line == "LOOP")
        return PrepFileSpace::kSectionLoop;
    else if (line == "IMPROPER")
        return PrepFileSpace::kSectionImproper;
    else if (line == "DONE")
        return PrepFileSpace::kSectionDone;
    return PrepFileSpace::kSectionOther;
}

/// Parse the loop section of each residue section and return a loop map
PrepFileResidue::Loop PrepFileResidue::ExtractLoops(std::ifstream &in_file)
{
    Loop loops;
    std::string line;
    std::stringstream ss;

    getline(in_file, line);
    while (!gmml::Trim(line).empty())                         /// Read file until blank line which determines the end of the section
    {
        ss.clear();
        ss.str(line);                                   /// Create a stream from the read line
        std::string atom_names[2];
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

/// Parse the improper dihedral section of each residue section and return a std::vector of improper dihedrals
std::vector<PrepFileResidue::Dihedral> PrepFileResidue::ExtractImproperDihedral(std::ifstream &in_file)
{
    std::string line;
    std::stringstream ss;
    std::vector<Dihedral> dihedrals;
    getline(in_file, line);
    while (!gmml::Trim(line).empty())                         /// Read file until blank line which determines the end of the section
    {
        std::string atom_names[4];
        Dihedral dihedral;
        ss.clear();
        ss.str(line);
        /// Extract improper atom types involving in a dihedral from each line
        ss >> atom_names[0]
           >> atom_names[1]
           >> atom_names[2]
           >> atom_names[3];
        /// Push all atoms into a std::vector of atom types
        for(int i = 0; i < 4; i++)
        {
            dihedral.push_back(atom_names[i]);
        }
        dihedrals.push_back(dihedral);                  /// Create a new dihedral into the std::vector of dihedrals
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
    for(PrepFileSpace::PrepFileAtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        PrepFileSpace::PrepFileAtom* atom = (*it);
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
    out << "Title: " << title_ << std::endl;
    out << std::setw(10) << "ResName"
        << std::setw(10) << "CrdType"
        << std::setw(10) << "Output"
        << std::setw(10) << "GeoType"
        << std::setw(15) << "DummyOmission"
        << std::setw(12) << "DummyType"
        << std::setw(12) << "DummyPos"
        << std::setw(10) << "Charge" << std::endl;
    out << std::setw(10) << name_;

    if(coordinate_type_ == PrepFileSpace::kINT)
        out << std::setw(10) << "INT";
    else if(coordinate_type_ == PrepFileSpace::kXYZ)
        out << std::setw(10) << "XYZ";
    else
        out << std::setw(10) << "--";

    if(output_format_ == PrepFileSpace::kBinary)
        out << std::setw(10) << "Binary";
    else if(output_format_ == PrepFileSpace::kFormatted)
        out << std::setw(10) << "NBinary";
    else
        out << std::setw(10) << "--";

    if(geometry_type_ == PrepFileSpace::kGeometryCorrect)
        out << std::setw(10) << "Correct";
    else if(geometry_type_ == PrepFileSpace::kGeometryChange)
        out << std::setw(10) << "Change";
    else
        out << std::setw(10) << "--";

    if(dummy_atom_omission_ == PrepFileSpace::kOmit)
        out << std::setw(15) << "YES";
    else if(dummy_atom_omission_ == PrepFileSpace::kNomit)
        out << std::setw(15) << "NO";
    else
        out << std::setw(15) << "--";

    out << std::setw(12) << dummy_atom_type_;

    if(dummy_atom_position_ == PrepFileSpace::kPositionAll)
        out << std::setw(12) << "ALL";
    else if (dummy_atom_position_ == PrepFileSpace::kPositionBeg)
        out << std::setw(12) << "BEG";
    else
        out << std::setw(12) << "--";

    out << std::setw(10) << charge_ << std::endl << std::endl;

    out << std::setw(3) << "#"
        << std::setw(6) << "Name"
        << std::setw(6) << "Type"
        << std::setw(3) << "TT"
        << std::setw(4) << "B#"
        << std::setw(4) << "A#"
        << std::setw(4) << "D#"
        << std::setw(10) << "Bond"
        << std::setw(10) << "Angle"
        << std::setw(10) << "Dihedral"
        << std::setw(10) << "Charge"
        << std::setw(10) << "Bonded"
        << std::endl;

    for(std::vector<PrepFileSpace::PrepFileAtom*>::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        (*it)->Print(out);
        std::vector<int> bonded_atoms = bonded_atoms_map[(*it)->GetIndex()];
        out << "\t";
        for(unsigned int i = 0; i < bonded_atoms.size(); i++)
        {
            if(i != bonded_atoms.size() - 1)
                out << this->GetAtomNameByIndex(bonded_atoms.at(i)) << ", ";
            else
                out << this->GetAtomNameByIndex(bonded_atoms.at(i));
        }
        out << std::endl;

    }

    out << std::endl << "Improper dihedrals" << std::endl;
    for(std::vector<Dihedral>::iterator it = improper_dihedrals_.begin(); it != improper_dihedrals_.end(); it++)
    {
        for(Dihedral::iterator it1 = it->begin(); it1 != it->end(); it1++)
        {
            out << std::setw(6) << (*it1);
        }
        out << std::endl;
    }

    out << std::endl << "Loops" << std::endl;
    for(Loop::iterator it = loops_.begin(); it != loops_.end(); it++)
    {
        out << std::setw(6) << this->GetAtomNameByIndex(it->first) << std::setw(6) << this->GetAtomNameByIndex(it->second) << std::endl;
    }

    out << std::endl;
}
