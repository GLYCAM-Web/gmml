#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ios>

#include "includes/ParameterSet/PrepFile/prepResidue.hpp"
#include "includes/ParameterSet/PrepFile/prepAtom.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/utils.hpp" //Trim

using prep::PrepResidue;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////

PrepResidue::PrepResidue(std::ifstream& in_file, std::string &line)
{
    std::string name, dummy_atom_type;
    std::istringstream ss;
    this->SetTitle(line);             /// Set title of the residue
    getline(in_file, line);             /// Blank line, skip
    getline(in_file, line);             /// Read the next line
    // I guess you can extract with spaces as the default delimiter with >> from a stringstream
    ss.str(line);                       /// Create an stream from the read line
    this->ExtractResidueName(ss);
    this->ExtractResidueCoordinateType(ss);
    this->ExtractResidueOutputFormat(ss);
    ss.clear();
    getline(in_file, line);             /// Read the next line
    this->ExtractResidueGeometryType(ss);
    this->ExtractResidueDummyAtomOmission(ss);
    ss >> dummy_atom_type;
    this->SetDummyAtomType(dummy_atom_type);
    this->ExtractResidueDummyAtomPosition(ss);
    getline(in_file, line);             /// Read the next line
    this->SetCharge(gmml::ConvertString<double>(line));
    /// Process atoms of the residue
    while (getline(in_file, line) && !gmml::Trim(line).empty())
    {
        this->addAtom(std::make_unique<PrepAtom>(line));
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
            case prep::kSectionLoop:
                this->ExtractLoops(in_file);
                break;
            case prep::kSectionImproper:
                this->ExtractImproperDihedral(in_file);
                break;
            case prep::kSectionDone:
                done = true;
                break;
            case prep::kSectionOther:
                gmml::log(__LINE__, __FILE__,  gmml::WAR, "Unrecognized section in prep file" );
                break;
        }
    }
}
//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
std::string PrepResidue::GetTitle() const
{
    return title_;
}

prep::CoordinateType PrepResidue::GetCoordinateType() const
{
    return coordinate_type_;
}

prep::OutputFormat PrepResidue::GetOutputFormat() const
{
    return output_format_;
}

prep::GeometryType PrepResidue::GetGeometryType() const
{
    return geometry_type_;
}

prep::DummyAtomOmission PrepResidue::GetDummyAtomOmission() const
{
    return dummy_atom_omission_;
}

std::string PrepResidue::GetDummyAtomType() const
{
    return dummy_atom_type_;
}

prep::DummyAtomPosition PrepResidue::GetDummyAtomPosition() const
{
    return dummy_atom_position_;
}

double PrepResidue::GetCharge() const
{
    return charge_;
}

PrepResidue::DihedralVector PrepResidue::GetImproperDihedrals() const
{
    return improper_dihedrals_;
}

std::vector<std::pair<std::string, std::string>> PrepResidue::GetLoops() const
{
	return loops_;
}

std::string PrepResidue::GetStringFormatOfCoordinateType(prep::CoordinateType coordinate_type) const
{
    switch(coordinate_type)
    {
        case prep::kINT:
            return "INT";
        case prep::kXYZ:
            return "XYZ";
        default:
            return "";
    }
}

std::string PrepResidue::GetStringFormatOfOutputFormat(prep::OutputFormat output_format) const
{
    switch(output_format)
    {
        case prep::kFormatted:
            return "FRM";
        case prep::kBinary:
            return "BIN";
        default:
            return "";
    }
}

std::string PrepResidue::GetStringFormatOfGeometryType(prep::GeometryType geometry_type) const
{
    switch(geometry_type)
    {
        case prep::kGeometryCorrect:
            return "CORRECT";
        case prep::kGeometryChange:
            return "CHANGE";
        default:
            return "";
    }
}

std::string PrepResidue::GetStringFormatOfDummyAtomPosition(DummyAtomPosition dummy_atom_position) const
{
    switch(dummy_atom_position)
    {
        case prep::kPositionAll:
            return "ALL";
        case prep::kPositionBeg:
            return "BEG";
        default:
            return "";
    }
}

std::string PrepResidue::GetStringFormatOfDummyAtomOmission(prep::DummyAtomOmission dummy_atom_omission) const
{
    switch(dummy_atom_omission)
    {
        case prep::kOmit:
            return "OMIT";
        case prep::kNomit:
            return "NOMIT";
        default:
            return "";
    }
}

std::string PrepResidue::GetStringFormatOfSectionType(SectionType section_type) const
{
    switch(section_type)
    {
        case prep::kSectionLoop:
            return "SectionLoop";
        case prep::kSectionImproper:
            return "SectionImproper";
        case prep::kSectionDone:
            return "prep::kSectionDone";
        case prep::kSectionOther:
            return "SectionOther";
        default:
            return "";
    }
}

prep::CoordinateType PrepResidue::GetCoordinateTypeFromString(std::string coordinate_type) const
{
    if(coordinate_type.compare("INT") == 0)
        return prep::kINT;
    if(coordinate_type.compare("XYZ") == 0)
        return prep::kXYZ;
    else
        return prep::kINT;
}

prep::OutputFormat PrepResidue::GetOutputFormatFromString(std::string output_format) const
{
    if(output_format.compare("Formatted") == 0)
        return prep::kFormatted;
    if(output_format.compare("Binary") == 0)
        return prep::kBinary;
    else
        return prep::kBinary;
}

prep::GeometryType PrepResidue::GetGeometryTypeFromString(std::string geometry_type) const
{
    if(geometry_type.compare("GeometryCorrect") == 0)
        return prep::kGeometryCorrect;
    if(geometry_type.compare("GeometryChange") == 0)
        return prep::kGeometryChange;
    else
        return prep::kGeometryCorrect;
}

prep::DummyAtomPosition PrepResidue::GetDummyAtomPositionFromString(std::string dummy_atom_position) const
{
    if(dummy_atom_position.compare("PositionAll") == 0)
        return prep::kPositionAll;
    if(dummy_atom_position.compare("PositionBeg") == 0)
        return prep::kPositionBeg;
    else
        return prep::kPositionBeg;
}

prep::DummyAtomOmission PrepResidue::GetDummyAtomOmissionFromString(std::string dummy_atom_omission) const
{
    if(dummy_atom_omission.compare("Omit") == 0)
        return prep::kOmit;
    if(dummy_atom_omission.compare("Nomit") == 0)
        return prep::kNomit;
    else
        return prep::kOmit;
}

prep::SectionType PrepResidue::GetSectionTypeFromString(std::string section_type) const
{
    if(section_type.compare("SectionLoop") == 0)
        return prep::kSectionLoop;
    if(section_type.compare("SectionImproper") == 0)
        return prep::kSectionImproper;
    if(section_type.compare("SectionDone") == 0)
        return prep::kSectionDone;
    if(section_type.compare("SectionOther") == 0)
        return prep::kSectionOther;
    else
        return prep::kSectionOther;
}

std::vector<std::string> PrepResidue::GetAtomNames() const
{
    std::vector<std::string> foundAtoms;
    for (auto &prepAtom : this->getAtoms())
    {
        foundAtoms.push_back(prepAtom->getName());
    }
    return foundAtoms;
}

std::vector<std::string> PrepResidue::GetHeavyAtomNames() const
{
    std::vector<std::string> foundAtoms;
    for (auto &prepAtom : this->getAtoms())
    {
        if(prepAtom->getName().at(0) != 'H' && prepAtom->getName() != "DUMM")
        {
            foundAtoms.push_back(prepAtom->getName());
        }
    }
    return foundAtoms;
}
//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////
void PrepResidue::SetTitle(const std::string title){
    title_ = title;
}

void PrepResidue::SetCoordinateType(prep::CoordinateType coordinate_type){
    coordinate_type_ = coordinate_type;
}

void PrepResidue::SetOutputFormat(prep::OutputFormat output_format){
    output_format_ = output_format;
}

void PrepResidue::SetGeometryType(prep::GeometryType geometry_type){
    geometry_type_ = geometry_type;
}

void PrepResidue::SetDummyAtomOmission(prep::DummyAtomOmission dummy_atom_omission){
    dummy_atom_omission_ = dummy_atom_omission;
}

void PrepResidue::SetDummyAtomType(const std::string dummy_atom_type){
    dummy_atom_type_ = dummy_atom_type;
}

void PrepResidue::SetDummyAtomPosition(DummyAtomPosition dummy_atom_position){
    dummy_atom_position_ = dummy_atom_position;
}

void PrepResidue::SetCharge(double charge){
    charge_ = charge;
}

void PrepResidue::SetImproperDihedrals(DihedralVector improper_dihedrals){
    improper_dihedrals_.clear();
    for(DihedralVector::iterator it = improper_dihedrals.begin(); it != improper_dihedrals.end(); it++)
    {
        improper_dihedrals_.push_back(*it);
    }
}

void PrepResidue::AddImproperDihedral(Dihedral improper_dihedral){
    improper_dihedrals_.push_back(improper_dihedral);
}

void PrepResidue::AddLoop(std::pair<std::string, std::string> loop)
{
	loops_.push_back(loop);
	return;
}
//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////

// can recursively look up incoming connections to find bond, angle and dihedral atoms. Just skip dummies.
//currentAtom->DeterminePosition(); // could do this once connectivity is set

void PrepResidue::SetConnectivities()
{
	//std::cout << "Attempting to set connectivities in prepResidue: " << this->getName() << std::endl;
	//std::cout << "Number of atoms: " << this->getAtoms().size() << std::endl;
	//std::cout << "First atom is " << this->getAtoms().front()->getName() << std::endl;
	std::vector<PrepAtom*> connectionPointStack;
	connectionPointStack.push_back(this->getAtoms().front());
	//while(currentAtom != this->getAtoms().end())
	for(auto &currentAtom : this->getAtoms())
	{
		connectionPointStack.back()->addBond(currentAtom);
		//std::cout << "Bonded " << connectionPointStack.back()->getName() << " to " << currentAtom->getName() << std::endl;;
		connectionPointStack.back()->visit();
		if (connectionPointStack.back()->GetVisits() >= connectionPointStack.back()->GetTopologicalType())
		{
			connectionPointStack.pop_back();
		}
		if(currentAtom->GetTopologicalType() > kTopTypeE)
		{
			connectionPointStack.push_back(currentAtom);
		}
		++currentAtom;
	}
}

/// Return residue name from a stream line which is the first column of the 3rd line in each residue section
void PrepResidue::ExtractResidueName(std::istream& ss)
{
    std::string name;
    ss >> name;
    this->setName(name);
    return;
}

/// Return coordinate type from a stream line which is the 2nd column of the 3rd line in each residue section
void PrepResidue::ExtractResidueCoordinateType(std::istream &ss)
{
    std::string s;
    ss >> s;
    if (s == "XYZ")
    {
    	this->SetCoordinateType(prep::kXYZ);
    }
    else
    {
    	this->SetCoordinateType(prep::kINT);
    }
    return;
}

/// Return output format from a stream line which is the 3rd column of the 3rd line in each residue section
void PrepResidue::ExtractResidueOutputFormat(std::istream& ss)
{
    int val;
    ss >> val;
    if (val == 1)
    	this->SetOutputFormat(prep::kBinary);
    else
    	this->SetOutputFormat(prep::kFormatted);
    return;
}

/// Return geometry type from a stream line which is the first column of the 4th line in each residue section
void PrepResidue::ExtractResidueGeometryType(std::istream& ss)
{
    std::string s;
    ss >> s;
    if (s == "CHANGE")
    	this->SetGeometryType(prep::kGeometryChange);
    else
    	this->SetGeometryType(prep::kGeometryCorrect);
    return;
}

/// Return dummy atom omission from a stream line which is the 2nd column of the 4th line in each residue section
void PrepResidue::ExtractResidueDummyAtomOmission(std::istream &ss)
{
    std::string s;
    ss >> s;
    if (s == "NOMIT")
    	this->SetDummyAtomOmission(prep::kNomit);
    else
    	this->SetDummyAtomOmission(prep::kOmit);
    return;
}

/// Return dummy atom position from a stream line which is the 4th column of the 4th line in each residue section
void PrepResidue::ExtractResidueDummyAtomPosition(std::istream &ss)
{
    std::string s;
    ss >> s;
    if (s == "ALL")
    	this->SetDummyAtomPosition(prep::kPositionAll);
    else
        this->SetDummyAtomPosition(prep::kPositionBeg);
    return;
}

/// Return a corresponding title from a stream line which may appear in each residue section
prep::SectionType PrepResidue::ExtractSectionType(std::string &line)
{
    if (line == "LOOP")
        return prep::kSectionLoop;
    else if (line == "IMPROPER")
        return prep::kSectionImproper;
    else if (line == "DONE")
        return prep::kSectionDone;
    return prep::kSectionOther;
}

/// Parse the loop section of each residue section and return a loop map
void PrepResidue::ExtractLoops(std::ifstream &in_file)
{
    std::string line;
    std::stringstream ss;
    getline(in_file, line);
    while (!gmml::Trim(line).empty())                         /// Read file until blank line which determines the end of the section
    {
        ss.clear();
        ss.str(line);                                   /// Create a stream from the read line
        std::string atom_names[2];
        ss >> atom_names[0] >> atom_names[1];           /// Extract atom names from the stream
        this->AddLoop(std::make_pair(atom_names[0], atom_names[1])); /// Add a new entry into the loop map of the residue
        getline(in_file, line);
    }
    return;
}

/// Parse the improper dihedral section of each residue section and return a std::vector of improper dihedrals
void PrepResidue::ExtractImproperDihedral(std::ifstream &in_file)
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
    this->SetImproperDihedrals(dihedrals);
    return;
}

//////////////////////////////////////////////////////////
//                     FUNCTIONS                        //
//////////////////////////////////////////////////////////
 double PrepResidue::CalculatePrepResidueCharge()
 {
    double residue_charge = 0.0;
    for(auto &atom : this->getAtoms())
    {
        residue_charge += atom->GetCharge();
    }
    return residue_charge;
 }

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void PrepResidue::Print(std::ostream& out)
{
//    BondedAtomIndexMap bonded_atoms_map = this->GetBondingsOfResidue();
    out << "Title: " << title_ << std::endl;
    out << std::setw(10) << "ResName"
        << std::setw(10) << "CrdType"
        << std::setw(10) << "Output"
        << std::setw(10) << "GeoType"
        << std::setw(15) << "DummyOmission"
        << std::setw(12) << "DummyType"
        << std::setw(12) << "DummyPos"
        << std::setw(10) << "Charge" << std::endl;
    out << std::setw(10) << this->getName();

    if(coordinate_type_ == prep::kINT)
        out << std::setw(10) << "INT";
    else if(coordinate_type_ == prep::kXYZ)
        out << std::setw(10) << "XYZ";
    else
        out << std::setw(10) << "--";

    if(output_format_ == prep::kBinary)
        out << std::setw(10) << "Binary";
    else if(output_format_ == prep::kFormatted)
        out << std::setw(10) << "NBinary";
    else
        out << std::setw(10) << "--";

    if(geometry_type_ == prep::kGeometryCorrect)
        out << std::setw(10) << "Correct";
    else if(geometry_type_ == prep::kGeometryChange)
        out << std::setw(10) << "Change";
    else
        out << std::setw(10) << "--";

    if(dummy_atom_omission_ == prep::kOmit)
        out << std::setw(15) << "YES";
    else if(dummy_atom_omission_ == prep::kNomit)
        out << std::setw(15) << "NO";
    else
        out << std::setw(15) << "--";

    out << std::setw(12) << dummy_atom_type_;

    if(dummy_atom_position_ == prep::kPositionAll)
        out << std::setw(12) << "ALL";
    else if (dummy_atom_position_ == prep::kPositionBeg)
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

//    for(std::vector<prep::PrepAtom*>::iterator it = atoms_.begin(); it != atoms_.end(); it++)
//    {
//        (*it)->Print(out);
//        std::vector<int> bonded_atoms = bonded_atoms_map[(*it)->GetIndex()];
//        out << "\t";
//        for(unsigned int i = 0; i < bonded_atoms.size(); i++)
//        {
//            if(i != bonded_atoms.size() - 1)
//                out << this->GetAtomNameByIndex(bonded_atoms.at(i)) << ", ";
//            else
//                out << this->GetAtomNameByIndex(bonded_atoms.at(i));
//        }
//        out << std::endl;
//
//    }

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
//    for(Loop::iterator it = loops_.begin(); it != loops_.end(); it++)
//    {
//        out << std::setw(6) << this->GetAtomNameByIndex(it->first) << std::setw(6) << this->GetAtomNameByIndex(it->second) << std::endl;
//    }

    out << std::endl;
}

void PrepResidue::Write(std::ostream &stream)
{
	stream << this->GetTitle() << std::endl << std::endl
			<< std::left << std::setw(4) << this->getName() << " " << std::right << std::setw(3) << this->GetStringFormatOfCoordinateType(this->GetCoordinateType()) << " "
			<< std::setw(1) << this->GetOutputFormat() << std::endl
			<< this->GetStringFormatOfGeometryType(this->GetGeometryType()) << " " << this->GetStringFormatOfDummyAtomOmission(this->GetDummyAtomOmission()) << " "
			<< this->GetDummyAtomType() << " " << this->GetStringFormatOfDummyAtomPosition(this->GetDummyAtomPosition()) << std::endl
			<< std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->GetCharge() << std::endl;
	for(auto &atom : this->getAtoms())
	{
		atom->Write(stream);
	}
	stream << std::endl;
	if (this->GetImproperDihedrals().size() > 0)
	{
		stream << "IMPROPER" << std::endl;
		for(auto &dihedral : this->GetImproperDihedrals())
		{
			stream << std::left << std::setw(4) << dihedral.at(0) << std::left << std::setw(4) << dihedral.at(1) << std::left << std::setw(4) << dihedral.at(2) << std::left << std::setw(4) << dihedral.at(3) << std::endl;
		}
		stream << std::endl;
	}
	if(this->GetLoops().size() > 0)
	{
		stream << "LOOP" << std::endl;
		for(auto &loop : this->GetLoops())
		{
			stream << loop.first << " " << loop.second << std::endl;
		}
	}
	stream << std::endl
			<< "DONE" << std::endl;
}

