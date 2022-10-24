#ifndef INCLUDES_PARAMETERSET_PREPFILE_PREPRESIDUE_HPP_
#define INCLUDES_PARAMETERSET_PREPFILE_PREPRESIDUE_HPP_

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/ParameterSet/PrepFile/prepAtom.hpp"

#include <string>
#include <map>
#include <vector>
#include <iostream>

namespace prep
{
enum CoordinateType { kINT, kXYZ };
enum OutputFormat { kFormatted = 0, kBinary = 1 };
enum GeometryType { kGeometryCorrect, kGeometryChange };
enum DummyAtomPosition { kPositionAll, kPositionBeg };
enum DummyAtomOmission { kOmit, kNomit };
enum SectionType { kSectionLoop, kSectionImproper, kSectionDone, kSectionOther };
class PrepResidue : public cds::Residue
{
public:
	//////////////////////////////////////////////////////////
	//                     TYPE DEFINITION                  //
	//////////////////////////////////////////////////////////
	typedef std::vector<std::string> Dihedral; // This looks poorly named
	typedef std::vector<Dihedral> DihedralVector;
	//////////////////////////////////////////////////////////
	//                       Constructor                    //
	//////////////////////////////////////////////////////////
	PrepResidue(std::ifstream& in_file, std::string& line);
	~PrepResidue() {std::cout << "PrepResidue dtor for " << this->getName() << ", ";}
	//////////////////////////////////////////////////////////
	//                       ACCESSOR                       //
	//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////
	//                       MUTATOR                        //
	//////////////////////////////////////////////////////////
	void SetTitle(std::string title);
	void SetCoordinateType(CoordinateType coordinate_type);
	void SetOutputFormat(OutputFormat output_format);
	void SetGeometryType(GeometryType geometry_type);
	void SetDummyAtomOmission(DummyAtomOmission dummy_atom_omission);
	void SetDummyAtomType(std::string dummy_atom_type);
	void SetDummyAtomPosition(DummyAtomPosition dummy_atom_position);
	void SetCharge(double charge);
	void SetImproperDihedrals(DihedralVector improper_dihedrals);
	void AddImproperDihedral(Dihedral improper_dihedral);
	void AddLoop(std::pair<std::string, std::string> loop);
	//////////////////////////////////////////////////////////
	//                     FUNCTIONS                        //
	//////////////////////////////////////////////////////////
	void Generate3dStructure();
	void DeleteDummyAtoms();
	void SetConnectivities();
	//void RecursivelySetConnectivities(std::vector<PrepAtom*>::iterator& currentAtom, std::vector<PrepAtom*>::iterator connectionPoint, std::vector<PrepAtom*>::iterator& lastAtom);
	std::vector<std::string> GetAtomNames() const;
	std::vector<std::string> GetHeavyAtomNames() const;
	double CalculatePrepResidueCharge();
	std::string Print() const;
	void Write(std::ostream &stream);
private:
	//////////////////////////////////////////////////////////
	//                         FUNCTIONS                    //
	//////////////////////////////////////////////////////////
	void ExtractResidueName(std::istream& ss);
	void ExtractResidueCoordinateType(std::istream& ss);
	void ExtractResidueOutputFormat(std::istream& ss);
	void ExtractResidueGeometryType(std::istream& ss);
	void ExtractResidueDummyAtomOmission(std::istream& ss);
	void ExtractResidueDummyAtomPosition(std::istream& ss);
	prep::SectionType ExtractSectionType(std::string& line);
	void ExtractLoops(std::ifstream& in_file);
	void ExtractImproperDihedral(std::ifstream& in_file);
	//////////////////////////////////////////////////////////
	//                         ACCESSOR                     //
	//////////////////////////////////////////////////////////
	std::string GetTitle() const;
	CoordinateType GetCoordinateType() const;
	OutputFormat GetOutputFormat() const;
	GeometryType GetGeometryType() const;
	DummyAtomOmission GetDummyAtomOmission() const;
	std::string GetDummyAtomType() const;
	DummyAtomPosition GetDummyAtomPosition() const;
	double GetCharge() const;
	DihedralVector GetImproperDihedrals() const;
	std::vector<std::pair<std::string, std::string>> GetLoops() const;
	std::string GetStringFormatOfCoordinateType(CoordinateType coordinate_type) const;
	std::string GetStringFormatOfOutputFormat(OutputFormat output_format) const;
	std::string GetStringFormatOfGeometryType(GeometryType geometry_type) const;
	std::string GetStringFormatOfDummyAtomPosition(DummyAtomPosition dummy_atom_position) const;
	std::string GetStringFormatOfDummyAtomOmission(DummyAtomOmission dummy_atom_omission) const;
	std::string GetStringFormatOfSectionType(SectionType section_type) const;
	CoordinateType GetCoordinateTypeFromString(std::string coordinate_type) const;
	OutputFormat GetOutputFormatFromString(std::string output_format) const;
	GeometryType GetGeometryTypeFromString(std::string geometry_type) const;
	DummyAtomPosition GetDummyAtomPositionFromString(std::string dummy_atom_position) const;
	DummyAtomOmission GetDummyAtomOmissionFromString(std::string dummy_atom_omission) const;
	SectionType GetSectionTypeFromString(std::string section_type) const;
	//////////////////////////////////////////////////////////
	//                         ATTRIBUTES                   //
	//////////////////////////////////////////////////////////
	std::string title_ = "";                         //!< Residue title; fill by the first line of each residue section of the file
	CoordinateType coordinate_type_ = prep::kINT;            //!< Coordinate type(INT, XYZ); fill by the 2nd column of the third line of each residue section of the file
	OutputFormat output_format_ = prep::kFormatted;                //!< Output format(Binary=1,Formatted=1); fill by the third column of the 3rd line of each residue section of the file
	GeometryType geometry_type_ = prep::kGeometryCorrect;                //!< Geometry type(CORRECT, CHANGE); fill by the first column of the 4th line of each residue section of the file
	DummyAtomOmission dummy_atom_omission_ = prep::kOmit;     //!< Dummy atom omission(OMIT, NOMIT); fill by the 3rd column of the 4th line of each residue section of the file
	std::string dummy_atom_type_ = "DU";               //!< Dummy atom type; fill by the 4th column of the 4th line of each residue section of the file
	DummyAtomPosition dummy_atom_position_ = prep::kPositionBeg;     //!< Dummy atom position(ALL, BEG); fill by the 5th column of the 4th line of each residue section of the file
	double charge_ = 0.0;                             //!< Total charge of the residue; fill by the 5th line of each residue section of the file
	DihedralVector improper_dihedrals_;  //!< Improper dihedrals; fill by all lines between IMPROPER title in each residue section of the file and a blank line in that section
	std::vector<std::pair<std::string, std::string>> loops_;                                //!< Loops; fill by all lines between LOOP title in each residue section of the file and a blank line in that section
	//	!< End of each residue section gets marked by DONE
	//	! \example
	//	 * A Sample of residue section in a prep file:
	//
	//            ROH for aglycon
	//
	//            ROH    INT 0
	//            CORRECT OMIT DU BEG
	//            -0.194
	//             1 DUMM DU  M  0 -1 -2  0.000     0.0       0.0     0.0
	//             2 DUMM DU  M  1  0 -1  1.000     0.0       0.0     0.0
	//             3 DUMM DU  M  2  1  0  1.000    90.0       0.0     0.0
	//             4 HO1  HO  M  3  2  1  1.000    90.0     180.0     0.445
	//             5 O1   OH  M  4  3  2  0.960   107.0     180.0    -0.639
	//
	//            DONE

};
} // namespace
#endif // INCLUDES_PARAMETERSET_PREPFILE_PREPRESIDUE_HPP_
