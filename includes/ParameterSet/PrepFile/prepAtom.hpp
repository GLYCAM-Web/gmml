#ifndef INCLUDES_PARAMETERSET_PREPFILE_PREPFILEATOM_HPP
#define INCLUDES_PARAMETERSET_PREPFILE_PREPFILEATOM_HPP

#include <string>
#include <iostream>
#include <iostream>
#include "includes/common.hpp"
#include "includes/CentralDataStructure/cdsAtom.hpp"

namespace prep
{
class PrepAtom : public cds::cdsAtom
{
public:
	//////////////////////////////////////////////////////////
	//                       Constructor                    //
	//////////////////////////////////////////////////////////
	PrepAtom(std::string& line);
	//////////////////////////////////////////////////////////
	//                         FUNCTIONS                    //
	//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////
	//                     DISPLAY FUNCTIONS                //
	//////////////////////////////////////////////////////////
	void Print(std::ostream& out = std::cerr);
	void Write(std::ostream &stream) const;
	//////////////////////////////////////////////////////////
	//                           ACCESSOR                   //
	//////////////////////////////////////////////////////////
	int GetIndex() const;
	std::string GetType() const;
	gmml::TopologicalType GetTopologicalType() const;
	int GetBondIndex() const;
	int GetAngleIndex() const;
	int GetDihedralIndex() const;
	double GetBondLength() const;
	double GetAngle() const;
	double GetDihedral() const;
	double GetCharge() const;
	//////////////////////////////////////////////////////////
	//                           MUTATOR                    //
	//////////////////////////////////////////////////////////
	void SetIndex(int index);
	void SetType(std::string type);
	void SetTopologicalType(gmml::TopologicalType topological_type);
	void SetBondIndex(int bond_index);
	void SetAngleIndex(int angle_index);
	void SetDihedralIndex(int dihedral_index);
	void SetBondLength(double bond_length);
	void SetAngle(double angle);
	void SetDihedral(double dihedral);
	void SetCharge(double charge);
private:
	//////////////////////////////////////////////////////////
	//                         FUNCTIONS                    //
	//////////////////////////////////////////////////////////
	gmml::TopologicalType ExtractAtomTopologicalType(std::istream& ss);
	std::string GetStringFormatOfTopologicalType(gmml::TopologicalType topological_type) const;
	std::string GetStringFormatOfTopologicalType() const;
	gmml::TopologicalType GetTopologicalTypeFromString(std::string topological_type) const;
	//////////////////////////////////////////////////////////
	//                         ATTRIBUTES                   //
	//////////////////////////////////////////////////////////
	int index_ = 0;                                 /*!< Atom index; fill by the first column of the residue section of the file */
	std::string type_ = "";                          /*!< Atom type; fill by the third column of the residue section of the file */
	gmml::TopologicalType topological_type_ = gmml::kTopTypeM;          /*!< Topological type (for chain extraction of the residue); fill by th 4th column of the residue section of the file */
	int bond_index_ = 0;;                            /*!< Bond index; fill by the 5th column of the residue section of the file */
	int angle_index_ = 0;;                           /*!< Angle index; fill by the 6th column of the residue section of the file */
	int dihedral_index_ = 0;                        /*!< Dihedral index; fill by the 7th column of the residue section of the file */
	double bond_length_ = gmml::dNotSet;                        /*!< Bond length; fill by the 8th column of the residue section of the file */
	double angle_ = gmml::dNotSet;                              /*!< Angle; fill by the 9th column of the residue section of the file */
	double dihedral_ = gmml::dNotSet;                           /*!< Dihedral; fill by the 10th column of the residue section of the file */
	double charge_ = gmml::dNotSet;                             /*!< Charge; fill by the 11th column of the residue section of the file */
	/*!< Sample line of the atom section of a prep file: 4 H1   H1  M  3  2  1  1.000    90.0     180.0     0.0 */
};
}

#endif // INCLUDES_PARAMETERSET_PREPFILE_PREPFILEATOM_HPP
