#ifndef PREPFILERESIDUE_HPP
#define PREPFILERESIDUE_HPP

#include <string>
#include <vector>
#include <map>

namespace PrepFileSpace
{
    enum CoordinateType { kINT, kXYZ };
    enum OutputFormat { kFormatted, kBinary };
    enum GeometryType { kGeometryCorrect, kGeometryChange };
    enum DummyAtomPosition { kPositionAll, kPositionBeg };
    enum DummyAtomOmission { kOmit, kNomit };
    enum SectionType { kSectionLoop, kSectionImproper, kSectionDone, kSectionOther };
    class PrepFileAtom;

    class PrepFileResidue
    {
        public:
            ///////////////////////////// TYPE DEFINITION /////////////////////////////
            // A mapping between two indices of atoms indicates a loop
            typedef std::map<int, int> Loop;
            // A collection of four atom types
            typedef std::vector<std::string> Dihedral;

            /////////////////////////////// CONSTRUCTOR ///////////////////////////////
            PrepFileResidue();

            ///////////////////////////// ACCESSOR ////////////////////////////////////
            int GetAtomIndexByName(const std::string& name);

            ///////////////////////////// FUNCTION ////////////////////////////////////
            PrepFileResidue* LoadFromStream(std::ifstream& in_file);
            std::string ExtractResidueName(std::istream& ss);
            CoordinateType ExtractResidueCoordinateType(std::istream& ss);
            OutputFormat ExtractResidueOutputFormat(std::istream& ss);
            GeometryType ExtractResidueGeometryType(std::istream& ss);
            DummyAtomOmission ExtractResidueDummyAtomOmission(std::istream& ss);
            DummyAtomPosition ExtractResidueDummyAtomPosition(std::istream& ss);
            SectionType ExtractSectionType(std::string& line);
            Loop ExtractLoops(std::ifstream& in_file);
            std::vector<Dihedral> ExtractImproperDihedral(std::ifstream& in_file);

            ////////////////////////// DISPLAY FUNCTION ///////////////////////////////
            void Print(std::ostream& out);

            ///////////////////////////// ATTRIBUTES //////////////////////////////////
            std::string title_;                         // Residue title; fill by the first line of each residue section of the file
            std::string name_;                          // Residue name; fill by the first column of the 3rd line of each residue section of the file
            CoordinateType coordinate_type_;            // Coordinate type(INT, XYZ); fill by the 2nd column of the third line of each residue section of the file
            OutputFormat output_format_;                // Output format(Binary=1,Formatted=1); fill by the third column of the 3rd line of each residue section of the file
            GeometryType geometry_type_;                // Geometry type(CORRECT, CHANGE); fill by the first column of the 4th line of each residue section of the file
            DummyAtomOmission dummy_atom_omission_;     // Dummy atom omission(OMIT, NOMIT); fill by the 3rd column of the 4th line of each residue section of the file
            std::string dummy_atom_type_;               // Dummy atom type; fill by the 4th column of the 4th line of each residue section of the file
            DummyAtomPosition dummy_atom_position_;     // Dummy atom position(ALL, BEG); fill by the 5th column of the 4th line of each residue section of the file
            double charge_;                             // Total charge of the residue; fill by the 5th line of each residue section of the file
            std::vector<PrepFileAtom*> atoms_;          // Atoms in the resisue; fill by all lines between 6th line of each residue section of the file and a blank line in that section
            std::vector<Dihedral> improper_dihedrals_;  // Improper dihedrals; fill by all lines between IMPROPER title in each residue section of the file and a blank line in that section
            Loop loops_;                                // Loops; fill by all lines between LOOP title in each residue section of the file and a blank line in that section
            // End of each residue section gets marked by DONE
            // A Sample of residue section in a prep file:
            /*
            ROH for aglycon

            ROH    INT 0
            CORRECT OMIT DU BEG
            -0.194
             1 DUMM DU  M  0 -1 -2  0.000     0.0       0.0     0.0
             2 DUMM DU  M  1  0 -1  1.000     0.0       0.0     0.0
             3 DUMM DU  M  2  1  0  1.000    90.0       0.0     0.0
             4 HO1  HO  M  3  2  1  1.000    90.0     180.0     0.445
             5 O1   OH  M  4  3  2  0.960   107.0     180.0    -0.639

            DONE
            */
    };
}

#endif // PREPFILERESIDUE_HPP
