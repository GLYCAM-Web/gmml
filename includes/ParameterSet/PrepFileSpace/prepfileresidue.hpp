#ifndef PREPFILERESIDUE_HPP
#define PREPFILERESIDUE_HPP

#include <string>
#include <vector>
#include <map>

namespace PrepFileSpace
{
    /*! \enum
      * Coordinate type enumerator
      */
    enum CoordinateType { kINT, kXYZ };
    /*! \enum
      * Output format enumerator
      */
    enum OutputFormat { kFormatted, kBinary };
    /*! \enum
      * Geometry type enumerator
      */
    enum GeometryType { kGeometryCorrect, kGeometryChange };
    /*! \enum
      * Dummy atom position enumerator
      */
    enum DummyAtomPosition { kPositionAll, kPositionBeg };
    /*! \enum
      * Dummy atom omission enumerator
      */
    enum DummyAtomOmission { kOmit, kNomit };
    /*! \enum
      * Section type enumerator
      */
    enum SectionType { kSectionLoop, kSectionImproper, kSectionDone, kSectionOther };
    class PrepFileAtom;

    class PrepFileResidue
    {
        public:
            //////////////////////////////////////////////////////////
            //                     TYPE DEFINITION                  //
            //////////////////////////////////////////////////////////
            /*! \def
              * A mapping between two indices of atoms indicates a loop
              */
            typedef std::map<int, int> Loop;
            /*! \def
              * A collection of four atom types
              */
            typedef std::vector<std::string> Dihedral;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PrepFileResidue();

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to an atom in the residue by its name
              * @return An atom in atmos_ attribute of the current object of this class that has the given name
              */
            int GetAtomIndexByName(const std::string& name);

            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to load information of the residue by the given stream
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a prep file
              * @return A residue that is contained in the given prep file
              */
            PrepFileResidue* LoadFromStream(std::ifstream& in_file);
            /*! \fn
              * A function that extracts the name of the residue from the given string stream
              * @param ss A string stream that contains the name of the residue
              * @return Name of the resdiue
              */
            std::string ExtractResidueName(std::istream& ss);
            /*! \fn
              * A function to parse and extract coordinate type of the residue from a given string stream
              * @param A string stream that has the information about the coordinate type of the residue
              * @return Coordinate type of the residue
              */
            CoordinateType ExtractResidueCoordinateType(std::istream& ss);
            /*! \fn
              * A function to parse and extract output format of the residue from a given string stream
              * @param A string stream that has the information about the output format of the residue
              * @return Output format of the residue
              */
            OutputFormat ExtractResidueOutputFormat(std::istream& ss);
            /*! \fn
              * A function to parse and extract geometry type of the residue from a given string stream
              * @param A string stream that has the information about the geometry type of the residue
              * @return Geometry type of the residue
              */
            GeometryType ExtractResidueGeometryType(std::istream& ss);
            /*! \fn
              * A function to parse and extract dummy atom omission of the residue from a given string stream
              * @param A string stream that has the information about the dummy atom omission of the residue
              * @return Dummy atom omission of the residue
              */
            DummyAtomOmission ExtractResidueDummyAtomOmission(std::istream& ss);
            /*! \fn
              * A function to parse and extract dummy atom position of the residue from a given string stream
              * @param A string stream that has the information about the dummy atom position of the residue
              * @return Dummy atom position of the residue
              */
            DummyAtomPosition ExtractResidueDummyAtomPosition(std::istream& ss);
            /*! \fn
              * A function to parse and extract section type of the residue from a given string stream
              * @param A string stream that has the information about the section type of the residue
              * @return Section type of the residue
              */
            SectionType ExtractSectionType(std::string& line);
            /*! \fn
              * A function to parse loop section and extraxt loop information of the residue from a given file stream
              * @param A file stream that contains loop section of the residue
              * @return A loop object that is a mapping between source and destination atoms in the residue
              */
            Loop ExtractLoops(std::ifstream& in_file);
            /*! \fn
              * A function to parse improper dihedral section and extraxt improper dihedral information of the residue from a given file stream
              * @param A file stream that contains improper dihedral section of the residue
              * @return A vector of improper dihedrals in the residue
              */
            std::vector<Dihedral> ExtractImproperDihedral(std::ifstream& in_file);

            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            std::string title_;                         /*!< Residue title; fill by the first line of each residue section of the file */
            std::string name_;                          /*!< Residue name; fill by the first column of the 3rd line of each residue section of the file */
            CoordinateType coordinate_type_;            /*!< Coordinate type(INT, XYZ); fill by the 2nd column of the third line of each residue section of the file */
            OutputFormat output_format_;                /*!< Output format(Binary=1,Formatted=1); fill by the third column of the 3rd line of each residue section of the file */
            GeometryType geometry_type_;                /*!< Geometry type(CORRECT, CHANGE); fill by the first column of the 4th line of each residue section of the file */
            DummyAtomOmission dummy_atom_omission_;     /*!< Dummy atom omission(OMIT, NOMIT); fill by the 3rd column of the 4th line of each residue section of the file */
            std::string dummy_atom_type_;               /*!< Dummy atom type; fill by the 4th column of the 4th line of each residue section of the file */
            DummyAtomPosition dummy_atom_position_;     /*!< Dummy atom position(ALL, BEG); fill by the 5th column of the 4th line of each residue section of the file */
            double charge_;                             /*!< Total charge of the residue; fill by the 5th line of each residue section of the file */
            std::vector<PrepFileAtom*> atoms_;          /*!< Atoms in the resisue; fill by all lines between 6th line of each residue section of the file and a blank line in that section */
            std::vector<Dihedral> improper_dihedrals_;  /*!< Improper dihedrals; fill by all lines between IMPROPER title in each residue section of the file and a blank line in that section */
            Loop loops_;                                /*!< Loops; fill by all lines between LOOP title in each residue section of the file and a blank line in that section */
            /*!< End of each residue section gets marked by DONE */
            /*! \example
            * A Sample of residue section in a prep file:

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
