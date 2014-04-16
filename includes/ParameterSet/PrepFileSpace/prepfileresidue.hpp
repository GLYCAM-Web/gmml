#ifndef PREPFILERESIDUE_HPP
#define PREPFILERESIDUE_HPP

#include <string>
#include <vector>
#include <map>
#include <iostream>

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
            typedef std::vector<PrepFileAtom*> PrepFileAtomVector;
            typedef std::vector<Dihedral> DihedralVector;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PrepFileResidue();




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
            DihedralVector ExtractImproperDihedral(std::ifstream& in_file);

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to an atom in the residue by its name
              * @return An atom in atmos_ attribute of the current object of this class that has the given name
              */
              int GetAtomIndexByName(const std::string& name);

///Delaram
            /*! \fn
              * An accessor function in order to access to the title attribute of the current object
              * The attribute is set by the contents of the given file
              * @return title_ of the current object of this class
              */
            std::string GetTitle();
            /*! \fn
              * An accessor function in order to access to name attribute of the current object
              * The attribute is set by the contents of the given file
              * @return name_ of the current object of this class
              */
            std::string GetName();
            /*! \fn
              * An accessor function in order to access to coordiante type of the current object
              * @return coordinate_type_ attribute of the current object of this class
              */
            CoordinateType GetCoordinateType();
            /*! \fn
              * An accessor function in order to access to output format of the current object
              * The attribute is set by the contents of the given file
              * @return output_format_ of the current object of this class
              */
            OutputFormat GetOutputFormat();
            /*! \fn
              * An accessor function in order to access to geometry type of the current object
              * The attribute is set by the contents of the given file
              * @return geometry_type_ of the current object of this class
              */
            GeometryType GetGeometryType();
            /*! \fn
              * An accessor function in order to access to dummy atom omission of the current object
              * The attribute is set by the contents of the given file
              * @return dummy_atom_omission_ of the current object of this class
              */
            DummyAtomOmission GetDummyAtomOmission();
            /*! \fn
              * An accessor function in order to access to dummy atom type of the current object
              * The attribute is set by the contents of the given file
              * @return dummy_atom_type_ of the current object of this class
              */
            std::string GetDummyAtomType();
            /*! \fn
              * An accessor function in order to access to angle attribute of the current object
              * The attribute is set by the contents of the given file
              * @return dummy_atom_position of the current object of this class
              */
            DummyAtomPosition GetDummyAtomPosition();
            /*! \fn
              * An accessor function in order to access to charge attribute of the current object
              * The attribute is set by the contents of the given file
              * @return charge_ of the current object of this class
              */
            double GetCharge();
            /*! \fn
              * An accessor function in order to access to atoms of the current object
              * The attribute is set by the contents of the given file
              * @return atoms_ of the current object of this class
              */
            PrepFileAtomVector GetAtoms();
            /*! \fn
              * An accessor function in order to access to improper dihedral of the current object
              * The attribute is set by the contents of the given file
              * @return improper_dihedral_ of the current object of this class
              */
            DihedralVector GetImproperDihedrals();
            /*! \fn
              * An accessor function in order to access to loops of the current object
              * The attribute is set by the contents of the given file
              * @return loops_ of the current object of this class
              */
            Loop GetLoops();

            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the title of the current object
              * Set the title_ attribute of the current object
              * @param title The index attribute of the current object
              */
            void SetTitle(std::string title);
            /*! \fn
              * A mutator function in order to set the name of the current object
              * Set the name_ attribute of the current object
              * @param name The force name of the current object
              */
            void SetName(std::string name);
            /*! \fn
              * A mutator function in order to set the coordinate_type of the current object
              * Set the coordinate_type_ attribute of the current object
              * @param coordinate_type The coordinate_type of the current object
              */
            void SetCoordinateType(CoordinateType coordinate_type);
            /*! \fn
              * A mutator function in order to set the output format of the current object
              * Set the output_format_ attribute of the current object
              * @param output_format The topological_type attribute of the current object
              */
            void SetOutputFormat(OutputFormat output_format);
            /*! \fn
              * A mutator function in order to set the geometry type of the current object
              * Set the geometry_type_ attribute of the current object
              * @param geometry_type The geometry_type of the current object
              */
            void SetBondIndex(GeometryType geometry_type);
            /*! \fn
              * A mutator function in order to set the dummy atom omission of the current object
              * Set the dummy_atom_omission_ attribute of the current object
              * @param dummy_atom_omission The dummy_atom_omission of the current object
              */
            void GetDummyAtomOmission(DummyAtomOmission dummy_atom_omission);
            /*! \fn
              * A mutator function in order to set the dummy atom type of the current object
              * Set the dummy_atom_type_ attribute of the current atom
              * @param dummy_atom_type The dummy_atom_type attribute of the current object
              */
            void SetDummyAtomType(std::string dummy_atom_type);
            /*! \fn
              * A mutator function in order to set the dummy atom position of the current object
              * Set the dummy_atom_position_ attribute of the current object
              * @param dummy_atom_position The dummy_atom_position of the current object
              */
            void SetDummyAtomPosition(DummyAtomPosition dummy_atom_position);
            /*! \fn
              * A mutator function in order to set the charge of the current object
              * Set the charge_ attribute of the current object
              * @param charge The charge of the current object
              */
            void SetCharge(double charge);
            /*! \fn
              * A mutator function in order to set the atoms of the current object
              * Set the atoms_ attribute of the current object
              * @param atoms The atoms of the current object
              */
            void SetAtoms(PrepFileAtomVector atoms);
            /*! \fn
              * A mutator function in order to add the atom of the current object
              * Add the atom_ attribute of the current object
              * @param atom The atom of the current object
              */
            void AddAtom(PrepFileAtom* atom);
            /*! \fn
              * A mutator function in order to set the improper dihedrals of the current object
              * Set the improper_dihedrals_ attribute of the current object
              * @param improper_dihedrals The improper_dihedrals of the current object
              */
            void SetImproperDihedrals(DihedralVector improper_dihedrals);
            /*! \fn
              * A mutator function in order to add the improper dihedral of the current object
              * Add the improper_dihedral_ attribute of the current object
              * @param improper_dihedral The improper_dihedral of the current object
              */
            void AddImproperDihedral(Dihedral improper_dihedral);
            /*! \fn
              * A mutator function in order to set the Loops of the current object
              * Set the loops_ attribute of the current object
              * @param loops The loops of the current object
              */
            void SetLoops(Loop loops);

///Delaram

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
            PrepFileAtomVector atoms_;          /*!< Atoms in the resisue; fill by all lines between 6th line of each residue section of the file and a blank line in that section */
            DihedralVector improper_dihedrals_;  /*!< Improper dihedrals; fill by all lines between IMPROPER title in each residue section of the file and a blank line in that section */
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
