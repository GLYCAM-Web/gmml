#ifndef PARAMETERFILE_HPP
#define PARAMETERFILE_HPP

#include <string>
#include <vector>
#include <map>


///////////////////// FORWARD DECLARATION //////////////////
namespace ParameterFileSpace
{
    class ParameterFileAtom;
    class ParameterFileBond;
    class ParameterFileAngle;
    class ParameterFileDihedral;
    class ParameterFile
    {
        public:
            ///////////////////////////////// TYPE DEFINITION ///////////////////////////////////////
            // A mapping between a atom type and its atom object
            typedef std::map<std::string, ParameterFileAtom*> AtomTypeMap;
            // A mapping between a connection of two atom types and its bond object
            typedef std::map<std::vector<std::string>, ParameterFileBond*> BondMap;
            // A mapping between a connection of three atom types and its angle object
            typedef std::map<std::vector<std::string>, ParameterFileAngle*> AngleMap;
            // A mapping between a connection of four atom types and its dihedral object
            typedef std::map<std::vector<std::string>, ParameterFileDihedral*> DihedralMap;

            //////////////////////////////////// CONSTRUCTOR ///////////////////////////////////////
            ParameterFile(const std::string& param_file);

            ///////////////////////////////////// ACCESSOR /////////////////////////////////////////
            const std::string& GetFilePath() const;
            const std::string& GetTitle() const;
            const AtomTypeMap& GetAtomTypes() const ;
            const BondMap& GetBonds() const;
            const AngleMap& GetAngles() const;
            const DihedralMap& GetDihedrals() const;

            ///////////////////////////////// FUNCTIONS ///////////////////////////////////////////
            void Read(std::ifstream& in_file);
            void ProcessAtomType(const std::string& line);
            void ProcessHydrophilicAtomType(const std::string& line);
            void ProcessBond(const std::string& line);
            void ProcessAngle(const std::string& line);
            void ProcessDihedral(std::string& line, int& line_number, std::ifstream& in_file);
            void ProcessImproperDihedral(std::string& line, int& line_number, std::ifstream& in_file);
            void ProcessHydrogenBond(const std::string& line);
            void ProcessEquivalentSymbols(const std::string& line);
            void ProcessPotentialParameter(const std::string& line);
            double ProcessDoubleDihedralDescription(const std::string& dscr, const std::string& key);

            ///////////////////////////////// DISPLAY FUNCTION ////////////////////////////////////
            void Print(std::ostream& out);

        private:
            //////////////////////////////////////// ATTRIBUTES ///////////////////////////////////
            std::string path_;
            std::string title_;          // Title of the parameter file
            AtomTypeMap atom_types_;     // A collection of mapping between atom type and its attributes
            BondMap bonds_;              // A collection of mapping between bond (double atom types) and its attributes
            AngleMap angles_;            // A collection of mapping between angle (tripple atom types) and its attributes
            DihedralMap dihedrals_;      // A collection of mapping between dihedral (quad atom types) and its attributes

    };
}

#endif // PARAMETERFILE_HPP
