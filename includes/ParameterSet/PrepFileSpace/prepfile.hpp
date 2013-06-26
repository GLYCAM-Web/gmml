#ifndef PREPFILE_HPP
#define PREPFILE_HPP

#include <map>
#include <string>

namespace PrepFileSpace
{
    class PrepFileResidue;
    class PrepFile
    {
        public:
            ///////////////////////////////// TYPE DEFINITION ///////////////////////////////////////
            typedef std::map< std::string, PrepFileResidue* > ResidueMap;

            //////////////////////////////////// CONSTRUCTOR ///////////////////////////////////////
            PrepFile(const std::string& prep_file);

            ///////////////////////////////////// ACCESSOR /////////////////////////////////////////
            ResidueMap& GetResidues();

            ///////////////////////////////// FUNCTIONS ////////////////////////////////////////////
            void Read(std::ifstream& in_file);
            PrepFileResidue* ProcessResidueSection(std::ifstream& in_file);

            ////////////////////////////////// DISPLAY FUNCTION ////////////////////////////////////
            void Print(std::ostream& out);

        private:
            //////////////////////////////////////// ATTRIBUTES ////////////////////////////////////
            std::string path_;              // Actual path of the given prep file
            ResidueMap residues_;           // Fill by PrepFileResidues
            //End of a prep file gets marked by STOP
    };
}

#endif // PREPFILE_HPP
