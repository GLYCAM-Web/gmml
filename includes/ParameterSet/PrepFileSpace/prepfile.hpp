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


        private:
            //////////////////////////////////////// ATTRIBUTES ////////////////////////////////////
            std::string path_;
            ResidueMap residues_;
    };
}

#endif // PREPFILE_HPP
