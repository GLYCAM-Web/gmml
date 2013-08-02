#ifndef LIBRARYFILE_HPP
#define LIBRARYFILE_HPP

#include <string>
#include <map>

namespace LibraryFileSpace
{
    class LibraryFileResidue;
    class LibraryFileAtom;
    class LibraryFile
    {
        public:
            ///////////////////////////////// TYPE DEFINITION ///////////////////////////////////////
            typedef std::map<std::string, LibraryFileResidue*> ResidueMap;

            //////////////////////////////////// CONSTRUCTOR ///////////////////////////////////////
            LibraryFile(const std::string& lib_file);

            ///////////////////////////////////// ACCESSOR /////////////////////////////////////////
            const std::string& GetFilePath() const;
            const ResidueMap& GetResidues() const;

            ///////////////////////////////// FUNCTIONS ///////////////////////////////////////////
            void Read(std::ifstream& in_file);
            LibraryFileAtom* ProcessAtom(std::string& line);

            ///////////////////////////////// DISPLAY FUNCTION ////////////////////////////////////
            void Print(std::ostream& out);

        private:
            std::string path_;          // Path of the given library file
            ResidueMap residues_;       // Map of residues including in the given file mapped to their names
    };
}

#endif // LIBRARYFILE_HPP
