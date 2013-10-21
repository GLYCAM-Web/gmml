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
            //////////////////////////////////////////////////////////
            //                     TYPE DEFINITION                  //
            //////////////////////////////////////////////////////////
            /*! \def
              * A mapping between a name of a residue and its residue object
              */
            typedef std::map<std::string, LibraryFileResidue*> ResidueMap;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Constructor
              * @param lib_file An existing library file path to be read
              */
            LibraryFile(const std::string& lib_file);

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to library file path of the current object
              * @return path_ attribute of the current object of this class
              */
            const std::string& GetFilePath() const;
            /*! \fn
              * An accessor function in order to access to residue map created based on the content of the given file
              * @return residues_ attribute of the current object of this class
              */
            const ResidueMap& GetResidues() const;

            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a library file
              */
            void Read(std::ifstream& in_file);
            /*! \fn
              * Process a line of the atom section of a library file and create a new atom object
              * @param line A line from the atom section of a library file
              * @return A new LibraryFileAtom instance created by the information of the given line
              */
            LibraryFileAtom* ProcessAtom(std::string& line);

            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::out);

        private:
            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            std::string path_;          /*!< Path of the given library file */
            ResidueMap residues_;       /*!< Map of residues included in the given file mapped to their names */
    };
}

#endif // LIBRARYFILE_HPP
