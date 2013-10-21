#ifndef PREPFILE_HPP
#define PREPFILE_HPP

#include <map>
#include <string>
#include <iostream>

namespace PrepFileSpace
{
    class PrepFileResidue;
    class PrepFile
    {
        public:
            //////////////////////////////////////////////////////////
            //                     TYPE DEFINITION                  //
            //////////////////////////////////////////////////////////
            /*! \def
              * A mapping between a residue name and its residue object
              */
            typedef std::map< std::string, PrepFileResidue* > ResidueMap;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Constructor
              * @param prep_file An existing prep file path to be read
              */
            PrepFile(const std::string& prep_file);

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to prep file path of the current object
              * @return path_ attribute of the current object of this class
              */
            ResidueMap& GetResidues();

            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a prep file
              */
            void Read(std::ifstream& in_file);
            /*! \fn
              * A function to residue section of the current prep file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a prep file
              */
            PrepFileResidue* ProcessResidueSection(std::ifstream& in_file);

            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            std::string path_;              /*!< Actual path of the given prep file */
            ResidueMap residues_;           /*!< Fill by PrepFileResidues */
            /*!< End of a prep file gets marked by STOP */
    };
}

#endif // PREPFILE_HPP
