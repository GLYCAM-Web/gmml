#ifndef PDBQTFILE_HPP
#define PDBQTFILE_HPP

#include <string>
#include <vector>
#include <map>

namespace PdbqtFileSpace{

    class PdbqtAtom;
    class PdbqtModelCard;
    class PdbqtAtomCard;
    class PdbqtBranchCard;
    class PdbqtFile
    {
        public:

            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of pdb atom
              */
            typedef std::vector<PdbqtAtom*> PdbqtAtomVector;
            /*! \typedef
              * A mapping between a
              */
            typedef std::map<std::string, PdbqtAtomVector* > PdbqtResidueAtomsMap;
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbqtFile();
            /*! \fn
              * Constructor
              * @param pdb_file An existing pdbqt file path to be read
              */
            PdbqtFile(const std::string& pdbqt_file);
            /*! \fn
              * Load PDB file
              */
            PdbqtFile* LoadPdbqtFile();
            /*! \fn
              * @param pdb_file An existing pdbqt file path to be read
              */
            PdbqtFile* LoadPdbqtFile(const std::string& pdbqt_file);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to pdbqt file path of the current object
              * @return path_ attribute of the current object of this class
              */
            std::string GetPath();
            /*! \fn
              * An accessor function in order to access to the models attribute of the current object
              * @return models_ attribute of the current object of this class
              */
            PdbqtModelCard* GetModels();
            /*! \fn
              * An accessor function in order to access to all atoms of all residues of the current object
              * @return atoms The vector of all atoms in a pdbqt file in the order that they appear in the file
              */
            PdbqtResidueAtomsMap GetAllAtomsInOrder(std::vector<std::string>& key_order);

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function in order to set the path of the pdbqt file
              * @param pdb_path Pdb file path
              */
            void SetPath(std::string pdbqt_path);
            /*! \fn
              * A mutator function in order to set the models card of the current object
              * Set the models_ attribute of the current pdb file
              * @param models The models attribute of the current object
              */
            void SetModels(PdbqtModelCard* models);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a pdbqt file
              */
            bool Read(std::ifstream& in_file);
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * @param in_stream A stream contains whole contents of a pdbqt file
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseCards(std::ifstream& in_stream);
            /*! \fn
              * A function to parse the model crad that has been given as a stream
              * @param stream A stream contains model card of a pdbqt file
              * @param line Current line in the stream
              * @return Boolean value that indicates parsing has been done successfully or not
              */
            bool ParseModelCard(std::ifstream& stream, std::string& line);

            /*! \fn
              * A function to create an output pdbqt file with the given name
              * @param pdbqt_file Output pdb file name
              */
            void Write(const std::string& pdbqt_file);
            /*! \fn
              * A function to write back model card of the pdbqt file into an output stream
              * @param stream Intermediate output stream in order to write model card
              */
            void ResolveModelCard(std::ofstream& stream);            
            /*! \fn
              * A function to write back root card of the pdbqt file into an output stream
              * @param stream Intermediate output stream in order to write root card
              */
            void ResolveRootCard(std::ofstream& stream, PdbqtAtomCard* atom_card);
            /*! \fn
              * A function to write back branch cards of the pdbqt file into an output stream
              * @param stream Intermediate output stream in order to write branch cards
              */
            void ResolveBranchCards(std::ofstream& stream, std::vector<PdbqtBranchCard*> branches);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdbqt file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string path_;
            PdbqtModelCard* models_;

    };
}
#endif // PDBQTFILE_HPP
