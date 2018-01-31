//#ifndef RESIDUE3DTEMPLATELIBRARY_HPP
//#define RESIDUE3DTEMPLATELIBRARY_HPP

//#include <string>
//#include <iostream>
//#include <vector>

//#include "./residueproperties.hpp"
//#include "./residue.hpp"
//#include "./gmml.hpp"

//namespace MolecularModeling
//{
//    class Residue3DTemplateLibrary
//{

//public:

//    //////////////////////////////////////////////////////////
//    //                    TYPE DEFINITION                   //
//    //////////////////////////////////////////////////////////
//    typedef std::vector<Residue*> ResidueVector;


//    //////////////////////////////////////////////////////////
//    //                       CONSTRUCTOR                    //
//    //////////////////////////////////////////////////////////
//    /*! \fn
//  * Default constructor
//  */
//    Residue3DTemplateLibrary();

//    /*! \fn
//  * Constructor to build a structure from a given input file
//  * @param file_path A file path that are required to build a structure
//  * @param type Type of the input which is selected from InputFileType enumerator
//  */
//    Residue3DTemplateLibrary(std::string file_path, gmml::InputFileType type);

//    /*! \fn
//  * Constructor to build a structure from multiple file types
//  * @param file_paths Set of set of file paths that are required to build a structure
//  * @param types Set of input file types of the inputs which are selected from InputFileType enumerator
//  */
//    Residue3DTemplateLibrary(std::vector<std::string> file_paths, std::vector<gmml::InputFileType> types);

//    //////////////////////////////////////////////////////////
//    //                       ACCESSOR                       //
//    //////////////////////////////////////////////////////////

//    /*! \fn
//  * Function to get residues from the library
//  */
//    ResidueVector GetAllResiduesFromLibrary();
//    /* BLF
//     * Is this a residuevector like what is in an assembly?  If
//     * so, this can't work.  This library will not have data in
//     * that sort of format.  It will be in all sorts of formats.
//     *
//     * So, maybe return a vector of strings?
//     */

//    /*! \fn
//  * Function to get residues from a type of the library
//  */
//    ResidueVector GetAllResiduesFromFileType(string typeoffile);
//    /* BLF
//     * Did you have some other function in mind for this?
//     *
//     * see comment above about residuevectors
//     */

//    /*! \fn
//  * Function to get residues based on type like amino
//  */
//    ResidueVector GetAllResiduesOfType(string typeofresidue);
//    /*
//    see comment above about residuevectors
//    */

//    /*! \fn
//  * Function to get all the file types in the library
//  */
//    std::string vector?  GetAllFileTypeFromLibrary();
//    /* would a vector of strings be a better return?
//     *
//     * Same applies to most other functions below here
//     * in this section
//     */

//    /*! \fn  BLF
//  * Function to list all the files in the library
//  */
//    std::string vector?  GetAllFilesFromLibrary();

//    /*! \fn  BLF
//  * Function to list all the files/types pairs in the library
//  */
//    std::string  vector? GetAllFilesTypesPairsFromLibrary();

//    /*! \fn  BLF
//  * Function to list files in the library that contain certain residue names
//  *
//  * This needs to be able to return more than one file per name
//  */
//    std::string vector?  GetFilesContainingResidueName(string vector of residue names);

//    /*! \fn  BLF
//  * Function to return the number of times each residue in a list is present
//  */


//    namespace 3DTemplateLibrary
//    {

//    vector of positive integers  GetNumberOfResidueNameOccurrences(string vector of residue names);



//    //////////////////////////////////////////////////////////
//    //                       MUTATOR                        //
//    //////////////////////////////////////////////////////////

//    BLF
//    Not sure it should be easy to change anything in the library.
//    So, maybe this section should stay empty.

//    The only exception I can think of would be to set notes or comments
//        about the library or its components, but I'm not sure that's necessary.


//    //////////////////////////////////////////////////////////
//    //                       FUNCTIONS                      //
//    //////////////////////////////////////////////////////////

//    /*! \fn
//  * Function to add new files to library
//  *
//  * BLF:  make this use pairs like the others do
//  */
//    void AddFile(std::vector<std::string> file_paths, gmml::InputFileType type);

//    /*! \fn  BLF
//  * Function to remove from library
//  */
//    void RemoveFile(std::string filepath);


//    /* BLF
//     * These need to go in the assembly class.  They don't go here.
//     *
//     * But, I'm including notes about them here.
//     */
//    /*! \fn
//  * Function to Build assembly from the file type library
//  *
//  * Important note:  this isn't a library of file types.  It is a library
//  * of files.  The files might be different types.
//  *
//  * What we most need is to be able to build an assembly of residue templates
//  * from the library.  Note that each residue might come from a different
//  * file and from a different type of file
//  */
//    Assembly BuildAssemblyOfResidueTemplates(std::string vector residue names, residuethreedtemplatelibrary residuelibrary);
//    /*! \fn
//  * Function to build subset assembly from the file type library
//  *
//  * This probably isn't necessary here.
//  */
//    Assembly BuildSubsetAssemblyFromFileType(std::string filetype);
//    /* But, someone might want (tho it would be a memory hog):
//     */
//    Assembly BuildAssemblyOfAllLibraryResidues(residuethreedtemplatelibrary residuelibrary);



//private:
//    //////////////////////////////////////////////////////////
//    //                       ATTRIBUTES                     //
//    //////////////////////////////////////////////////////////


//};
//}

//#endif // RESIDUE3DTEMPLATELIBRARY_HPP
