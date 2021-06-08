// Created by: Dave Montgomery

#ifndef PDBDATABASEREFERENCESECTION_HPP
#define PDBDATABASEREFERENCESECTION_HPP

#include <string>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    class PdbDatabaseReference;

    class PdbDatabaseReferenceSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of database references
              */
            typedef std::vector<PdbDatabaseReference*> DatabaseReferenceVector;


            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbDatabaseReferenceSection();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbDatabaseReferenceSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the database_reference in a database_reference card
              * @return database_reference_ attribute of the current object of this class
              */
            DatabaseReferenceVector GetDatabaseReferences();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the database_reference attribute of the current object
              * Set the database_reference_ attribute of the current database_reference card
              * @param database_reference The database_reference attribute of the current object
              */
            void SetDatabaseReferences(DatabaseReferenceVector database_reference);
            /*! \fn
              * A function in order to add the database_reference attribute to the current object
              * Set the database_reference_ attribute of the current database_reference card
              * @param database_reference The database_reference attribute of the current object
              */
            void AddDatabaseReferences(PdbDatabaseReference *database_reference);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////
            std::string GetUniprotIDs();

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the database_reference card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            DatabaseReferenceVector database_reference_;      /*!< List of database_references >*/

    };
}

#endif // PDBDATABASEREFERENCESECTION_HPP
