// Created by: Dave Montgomery

#ifndef PDBSEQUENCEADVANCEDCARD_HPP
#define PDBSEQUENCEADVANCEDCARD_HPP

#include <string>
#include <iostream>


namespace PdbFileSpace
{
    class PdbSequenceAdvancedCard
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbSequenceAdvancedCard();
            /*! \fn
              * Constructor with required parameters
              * @param line
              */
            PdbSequenceAdvancedCard(std::string& line);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a sequence advanced card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the ID Code attribute in a sequence advanced card
              * @return id_code_ attribute of the current object of this class
              */
            std::string GetIdentifierCode();
            /*! \fn
              * An accessor function in order to access to the Residue Name attribute in a sequence advanced card
              * @return residue_name_ attribute of the current object of this class
              */
            std::string GetResidueName();
            /*! \fn
              * An accessor function in order to access to the Chain ID attribute in a sequence advanced card
              * @return chain_id_ attribute of the current object of this class
              */
            std::string GetChainId();
            /*! \fn
              * An accessor function in order to access to the Sequence Number attribute in a sequence advanced card
              * @return seq_num_ attribute of the current object of this class
              */
            int GetSequenceNumber();
            /*! \fn
              * An accessor function in order to access to the InsertionCode attribute in a sequence advanced card
              * @return i_code_ attribute of the current object of this class
              */
            std::string GetInsertionCode();
            /*! \fn
              * An accessor function in order to access to the Database attribute in a sequence advanced card
              * @return database_ attribute of the current object of this class
              */
            std::string GetDatabase();
            /*! \fn
              * An accessor function in order to access to the Database Accession attribute in a sequence advanced card
              * @return db_accession_ attribute of the current object of this class
              */
            std::string GetDatabaseAccession();
            /*! \fn
              * An accessor function in order to access to the Database Residue attribute in a sequence advanced card
              * @return db_res_ attribute of the current object of this class
              */
            std::string GetDatabaseResidue();
            /*! \fn
              * An accessor function in order to access to the Database Sequence attribute in a sequence advanced card
              * @return db_seq_ attribute of the current object of this class
              */
            int GetDatabaseSequence();
            /*! \fn
              * An accessor function in order to access to the Conflict attribute in a sequence advanced card
              * @return conflict_ attribute of the current object of this class
              */
            std::string GetConflict();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current sequence advanced card
              * @param record_name The record name attribute of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the Identifier Code attribute of the current object
              * Set the id_code_ attribute of the current sequence advanced card
              * @param id_code The Identifier Code attribute of the current object
              */
            void SetIdentifierCode(std::string id_code);
            /*! \fn
              * A mutator function in order to set the Residue Name attribute of the current object
              * Set the residue_name_ attribute of the current sequence advanced card
              * @param residue_name The Residue Name attribute of the current object
              */
            void SetResidueName(std::string residue_name);
            /*! \fn
              * A mutator function in order to set the Chain Id attribute of the current object
              * Set the chain_id_ attribute of the current sequence advanced card
              * @param chain_id The Chain Id attribute of the current object
              */
            void SetChainId(std::string chain_id);
            /*! \fn
              * A mutator function in order to set the Sequence Number attribute of the current object
              * Set the seq_num_ attribute of the current sequence advanced card
              * @param seq_num The Sequence Number attribute of the current object
              */
            void SetSequenceNumber(int seq_num);
            /*! \fn
              * A mutator function in order to set the Insertion Code attribute of the current object
              * Set the i_code_ attribute of the current sequence advanced card
              * @param i_code The Insertion Code attribute of the current object
              */
            void SetInsertionCode(std::string i_code);
            /*! \fn
              * A mutator function in order to set the Database attribute of the current object
              * Set the database_ attribute of the current sequence advanced card
              * @param database The Database attribute of the current object
              */
            void SetDatabase(std::string database);
            /*! \fn
              * A mutator function in order to set the Database Accession attribute of the current object
              * Set the db_accession_ attribute of the current sequence advanced card
              * @param db_accession The Database Accession attribute of the current object
              */
            void SetDatabaseAccession(std::string db_accession);
            /*! \fn
              * A mutator function in order to set the Database Residue attribute of the current object
              * Set the db_res_ attribute of the current sequence advanced card
              * @param db_res The Database Residue attribute of the current object
              */
            void SetDatabaseResidue(std::string db_res);
            /*! \fn
              * A mutator function in order to set the Database Sequence attribute of the current object
              * Set the db_seq_ attribute of the current sequence advanced card
              * @param db_seq The Database Sequence attribute of the current object
              */
            void SetDatabaseSequence(int db_seq);
            /*! \fn
              * A mutator function in order to set the Conflict attribute of the current object
              * Set the conflict_ attribute of the current sequence advanced card
              * @param conflict The Conflict attribute of the current object
              */
            void SetConflict(std::string conflict);

/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the sequence advanced card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;                       /*!< Record name of sequence advanced card which is the first column of each line of the card >*/
            std::string id_code_;                           /*!< ID code of this PDB entry >*/
            std::string residue_name_;                      /*!< Name of the PDB residue in conflict >*/
            std::string chain_id_;                          /*!< PDB chain identifier >*/
            int seq_num_;                                   /*!< PDB sequence number >*/
            std::string i_code_;                            /*!< PDB insertion code >*/
            std::string database_;                          /*!< Database of the sequence advanced card >*/
            std::string db_accession_;                      /*!< Sequence database accession number >*/
            std::string db_res_;                            /*!< Sequence database residue name >*/
            int db_seq_;                                    /*!< Sequence database sequence number >*/
            std::string conflict_;                          /*!< Conflict comment >*/
    };
}

#endif // PDBSEQUENCEADVANCEDCARD_HPP
