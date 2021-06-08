// Created by: Dave Montgomery

#ifndef PDBDATABASEREFERENCE_HPP
#define PDBDATABASEREFERENCE_HPP

#include <string>
#include <iostream>


namespace PdbFileSpace
{
    class PdbDatabaseReference
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbDatabaseReference();
            /*! \fn
              * Constructor with required parameters
              * @param line
              */
            PdbDatabaseReference(std::string& line);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a database reference
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the id_code_ attribute in a database reference
              * @return id_code_ attribute of the current object of this class
              */
            std::string GetIDCode();
            /*! \fn
              * An accessor function in order to access to the Chain ID attribute in a database reference
              * @return chain_id_ attribute of the current object of this class
              */
            std::string GetChainID();
            /*! \fn
              * An accessor function in order to access to the Sequence Begin attribute in a database reference
              * @return seq_begin_ attribute of the current object of this class
              */
            int GetSeqBegin();
            /*! \fn
              * An accessor function in order to access to the Insert Begin attribute in a database reference
              * @return insert_begin_ attribute of the current object of this class
              */
            std::string GetInsertBegin();
            /*! \fn
              * An accessor function in order to access to the Sequence End attribute in a database reference
              * @return seq_end_ attribute of the current object of this class
              */
            int GetSeqEnd();
            /*! \fn
              * An accessor function in order to access to the Insert End attribute in a database reference
              * @return insert_end_ attribute of the current object of this class
              */
            std::string GetInsertEnd();
            /*! \fn
              * An accessor function in order to access to the Database attribute in a database reference
              * @return database_ attribute of the current object of this class
              */
            std::string GetDatabase();
            /*! \fn
              * An accessor function in order to access to the Database Accession attribute in a database reference
              * @return db_accession_ attribute of the current object of this class
              */
            std::string GetDatabaseAccession();
            /*! \fn
              * An accessor function in order to access to the Database ID Code attribute in a database reference
              * @return db_id_code_ attribute of the current object of this class
              */
            std::string GetDatabaseIDCode();
            /*! \fn
              * An accessor function in order to access to the Database Sequence Begin attribute in a database reference
              * @return db_seq_begin_ attribute of the current object of this class
              */
            int GetDatabaseSeqBegin();
            /*! \fn
              * An accessor function in order to access to the Database Insert Begin attribute in a database reference
              * @return db_ins_beg_ attribute of the current object of this class
              */
            std::string GetDatabaseInsBegin();
            /*! \fn
              * An accessor function in order to access to the Database Sequence End attribute in a database reference
              * @return db_seq_end_ attribute of the current object of this class
              */
            int GetDatabaseSeqEnd();
            /*! \fn
              * An accessor function in order to access to the Database Insert End attribute in a database reference
              * @return db_ins_end_ attribute of the current object of this class
              */
            std::string GetDatabaseInsEnd();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current database reference
              * @param record_name The record name attribute of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the ID Code attribute of the current object
              * Set the id_code_ attribute of the current database reference
              * @param id_code The ID Code attribute of the current object
              */
            void SetIDCode(const std::string id_code);
            /*! \fn
              * A mutator function in order to set the Chain ID attribute of the current object
              * Set the chain_id_ attribute of the current database reference
              * @param chain_id The Chain ID attribute of the current object
              */
            void SetChainId(const std::string chain_id);
            /*! \fn
              * A mutator function in order to set the Sequence Begin attribute of the current object
              * Set the seq_begin_ attribute of the current database reference
              * @param seq_begin The Sequence Begin attribute of the current object
              */
            void SetSeqBegin(const int seq_begin);
            /*! \fn
              * A mutator function in order to set the Insert Begin attribute of the current object
              * Set the insert_begin_ attribute of the current database reference
              * @param insert_begin The Insert Begin attribute of the current object
              */
            void SetInsertBegin(const std::string insert_begin);
            /*! \fn
              * A mutator function in order to set the Sequence End attribute of the current object
              * Set the seq_end_ attribute of the current database reference
              * @param seq_end The Sequence End attribute of the current object
              */
            void SetSeqEnd(const int seq_end);
            /*! \fn
              * A mutator function in order to set the Insert End attribute of the current object
              * Set the insert_end_ attribute of the current database reference
              * @param insert_end The Insert End attribute of the current object
              */
            void SetInsertEnd(const std::string insert_end);
            /*! \fn
              * A mutator function in order to set the Database attribute of the current object
              * Set the database_ attribute of the current database reference
              * @param database The Database attribute of the current object
              */
            void SetDatabase(const std::string database);
            /*! \fn
              * A mutator function in order to set the Database Accession attribute of the current object
              * Set the db_accession_ attribute of the current database reference
              * @param db_accession The Database Accession attribute of the current object
              */
            void SetDatabaseAccession(const std::string db_accession);
            /*! \fn
              * A mutator function in order to set the Database ID Code attribute of the current object
              * Set the db_id_code_ attribute of the current database reference
              * @param db_id_code The Database ID Code attribute of the current object
              */
            void SetDatabaseIDCode(const std::string db_id_code);
            /*! \fn
              * A mutator function in order to set the Database Sequence Begin attribute of the current object
              * Set the db_seq_begin_ attribute of the current database reference
              * @param db_seq_begin The Database Sequence Begin attribute of the current object
              */
            void SetDatabaseSeqBegin(const int db_seq_begin);
            /*! \fn
              * A mutator function in order to set the Database Insert Begin attribute of the current object
              * Set the db_ins_begin_ attribute of the current database reference
              * @param db_ins_begin The Database Insert Begin attribute of the current object
              */
            void SetDatabaseInsBegin(const std::string db_ins_begin);
            /*! \fn
              * A mutator function in order to set the Database Sequence End attribute of the current object
              * Set the db_seq_end_ attribute of the current database reference
              * @param db_seq_end The Database Sequence End attribute of the current object
              */
            void SetDatabaseSeqEnd(const int db_seq_end);
            /*! \fn
              * A mutator function in order to set the Database Insert End attribute of the current object
              * Set the db_ins_end_ attribute of the current database reference
              * @param db_ins_end The Database Insert End attribute of the current object
              */
            void SetDatabaseInsEnd(const std::string db_ins_end);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the database reference contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;                       /*!< Record name of database reference which is the first column of each line of the card >*/
            std::string id_code_;                           /*!< ID code of this entry >*/
            std::string chain_id_;                          /*!< Chain ID of the database reference >*/
            int seq_begin_;                                 /*!< Initial sequence number of the PDB sequence segment, right justified >*/
            std::string insert_begin_;                      /*!< Initial insertion code of the PDB sequence segment >*/
            int seq_end_;                                   /*!< Ending sequence number of the PDB sequence segment, right justified >*/
            std::string insert_end_;                        /*!< Ending insertion code of the PDB sequence segment >*/
            std::string database_;                          /*!< Sequence database name >*/
            std::string db_accession_;                      /*!< Sequence database accession code >*/
            std::string db_id_code_;                        /*!< Sequence  database identification code >*/
            int db_seq_begin_;                              /*!< Initial sequence number of the database segment >*/
            std::string db_ins_beg_;                        /*!< Insertion code of initial residue of the segment, if PDB is the reference >*/
            int db_seq_end_;                                /*!< Ending sequence number of the database segment >*/
            std::string db_ins_end_;                        /*!< Insertion code of the ending residue of the segment, if PDB is the reference >*/

    };
}

#endif // PDBDATABASEREFERENCE_HPP
