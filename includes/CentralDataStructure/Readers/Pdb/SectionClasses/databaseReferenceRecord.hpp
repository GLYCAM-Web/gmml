#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_DATABASEREFERENCERECORD_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_DATABASEREFERENCERECORD_HPP

#include <string>
#include <iostream>

namespace pdb
{
    class DatabaseReference
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            DatabaseReference(std::string& line);
            std::string GetUniprotID() const;
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cerr) const;
            void Write(std::ostream& stream) const;
        private:
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetRecordName() const;
            std::string GetIDCode() const;
            std::string GetChainID() const;
            int GetSeqBegin() const;
            std::string GetInsertBegin() const;
            int GetSeqEnd() const;
            std::string GetInsertEnd() const;
            std::string GetDatabase() const;
            std::string GetDatabaseAccession() const;
            std::string GetDatabaseIDCode() const;
            int GetDatabaseSeqBegin() const;
            std::string GetDatabaseInsBegin() const;
            int GetDatabaseSeqEnd() const;
            std::string GetDatabaseInsEnd() const;
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetRecordName(const std::string record_name);
            void SetIDCode(const std::string id_code);
            void SetChainId(const std::string chain_id);
            void SetSeqBegin(const int seq_begin);
            void SetInsertBegin(const std::string insert_begin);
            void SetSeqEnd(const int seq_end);
            void SetInsertEnd(const std::string insert_end);
            void SetDatabase(const std::string database);
            void SetDatabaseAccession(const std::string db_accession);
            void SetDatabaseIDCode(const std::string db_id_code);
            void SetDatabaseSeqBegin(const int db_seq_begin);
            void SetDatabaseInsBegin(const std::string db_ins_begin);
            void SetDatabaseSeqEnd(const int db_seq_end);
            void SetDatabaseInsEnd(const std::string db_ins_end);
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

#endif // GMML_INCLUDES_INPUTSET_PDBFILE_DATABASEREFERENCERECORD_HPP
