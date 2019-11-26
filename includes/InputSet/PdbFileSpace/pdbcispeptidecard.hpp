// Created by: Dave Montgomery

#ifndef PDBCISPEPTIDECARD_HPP
#define PDBCISPEPTIDECARD_HPP

#include <string>
#include <iostream>


namespace PdbFileSpace
{
    class PdbCISPeptideCard
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbCISPeptideCard();
            /*! \fn
              * Constructor with required parameters
              * @param line
              */
            PdbCISPeptideCard(std::string& line);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a CIS peptide card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the Serial Number attribute in a CIS peptide card
              * @return serial_number_ attribute of the current object of this class
              */
            int GetSerialNumber();
            /*! \fn
              * An accessor function in order to access to the Peptide 1 Residue Name attribute in a CIS peptide card
              * @return pep_1_ attribute of the current object of this class
              */
            std::string GetPeptide1ResidueName();
            /*! \fn
              * An accessor function in order to access to the Peptide 1 Chain ID attribute in a CIS peptide card
              * @return chain_id_1_ attribute of the current object of this class
              */
            std::string GetPeptide1ChainId();
            /*! \fn
              * An accessor function in order to access to the Peptide 1 Sequence Number attribute in a CIS peptide card
              * @return seq_num_1_ attribute of the current object of this class
              */
            int GetPeptide1SequenceNumber();
            /*! \fn
              * An accessor function in order to access to the Peptide 1 Insertion Code attribute in a CIS peptide card
              * @return i_code_1_ attribute of the current object of this class
              */
            std::string GetPeptide1InsertionCode();
            /*! \fn
              * An accessor function in order to access to the Peptide 2 Residue Name attribute in a CIS peptide card
              * @return pep_2_ attribute of the current object of this class
              */
            std::string GetPeptide2ResidueName();
            /*! \fn
              * An accessor function in order to access to the Peptide 2 Chain ID attribute in a CIS peptide card
              * @return chain_id_2_ attribute of the current object of this class
              */
            std::string GetPeptide2ChainId();
            /*! \fn
              * An accessor function in order to access to the Peptide 2 Sequence Number attribute in a CIS peptide card
              * @return seq_num_2_ attribute of the current object of this class
              */
            int GetPeptide2SequenceNumber();
            /*! \fn
              * An accessor function in order to access to the Peptide 2 Insertion Code attribute in a CIS peptide card
              * @return i_code_2_ attribute of the current object of this class
              */
            std::string GetPeptide2InsertionCode();
            /*! \fn
              * An accessor function in order to access to the Model Number attribute in a CIS peptide card
              * @return mod_num_ attribute of the current object of this class
              */
            int GetModelNumber();
            /*! \fn
              * An accessor function in order to access to the Measure attribute in a CIS peptide card
              * @return measure_ attribute of the current object of this class
              */
            float GetMeasure();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current CIS peptide card
              * @param record_name The record name attribute of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the Serial Number attribute of the current object
              * Set the serial_number_ attribute of the current CIS peptide card
              * @param serial_number The Serial Number attribute of the current object
              */
            void SetSerialNumber(int serial_number);
            /*! \fn
              * A mutator function in order to set the Peptide 1 Residue Name attribute of the current object
              * Set the pep_1_ attribute of the current CIS peptide card
              * @param pep_1 The Peptide 1 Residue Name attribute of the current object
              */
            void SetPeptide1ResidueName(std::string pep_1);
            /*! \fn
              * A mutator function in order to set the Peptide 1 Chain Id attribute of the current object
              * Set the chain_id_1_ attribute of the current CIS peptide card
              * @param chain_id_1 The Peptide 1 Chain Id attribute of the current object
              */
            void SetPeptide1ChainId(std::string chain_id_1);
            /*! \fn
              * A mutator function in order to set the Peptide 1 Sequence Number attribute of the current object
              * Set the seq_num_1_ attribute of the current CIS peptide card
              * @param seq_num_1 The Peptide 1 Sequence Number attribute of the current object
              */
            void SetPeptide1SequenceNumber(int seq_num_1);
            /*! \fn
              * A mutator function in order to set the Peptide 1 Insertion Code attribute of the current object
              * Set the i_code_1_ attribute of the current CIS peptide card
              * @param i_code_1 The Peptide 1 Insertion Code attribute of the current object
              */
            void SetPeptide1InsertionCode(std::string i_code_1);
            /*! \fn
              * A mutator function in order to set the Peptide 2 Residue Name attribute of the current object
              * Set the pep_2_ attribute of the current CIS peptide card
              * @param pep_2 The Peptide 2 Residue Name attribute of the current object
              */
            void SetPeptide2ResidueName(std::string pep_2);
            /*! \fn
              * A mutator function in order to set the Peptide 2 Chain Id attribute of the current object
              * Set the chain_id_2_ attribute of the current CIS peptide card
              * @param chain_id_2 The Peptide 2 Chain Id attribute of the current object
              */
            void SetPeptide2ChainId(std::string chain_id_2);
            /*! \fn
              * A mutator function in order to set the Peptide 2 Sequence Number attribute of the current object
              * Set the seq_num_2_ attribute of the current CIS peptide card
              * @param seq_num_2 The Peptide 2 Sequence Number attribute of the current object
              */
            void SetPeptide2SequenceNumber(int seq_num_2);
            /*! \fn
              * A mutator function in order to set the Peptide 2 Insertion Code attribute of the current object
              * Set the i_code_2_ attribute of the current CIS peptide card
              * @param i_code_2 The Peptide 2 Insertion Code attribute of the current object
              */
            void SetPeptide2InsertionCode(std::string i_code_2);
            /*! \fn
              * A mutator function in order to set the Model Number attribute of the current object
              * Set the mod_num_ attribute of the current CIS peptide card
              * @param mod_num The Model Number attribute of the current object
              */
            void SetModelNumber(int mod_num);
            /*! \fn
              * A mutator function in order to set the Measure attribute of the current object
              * Set the measure_ attribute of the current CIS peptide card
              * @param measure The Measure attribute of the current object
              */
            void SetMeasure(float measure);

/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the CIS peptide card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;                       /*!< Record name of CIS peptide card which is the first column of each line of the card >*/
            int serial_number_;                             /*!< Record Serial Number >*/
            std::string pep_1_;                             /*!< Residue Name >*/
            std::string chain_id_1_;                        /*!< Chain identifier >*/
            int seq_num_1_;                                 /*!< Residue sequence number >*/
            std::string i_code_1_;                          /*!< Insertion code >*/
            std::string pep_2_;                             /*!< Residue Name >*/
            std::string chain_id_2_;                        /*!< Chain identifier >*/
            int seq_num_2_;                                 /*!< Residue sequence number >*/
            std::string i_code_2_;                          /*!< Insertion code >*/
            int mod_num_;                                   /*!< Identifies the specific Model >*/
            float measure_;                                 /*!< Angle measurement in degrees >*/
    };
}

#endif // PDBCISPEPTIDECARD_HPP
