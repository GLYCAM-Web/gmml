// Created by: Dave Montgomery

#ifndef PDBJOURNALSECTION_HPP
#define PDBJOURNALSECTION_HPP

#include <string>
#include <sstream>
#include <iostream>
#include <vector>

namespace PdbFileSpace
{
    class PdbJournalSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbJournalSection();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbJournalSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the record name in a journal section
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the authors in a journal section
              * @return authors_ attribute of the current object of this class
              */
            std::vector<std::string> GetAuthors();
            /*! \fn
              * An accessor function in order to access to the title in a journal section
              * @return title_ attribute of the current object of this class
              */
            std::string GetTitle();
            /*! \fn
              * An accessor function in order to access to the editors in a journal section
              * @return editors_ attribute of the current object of this class
              */
            std::vector<std::string> GetEditors();
            /*! \fn
              * An accessor function in order to access to the reference in a journal section
              * @return reference_ attribute of the current object of this class
              */
            std::string GetReference();
            /*! \fn
              * An accessor function in order to access to the publisher in a journal section
              * @return publisher_ attribute of the current object of this class
              */
            std::string GetPublisher();
            /*! \fn
              * An accessor function in order to access to the Reference Numbers in a journal section
              * @return reference_nums_ attribute of the current object of this class
              */
            std::vector<std::string> GetReferenceNumbers();
            /*! \fn
              * An accessor function in order to access to the PMID in a journal section
              * @return pmid_ attribute of the current object of this class
              */
            std::string GetPMID();
            /*! \fn
              * An accessor function in order to access to the DOI in a journal section
              * @return doi_ attribute of the current object of this class
              */
            std::string GetDOI();
            /*! \fn
              * An accessor function in order to access to the text in a journal section
              * @return text_ attribute of the current object of this class
              */
            std::string GetText();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current journal section
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the authors of the current object
              * Set the authors_ attribute of the current journal section
              * @param authors The authors of the current object
              */
            void SetAuthors(std::vector<std::string> authors);
            /*! \fn
              * A mutator function in order to set the title of the current object
              * Set the title_ attribute of the current journal section
              * @param title The title of the current object
              */
            void SetTitle(const std::string title);
            /*! \fn
              * A mutator function in order to set the editors of the current object
              * Set the editors_ attribute of the current journal section
              * @param editors The editors of the current object
              */
            void SetEditors(std::vector<std::string> editors);
            /*! \fn
              * A mutator function in order to set the reference of the current object
              * Set the reference_ attribute of the current journal section
              * @param reference The reference of the current object
              */
            void SetReference(const std::string reference);
            /*! \fn
              * A mutator function in order to set the publisher of the current object
              * Set the publisher_ attribute of the current journal section
              * @param publisher The publisher of the current object
              */
            void SetPublisher(const std::string publisher);
            /*! \fn
              * A mutator function in order to set the reference_nums of the current object
              * Set the reference_nums_ attribute of the current journal section
              * @param reference_nums The reference_nums of the current object
              */
            void SetReferenceNumbers(std::vector<std::string> reference_nums);
            /*! \fn
              * A mutator function in order to set the pmid of the current object
              * Set the pmid_ attribute of the current journal section
              * @param pmid The pmid of the current object
              */
            void SetPMID(const std::string pmid);
            /*! \fn
              * A mutator function in order to set the doi of the current object
              * Set the doi_ attribute of the current journal section
              * @param doi The doi of the current object
              */
            void SetDOI(const std::string doi);
            /*! \fn
              * A mutator function in order to set the text of the current object
              * Set the text_ attribute of the current journal section
              * @param text The text of the current object
              */
            void SetText(const std::string text);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb journal section contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;                 /*!< Record name which appears in the first column of each line in a pdb file >*/
            std::vector<std::string> authors_;        /*!< List of authors that appear in Journal record of a pdb file >*/
            std::string title_;                       /*!< Title that appears in Journal record of a pdb file >*/
            std::vector<std::string> editors_;        /*!< List of editors that appear in Journal record of a pdb file >*/
            std::string reference_;                   /*!< Reference that appears in Journal record of a pdb file >*/
            std::string publisher_;                   /*!< Publisher that appears in Journal record of a pdb file >*/
            std::vector<std::string> reference_nums_; /*!< List of reference numbers that appear in Journal record of a pdb file >*/
            std::string pmid_;                        /*!< Pub Med ID number that appears in Journal record of a pdb file >*/
            std::string doi_;                         /*!< DOI number that appears in Journal record of a pdb file >*/
            std::string text_;                        /*!< Text in a Journal Section >*/


    };
}

#endif // PDBJOURNALSECTION_HPP
