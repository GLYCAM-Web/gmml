#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_JOURNALRECORD_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_JOURNALRECORD_HPP

#include <string>
#include <iostream>
#include <vector>

namespace pdb
{
    class JournalRecord
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        JournalRecord();
        JournalRecord(std::stringstream& stream_block);

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline const std::string& GetRecordName() const
        {
            return record_name_;
        }

        inline const std::vector<std::string>& GetAuthors() const
        {
            return authors_;
        }

        inline const std::string& GetTitle() const
        {
            return title_;
        }

        inline const std::vector<std::string>& GetEditors() const
        {
            return editors_;
        }

        inline const std::string& GetReference() const
        {
            return reference_;
        }

        inline const std::string& GetPublisher() const
        {
            return publisher_;
        }

        inline const std::vector<std::string>& GetReferenceNumbers() const
        {
            return reference_nums_;
        }

        inline const std::string& GetPMID() const
        {
            return pmid_;
        }

        inline const std::string& GetDOI() const
        {
            return doi_;
        }

        inline const std::string& GetText() const
        {
            return text_;
        }

        //////////////////////////////////////////////////////////
        //                       DISPLAY FUNCTION               //
        //////////////////////////////////////////////////////////
        void Print(std::ostream& out = std::cerr) const;
        void Write(std::ostream& stream) const;

      private:
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        void SetRecordName(const std::string record_name);
        void SetAuthors(std::vector<std::string> authors);
        void SetTitle(const std::string title);
        void SetEditors(std::vector<std::string> editors);
        void SetReference(const std::string reference);
        void SetPublisher(const std::string publisher);
        void SetReferenceNumbers(std::vector<std::string> reference_nums);
        void SetPMID(const std::string pmid);
        void SetDOI(const std::string doi);
        void SetText(const std::string text);
        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        std::string record_name_ = ""; /*!< Record name which appears in the first column of each line in a pdb file >*/
        std::vector<std::string> authors_; /*!< List of authors that appear in Journal record of a pdb file >*/
        std::string title_ = "";           /*!< Title that appears in Journal record of a pdb file >*/
        std::vector<std::string> editors_; /*!< List of editors that appear in Journal record of a pdb file >*/
        std::string reference_ = "";       /*!< Reference that appears in Journal record of a pdb file >*/
        std::string publisher_ = "";       /*!< Publisher that appears in Journal record of a pdb file >*/
        std::vector<std::string>
            reference_nums_;    /*!< List of reference numbers that appear in Journal record of a pdb file >*/
        std::string pmid_ = ""; /*!< Pub Med ID number that appears in Journal record of a pdb file >*/
        std::string doi_  = ""; /*!< DOI number that appears in Journal record of a pdb file >*/
        std::string text_ = ""; /*!< Text in a Journal Section >*/
    };
} // namespace pdb
#endif // GMML_INCLUDES_INPUTSET_PDBFILE_JOURNALRECORD_HPP
