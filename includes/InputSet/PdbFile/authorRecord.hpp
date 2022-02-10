#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_AUTHORRECORD_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_AUTHORRECORD_HPP

#include <string>
#include <sstream>
#include <iostream>

namespace pdb
{
    class AuthorRecord
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            AuthorRecord();
            AuthorRecord(std::string record_name, std::string author);
            AuthorRecord(std::stringstream& stream_block);
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            inline const std::string& GetRecordName() const {return record_name_;}
            inline const std::string& GetAuthor() const {return author_;}
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
            void SetAuthor(const std::string author);
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Record name which appears in the first column of each line in a pdb file >*/
            std::string author_;              /*!< Author that appears in KEYWORD record of a pdb file >*/
    };
}

#endif // GMML_INCLUDES_INPUTSET_PDBFILE_AUTHORRECORD_HPP
