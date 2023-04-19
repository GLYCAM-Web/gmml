#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_HEADERRECORD_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_HEADERRECORD_HPP

#include <string>
#include <sstream>
#include <iostream>

namespace pdb
{
    class HeaderRecord
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            HeaderRecord();
            HeaderRecord(const std::string& record_name, const std::string& classification, const std::string& deposition_date, const std::string& identifier_code);
            HeaderRecord(std::stringstream& stream_block);
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetRecordName() const;
            std::string GetClassification() const;
            std::string GetDepositionDate() const;
            std::string GetIdentifierCode() const;
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////            
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
            void SetClassification(const std::string classification);
            void SetDepositionDate(const std::string deposition_date);
            void SetIdentificationCode(const std::string identifier_code);
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Record name of headr card in a pdb file >*/
            std::string classification_;        /*!< Classification of the pdb file >*/
            std::string deposition_date_;       /*!< Date of deposition >*/
            std::string identifier_code_;       /*!< Identifier code of the pdb file >*/
    };
}
#endif // GMML_INCLUDES_INPUTSET_PDBFILE_HEADERRECORD_HPP
