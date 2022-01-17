#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_REMARKRECORD_HPP
#define GMML_INCLUDES_INPUTSET_PDBFILE_REMARKRECORD_HPP

#include <string>
#include <iostream>

namespace pdb
{
    class RemarkRecord
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            RemarkRecord();
            RemarkRecord(std::stringstream& stream_block);
            RemarkRecord(const std::string record_name, const std::string remark_cards);
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            inline const std::string& GetRecordName() const {return record_name_;}
            inline const std::string& GetRemarks() const {return remark_cards_;}
            inline const float& GetResolution() const {return resolution_;}
            inline const float& GetBFactor() const {return b_factor_;}
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cerr) const;
        private:
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetRecordName(const std::string record_name);
            void SetRemarks(const std::string remark_cards);
            void SetResolution(const float resolution);
            void SetBFactor(const float b_factor);
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Name of remark card record which is in the first column of each line of a pdb file >*/
            std::string remark_cards_;          /*!< Remarks that are in remark card of a pdb file >*/
            float resolution_;            /*!< Resolution of PDB >*/
            float b_factor_;              /*!< B Factor of PDB >*/
    };
}

#endif // GMML_INCLUDES_INPUTSET_PDBFILE_REMARKRECORD_HPP
