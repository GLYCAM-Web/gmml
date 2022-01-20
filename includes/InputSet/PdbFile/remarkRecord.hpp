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
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
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
            void SetResolution(const float resolution);
            void SetBFactor(const float b_factor);
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            float resolution_;            /*!< Resolution of PDB >*/
            float b_factor_;              /*!< B Factor of PDB >*/
    };
}

#endif // GMML_INCLUDES_INPUTSET_PDBFILE_REMARKRECORD_HPP
