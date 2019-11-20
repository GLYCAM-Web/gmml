// Created by: Dave Montgomery

#ifndef PDBMASTERCARD_HPP
#define PDBMASTERCARD_HPP

#include <string>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbMasterCard
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbMasterCard();
            /*! \fn
              * A constructor that get a stream block of master card and parse the whole block to fill the related fields
              * @param stream_block A whole block of master card in a pdb file
              */
            PdbMasterCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a master card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the NumRemark in a master card
              * @return num_remark_ attribute of the current object of this class
              */
            int GetNumRemark();
            /*! \fn
              * An accessor function in order to access to the NumHet in a master card
              * @return num_het_ attribute of the current object of this class
              */
            int GetNumHet();
            /*! \fn
              * An accessor function in order to access to the NumHelix in a master card
              * @return num_helix_ attribute of the current object of this class
              */
            int GetNumHelix();
            /*! \fn
              * An accessor function in order to access to the NumSheet in a master card
              * @return num_sheet_ attribute of the current object of this class
              */
            int GetNumSheet();
            /*! \fn
              * An accessor function in order to access to the NumSite in a master card
              * @return num_site_ attribute of the current object of this class
              */
            int GetNumSite();
            /*! \fn
              * An accessor function in order to access to the NumXForm in a master card
              * @return num_x_form_ attribute of the current object of this class
              */
            int GetNumXForm();
            /*! \fn
              * An accessor function in order to access to the NumCoord in a master card
              * @return num_coord_ attribute of the current object of this class
              */
            int GetNumCoord();
            /*! \fn
              * An accessor function in order to access to the NumTer in a master card
              * @return num_ter_ attribute of the current object of this class
              */
            int GetNumTer();
            /*! \fn
              * An accessor function in order to access to the NumConnect in a master card
              * @return num_connect_ attribute of the current object of this class
              */
            int GetNumConnect();
            /*! \fn
              * An accessor function in order to access to the NumSeq in a master card
              * @return num_seq_ attribute of the current object of this class
              */
            int GetNumSeq();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current master card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the NumRemark attribute of the current object
              * Set the num_remark_ attribute of the current master card
              * @param num_remark The NumRemark attribute of the current object
              */
            void SetNumRemark(int num_remark);
            /*! \fn
              * A mutator function in order to set the NumHet of the current object
              * Set the num_het_ attribute of the current master card
              * @param num_het The NumHet attribute of the current object
              */
            void SetNumHet(int num_het);
            /*! \fn
              * A mutator function in order to set the NumHelix attribute of the current object
              * Set the num_helix_ attribute of the current master card
              * @param num_helix The NumHelix attribute of the current object
              */
            void SetNumHelix(int num_helix);
            /*! \fn
              * A mutator function in order to set the NumSheet of the current object
              * Set the num_sheet_ attribute of the current master card
              * @param num_sheet The NumSheet attribute of the current object
              */
            void SetNumSheet(int num_sheet);
            /*! \fn
              * A mutator function in order to set the NumSite of the current object
              * Set the num_site_ attribute of the current master card
              * @param num_site The NumSite attribute of the current object
              */
            void SetNumSite(int num_site);
            /*! \fn
              * A mutator function in order to set the NumXForm of the current object
              * Set the num_x_form_ attribute of the current master card
              * @param num_x_form The NumXForm attribute of the current object
              */
            void SetNumXForm(int num_x_form);
            /*! \fn
              * A mutator function in order to set the NumCoord of the current object
              * Set the num_coord_ attribute of the current master card
              * @param num_coord The NumCoord of the current object
              */
            void SetNumCoord(int num_coord);
            /*! \fn
              * A mutator function in order to set the NumTer of the current object
              * Set the num_ter_ attribute of the current master card
              * @param num_ter The NumTer of the current object
              */
            void SetNumTer(int num_ter);
            /*! \fn
              * A mutator function in order to set the NumConnect of the current object
              * Set the num_connect_ attribute of the current master card
              * @param num_connect The NumConnect of the current object
              */
            void SetNumConnect(int num_connect);
            /*! \fn
              * A mutator function in order to set the NumSeq of the current object
              * Set the num_seq_ attribute of the current master card
              * @param num_seq The NumSeq of the current object
              */
            void SetNumSeq(int num_seq);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the master card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Record name of Master card in a pdb file: "MASTER" */
            int num_remark_;                    /*!< Number of REMARK Cards >*/
            int num_het_;                       /*!< Number of HET Cards >*/
            int num_helix_;                     /*!< Number of HELIX Cards >*/
            int num_sheet_;                     /*!< Number of SHEET Cards >*/
            int num_site_;                      /*!< Number of SITE Cards >*/
            int num_x_form_;                     /*!< Number of ORIGX+SCALE+MTRIX Cards >*/
            int num_coord_;                     /*!< Number of ATOM+HETATM Cards >*/
            int num_ter_;                        /*!< Number of TER Cards >*/
            int num_connect_;                    /*!< Number of CONECT Cards >*/
            int num_seq_;                        /*!< Number of SEQRES Cards >*/

    };
}


#endif // PDBMASTERCARD_HPP
