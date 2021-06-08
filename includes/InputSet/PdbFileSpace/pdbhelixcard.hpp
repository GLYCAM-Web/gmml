// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBHELIXCARD_HPP
#define PDBHELIXCARD_HPP

#include <string>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
/*! \enum
  * Pdb helix class enumerator
  */
    enum PdbHelixClass
    {
        RIGHT_HANDED_ALPHA = 1,
        RIGHT_HANDED_OMEGA = 2,
        RIGHT_HANDED_PI = 3,
        RIGHT_HANDED_GAMMA = 4,
        RIGHT_HANDED_310 = 5,
        LEFT_HANDED_ALPHA = 6,
        LEFT_HANDED_OMEGA_ = 7,
        LEFT_HANDED_GAMMA_ = 8,
        RIBBON_27 = 9,
        POLYPROLINE = 10,
        UnknownHelix
    };

    class PdbHelixResidue;
    class PdbHelixCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of helix residues
              */
            typedef std::vector<PdbHelixResidue*> HelixResidueVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbHelixCard();
            /*! \fn
              * Constructor with required parameters
              * @param helix_id Helix identifier
              * @param helix_serial_number Serial number of the helix
              * @param helix_residues List of residues in the helix
              * @param helix_class Classification of the helix
              * @param comment Comment
              * @param helix_length Length of the helix
              */
            PdbHelixCard(const std::string& helix_id, int helix_serial_number, HelixResidueVector helix_residues,
                     PdbHelixClass helix_class, const std::string& comment, double helix_length);
            /*! \fnhelix_cards_
              * Constructor with required parameters
              * @param stream_block
              */
            PdbHelixCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the helix id in a helix class
              * @return helix_id_ attribute of the current object of this class
              */
            std::string GetHelixId();
            /*! \fn
              * An accessor function in order to access to the helix serial number in a helix class
              * @return helix_serial_number_ attribute of the current object of this class
              */
            int GetHelixSerialNumber();
            /*! \fn
              * An accessor function in order to access to the helix residues in a helix class
              * @return helix_residues_ attribute of the current object of this class
              */
            HelixResidueVector GetHelixResidues();
            /*! \fn
              * An accessor function in order to access to the helix class in a helix class
              * @return helix_class_ attribute of the current object of this class
              */
            PdbHelixClass GetHelixClass();
            /*! \fn
              * An accessor function in order to access to the comment in a helix class
              * @return comment_ attribute of the current object of this class
              */
            std::string GetComment();
            /*! \fn
              * An accessor function in order to access to the helix length in a helix class
              * @return helix_length_ attribute of the current object of this class
              */
            double GetHelixLength();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the helix id of the current object
              * Set the helix_id_ attribute of the current helix
              * @param helix_id The helix id of the current object
              */
            void SetHelixId(const std::string helix_id);
            /*! \fn
              * A mutator function in order to set the helix serial number of the current object
              * Set the helix_serial_number_ attribute of the current helix
              * @param helix_serial_number The helix serial number of the current object
              */
            void SetHelixSerialNumber(int helix_serial_number);
            /*! \fn
              * A mutator function in order to set the helix residues of the current object
              * Set the helix_residues_ attribute of the current helix
              * @param helix_residues The helix residues of the current object
              */
            void SetHelixResidues(const HelixResidueVector helix_residues);
            /*! \fn
              * A function in order to add the helix residue to the current object
              * Set the helix_residue_ attribute of the current helix
              * @param helix_residue The helix residue of the current object
              */
            void AddHelixResidue(PdbHelixResidue* helix_residue);
            /*! \fnhelix_cards_
              * A mutator function in order to set the helix class of the current object
              * Set the helix_class_ attribute of the current helix
              * @param helix_class The helix class of the current object
              */
            void SetHelixClass(PdbHelixClass helix_class);
            /*! \fn
              * A mutator function in order to sehelix_cards_t the comment of the current object
              * Set the comment_ attribute of the current helix
              * @param comment The commentd of the current object
              */
            void SetComment(const std::string& comment);
            /*! \fn
              * A mutator function in order to set the helix length of the current object
              * Set the helix_length_ attribute of the current helix
              * @param helix_length The helix length of the current object
              */
            void SetHelixLength(double helix_length);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the helix contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string helix_id_;                  /*!< Helix identifier >*/
            int helix_serial_number_;               /*!< Helix serial number >*/
            HelixResidueVector helix_residues_;     /*!< List of helix residues >*/
            PdbHelixClass helix_class_;             /*!< Classification of the helix >*/
            std::string comment_;                   /*!< Comment >*/
            double helix_length_;                   /*!< Length of the helix >*/
    };
}

#endif // PDBHELIXCARD_HPP
