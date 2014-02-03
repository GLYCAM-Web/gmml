// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBHELIX_HPP
#define PDBHELIX_HPP

#include <string>
#include <vector>

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
        POLYPROLINE = 10
    };

    class PdbHelixResidue;
    class PdbHelix
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<PdbHelixResidue*> HelixResidueVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbHelix();
            /*! \fn
              * Constructor with required parameters
              * @param helix_id
              * @param helix_serial_number
              * @param helix_residues
              * @param helix_class
              * @param comment
              * @param helix_length
              */
            PdbHelix(const std::string& helix_id, int helix_serial_number, HelixResidueVector helix_residues,
                     PdbHelixClass helix_class, const std::string& comment, double helix_length);
            PdbHelix(std::stringstream& specification_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
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

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
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
            /*! \fn
              * A mutator function in order to set the helix class of the current object
              * Set the helix_class_ attribute of the current helix
              * @param helix_class The helix class of the current object
              */
            void SetHelixClass(PdbHelixClass helix_class);
            /*! \fn
              * A mutator function in order to set the comment of the current object
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

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string helix_id_;
            int helix_serial_number_;
            HelixResidueVector helix_residues_;
            PdbHelixClass helix_class_;
            std::string comment_;
            double helix_length_;
    };
}

#endif // PDBHELIX_HPP
