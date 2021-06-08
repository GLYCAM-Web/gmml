// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBCRYSTALLOGRAPHICCARD_HPP
#define PDBCRYSTALLOGRAPHICCARD_HPP

#include <string>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbCrystallographicCard
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbCrystallographicCard();
            /*! \fn
              * A constructor that get a stream block of crystallography card and parse the whole block to fill the related fields
              * @param stream_block A whole block of crystallography card in a pdb file
              */
            PdbCrystallographicCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a crystallographic card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the A in a crystallographic card
              * @return a_ attribute of the current object of this class
              */
            double GetA();
            /*! \fn
              * An accessor function in order to access to the B in a crystallographic card
              * @return b_ attribute of the current object of this class
              */
            double GetB();
            /*! \fn
              * An accessor function in order to access to the C in a crystallographic card
              * @return c_ attribute of the current object of this class
              */
            double GetC();
            /*! \fn
              * An accessor function in order to access to the alpha in a crystallographic card
              * @return alpha_ attribute of the current object of this class
              */
            double GetAlpha();
            /*! \fn
              * An accessor function in order to access to the beta in a crystallographic card
              * @return beta_ attribute of the current object of this class
              */
            double GetBeta();
            /*! \fn
              * An accessor function in order to access to the gamma in a crystallographic card
              * @return gamma_ attribute of the current object of this class
              */
            double GetGamma();
            /*! \fn
              * An accessor function in order to access to the space group in a crystallographic card
              * @return space_group_ attribute of the current object of this class
              */
            std::string GetSpaceGroup();
            /*! \fn
              * An accessor function in order to access to the Z value in a crystallographic card
              * @return z_value_ attribute of the current object of this class
              */
            int GetZValue();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current crystallographic card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the A attribute of the current object
              * Set the a_ attribute of the current crystallographic card
              * @param a The A attribute of the current object
              */
            void SetA(double a);        
            /*! \fn
              * A mutator function in order to set the B of the current object
              * Set the b_ attribute of the current crystallographic card
              * @param b The B attribute of the current object
              */
            void SetB(double b);
            /*! \fn
              * A mutator function in order to set the C attribute of the current object
              * Set the c_ attribute of the current crystallographic card
              * @param c The C attribute of the current object
              */
            void SetC(double c);
            /*! \fn
              * A mutator function in order to set the alpha of the current object
              * Set the alpha_ attribute of the current crystallographic card
              * @param alpha The alpha attribute of the current object
              */
            void SetAlpha(double alpha);
            /*! \fn
              * A mutator function in order to set the beta of the current object
              * Set the beta_ attribute of the current crystallographic card
              * @param beta The beta attribute of the current object
              */
            void SetBeta(double beta);
            /*! \fn
              * A mutator function in order to set the gamma of the current object
              * Set the gamma_ attribute of the current crystallographic card
              * @param gamma The gamma attribute of the current object
              */
            void SetGamma(double gamma);
            /*! \fn
              * A mutator function in order to set the space group of the current object
              * Set the space_group_ attribute of the current crystallographic card
              * @param space_group The space group of the current object
              */
            void SetSpaceGroup(const std::string space_group);
            /*! \fn
              * A mutator function in order to set the z_value of the current object
              * Set the z_value_ attribute of the current crystallographic card
              * @param z_value The z value of the current object
              */
            void SetZValue(int z_value);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the crystallography card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Record name of crystallography card in a pdb file: "CRYST1" */
            double a_;                          /*!< a >*/
            double b_;                          /*!< b >*/
            double c_;                          /*!< c >*/
            double alpha_;                      /*!< alpha >*/
            double beta_;                       /*!< beta >*/
            double gamma_;                      /*!< gamma >*/
            std::string space_group_;           /*!< space group >*/
            int z_value_;                       /*!< z value >*/

    };
}


#endif // PDBCRYSTALLOGRAPHICCARD_HPP
