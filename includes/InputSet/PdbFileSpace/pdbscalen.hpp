// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBSCALEN_HPP
#define PDBSCALEN_HPP

#include <string>
#include <sstream>
#include <iostream>

#include "../../GeometryTopology/coordinate.hpp"

namespace PdbFileSpace
{
    class PdbScaleN
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbScaleN();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbScaleN(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the record name in a scale n
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the N attribute in a scale n
              * @return n_ attribute of the current object of this class
              */
            int GetN();
            /*! \fn
              * An accessor function in order to access to the sclae vector in a scale n
              * @return scale_vector_ attribute of the current object of this class
              */
            GeometryTopology::Coordinate GetScaleVector();
            /*! \fn
              * An accessor function in order to access to the U attribute in a scale n
              * @return u_ attribute of the current object of this class
              */
            double GetU();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the record name attribute of the current object
              * Set the record_name_ attribute of the current scale n
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the scale N attribute of the current object
              * Set the n_ attribute of the current scale n card
              * @param n The N attribute of the current object
              */
            void SetN(int n);
            /*! \fn
              * A mutator function in order to set the scale vector of the current object
              * Set the scale_vector_ attribute of the current scale n card
              * @param scale_vector The scale vector attribute of the current object
              */
            void SetScaleVector(GeometryTopology::Coordinate scale_vector);
            /*! \fn
              * A mutator function in order to set the U attribute of the current object
              * Set the u_ attribute of the current scale n card
              * @param u The U attribute of the current object
              */
            void SetU(double u);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb scale n contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;               /*!< Name of scale record that appears in the first column of each line of a pdb file >*/
            int n_;                                 /*!< n >*/
            GeometryTopology::Coordinate scale_vector_;     /*!< Scale vector which is displayed in a coordinate >*/
            double u_;                              /*!< u value of a scale >*/
    };
}

#endif // PDBSCALEN_HPP
