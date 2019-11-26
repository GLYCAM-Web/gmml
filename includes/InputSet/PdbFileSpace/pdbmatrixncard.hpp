// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBMATRIXNCARD_HPP
#define PDBMATRIXNCARD_HPP

#include <string>
#include <iostream>

#include "../../GeometryTopology/coordinate.hpp"

namespace PdbFileSpace
{
    class PdbMatrixNCard
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbMatrixNCard();
            /*! \fn
              * Constructor with required parameters
              * @param line
              */
            PdbMatrixNCard(std::string& line);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a matrix n
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the N attribute in a matrix n
              * @return n_ attribute of the current object of this class
              */
            int GetN();
            /*! \fn
              * An accessor function in order to access to the serial number in a matrix n
              * @return serial_number_ attribute of the current object of this class
              */
            int GetSerialNumber();
            /*! \fn
              * An accessor function in order to access to the transformation vector in a matrix n
              * @return transformation_vector_ attribute of the current object of this class
              */
            GeometryTopology::Coordinate GetTransformationVector();
            /*! \fn
              * An accessor function in order to access to the V attribute in a matrix n
              * @return v_ attribute of the current object of this class
              */
            double GetV();
            /*! \fn
              * An accessor function in order to access to the i given attribute in a matrix n
              * @return i_given_ attribute of the current object of this class
              */
            int GetIGiven();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current matrix n
              * @param record_name The record name attribute of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the N attribute of the current object
              * Set the n_ attribute of the current matrix n
              * @param n The N attribute of the current object
              */
            void SetN(int n);
            /*! \fn
              * A mutator function in order to set the serial number of the current object
              * Set the serial_number_ attribute of the current matrix n
              * @param serial_number The serial number attribute of the current object
              */
            void SetSerialNumber(int serial_number);
            /*! \fn
              * A mutator function in order to set the transfomration vector of the current object
              * Set the transfomration_vector_ attribute of the current matrix n
              * @param transfomration_vector The transfomration vector attribute of the current object
              */
            void SetTransformationVector(GeometryTopology::Coordinate transfomration_vector);
            /*! \fn
              * A mutator function in order to set the V attribute of the current object
              * Set the v_ attribute of the current matrix n
              * @param v The V attribute of the current object
              */
            void SetV(double v);
            /*! \fn
              * A mutator function in order to set the i given of the current object
              * Set the i_given_ attribute of the current matrix n
              * @param i_given The i given attribute of the current object
              */
            void SetIGiven(int i_given);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the matrix n contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;                       /*!< Record name of matrix card which is the first column of each line of the card >*/
            int n_;                                         /*!< n >*/
            int serial_number_;                             /*!< Model serial number >*/
            GeometryTopology::Coordinate transfomration_vector_;    /*!< Transformation vector >*/
            double v_;                                      /*!< velocity >*/
            int i_given_;                                   /*!< i_given >*/
    };
}

#endif // PDBMATRIXNCARD_HPP
