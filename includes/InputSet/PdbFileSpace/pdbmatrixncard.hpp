// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBMATRIXNCARD_HPP
#define PDBMATRIXNCARD_HPP

#include <string>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    class PdbMatrixN;

    class PdbMatrixNCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of matrices
              */
            typedef std::vector<PdbMatrixN*> MatrixNVector;
            /*! \typedef
              * List of three lists of matrices
              */
            typedef std::vector<MatrixNVector> MatrixNVectorVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbMatrixNCard();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbMatrixNCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the matrix n in a matrix n card
              * @return matrix_n_ attribute of the current object of this class
              */
            MatrixNVectorVector GetMatrixN();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the matrix n attribute of the current object
              * Set the matrix_n_ attribute of the current matrix n card
              * @param matrix_n The matrix_n attribute of the current object
              */
            void SetMatrixN(MatrixNVector matrix_n);
            /*! \fn
              * A function in order to add the matrix n attribute to the current object
              * Set the matrix_ attribute of the current matrix n card
              * @param matrix The matrix attribute of the current object
              */
            void AddMatrixN(MatrixNVector matrix);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the matrix n card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            MatrixNVectorVector matrix_n_;      /*!< List of three lists of matrices >*/

    };
}

#endif // PDBMATRIXNCARD_HPP
