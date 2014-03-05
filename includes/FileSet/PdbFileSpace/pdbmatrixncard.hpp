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
            typedef std::vector< std::vector <PdbMatrixN*> > MatrixNVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbMatrixNCard();
            PdbMatrixNCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the matrix n in a matrix n card
              * @return matrix_n_ attribute of the current object of this class
              */
            MatrixNVector GetMatrixN();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the matrix n attribute of the current object
              * Set the matrix_n_ attribute of the current matrix n card
              * @param a The matrix_n attribute of the current object
              */
            void SetMatrixN(MatrixNVector matrix_n);
            /*! \fn
              * A function in order to add the matrix n attribute to the current object
              * Set the matrix_ attribute of the current matrix n card
              * @param a The matrix attribute of the current object
              */
            void AddMatrixN(std::vector <PdbMatrixN*> matrix);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            MatrixNVector matrix_n_;

    };
}

#endif // PDBMATRIXNCARD_HPP
