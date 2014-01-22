#ifndef PDBMATRIXNCARD_HPP
#define PDBMATRIXNCARD_HPP

#include <string>
#include <vector>

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
            PdbMatrixNCard();
            PdbMatrixNCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            MatrixNVector GetMatrixN();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetMatrixN(const MatrixNVector matrix_n);
            void AddMatrixN(std::vector <PdbMatrixN*> matrix);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            MatrixNVector matrix_n_;

    };
}

#endif // PDBMATRIXNCARD_HPP
