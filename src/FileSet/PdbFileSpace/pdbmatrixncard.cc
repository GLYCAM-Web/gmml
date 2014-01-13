#include "../../../includes/FileSet/PdbFileSpace/pdbmatrixn.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbmatrixncard.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbMatrixNCard::PdbMatrixNCard() {}
PdbMatrixNCard::PdbMatrixNCard(stringstream &stream_block)
{
    string line;
    stringstream ss;
    getline(stream_block, line);
    line = Trim(line);
    while (!Trim(line).empty())
    {
        ss << line;
        PdbMatrixN* matrix = new PdbMatrixN(ss);
        AddMatrixN(matrix);
        getline(stream_block, line);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbMatrixNCard::MatrixNVector PdbMatrixNCard::GetMatrixN(){
    return matrix_n_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbMatrixNCard::SetMatrixN(const MatrixNVector matrix_n){
    matrix_n_ = matrix_n;
}

void PdbMatrixNCard::AddMatrixN(PdbMatrixN *matrix)
{
    matrix_n_.push_back(matrix);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////



