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
    getline(stream_block, line);
    line = Trim(line);
    vector <PdbMatrixN*> matrix_1, matrix_2, matrix_3;
    while (!Trim(line).empty())
    {
        int index = ConvertString<int>(line.substr(5, 1));
        PdbMatrixN* matrix = new PdbMatrixN(line);
        switch (index) {
        case 1 :
           matrix_1.push_back(matrix);
           break;
        case 2:
            matrix_2.push_back(matrix);
            break;
        case 3:
            matrix_3.push_back(matrix);
            break;
        }
        getline(stream_block, line);
    }
    this->AddMatrixN(matrix_1);
    this->AddMatrixN(matrix_2);
    this->AddMatrixN(matrix_3);
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

void PdbMatrixNCard::AddMatrixN(vector <PdbMatrixN*> matrix)
{
    matrix_n_.push_back(matrix);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////



