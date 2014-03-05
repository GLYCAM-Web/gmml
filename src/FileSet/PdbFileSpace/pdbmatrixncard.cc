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
    string temp = line;
    vector <PdbMatrixN*> matrix_1, matrix_2, matrix_3;
    while (!Trim(temp).empty())
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
        temp = line;
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
void PdbMatrixNCard::SetMatrixN(MatrixNVector matrix_n){
    matrix_n_.clear();
    for(MatrixNVector::iterator it = matrix_n.begin(); it != matrix_n.end(); it++)
    {
        matrix_n.push_back(*it);
    }
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
void PdbMatrixNCard::Print(ostream &out)
{
    for(unsigned int i = 0; i < matrix_n_.size(); i++)
    {
        for(vector<PdbMatrixN*>::iterator it = matrix_n_.at(i).begin(); it != matrix_n_.at(i).end(); it++)
            (*it)->Print(out);
    }
    out << endl;
}
