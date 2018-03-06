#include "../../../includes/InputSet/PdbFileSpace/pdbmatrixnsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmatrixncard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbMatrixNSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbMatrixNSection::PdbMatrixNSection() {}
PdbMatrixNSection::PdbMatrixNSection(std::stringstream &stream_block)
{
    std::string line;
    getline(stream_block, line);
    std::string temp = line;
    MatrixNVector matrix_1, matrix_2, matrix_3;
    while (!gmml::Trim(temp).empty())
    {
        int index;
        if(line.substr(5,1) == " ")
            index = gmml::iNotSet;
        else
            index = gmml::ConvertString<int>(line.substr(5, 1));
        PdbMatrixNCard* matrix = new PdbMatrixNCard(line);
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
            case gmml::iNotSet:
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

PdbMatrixNSection::MatrixNVectorVector PdbMatrixNSection::GetMatrixN(){
    return matrix_n_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbMatrixNSection::SetMatrixN(MatrixNVector matrix_n){
    matrix_n_.clear();
    for(MatrixNVector::iterator it = matrix_n.begin(); it != matrix_n.end(); it++)
    {
        matrix_n.push_back(*it);
    }
}

void PdbMatrixNSection::AddMatrixN(MatrixNVector matrix)
{
    matrix_n_.push_back(matrix);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbMatrixNSection::Print(std::ostream &out)
{
    for(unsigned int i = 0; i < matrix_n_.size(); i++)
    {
        for(MatrixNVector::iterator it = matrix_n_.at(i).begin(); it != matrix_n_.at(i).end(); it++)
            (*it)->Print(out);
    }
    out << std::endl;
}
