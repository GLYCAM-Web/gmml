#include <iostream>
#include "includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"

using namespace ParameterFileSpace;


int main()
{
    ParameterFile* temp = new ParameterFile("dat/parm99.dat");
    //std::cout << "HERE";
    temp->Print(std::cout);
    return 0;
}
