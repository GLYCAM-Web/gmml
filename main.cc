#include <iostream>
#include "includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"

using namespace ParameterFileSpace;


int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        ParameterFile* temp = new ParameterFile(argv[1]);
        //std::cout << "HERE";
        temp->Print(std::cout);
    }
    return 0;
}
