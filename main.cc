#include <iostream>
#include "includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfile.hpp"

#define PrepFileTest

using namespace ParameterFileSpace;
using namespace PrepFileSpace;

#ifdef ParameterFileTest
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
#endif

#ifdef PrepFileTest
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        PrepFile* temp = new PrepFile(argv[1]);
        temp->Print(std::cout);
    }
    return 0;
}
#endif
