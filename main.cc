#include <iostream>

#define LibFileTest

#ifdef ParameterFileTest
#include "includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
using namespace ParameterFileSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        ParameterFile* temp = new ParameterFile(argv[1]);
        temp->Print(std::cout);
    }
    return 0;
}
#endif

#ifdef PrepFileTest
#include "includes/ParameterSet/PrepFileSpace/prepfile.hpp"
using namespace PrepFileSpace;
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

#ifdef LibFileTest
#include "includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
using namespace LibraryFileSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        LibraryFile* temp = new LibraryFile(argv[1]);
        temp->Print(std::cout);
    }
    return 0;
}
#endif

#ifdef CrdFileTest
#include "../../../includes/FileSet/CoordinateFileSpace/coordinatefile.hpp"
using namespace CoordinateFileSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        CoordinateFile* temp = new CoordinateFile(argv[1]);
        temp->Print(std::cout);
    }
    return 0;
}
#endif
