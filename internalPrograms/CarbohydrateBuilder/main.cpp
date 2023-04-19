#include "includes/CodeUtils/strings.hpp"
#include "includes/CentralDataStructure/InternalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h> // stat

int main(int argc, char** argv)
{
    if (argc != 5)
    {
        std::cerr << "Usage:   ./buildOligosaccharideLibrary <List IDs and sequences> <Char delimiter used in list> <Folder to put outputs> <Prepfile>\n";
        std::cerr << "Example: ./buildOligosaccharideLibrary exampleLibrary.txt _ outputs/ ../../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep\n";
        std::cerr << "Don't use a delimiter that appears in glycan sequences or ids. Like - or , or [] etc\n";
        std::exit(EXIT_FAILURE); 
    }
    //Convert command line inputs to legible variables
    std::ifstream infile(argv[1]);
    char delimiter = argv[2][0]; // The second [0] gets me the first element of the argv which is type char**
    std::string outputFolderName = argv[3];
    std::string prepFile = argv[4];
    struct stat info;
    if( stat(argv[3], &info) != 0)
    {
        std::cerr << "Folder " << outputFolderName << "/ does not exist and it isn't my job to make it.\n";
        std::exit(EXIT_FAILURE);  
    }
    std::string line;
    while (std::getline(infile, line))
    {
        std::vector<std::string> splitLine = codeUtils::split(line, delimiter);
        std::string inputSequence = splitLine.at(1);
        std::cout << "\n*********************\nBuilding " << inputSequence << "\n*********************\n";
        try // ToDo Ok fix this: CondensedSequence::carbohydrateBuilder the sequenceManipulator constructor can throw.
        {
            cdsCondensedSequence::Carbohydrate carbBuilder(inputSequence, prepFile);
            if (carbBuilder.IsStatusOk())
            {
                //carbBuilder.Print();
                std::string inputGlycanID = splitLine.at(0);
                //carbBuilder.GenerateSingle3DStructureDefaultFiles(outputFolderName, inputGlycanID);
                //CondensedSequence::carbohydrateBuilder carbBuilder(inputSequence, prepFile);

                carbBuilder.Generate3DStructureFiles(outputFolderName, inputGlycanID);
                if (!carbBuilder.IsStatusOk()) // This is bad. Fix me once gems can catch what the carbBuilder throws.
                {
                    std::cerr << "Error thrown by the carbohydrateBuilder in gmml during 3D structure generation was: " << carbBuilder.GetStatusMessage() << std::endl;
                }
            }
            else
            {
                std::cerr << "Error thrown by the carbohydrateBuilder in gmml during construction was: " << carbBuilder.GetStatusMessage() << std::endl;
            }
        }
        catch(const std::string &exceptionMessage)
        {
            gmml::log(__LINE__, __FILE__, gmml::ERR, "carbohydrateBuilder class caught this exception message: " + exceptionMessage);
            std::cerr << "Error thrown by the carbohydrateBuilder in gmml during construction was: " << exceptionMessage << std::endl;
        }
        catch (const std::runtime_error &error)
        {
            gmml::log(__LINE__, __FILE__, gmml::ERR, error.what());
            std::cerr << "Error thrown by the carbohydrateBuilder in gmml during construction was: " << error.what() << std::endl;
        }
        catch (...)
        {
            gmml::log(__LINE__, __FILE__, gmml::ERR, "carbohydrateBuilder class caught a throw that was not anticipated. Curious. Death cometh?");
            std::cerr << "ERROR carbohydrateBuilder caught a throw type that was not anticipated. Pretty please report how you got to this to glycam@gmail.com.";
        }
    }
}
