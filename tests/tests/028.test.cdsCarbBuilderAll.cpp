#include "includes/CodeUtils/strings.hpp"
#include "includes/CentralDataStructure/InternalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"
#include "includes/CentralDataStructure/InternalPrograms/Sequence/sequence.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h> // stat

int main(int argc, char** argv)
{
    if (argc != 5)
    {
        std::cerr << "Usage:  " << argv[0]
                  << " <List IDs and sequences> <Char delimiter used in list> "
                     "<Folder to put outputs> <Prepfile>\n";
        std::cerr << "Example: " << argv[0]
                  << " exampleLibrary.txt _ outputs/ "
                     "../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep\n";
        std::cerr << "Don't use a delimiter that appears in glycan sequences or ids. Like - or , or [] etc\n";
        std::exit(EXIT_FAILURE);
    }
    // Convert command line inputs to legible variables
    std::ifstream infile(argv[1]);
    char delimiter = argv[2][0]; // The second [0] gets me the first element of the argv which is type char**
    std::string outputFolderName = argv[3];
    std::string prepFile         = argv[4];
    struct stat info;
    if (stat(argv[3], &info) != 0)
    {
        std::cerr << "Folder " << outputFolderName << "/ does not exist and it isn't my job to make it.\n";
        std::exit(EXIT_FAILURE);
    }
    std::string s1 = "DTalp[2S,3Me]a1-6DManpa1-6[DAllpb1-3][DNeup5Aca2-6DGalpb1-4DGlcp[3S]b1-2DAltpa1-4]DManpb1-4DGulp["
                     "6Me]b1-4DGlcpNAcb1-OH";
    std::string s2 = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-"
                     "4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
    std::string s3 = "DGlcpNAcb1-4DGlcpAb1-4DGlcpAb1-3DGalpb1-3DGalpb1-4DXylpb1-OH";
    std::string s4 = "dUA[2S]1-4DGlcpNAc[3S,6S]a1-4LIdopA(2SO)[2S]a1-4DGlcpNSa1-4DGlcpA[2S]b1-4DGlcpAb1-3DGalpb1-"
                     "3DGalpb1-4DXylpb1-OH";
    std::string s5 = "dUA[2S]a1-4DGlcpNSa1-4LIdopA[2S]a1-4DGlcpNSa1-4LIdopA(4C1)a1-4DGlcpNS[6S]a1-OH";
    std::string s6 = "DGlcpa1-2DFrufb";
    std::string s7 = "DFrufb2-1DGlcpa";
    std::string s8 =
        "DNeup5Ac&Label=residue-9;a2-6&Label=link-7;DGalp&Label=residue-8;b1-4&Label=link-6;DGlcpNAc&Label=residue-6;["
        "3&Label=link-5;S&Label=residue-7;]b1-2&Label=link-4;DManp&Label=residue-5;a1-3&Label=link-3;[DGlcpNAc&Label="
        "residue-10;b1-4&Label=link-8;][DManp&Label=residue-12;[2&Label=link-11;S&Label=residue-13;,3&Label=link-12;Me&"
        "Label=residue-14;]a1-6&Label=link-10;DManp&Label=residue-11;a1-6&Label=link-9;]DManp&Label=residue-4;b1-4&"
        "Label=link-2;DGlcpNAc&Label=residue-3;[6&Label=link-13;Me&Label=residue-15;]b1-4&Label=link-1;DGlcpNAc&Label="
        "residue-2;b1-1&Label=link-0;-OH&Label=residue-1;";
    std::string s9  = "";
    std::string s10 = "There will be cake.";
    std::string s11 = "DGlcpNAcb1-4DGlcpAb1-4DGlcpAb1-3DGalpb1-3DGalpb1-4DXylpb1-OH ";
    std::string s12 = "DGlcpNAcb1-4DGlcpAb1-4DGlcpAb1-3DGalpb1-3]DGalpb1-4DXylpb1-OH";
    std::string s13 = "DGlcpNAcb1-4DGlcpAb1-4DGlcpAb1-3DGalp[Boo]b1-3DGalpb1-4DXylpb1-OH";
    std::string s14 = "dUA[2S]1-4DGlcpNAc[3S,6S]a1-4LIdopA(2SO)[2S]a1-4LIdopA(2SO)a1-4DGlcpNSa1-4DGlcpA[2S]b1-OH";
    std::string s15 = "DGlNAcb1-OH";
    std::string s16 = "DManpa1-4DManpa1-4DManpa1-4DManpa1-4DManpa1-4DManp[6D]a1-4DManp[2S,6S]a1-4DManpa1-OME";
    std::string s17 = "DManpa[6S,2S]1-OME";
    std::string s18 = "DNeup5Aca2-3DManp[2S,6Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-"
                      "2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
    std::string line;
    while (std::getline(infile, line))
    {
        if (line.empty())
        {
            continue;
        }
        std::vector<std::string> splitLine = codeUtils::split(line, delimiter);
        if (splitLine.size() != 2)
        {
            std::cerr << "Encountered problem when splitting this line >>>" << line << "<<< from file >>>" << argv[1]
                      << "<<< into an ID and carb string separated by your specified delimiter: >>>" << argv[2]
                      << "<<<\n";
            std::exit(EXIT_FAILURE);
        }
        std::string inputGlycanID = splitLine.at(0);
        std::string inputSequence = splitLine.at(1);
        std::cout << "\n*********************\nBuilding id " << inputGlycanID << ":" << inputSequence
                  << "\n*********************\n";
        try
        { // I want to first build and reorder the sequence like gems does, and pass that to the carbohydrateBuilder
            CondensedSequence::Sequence sequence(inputSequence);
            std::cout << "Input:\n" << sequence.getInputSequence() << "\n";
            std::cout << "Interpreted:\n" << sequence.getInterpretedSequence() << "\n";
            std::cout << "IndexOrdered:\n" << sequence.getIndexOrdered() << "\n";
            std::cout << "IndexOrderedLabeled:\n" << sequence.getIndexLabeled() << "\n";
            std::cout << "Parsed and labelled " << inputGlycanID << " with no exceptions thrown.\n\n";
            cdsCondensedSequence::carbohydrateBuilder carbBuilder(sequence.getIndexOrdered(), prepFile);
            bool likelyShapesOnly = true; // Not sure I like this. Two functions probably more readable.
            std::cout << "Number of residues for this sequence is " << carbBuilder.GetCarbohydrate().GetResidueCount()
                      << "\n";
            std::cout << "Number of likely shapes is " << carbBuilder.GetNumberOfShapes(likelyShapesOnly) << "\n";
            std::cout << "Number of possible shapes is " << carbBuilder.GetNumberOfShapes() << "\n";
            for (auto& residue : carbBuilder.GetCarbohydrate().getResidues())
            {
                std::cout << residue->getStringId() << "\n";
            }
            for (auto& linkageInfo : carbBuilder.GenerateUserOptionsDataStruct())
            {
                std::cout << "Name: " << linkageInfo.linkageName_
                          << ", LinkageIndex: " << linkageInfo.indexOrderedLabel_
                          << ", Res1: " << linkageInfo.firstResidueNumber_
                          << ", Res2: " << linkageInfo.secondResidueNumber_ << "\n";
                for (auto& dihedralInfo : linkageInfo.likelyRotamers_)
                {
                    std::cout << "    LikelyRotamers: " << dihedralInfo.dihedralName_;
                    for (auto& rotamer : dihedralInfo.rotamers_)
                    {
                        std::cout << ", " << rotamer;
                    }
                    std::cout << "\n";
                }
            }
            carbBuilder.GenerateSingle3DStructureDefaultFiles(outputFolderName, inputGlycanID);
        }
        catch (const std::runtime_error& error)
        {
            gmml::log(__LINE__, __FILE__, gmml::ERR, error.what());
            std::cerr << "Error thrown by the carbohydrateBuilder in gmml during construction was: " << error.what()
                      << std::endl;
        }
        catch (...)
        {
            gmml::log(__LINE__, __FILE__, gmml::ERR,
                      "carbohydrateBuilder class caught a throw that was not anticipated. Curious. Death cometh?");
            std::cerr << "ERROR carbohydrateBuilder caught a throw type that was not anticipated. Pretty please report "
                         "how you got to this to glycam@gmail.com.";
        }
    }
}
