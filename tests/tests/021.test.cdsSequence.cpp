#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/CentralDataStructure/Writers/cdsOffWriter.hpp"
#include "includes/CodeUtils/logging.hpp"
#include <iostream>

int main ()
{	
    std::string s1 = "DTalp[2S,3Me]a1-6DManpa1-6[DAllpb1-3][DNeup5Aca2-6DGalpb1-4DGlcp[3S]b1-2DAltpa1-4]DManpb1-4DGulp[6Me]b1-4DGlcpNAcb1-OH";
    std::string s2 = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
    std::string s3 = "DGlcpNAcb1-4DGlcpAb1-4DGlcpAb1-3DGalpb1-3DGalpb1-4DXylpb1-OH";
    std::string s4 = "dUA[2S]1-4DGlcpNAc[3S,6S]a1-4LIdopA(2SO)[2S]a1-4DGlcpNSa1-4DGlcpA[2S]b1-4DGlcpAb1-3DGalpb1-3DGalpb1-4DXylpb1-OH";
    std::string s5 = "dUA[2S]a1-4DGlcpNSa1-4LIdopA[2S]a1-4DGlcpNSa1-4LIdopA(4C1)a1-4DGlcpNS[6S]a1-OH";
    std::string s6 = "DGlcpa1-2DFrufb";
    std::string s7 = "DFrufb2-1DGlcpa";
    std::string s8 = "DNeup5Ac&Label=residue-9;a2-6&Label=link-7;DGalp&Label=residue-8;b1-4&Label=link-6;DGlcpNAc&Label=residue-6;[3&Label=link-5;S&Label=residue-7;]b1-2&Label=link-4;DManp&Label=residue-5;a1-3&Label=link-3;[DGlcpNAc&Label=residue-10;b1-4&Label=link-8;][DManp&Label=residue-12;[2&Label=link-11;S&Label=residue-13;,3&Label=link-12;Me&Label=residue-14;]a1-6&Label=link-10;DManp&Label=residue-11;a1-6&Label=link-9;]DManp&Label=residue-4;b1-4&Label=link-2;DGlcpNAc&Label=residue-3;[6&Label=link-13;Me&Label=residue-15;]b1-4&Label=link-1;DGlcpNAc&Label=residue-2;b1-1&Label=link-0;-OH&Label=residue-1;";
    std::string s9 = "";
    std::string s10 = "There will be cake.";
    std::string s11 = "DGlcpNAcb1-4DGlcpAb1-4DGlcpAb1-3DGalpb1-3DGalpb1-4DXylpb1-OH ";
    std::string s12 = "DGlcpNAcb1-4DGlcpAb1-4DGlcpAb1-3DGalpb1-3]DGalpb1-4DXylpb1-OH";
    std::string s13 = "DGlcpNAcb1-4DGlcpAb1-4DGlcpAb1-3DGalp[Boo]b1-3DGalpb1-4DXylpb1-OH";
    std::string s14 = "dUA[2S]1-4DGlcpNAc[3S,6S]a1-4LIdopA(2SO)[2S]a1-4LIdopA(2SO)a1-4DGlcpNSa1-4DGlcpA[2S]b1-OH";
    std::string s15 = "DGlNAcb1-OH";
    std::string s16 = "DManpa1-4DManpa1-4DManpa1-4DManpa1-4DManpa1-4DManp[6D]a1-4DManp[2S,6S]a1-4DManpa1-OME";
    std::string s17 = "DManpa[6S,2S]1-OME";
    //std::vector<std::string> sequences {s1, s2, s3, s4, s5, s6, s7};
   // std::vector<std::string> sequences {s17};
    std::vector<std::string> sequences {s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17};
//    std::string prepFilePath = "/programs/gems/gmml/dat/prep/GLYCAM_06j-1_GAGS.prep";
    std::string prepFilePath = "../dat/prep/GLYCAM_06j-1_GAGS.prep";
    int loopCounter = 0;
    std::cout << "-----------------------------------------------------------------------------------------------------\n\n";
    for (auto &sequence : sequences)
    {	
        try
        {
            loopCounter++;
            std::cout << "Sequence " << loopCounter << ": " << sequence << "\n";
            //MolecularModeling::Assembly ass(sequence, prepFilePath); // WARNING. Just a test. 3D structures are not correct.
            cdsCondensedSequence::Carbohydrate theVanToMordor(sequence, prepFilePath);
            if (theVanToMordor.GetStatusType() != "OK")
            { // This is dumb, but I haven't quite finalized how to handle the error throwing in e.g. glycoprotein builder using the carb builder vs here.
                throw std::runtime_error(theVanToMordor.GetStatusMessage());
            }
            std::cout << "Structure created without throwing an exception for: " << sequence << "\n\n";
            for(auto &residue : theVanToMordor.getResidues())
            {
                std::cout << "Residue:\n";
                residue->Print(std::cout);
                std::cout << "\n";
            }
            std::cout << "*****\n\n\n";
            // Outputs for fun:

            theVanToMordor.setName("HiMyNameIs " + sequence);
            std::string fileName = "./021.sequenceAsPdbFile.pdb";
            std::ofstream outFileStream;
            try
            {
                std::cout << "Writing: " << fileName << "\n";
                outFileStream.open(fileName.c_str());
                theVanToMordor.WritePdb(outFileStream);
                outFileStream.close();
            }
            catch(...)
            {
                gmml::log(__LINE__,__FILE__,gmml::ERR, "Error when writing pdbFile class to file:\n" + fileName);
                throw std::runtime_error("Error when writing pdbFile class to file:\n" + fileName);
            }
            //OFF molecule
            fileName = "./021.sequenceAsOffFile.off";
            try
            {

                std::cout << "Writing: " << fileName << "\n";
                std::ofstream outFileStream;
                outFileStream.open(fileName.c_str());
                cds::WriteMoleculeToOffFile(theVanToMordor.getResidues(), outFileStream, theVanToMordor.getName());
                outFileStream.close();
            }
            catch(...)
            {
                gmml::log(__LINE__,__FILE__,gmml::ERR, "Error when writing to file:\n" + fileName);
                throw std::runtime_error("Error when writing to file:\n" + fileName);
            }
            std::cout << "Beep beep the van has arrived\n";
        }
        catch (const std::string &exception)
        {
            std::cerr << "Test level caught error: " << exception << std::endl;
        }
        catch (const std::runtime_error &error)
        {
            std::cerr << "Test level caught runtime error: " << error.what() << std::endl;
        }
        std::cout << "-----------------------------------------------------------------------------------------------------\n\n";
	}
	return 0;
}

