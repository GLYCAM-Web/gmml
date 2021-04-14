
#include <iostream>
#include "includes/InputSet/CondensedSequence/sequenceParser.hpp"
#include "includes/InputSet/CondensedSequence/sequenceManipulator.hpp"


int main ()
{	
    std::string s1 = "DTalp[2S,3Me]a1-6DManpa1-6[DAllpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcp[3S]b1-2DAltpa1-3]DManpb1-4DGulpNAc[6Me]b1-4DGlcpNAcb1-OH";
    std::string s2 = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH";
    std::string s3 = "DGlcpNAcb1-4DGlcpAb1-4DGlcpAb1-3DGalpb1-3DGalpb1-4DXylpb1-OH";
    std::string s4 = "dUA[2S]1-4DGlcpNAc[3S,6S]a1-4LIdopA(2SO)[2S]a1-4DGlcpNSa1-4DGlcpA[2S]b1-4DGlcpAb1-3DGalpb1-3DGalpb1-4DXylpb1-OH";
    std::string s5 = "dUA[2S]a1-4DGlcpNSa1-4LIdopA[2S]a1-4DGlcpNSa1-4LIdopA(4C1)a1-4DGlcpNS[6S]a1-OH";
    std::vector<std::string> sequences {s1, s2, s3, s4, s5};
    // Test sequenceParser
    for (auto &sequence : sequences)
    {
    	CondensedSequence::SequenceParser parser(sequence);
    	std::cout << parser.Print() << "\n\n";
    	std::cout << "About to go out of scope.\n";
	}
	// Now test sequenceManipulator
	for (auto &sequence : sequences)
    {
    	CondensedSequence::SequenceManipulator man(sequence);
    	std::cout << man.Print() << "\n\n";
    	man.ReorderSequence();
    	man.LabelSequence();
    	std::cout << "PrintLabelledSequence:" << std::endl;
    	man.PrintLabelledSequence();
    	std::cout << "About to go out of scope.\n";
	}
	return 0;
}

