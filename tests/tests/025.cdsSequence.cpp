#include "includes/CentralDataStructure/InternalPrograms/Sequence/sequence.hpp"
#include <iostream>

using CondensedSequence::Sequence;

int main(int argc, char* argv[])
{
    std::string inputSequence = argv[1];
    try
    {
        Sequence sequence(inputSequence);
        std::cout << "Input:\n" << sequence.getInputSequence() << "\n";
        std::cout << "Interpreted:\n" << sequence.getInterpretedSequence() << "\n";
        std::cout << "IndexOrdered:\n" << sequence.getIndexOrdered() << "\n";
        std::cout << "IndexOrderedLabeled:\n" << sequence.getIndexLabeled() << "\n";
    }
    catch (std::runtime_error& error)
    {
        std::cout << "Test level caught this: " << error.what();
    }
    catch (...)
    {
        std::cout << "Test level caught unknown error\n";
    }
    return 0;
}
