#include "../../../../includes/CentralDataStructure/InternalPrograms/Sequence/sequence.hpp"

#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/CodeUtils/logging.hpp"

using cdsCondensedSequence::SequenceManipulator;
using CondensedSequence::Sequence;

Sequence::Sequence(std::string condensedSequence)
{
    try
    {
        SequenceManipulator manipulator(condensedSequence);
        this->setInputSequence(condensedSequence);
        this->setInterpretedSequence(manipulator.Print());
        this->setIndexOrdered(manipulator.ReorderSequence());
        bool withLabels = true;
        this->setIndexLabeled(manipulator.Print(withLabels));
    }
    catch (const std::string& exceptionMessage)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR,
                  "Sequence class constructor caught this exception message: " + exceptionMessage);
        this->SetStatus("ERROR", exceptionMessage);
    }
    catch (const std::runtime_error& error)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, error.what());
        this->SetStatus("ERROR", error.what());
    }
    catch (...)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR,
                  "Sequence class constructor caught a throw that was not anticipated. Death?");
        this->SetStatus("ERROR", "carbohydrateBuilder constructor caught a throw type that was not anticipated. Please "
                                 "report how you got to this to glycam@gmail.com.");
    }
}
