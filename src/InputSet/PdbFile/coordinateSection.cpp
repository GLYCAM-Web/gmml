#include "includes/InputSet/PdbFile/coordinateSection.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"
using pdb::CoordinateSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
CoordinateSection::CoordinateSection() {}

CoordinateSection::CoordinateSection(std::stringstream &stream_block)
{
    int currentModelNumber = 1;
   // std::string currentResidueId = "?_?_?"; // This is just for organizing into residues.
   // std::string previousResidueId = "IJustNeedThisToBeDifferentFromCurrentResidueIdForTheinitialCheckOkDontJudgeMeImDoingMyBestOutHere";
    std::string line;
    while(getline(stream_block, line))
    {
        std::string recordName = codeUtils::RemoveWhiteSpace(line.substr(0,6));
        // MODEL cards
        if(recordName == "MODEL")
        {
            try
            {
                currentModelNumber = std::stoi(codeUtils::RemoveWhiteSpace(line.substr(10,4)));
            }
            catch (...)
            {
                gmml::log(__LINE__, __FILE__, gmml::ERR, "Model issue: this ain't an int: " + codeUtils::RemoveWhiteSpace(line.substr(10,4)));
                currentModelNumber = 1; // Seems like a reasonable default.
            }
        }
        // ATOM
        else if ( (recordName == "ATOM") || (recordName == "HETATM") )
        {
            atomRecords_.emplace_back(line, currentModelNumber);
            //AtomRecord& newAtomRecord = atomRecords_.emplace_back(line, currentModelNumber);
            // This next thing is to pre-organize into residues. Maybe dumb.
//            currentResidueId = newAtomRecord.GetResidueId();
//            if(currentResidueId == previousResidueId)
//            {
//                this->GetCurrentResidue().AddAtom(&newAtomRecord);
//                std::cout << "Added to current residue\n";
//            }
//            else
//            {
//                this->CreateNewResidue(&newAtomRecord);
//                std::cout << "New residue created!\n";
//            }
//            previousResidueId = currentResidueId;
        }
    }
}

//////////////////////////////////////////////////////////
//                      ACCESSOR                        //
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//                      MUTATOR                         //
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//                      FUNCTIONS                       //
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
////////////////////////////////////////////////////////
void CoordinateSection::Print(std::ostream &out) const
{
    out << "The atom records are: " << "\n";
    for (auto &atomRecord : this->GetAtomRecords())
    {
        atomRecord.Print(out);
    }
//    out << "The number of residues is: " << this->GetResidues().size() << "\n";
//    out << "The residues are: " << "\n";
//    for (auto &residue : this->GetResidues())
//    {
//        residue.Print();
//    }
}
