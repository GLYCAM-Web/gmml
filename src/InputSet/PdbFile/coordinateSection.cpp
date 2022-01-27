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
//    PdbResidue* currentResidue;
//    std::string currentResidueId = "?_?_?_?_?"; // model_resname_insertionCode_resnumber_chain_
//    std::string previousResidueId = "IJustNeedThisToBeDifferentFromCurrentResidueIdForTheinitialCheckOkDontJudgeMeImDoingMyBestOutHere";
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
 //           atomRecords_.emplace_back(line, currentModelNumber);
            atomRecords_.push_back(std::make_unique<AtomRecord>(line, currentModelNumber));
            // This next thing is to pre-organize into residues. Maybe dumb.
//            currentResidueId = newAtomRecord.GetResidueId();
//            if(currentResidueId == previousResidueId)
//            {
//                currentResidue->AddAtom(&newAtomRecord);
//                std::cout << "Added to current residue\n";
//            }
//            else
//            {
//                currentResidue = this->CreateNewResidue(&newAtomRecord);
//                std::cout << "New residue created!\n";
//            }
//            previousResidueId = currentResidueId;
        }
//        std::cerr << "Residues so far are:\n";
//        for (auto &residue : residues_)
//        {
//            residue.Print();
//        }
    }
}

std::vector<pdb::PdbResidue> CoordinateSection::GetResidues() const
{
    std::vector<pdb::PdbResidue> residues;
    std::string id = "?_?_?_?_?"; // model_resname_insertionCode_resnumber_chain_
    std::string previousId = "AStringThatIsUniqueFromIdForTheFirstLoop";
    for(auto &atomUniquePtr : atomRecords_)
    {
        id = atomUniquePtr->GetResidueId();
        std::cout << "Comparing " << id << " with " << previousId << "\n";
        if(id == previousId)
        {
            residues.back().AddAtom(atomUniquePtr.get());
        }
        else
        {
            std::cout << "Creating new residue\n";
            residues.emplace_back(atomUniquePtr.get());
        }
        previousId = id;
    }
    return residues;
}

std::vector<pdb::PdbResidue> CoordinateSection::FindResidues(const std::string selector)
{
    std::vector<pdb::PdbResidue> matchingResidues;
    for(auto &residue : this->GetResidues())
    {
        std::size_t found = residue.GetId().find(selector);
        if(found != std::string::npos)
        {
            matchingResidues.push_back(residue);
            std::cout << residue.GetId() << " contains " << selector << "\n";
        }
    }
    return matchingResidues;
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
//    out << "The atom records are: " << "\n";
//    for (auto &atomRecord : this->GetAtomRecords())
//    {
//        atomRecord.Print(out);
//    }
    out << "The number of residues is: " << this->GetResidues().size() << "\n";
    out << "The residues are: " << "\n";
    for (auto &residue : this->GetResidues())
    {
        residue.Print();
    }
}
