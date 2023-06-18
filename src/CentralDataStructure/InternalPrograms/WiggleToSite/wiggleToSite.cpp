#include "includes/CentralDataStructure/InternalPrograms/WiggleToSite/wiggleToSite.hpp"
#include "includes/CentralDataStructure/InternalPrograms/WiggleToSite/wiggleToSite.hpp"

#include "../../../../includes/CentralDataStructure/Selections/templatedSelections.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbSelections.hpp" //select
#include "includes/CentralDataStructure/Editors/superimposition.hpp"
#include "includes/CentralDataStructure/Shapers/atomToCoordinateInterface.hpp"
#include "includes/CodeUtils/metropolisCriterion.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp" // calculateGeometricCenter

//Prototype: Working and producing useful data in 1.5 days. Included fixing some things in the CDS.
using gmmlPrograms::WiggleToSite;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
WiggleToSite::WiggleToSite(WiggleToSiteInputs inputStruct) : substrate_(inputStruct.substrateFile_), carbohydrate_(inputStruct.carbohydrateSequence_)
{
    this->getCarbohydrate().Generate3DStructureFiles("./", "initial");
    std::vector<Coordinate*> carbohydrateCoordinates = cds::getCoordinatesFromAtoms(this->getCarbohydrate().getAtoms());
    //const Residue* superimpositionTarget = codeUtils::findElementWithNumber(this->getSubstrate().getResidues(), inputStruct.superimpositionTargetResidue_);
    const Residue* superimpositionTarget = pdb::residueSelector(this->getSubstrate(), inputStruct.superimpositionTargetResidue_, inputStruct.substrateModelNumber_);
    const Residue* wigglingTarget = pdb::residueSelector(this->getSubstrate(), inputStruct.wigglingTargetResidue_, inputStruct.substrateModelNumber_);
    if(superimpositionTarget == nullptr || wigglingTarget == nullptr)
    {
        std::stringstream ss;
        ss << "Selection for superimposition target: " << inputStruct.superimpositionTargetResidue_ << "\nOr selection for wiggling target: " << inputStruct.wigglingTargetResidue_ << " was not found\n";
        throw std::runtime_error(ss.str());
    }
//    std::cout << "Super target is " << superimpositionTarget->getStringId() << std::endl;
    Residue* superimposeMe = codeUtils::findElementWithNumber(this->getCarbohydrate().getResidues(), inputStruct.carbohydrateSuperimpositionResidue_);
    Residue* wiggleMe = codeUtils::findElementWithNumber(this->getCarbohydrate().getResidues(), inputStruct.carbohydrateWigglingResidue_);
//    std::cout << "Residue to wiggle is " << wiggleMe->getStringId();
    this->superimpose(carbohydrateCoordinates, superimpositionTarget, superimposeMe);
    this->getCarbohydrate().Generate3DStructureFiles("./", "superimposed");
    this->determineWiggleLinkages(superimposeMe, wiggleMe);
    std::vector<Atom*> substrateWithoutSuperimpositionAtoms = codeUtils::findElementsNotInVector(this->getSubstrate().getAtoms(), superimpositionTarget->getAtoms());
    std::vector<Atom*> substrateAtomsToAvoidOverlappingWith = codeUtils::findElementsNotInVector(substrateWithoutSuperimpositionAtoms, wigglingTarget->getAtoms());
    this->coordsToAvoids_ = cds::getCoordinatesFromAtoms(substrateAtomsToAvoidOverlappingWith);
    this->setCurrentOverlapCount(cds::CountOverlappingCoordinates(this->getCoordsToAvoid(), this->getCarbohydrate().getCoordinates()));
    // call below function.
//    std::cout << "Finished reading and ready to rock captain" << std::endl;
    this->wiggleMeCoordinates_ = {wiggleMe->FindAtom("C1")->getCoordinate(), wiggleMe->FindAtom("C3")->getCoordinate(), wiggleMe->FindAtom("C5")->getCoordinate()};
    this->wiggleTargetCoordinates_ = {wigglingTarget->FindAtom("C1")->getCoordinate(), wigglingTarget->FindAtom("C3")->getCoordinate(), wigglingTarget->FindAtom("C5")->getCoordinate()};
    if (wiggleMeCoordinates_.size() < 3 || wiggleTargetCoordinates_.size() < 3)
    {
        throw std::runtime_error("Did not find the cooordinates of the atoms required for wiggling\n");
    }
    this->setCurrentDistance(this->calculateDistance());
    int structureCount = this->minimizeDistance(inputStruct.persistCycles_, !inputStruct.isDeterministic_);
    this->minimizeDistance(inputStruct.persistCycles_, false, structureCount);
    this->getCarbohydrate().Generate3DStructureFiles("./", "finished");
}

int WiggleToSite::minimizeDistance(int persistCycles, bool useMonteCarlo, int structureCount)
{
//    std::cout << "Starting to wiggle!" << std::endl;
    int cycle = 0;
    while ( (cycle < persistCycles) )
    {
        ++cycle;
        std::stringstream ss;
        ss << "Cycle " << cycle << "/" << persistCycles << "\n";
        gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
        this->randomizeLinkageOrder();
        for(auto & linkage : this->getWiggleLinkages())
        {
            linkage.SetRandomShapeUsingMetadata(true);
            if ( this->acceptDistance(useMonteCarlo) && this->acceptOverlaps() )
            {
                cycle = 0; // reset when it improves
//                std::cout << "Improved distance: " << this->getCurrentDistance() << "\n";
//                std::stringstream fileName;
//                fileName << std::setfill('0') << std::setw(6) << ++structureCount << "_" << this->getCurrentDistance();
//                this->getCarbohydrate().Generate3DStructureFiles("./", "best" + fileName.str());
            }
            else
            {
                linkage.SetShapeToPrevious();
            }
        }
//        if (structureCount > 500)
//        {
//            return structureCount;
//        }
    }
//    std::cout << "ALL DONE HON\n";
    return structureCount;
}

//////////////////////////////////////////////////////////
//                  PRIVATE FUNCTIONS                   //
//////////////////////////////////////////////////////////
void WiggleToSite::superimpose(std::vector<Coordinate*>& carbohydrateCoordinates, const Residue* superimpositionTarget, Residue* superimposeMe)
{
    // Limiting the selection to just these atoms as sometimes hydrogens or an oxygen is missing from xtal. That's ok.
    std::vector<Coordinate*> superimposeMeCoordinates = {superimposeMe->FindAtom("C1")->getCoordinate(), superimposeMe->FindAtom("C3")->getCoordinate(), superimposeMe->FindAtom("C5")->getCoordinate()};
    std::vector<Coordinate*> superTargetCoordinates = {superimpositionTarget->FindAtom("C1")->getCoordinate(), superimpositionTarget->FindAtom("C3")->getCoordinate(), superimpositionTarget->FindAtom("C5")->getCoordinate()};
    cds::Superimpose(superimposeMeCoordinates, superTargetCoordinates, carbohydrateCoordinates); // "alsoMoving" are the carbohydrate Coordinates
    return;
}

std::vector<cds::ResidueLinkage>& WiggleToSite::determineWiggleLinkages(Residue* startResidue, Residue* endResidue)
{
       std::vector<Residue*> residuesInPath;
       bool targetFound = false;
       std::vector<Residue*> visitedResidues;
       std::cout << "Gonna find path between " << startResidue->getName() << " and " << endResidue->getName() << "\n";
       codeUtils::findPathBetweenElementsInGraph(startResidue, endResidue, visitedResidues, residuesInPath, targetFound);
       Residue* previousResidue = nullptr; // wanna skip the first iteration
       for(auto & residue : residuesInPath)
       {
           std::cout << residue->getName() << "_" << residue->getNumber() << ", ";
           if (previousResidue != nullptr)
           {
               wiggleLinkages_.emplace_back(cds::ResidueLinkage(previousResidue, residue));
           }
           previousResidue = residue;
       }
       std::cout << "\nLinkages I behold:\n" << std::endl;
       for(auto & linkage : this->getWiggleLinkages())
       {
           std::cout << linkage.GetName() << ": " << linkage.GetNumberOfShapes() << std::endl;
       }
       return this->getWiggleLinkages();
}

double WiggleToSite::calculateDistance()
{
    return wiggleTargetCoordinates_.at(0)->Distance(wiggleMeCoordinates_.at(0));
//    double totalDistance = 0.0;
//    totalDistance += wiggleTargetCoordinates_.at(0)->Distance(wiggleMeCoordinates_.at(0));
//    totalDistance += wiggleTargetCoordinates_.at(1)->Distance(wiggleMeCoordinates_.at(1));
//    totalDistance += wiggleTargetCoordinates_.at(2)->Distance(wiggleMeCoordinates_.at(2));
//    return totalDistance;
//    //return cds::calculateGeometricCenter(wiggleTargetCoordinates_).Distance(&cds::calculateGeometricCenter(wiggleMeCoordinates_));
}

bool WiggleToSite::acceptOverlaps()
{
    unsigned int overlapCount = cds::CountOverlappingCoordinates(this->getCoordsToAvoid(), this->getCarbohydrate().getCoordinates());
    if (overlapCount > this->getCurrentOverlapCount())
    {
        return false;
    }
    this->setCurrentOverlapCount(overlapCount);
    return true;
}

bool WiggleToSite::acceptDistance(bool useMonteCarlo)
{
    double distance = this->calculateDistance();
    if (distance < this->getCurrentDistance() || (useMonteCarlo && monte_carlo::accept_via_metropolis_criterion(distance - this->getCurrentDistance()))  )
    {
        this->setCurrentDistance(distance);
        return true;
    }
    return false;
}

