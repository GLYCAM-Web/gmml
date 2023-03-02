#include "includes/CentralDataStructure/InternalPrograms/WiggleToSite/wiggleToSite.hpp"
#include "includes/CentralDataStructure/InternalPrograms/WiggleToSite/wiggleToSite.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Editors/superimposition.hpp"
#include "includes/CentralDataStructure/Shapers/atomToCoordinateInterface.hpp"
//#include "includes/CentralDataStructure/Selections/shaperSelections.hpp"

#include "includes/CodeUtils/templatedSelections.hpp"

void gmmlPrograms::wiggleToSite(WiggleToSiteInputs inputStruct)
{
    cdsCondensedSequence::Carbohydrate carbohydrate(inputStruct.carbohydrateSequence_);
    pdb::PdbFile substrate(inputStruct.substrateFile_);
    cds::Residue* superimpositionTarget = codeUtils::findElementWithNumber(substrate.getResidues(), inputStruct.superimpositionTargetResidue_);
    cds::Residue* superimposeMe = codeUtils::findElementWithNumber(carbohydrate.getResidues(), inputStruct.carbohydrateSuperimpositionResidue_);


    cds::Residue* wigglingTarget = codeUtils::findElementWithNumber(substrate.getResidues(), inputStruct.wigglingTargetResidue_);
    cds::Residue* wiggleMe = codeUtils::findElementWithNumber(carbohydrate.getResidues(), inputStruct.carbohydrateWigglingResidue_);

    std::vector<Atom*> substrateWithoutSuperimpositionAtoms = codeUtils::findElementsNotInVector(substrate.getAtoms(), superimpositionTarget->getAtoms());
    std::vector<Atom*> substrateAtomsToAvoidOverlappingWith = codeUtils::findElementsNotInVector(substrateWithoutSuperimpositionAtoms, wigglingTarget->getAtoms());
    std::vector<cds::Coordinate*> coordsToAvoids = cds::getCoordinatesFromAtoms(substrateAtomsToAvoidOverlappingWith);
    // call below function.
    std::cout << "Finished reading and ready to rock captain" << std::endl;
    gmmlPrograms::wiggleToSite(coordsToAvoids, &carbohydrate, wigglingTarget, wiggleMe, superimpositionTarget, superimposeMe);
    carbohydrate.Generate3DStructureFiles("./", "bob");
}

void gmmlPrograms::wiggleToSite(const std::vector<Coordinate*> avoidOverlappingWithMe, cdsCondensedSequence::Carbohydrate* carbohydrate, const cds::Residue* wiggleTarget, cds::Residue* wiggleMe, const cds::Residue* superimpositionTarget, cds::Residue* superimposeMe)
{
    // Superimposition
    //std::cout << carbohydrate->getAtoms().front()->getCoordinate()->ToString() << std::endl;
    std::vector<Coordinate*> carbohydrateCoordinates = cds::getCoordinatesFromAtoms(carbohydrate->getAtoms());
        // Limiting the selection to just these atoms as sometimes hydrogens or an oxygen is missing from xtal. That's ok.
    std::vector<Coordinate*> superimposeMeCoordinates = {superimposeMe->FindAtom("C1")->getCoordinate(), superimposeMe->FindAtom("C3")->getCoordinate(), superimposeMe->FindAtom("C5")->getCoordinate()};
    std::vector<Coordinate*> superTargetCoordinates = {superimpositionTarget->FindAtom("C1")->getCoordinate(), superimpositionTarget->FindAtom("C3")->getCoordinate(), superimpositionTarget->FindAtom("C5")->getCoordinate()};
    cds::Superimpose(superimposeMeCoordinates, superTargetCoordinates, carbohydrateCoordinates); // alsoMoving is all the carbohydrateCoordinates
    // Wiggling to target/site
    std::vector<Coordinate*> wiggleMeCoordinates = {wiggleMe->FindAtom("C1")->getCoordinate(), wiggleMe->FindAtom("C3")->getCoordinate(), wiggleMe->FindAtom("C5")->getCoordinate()};
    std::vector<Coordinate*> wiggleTargetCoordinates = {wiggleTarget->FindAtom("C1")->getCoordinate(), wiggleTarget->FindAtom("C3")->getCoordinate(), wiggleTarget->FindAtom("C5")->getCoordinate()};
//    for(auto & coord : wiggleMeCoordinates)
//    {
//        std::cout << coord->ToString() << std::endl;
//    }

    // Explore from wiggleMe to superimposeMe, get all linkages inbetween
    std::vector<cds::Residue*> residuesInPath;
    bool targetFound = false;
    std::vector<cds::Residue*> visitedResidues;
    codeUtils::findPathBetweenElementsInGraph(superimposeMe, wiggleMe, visitedResidues, residuesInPath, targetFound);
    std::vector<cds::ResidueLinkage> wiggleLinkages;
    cds::Residue* previousResidue = nullptr; // wanna skip the first iteration
    for(auto & residue : residuesInPath)
    {
        if (previousResidue != nullptr)
        {
            wiggleLinkages.emplace_back(cds::ResidueLinkage(previousResidue, residue));
        }
        //std::cout << residue->getId() << ", ";
        previousResidue = residue;
    }

    std::cout << "Linkages I behold:\n" << std::endl;
    for(auto & linkage : wiggleLinkages)
    {
        std::cout << linkage.GetName() << ": " << linkage.GetNumberOfShapes() << std::endl;
    }
    int cycle = 0;
    bool stop = false;
    int persistCycles = 500;
    unsigned int totalStructureCounter = 0;
    double bestDistance = wiggleMeCoordinates.front()->Distance(wiggleTargetCoordinates.front()); //ToDo average all of them instead of just front();
    while ( (cycle < persistCycles) && (stop == false) )
    {
        ++cycle;
        std::stringstream ss;
        ss << "Cycle " << cycle << "/" << persistCycles << "\n";
        gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
        //std::random_shuffle(wiggleLinkages.begin(), wiggleLinkages.end());
        for(auto & linkage : wiggleLinkages)
        {
            linkage.SetRandomShapeUsingMetadata(true);
            double newDistance = wiggleMeCoordinates.front()->Distance( wiggleTargetCoordinates.front() );
            std::stringstream cycleCount;
            cycleCount << std::setfill('0') << std::setw(6) << ++totalStructureCounter;
//            if (bestDistance < newDistance)
//            {
//                linkage.SetShapeToPrevious();
//            }
//            else
//            {
//                bestDistance = newDistance;
//                cycle = 0; // reset when it improves
//                carbohydrate->Generate3DStructureFiles("./", "best" + cycleCount.str());
//            }
            carbohydrate->Generate3DStructureFiles("./", "all" + cycleCount.str());
        }
    }
//    }

    std::cout << "ALL DONE HON\n";
}







