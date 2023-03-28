#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinBuilder.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
//#include "includes/CentralDataStructure/Overlaps/beadResidues.hpp"
#include "includes/CodeUtils/metropolisCriterion.hpp"
//#include "includes/InternalPrograms/functionsForGMML.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/strings.hpp" // split
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/Writers/cdsOffWriter.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/cdsFunctions.hpp" // bondAtomsByDistance
#include "includes/CentralDataStructure/Selections/residueSelections.hpp" // selectResiduesByType
// ToDo Check for negative overlap in case the funk gets funky.
// ToDo The cout for accepting overlap changes doesn't match the values printed. Why? Is it making them high, not printing, accepting lower, then rejecting, etc?
using cds::Assembly;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
GlycoproteinBuilder::GlycoproteinBuilder(glycoprotein::GlycoproteinBuilderInputs inputStruct)
{
    try
    {
        this->SetIsDeterministic(inputStruct.isDeterministic_);
        this->SetNumberOfOutputStructures(inputStruct.number3DStructures_);
        this->SetPersistCycles(inputStruct.persistCycles_);
        this->SetOverlapTolerance(inputStruct.overlapTolerance_);
        pdb::PdbFile pdbFile(inputStruct.substrateFileName_);
        glycoprotein_ = std::move(*(pdbFile.getAssemblies().front()));
        cds::bondAtomsByDistance(glycoprotein_.getAtoms());
        gmml::log(__LINE__, __FILE__, gmml::INF, "Attaching Glycans To Glycosites.");
        this->CreateGlycosites(inputStruct.glycositesInputVector_, inputStruct.prepFileLocation_);
    }
    catch (const std::string &errorMessage)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
        throw std::runtime_error(errorMessage);
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Initialization of Glycoprotein builder complete!");
    return;
}
//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void GlycoproteinBuilder::WriteOutputFiles()
{
    // Pdb, original numbering.
    std::string fileName = "GlycoProtein_All_Resolved.pdb";
    std::ofstream outFileStream;
    outFileStream.open(fileName.c_str());
    cds::writeMoleculeToPdb(outFileStream, this->getGlycoprotein()->getResidues());
    outFileStream.close();
    // Off file, serializes.
    fileName = "GlycoProtein_All_Resolved.off";
	outFileStream.open(fileName.c_str());
	cds::WriteAssemblyToOffFile(this->getGlycoprotein(), outFileStream, "GLYCOPROTEINBUILDER");
	outFileStream.close();
	// Pdb, serialized numbering.
	fileName = "GlycoProtein_All_Resolved_Serialized.pdb";
	outFileStream.open(fileName.c_str());
	cds::writeMoleculeToPdb(outFileStream, this->getGlycoprotein()->getResidues());
//    this->DeleteSitesIterativelyWithAtomicOverlapAboveTolerance(this->GetGlycosites(), this->GetOverlapTolerance());
//	std::stringstream logss;	
//    logss << "Atomic overlap is " << this->CalculateOverlaps(ATOMIC) << "\n";
//    logss << "Global overlap is " << this->GetGlobalOverlap() << "\n";
//	  gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    this->PrintDihedralAnglesAndOverlapOfGlycosites();
//	  this->WritePDBFile(this->GetGlycoproteinAssembly(), this->GetWorkingDirectory(), "GlycoProtein_Resolved", false);
    return;
}

void GlycoproteinBuilder::ResolveOverlaps()
{
    std::cout << "In here" << std::endl;
	bool randomize = !this->GetIsDeterministic();
    std::cout << "In here" << std::endl;
	if (randomize)
	{ // First try a very fast/cheap approach
	    std::cout << "Trying to be dumb" << std::endl;
		if (this->DumbRandomWalk()) // returns true if it fully resolves overlaps.
		{
			return;
		}
	}
	bool useMonteCarlo = true;
	bool wiggleFirstLinkageOnly = true; // Only happens when passed to Wiggle function.
	std::stringstream logss;
    std::cout << "Need to wiggle" << std::endl;
	this->Wiggle(RESIDUE, this->GetPersistCycles(), wiggleFirstLinkageOnly);
	logss << "1. Post WiggleFirst Overlaps ResidueResolution: " << this->CalculateOverlaps(RESIDUE) << ". AtomicResolution: " << this->CalculateOverlaps(ATOMIC) << std::endl;
	std::cout << logss.str();
	if (randomize)
	{
		this->RandomDescent(RESIDUE, this->GetPersistCycles(), useMonteCarlo);
		logss << "2. Post Monte Carlo Overlaps ResidueResolution: " << this->CalculateOverlaps(RESIDUE) << ". AtomicResolution: " << this->CalculateOverlaps(ATOMIC) << std::endl;
	    std::cout << logss.str();
	}
	this->Wiggle(RESIDUE, this->GetPersistCycles());
	logss << "3. Post Wiggle Overlaps ResidueResolution: " << this->CalculateOverlaps(RESIDUE) << ". AtomicResolution: " << this->CalculateOverlaps(ATOMIC) << std::endl;
    std::cout << logss.str();
//	this->Wiggle(BEAD, this->GetPersistCycles(), wiggleFirstLinkageOnly);
//	logss << "4. Post WiggleFirst Overlaps Bead: " << this->CalculateOverlaps(BEAD) << ". Atomic: " << this->CalculateOverlaps(ATOMIC) << std::endl;
//	if (randomize)
//	{
//		this->RandomDescent(BEAD, this->GetPersistCycles(), useMonteCarlo);
//		logss << "5. Post Monte Carlo Overlaps Bead: " << this->CalculateOverlaps(BEAD) << ". Atomic: " << this->CalculateOverlaps(ATOMIC) << std::endl;
//	}
//	this->Wiggle(BEAD, this->GetPersistCycles());
//	logss << "6. Post Wiggle Overlaps Bead: " << this->CalculateOverlaps(BEAD) << ". Atomic: " << this->CalculateOverlaps(ATOMIC) << std::endl;
    std::cout << logss.str();
    std::cout << logss.str();
	this->Wiggle(ATOMIC, this->GetPersistCycles(), wiggleFirstLinkageOnly);
	logss << "7. Post WiggleFirst Overlaps Atomic: " << this->CalculateOverlaps(ATOMIC) << std::endl;
    std::cout << logss.str();
	if (randomize)
	{
		this->RandomDescent(ATOMIC, this->GetPersistCycles(), useMonteCarlo);
		logss << "8. Post Monte Carlo Overlaps Atomic: " << this->CalculateOverlaps(ATOMIC) << std::endl;
	}
    std::cout << logss.str();
	this->Wiggle(ATOMIC, this->GetPersistCycles());
    std::cout << logss.str();
	logss << "9. Post Wiggle Overlaps Atomic: " << this->CalculateOverlaps(ATOMIC) << std::endl;
	int interval = 1; // Default interval when wiggle is 5 degrees.
	this->Wiggle(ATOMIC, this->GetPersistCycles(), wiggleFirstLinkageOnly, interval);
	logss << "10. Post WiggleFirst Interval 1 Overlaps Atomic: " << this->CalculateOverlaps(ATOMIC) << std::endl;
	this->PrintDihedralAnglesAndOverlapOfGlycosites();
	this->CalculateOverlaps(ATOMIC, ALL);
	gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    return;
}

void GlycoproteinBuilder::RandomDescent(Resolution resolutionLevel, int persistCycles, bool monte_carlo)
{
	std::stringstream logss;
	logss << "Random Decent, persisting for " << persistCycles << " cycles and monte carlo is set as " << std::boolalpha << monte_carlo << ".\n";
	int cycle = 1;
	bool stop = false;
	double previous_glycan_overlap, new_glycan_overlap, previous_protein_overlap, new_protein_overlap;
	double lowest_global_overlap = this->CalculateOverlaps(resolutionLevel, ALL);
	double new_global_overlap;
	double overlap_difference = 0.0;
	std::vector<GlycosylationSite*> sites_with_overlaps = this->DetermineSitesWithOverlap(resolutionLevel, ALL);
	if (sites_with_overlaps.size() == 0)
	{
		logss << "Stopping with all overlaps resolved.\n";
		stop = true;
	}
	//logss << "Initial torsions and overlaps:\n";
	//this->PrintDihedralAnglesAndOverlapOfGlycosites();
	while ( (cycle < persistCycles) && (stop == false) )
	{
		logss << "Cycle " << cycle << "/" << persistCycles << "\n";
		++cycle;
		std::random_shuffle (sites_with_overlaps.begin(), sites_with_overlaps.end());
		for(auto &current_glycosite : sites_with_overlaps)
		{
			// logss << "Checking " << current_glycosite->GetResidue()->GetId() << "\n";
			previous_glycan_overlap = current_glycosite->CalculateOverlaps(resolutionLevel, GLYCAN);
			previous_protein_overlap = current_glycosite->CalculateOverlaps(resolutionLevel, PROTEIN);
			current_glycosite->SetRandomDihedralAnglesUsingMetadata();
			//logss << "Site: " << current_glycosite->GetResidueNumber() << "\n";
			new_glycan_overlap = current_glycosite->CalculateOverlaps(resolutionLevel, GLYCAN);
			new_protein_overlap = current_glycosite->CalculateOverlaps(resolutionLevel, PROTEIN);
			overlap_difference = (new_glycan_overlap + (new_protein_overlap*5)) - (previous_glycan_overlap + (previous_protein_overlap*5));
			if (overlap_difference >= 0.0) // if the change made it worse
			{
				current_glycosite->ResetDihedralAngles();
			}
			else if ( (monte_carlo) && (! monte_carlo::accept_via_metropolis_criterion(overlap_difference)) )
			{
				current_glycosite->ResetDihedralAngles();
			}
			else
			{
				logss << "RandomDescent accepted a change of " << overlap_difference << "\n";
			}
		}
		//logss << "Updating list of sites with overlaps." << std::endl;
		new_global_overlap = this->CalculateOverlaps(resolutionLevel, ALL);
		sites_with_overlaps = this->DetermineSitesWithOverlap(resolutionLevel, ALL); // Moved glycans may clash with other glycans. Need to check.
		if (sites_with_overlaps.size() == 0)
		{
			logss << "Stopping with all overlaps resolved and global overlap is " << new_global_overlap <<  "\n";
			stop = true;
		}
		//gmml::WritePDBFile(glycosites.at(0).GetGlycoprotein(), cycle, "current", new_global_overlap);
		//	        std::stringstream ss;
		//	        ss << "current_cycle_" << cycle << "_overlap" << new_global_overlap;
		//	        this->WriteOutputFile(ss.str());
		if ( lowest_global_overlap > new_global_overlap + 1 )
		{
			//   logss << "Lowest: " << lowest_global_overlap << ", Current: " << new_global_overlap << "\n";
			//std::stringstream ss;
			//ss << "bestRandom_overlap" << new_global_overlap;
			//gmml::WritePDBFile(this->GetGlycoproteinAssembly(), this->GetWorkingDirectory(), ss.str());
			lowest_global_overlap = new_global_overlap;
			logss << "Found a lower overlap. Resetting cycle count\n";
			cycle = 1;
		}
	}
	gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
	return;
}

void GlycoproteinBuilder::Wiggle(Resolution resolutionLevel, int persistCycles, bool firstLinkageOnly, int interval)
{
	std::stringstream logss;
    if (resolutionLevel == RESIDUE)
		logss << "GPWiggle with Resolution level RESIDUE ";
	if (resolutionLevel == ATOMIC)
		logss << "GPwiggle with Resolution level ATOMIC ";
	logss << persistCycles << " persistCycles, " << interval << " interval\n";
	std::cout << logss.str();
    std::vector<GlycosylationSite*> sites_with_overlaps = this->DetermineSitesWithOverlap(resolutionLevel, ALL);
    int cycle = 1;
    bool stop = false;
    int savedOverlap = this->CalculateOverlaps(resolutionLevel, ALL);
    if (sites_with_overlaps.size() == 0)
    {
        logss << "Stopping with all overlaps resolved.\n";
        stop = true;
    }
    while ( (cycle < persistCycles) && (stop == false) )
    {
        ++cycle;
        logss << "Cycle " << cycle << "/" << persistCycles << "\n";
        std::cout << logss.str();
        std::random_shuffle (sites_with_overlaps.begin(), sites_with_overlaps.end());
        for(auto &glycosite : sites_with_overlaps)
        {
            glycosite->Wiggle(firstLinkageOnly);
        }
        // Check which sites still have overlaps, stop if none.
        sites_with_overlaps = this->DetermineSitesWithOverlap(resolutionLevel, ALL); // Moved glycans may clash with other glycans. Need to check.
        if (sites_with_overlaps.size() == 0)
        {
            logss << "Stopping with all overlaps resolved.\n";
            stop = true;
        }
        // The above will trigger re-calculation of overlaps, so can do this:
        if (savedOverlap > (this->CalculateOverlaps(resolutionLevel, ALL)))
        {
//        	std::stringstream ss;
//        	ss << "best_Wiggle_overlap" << this->GetGlobalOverlap();
//        	gmml::WritePDBFile(this->GetGlycoproteinAssembly(), this->GetWorkingDirectory(), ss.str());
        	logss << "Resetting cycle count to zero, will persist for another " << persistCycles << " cycles." << std::endl;
        	savedOverlap = this->CalculateOverlaps(resolutionLevel, ALL);
        	cycle = 1;
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    return;
}

bool GlycoproteinBuilder::DumbRandomWalk(int maxCycles)
{
	std::stringstream logss;
	logss << "Starting DumbRandomWalk\n";
    std::cout << logss.str();
    int cycle = 1;
    std::vector<GlycosylationSite*> sites_with_overlaps = DetermineSitesWithOverlap();
    while (cycle < maxCycles)
    {
        ++cycle;
        logss << "Cycle " << cycle << " of " << maxCycles << std::endl;
        std::cout << logss.str();
        for(auto &currentGlycosite : sites_with_overlaps)
        {
        	currentGlycosite->SetRandomDihedralAnglesUsingMetadata();
        }
        logss << "Updating list of sites with overlaps." << std::endl;
        sites_with_overlaps = this->DetermineSitesWithOverlap(); // Moved glycans may clash with other glycans. Need to check.
        if (sites_with_overlaps.size() == 0)
        {
            logss << "DumbRandomWalk resolved the overlaps. Stopping\n";
            return true;
        }
        std::cout << logss.str();
    }
	logss << "DumbRandomWalk did not resolve the overlaps\n";
    gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    return false;
}

void GlycoproteinBuilder::CreateGlycosites(std::vector<glycoprotein::GlycositeInput> glycositesInputVector, std::string prepFileLocation)
{
    std::vector<Residue*> proteinResidues = this->getGlycoprotein()->getResidues(); // Before any glycans are added.
	for (auto &glycositeInput : glycositesInputVector)
	{
	    gmml::log(__LINE__, __FILE__, gmml::INF, "Creating glycosite on residue " + glycositeInput.proteinResidueId_ + " with glycan " + glycositeInput.glycanInput_ );
	    Carbohydrate* carb = static_cast<Carbohydrate*>(glycoprotein_.addMolecule(std::make_unique<Carbohydrate>(glycositeInput.glycanInput_, prepFileLocation)));
	    if (!carb->IsStatusOk()) // status checking needs to die. Still haven't figured out catching it in gems/Python smh.
	    {
	        throw std::runtime_error(carb->GetStatusMessage());
	    }
        Residue* glycositeResidue = this->SelectResidueFromInput(glycositeInput.proteinResidueId_);
        if (glycositeResidue == nullptr)
        {
            std::cout << "Did not find glycosite residue" << std::endl;
            throw std::runtime_error("Did not find a residue with id matching " + glycositeInput.proteinResidueId_);
        }
	    std::vector<Residue*> otherResidues = proteinResidues;
	    otherResidues.erase(std::remove(otherResidues.begin(), otherResidues.end(), glycositeResidue), otherResidues.end());
	    gmml::log(__LINE__, __FILE__, gmml::INF, "About to emplace_back to glycosites with: " + glycositeInput.proteinResidueId_ + " and glycan " + glycositeInput.glycanInput_);
        unsigned int highestResidueNumber = cdsSelections::findHighestResidueNumber(this->getGlycoprotein()->getResidues());
		glycosites_.emplace_back(glycositeResidue, carb, otherResidues, highestResidueNumber);
	    std::cout << "Done with glycan" << std::endl;
		gmml::log(__LINE__, __FILE__, gmml::INF, "Completed creating glycosite on residue " + glycositeInput.proteinResidueId_ + " with glycan " + glycositeInput.glycanInput_);
	}
    std::cout << "Done attaching all glycans" << std::endl;
    this->SetOtherGlycosites();
    return;
}

Residue* GlycoproteinBuilder::SelectResidueFromInput(const std::string userSelection)
{ // Chain_residueNumber_insertionCode* *optional.
    std::vector<std::string> splitUserSelection = codeUtils::split(userSelection, '_');
    if(splitUserSelection.size() < 2)
    {
        throw std::runtime_error("userSelection (" + userSelection + ") for residue to glycosylate is incorrect format.\nMust be chain_residueNumber_insertionCode.\nInsertionCode is optional. Chain can be ? if no chain numbers are in input.\nExamples: ?_24_? or ?_24 will use the first residue it encounters numbered 24. A_24_B is A chain, residue 24, insertion code B");
    }
    std::string userSelectedChain = "";
    if(splitUserSelection.at(0) != "?")
    {
        userSelectedChain = splitUserSelection.at(0);
    }
    std::string userSelectedResidue = splitUserSelection.at(1);
    if(splitUserSelection.size() == 3)
    {
        userSelectedResidue += splitUserSelection.at(2); // So will be 24A or just 24.
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "We working with " + userSelectedChain + "_" + userSelectedResidue);
    for (auto &residue : this->getGlycoprotein()->getResidues())
    {
        pdb::PdbResidue* pdbResidue = static_cast<pdb::PdbResidue*>(residue);
        //std::cout << pdbResidue->getChainId() << "_";
        std::cout << pdbResidue->getNumberAndInsertionCode() << "\n";
        if ( (pdbResidue->getChainId() == userSelectedChain) && (pdbResidue->getNumberAndInsertionCode() == userSelectedResidue) )
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Id of selected glycosite: " + pdbResidue->printId());
            return residue;
        }
    }
    return nullptr;
}

void GlycoproteinBuilder::SetOtherGlycosites()
{
    for (auto &glycosite : this->GetGlycosites())
    {
    	glycosite.SetOtherGlycosites(this->GetGlycosites());
    }
    return;
}

void GlycoproteinBuilder::PrintDihedralAnglesAndOverlapOfGlycosites()
{
    for(auto &glycosite : this->GetGlycosites())
    {
        glycosite.Print("All");
    }
    return;
}

void GlycoproteinBuilder::SetRandomDihedralAnglesUsingMetadata()
{
    for(auto &glycosite : this->GetGlycosites())
    {
        glycosite.SetRandomDihedralAnglesUsingMetadata();
    }
    return;
}

int GlycoproteinBuilder::CalculateOverlaps(Resolution resolutionLevel, MoleculeType moleculeType)
{
    int overlap = 0;
    for(auto &glycosite : this->GetGlycosites())
    {
        overlap += glycosite.CalculateOverlaps(resolutionLevel, moleculeType);
    }
    return overlap;
}

std::vector<GlycosylationSite*> GlycoproteinBuilder::DetermineSitesWithOverlap(Resolution resolutionLevel, MoleculeType moleculeType)
{
    std::vector<GlycosylationSite*> sites_to_return;
    int overlap = 0;
    for(auto & glycosite : this->GetGlycosites())
    {
        overlap = glycosite.CalculateOverlaps(resolutionLevel, moleculeType);
        if ( overlap > this->GetOverlapTolerance())
        {
            sites_to_return.push_back(&glycosite);
        }
    }
    return sites_to_return;
}

void GlycoproteinBuilder::DeleteSitesIterativelyWithAtomicOverlapAboveTolerance(std::vector<GlycosylationSite> &glycosites, int tolerance)
{
	std::stringstream logss;
    logss << "Atomic overlap before deleting sites is " << this->CalculateOverlaps(ATOMIC) << "\n";
    bool continue_deleting = true;
    // While overlap for any site is > tolerance delete site with highest overlap then re-calculate overlaps as there may be glycan-glycan overlap.
    while (continue_deleting)
    {
    	GlycosylationSite *worst_site = glycosites.data(); // Pointer to the first glycosite. Remember an erase/remove "advances"
    	for (auto &glycosite : this->GetGlycosites())
    	{
    	    if  (glycosite.CalculateOverlaps(ATOMIC) > worst_site->CalculateOverlaps(ATOMIC))
    	    {
    	        worst_site = &glycosite;
    	    }
    	}
        int worst_site_overlap = worst_site->CalculateOverlaps(ATOMIC);
        if ( worst_site_overlap > tolerance)
        {
            continue_deleting = true;
            //worst_site->Remove(this->GetGlycoproteinAssemblyPtr());
            worst_site->Rename_Protein_Residue_From_GLYCAM_To_Standard();
            logss << "Site " << worst_site->GetResidueId() << ": " << worst_site_overlap << " :" << "Removed\n";
            glycosites.erase(std::remove(glycosites.begin(), glycosites.end(), *worst_site), glycosites.end());
            this->SetOtherGlycosites(); // This needs to be updated after you delete each site, so the others no longer have a pointer to it. Great design Oliver well done you nailed it high five.
        }
        else
        {
            continue_deleting = false;
        }
        if(glycosites.empty()) // If we have deleted every site
        {
            continue_deleting = false;
        }
    }
   	gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    return;
}
