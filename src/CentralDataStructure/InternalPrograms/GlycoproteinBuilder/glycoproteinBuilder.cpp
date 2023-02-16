#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycoproteinBuilder.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/Overlaps/beadResidues.hpp"
#include "includes/CodeUtils/metropolisCriterion.hpp"
//#include "includes/InternalPrograms/functionsForGMML.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/directories.hpp"
#include "includes/CodeUtils/strings.hpp" // split
#include "includes/CentralDataStructure/Writers/pdbWriter.hpp"
#include "includes/CentralDataStructure/Writers/cdsOffWriter.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbFile.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidue.hpp"
#include "includes/CentralDataStructure/cdsFunctions.hpp" // bondAtomsByDistance
#include "includes/CentralDataStructure/Selections/residueSelections.hpp" // selectResiduesByType

//#include "includes/ParameterSet/OffFileSpace/offfile.hpp"

// ToDo Check for negative overlap in case the funk gets funky.
// ToDo The cout for accepting overlap changes doesn't match the values printed. Why? Is it making them high, not printing, accepting lower, then rejecting, etc?
using cds::Assembly;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
GlycoproteinBuilder::GlycoproteinBuilder(std::string inputFile)
{
    GlycoproteinBuilderInputs inputStruct = GPInputs::readGPInputFile(codeUtils::Find_Program_workingDirectory(), inputFile);
	this->InitializeGlycoproteinBuilder(inputStruct);
}
GlycoproteinBuilder::GlycoproteinBuilder(std::string inputFile, std::string workingDirectory)
{
    try
    {
        GlycoproteinBuilderInputs inputStruct = GPInputs::readGPInputFile(workingDirectory, inputFile);
        this->InitializeGlycoproteinBuilder(inputStruct);
    }
    catch (std::runtime_error const &error)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, error.what());
        this->SetStatus("ERROR", error.what());
    }
}
GlycoproteinBuilder::GlycoproteinBuilder(GlycoproteinBuilderInputs inputStruct)
{
    try
    {
        this->InitializeGlycoproteinBuilder(inputStruct);
    }
    catch (std::runtime_error const &error)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, error.what());
        this->SetStatus("ERROR", error.what());
    }
}
//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void GlycoproteinBuilder::SetWorkingDirectory(std::string workingDirectory)
{
	if (workingDirectory == "Default") // This value is set in the default constructor for the input struct
	{
		workingDirectory_ = codeUtils::Find_Program_workingDirectory();
	}
	else
	{
		workingDirectory_ = workingDirectory;
	}
    gmml::log(__LINE__, __FILE__, gmml::INF, "WorkingDirectory is set to  " + workingDirectory_);
	return;
}

void GlycoproteinBuilder::SetPrepFileLocation(std::string prepFileLocation)
{
	if (prepFileLocation == "Default") // This value is set in the default constructor for the input struct
	{
		prepFileLocation_ = codeUtils::Find_Program_Installation_Directory() + "/../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep";
	}
	else
	{
		prepFileLocation_ = prepFileLocation;
	}
	return;
}

void GlycoproteinBuilder::SetIsDeterministic(std::string isDeterministic)
{
	if (isDeterministic == "true")
	{
		isDeterministic_ = true;
	}
	else if (isDeterministic == "false")
	{
		isDeterministic_ = false;
	}
	else
	{
		throw std::runtime_error("isDeterminstic not set to either \"true\" or \"false\" in inputs, value is " + isDeterministic);
	}
	return;
}
//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

void GlycoproteinBuilder::InitializeGlycoproteinBuilder(GlycoproteinBuilderInputs inputStruct)
{
    try
    {
        this->SetWorkingDirectory(inputStruct.workingDirectory_);
        this->SetPrepFileLocation(inputStruct.prepFileLocation_);
        this->SetIsDeterministic(inputStruct.isDeterministic_);
        this->ConvertInputStructEntries(inputStruct);
        pdb::PdbFile pdbFile(this->GetWorkingDirectory() + inputStruct.substrateFileName_);
//        for (auto &assembly : pdbFile.getAssemblies())
//        { // Per assembly as depending on input pdb file, like an NMR structure with lots of models, ensemble.getAtoms() would be very wrong.
//            cds::bondAtomsByDistance(assembly->getAtoms());
//        }
        // Dumb but need to get it working again before I fix this. What should I do when there are lots of models?
        // Don't copy or you lose the ability to cast to pdbResidue.. :(
        glycoprotein_ = std::move(*(pdbFile.getAssemblies().front()));
        cds::bondAtomsByDistance(glycoprotein_.getAtoms());
        gmml::log(__LINE__, __FILE__, gmml::INF, "Attaching Glycans To Glycosites.");
        this->CreateGlycosites(inputStruct.glycositesInputVector_);
    }
    catch (const std::string &errorMessage)
    {
        gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
        throw std::runtime_error(errorMessage);
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Initialization of Glycoprotein builder complete!");
    return;
}

void GlycoproteinBuilder::ConvertInputStructEntries(GlycoproteinBuilderInputs inputStruct)
{
	try
	{
		this->SetNumberOfOutputStructures(std::stoi(inputStruct.number3DStructures_));
		this->SetPersistCycles(std::stoi(inputStruct.persistCycles_));
		this->SetOverlapTolerance(std::stod(inputStruct.overlapTolerance_));
	}
	catch (...)
	{
		std::string errorMessage = "Error converting types in GlycoproteinBuilder::ConvertInputStructEntries(). Check that the values you entered in the input are convertible to integers.";
        gmml::log(__LINE__, __FILE__, gmml::ERR, errorMessage);
		throw std::runtime_error(errorMessage);
	}
	return;
}

void GlycoproteinBuilder::WriteOutputFiles()
{
    // Pdb, original numbering.
    std::string fileName = this->GetWorkingDirectory() + "GlycoProtein_All_Resolved.pdb";
    std::ofstream outFileStream;
    outFileStream.open(fileName.c_str());
    cds::writeMoleculeToPdb(outFileStream, this->GetGlycoproteinAssembly().getResidues());
    outFileStream.close();
    // Off file, serializes.
    fileName = this->GetWorkingDirectory() + "GlycoProtein_All_Resolved.off";
	outFileStream.open(fileName.c_str());
	cds::WriteAssemblyToOffFile(&(this->GetGlycoproteinAssembly()), outFileStream, "GLYCOPROTEINBUILDER");
	outFileStream.close();
	// Pdb, serialized numbering.
	fileName = this->GetWorkingDirectory() + "GlycoProtein_All_Resolved_Serialized.pdb";
	outFileStream.open(fileName.c_str());
	cds::writeMoleculeToPdb(outFileStream, this->GetGlycoproteinAssembly().getResidues());
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
	this->Wiggle(BEAD, this->GetPersistCycles(), wiggleFirstLinkageOnly);
	logss << "1. Post WiggleFirst Overlaps Bead: " << this->CalculateOverlaps(BEAD) << ". Atomic: " << this->CalculateOverlaps(ATOMIC) << std::endl;
	std::cout << logss.str();
	if (randomize)
	{
		this->RandomDescent(BEAD, this->GetPersistCycles(), useMonteCarlo);
		logss << "2. Post Monte Carlo Overlaps Bead: " << this->CalculateOverlaps(BEAD) << ". Atomic: " << this->CalculateOverlaps(ATOMIC) << std::endl;
	    std::cout << logss.str();
	}
	this->Wiggle(BEAD, this->GetPersistCycles());
	logss << "3. Post Wiggle Overlaps Bead: " << this->CalculateOverlaps(BEAD) << ". Atomic: " << this->CalculateOverlaps(ATOMIC) << std::endl;
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
	logss << "Removing beads!\n";
    std::cout << logss.str();
	beads::Remove_Beads(this->GetGlycoproteinAssembly()); //Remove beads and write a final PDB & PRMTOP
	logss << "Beads removed!\n";
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
	this->CalculateOverlaps(ATOMIC, ALL, true, true);
	gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    return;
}

void GlycoproteinBuilder::RandomDescent(OverlapType overlapType, int persistCycles, bool monte_carlo)
{
	std::stringstream logss;
	logss << "Random Decent, persisting for " << persistCycles << " cycles and monte carlo is set as " << std::boolalpha << monte_carlo << ".\n";
	int cycle = 1;
	bool stop = false;
	bool record_overlap = false;
	double previous_glycan_overlap, new_glycan_overlap, previous_protein_overlap, new_protein_overlap;
	double lowest_global_overlap = this->GetGlobalOverlap();
	double new_global_overlap;
	double overlap_difference = 0.0;
	std::vector<GlycosylationSite*> sites_with_overlaps = this->DetermineSitesWithOverlap(this->GetOverlapTolerance(), overlapType);
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
			previous_glycan_overlap = current_glycosite->GetGlycanOverlap();
			previous_protein_overlap = current_glycosite->GetProteinOverlap();
			current_glycosite->SetRandomDihedralAnglesUsingMetadata();
			//logss << "Site: " << current_glycosite->GetResidueNumber() << "\n";
			new_glycan_overlap = current_glycosite->CalculateOverlaps(overlapType, GLYCAN, record_overlap);
			new_protein_overlap = current_glycosite->CalculateOverlaps(overlapType, PROTEIN, record_overlap);
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
		new_global_overlap = this->GetGlobalOverlap();
		sites_with_overlaps = this->DetermineSitesWithOverlap(this->GetOverlapTolerance(), overlapType); // Moved glycans may clash with other glycans. Need to check.
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

void GlycoproteinBuilder::Wiggle(OverlapType overlapType, int persistCycles, bool firstLinkageOnly, int interval)
{
	std::stringstream logss;
    if (overlapType == BEAD)
		logss << "GPWiggle BEAD ";
	if (overlapType == ATOMIC)
		logss << "GPwiggle ATOMIC ";
	logss << persistCycles << " persistCycles, " << interval << " interval\n";
	std::cout << logss.str();
    std::vector<GlycosylationSite*> sites_with_overlaps = this->DetermineSitesWithOverlap(this->GetOverlapTolerance(), overlapType);
    int cycle = 1;
    bool stop = false;
    double savedOverlap = this->GetGlobalOverlap();
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
            glycosite->Wiggle(overlapType, firstLinkageOnly, this->GetOverlapTolerance());
        }
        // Check which sites still have overlaps, stop if none.
        sites_with_overlaps = this->DetermineSitesWithOverlap(this->GetOverlapTolerance(), overlapType); // Moved glycans may clash with other glycans. Need to check.
        if (sites_with_overlaps.size() == 0)
        {
            logss << "Stopping with all overlaps resolved.\n";
            stop = true;
        }
        // The above will trigger re-calculation of overlaps, so can do this:
        if (savedOverlap > (this->GetGlobalOverlap() + 0.1)) // A significant improvement has happened.
        {
//        	std::stringstream ss;
//        	ss << "best_Wiggle_overlap" << this->GetGlobalOverlap();
//        	gmml::WritePDBFile(this->GetGlycoproteinAssembly(), this->GetWorkingDirectory(), ss.str());
        	logss << "Resetting cycle count to zero, will persist for another " << persistCycles << " cycles." << std::endl;
        	savedOverlap = this->GetGlobalOverlap();
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
    std::vector<GlycosylationSite*> sites_with_overlaps = DetermineSitesWithOverlap(this->GetOverlapTolerance());
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
        sites_with_overlaps = this->DetermineSitesWithOverlap(this->GetOverlapTolerance()); // Moved glycans may clash with other glycans. Need to check.
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

void GlycoproteinBuilder::CreateGlycosites(std::vector<GlycositeInput> glycositesInputVector)
{
	for (auto &glycositeInput : glycositesInputVector)
	{
	    gmml::log(__LINE__, __FILE__, gmml::INF, "Creating glycosite on residue " + glycositeInput.proteinResidueId_ + " with glycan " + glycositeInput.glycanInput_ );
	    std::cout << "Checking for residue" << std::endl;
	    Residue* glycositeResidue = this->SelectResidueFromInput(glycositeInput.proteinResidueId_);
	    if (glycositeResidue == nullptr)
	    {
	        std::cout << "Did not find ersidue" << std::endl;
	        throw std::runtime_error("Did not find a residue with id matching " + glycositeInput.proteinResidueId_);
	    }
	    std::vector<Residue*> otherResidues = this->GetGlycoproteinAssembly().getResidues();
	    otherResidues.erase(std::remove(otherResidues.begin(), otherResidues.end(), glycositeResidue), otherResidues.end());
		glycosites_.emplace_back(glycositeResidue, otherResidues, glycositeInput.glycanInput_, this->GetPrepFileLocation());
	    std::cout << "Done with glycan" << std::endl;
		gmml::log(__LINE__, __FILE__, gmml::INF, "Completed creating glycosite on residue " + glycositeInput.proteinResidueId_ + " with glycan " + glycositeInput.glycanInput_);
	}
    std::cout << "Done attaching all glycans" << std::endl;
    this->SetOtherGlycosites();
    gmml::log(__LINE__, __FILE__, gmml::INF, "Adding beads");
    this->Add_Beads(this->GetGlycoproteinAssembly(), this->GetGlycosites());
    this->UpdateAtomsThatMoveInLinkages(); // Must update to include beads.
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
    for (auto &residue : this->GetGlycoproteinAssemblyPtr()->getResidues())
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

void GlycoproteinBuilder::SetDefaultDihedralAnglesUsingMetadata()
{
    for(auto &glycosite : this->GetGlycosites())
    {
        glycosite.SetDefaultDihedralAnglesUsingMetadata();
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

double GlycoproteinBuilder::GetGlobalOverlap()
{
	//std::stringstream logss;
    //logss << "Calculating global overlap\n";
    double global_overlap = 0.0;
    for (std::vector<GlycosylationSite>::iterator current_glycosite = glycosites_.begin(); current_glycosite != glycosites_.end(); ++current_glycosite)
    {
        global_overlap += current_glycosite->GetOverlap();
        //logss << "Current site is " << current_glycosite->GetResidueNumber() << ", overlap: " << current_glycosite->GetOverlap() << ", making global: " << global_overlap << "\n";
    }
    //gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    return global_overlap;
}

//// random number generator; allows full range rotation
//double GlycoproteinBuilder::RandomAngle_360range()
//{
//    return (rand() % 360) + 1 - 180;
//}
/*double GlycoproteinBuilder::RandomAngle_range(int min, int max)
{
   // double angle = rand() % (max + 1 - min) + min;
    //std::cout << "Angle in range " << min << " - " << max << " is " << angle << "\n";
    // return angle;
    return rand() % (max + 1 - min) + min;
}*/
//// random number generator; specify a maximum step size relative to a start point
//double GlycoproteinBuilder::RandomAngle_PlusMinusX(double start_point, int max_step_size)
//{
//    return start_point + (rand() % (max_step_size * 2) + 1) - max_step_size;
//}
//double GlycoproteinBuilder::GetNewAngleScaledToDegreeOfOverlap(double current_angle, double overlap, int number_of_atoms)
//{
//    int max_step_size = 1 + std::round( 180 * ( overlap / number_of_atoms ) ); // Always allow at least 1 degrees of movement
//    return RandomAngle_PlusMinusX(current_angle, max_step_size);
//}



double GlycoproteinBuilder::CalculateOverlaps(OverlapType overlapType, MoleculeType moleculeType, bool recordOverlap, bool printOverlap)
{
    double overlap = 0.0;
    for(auto &glycosite : this->GetGlycosites())
    {
        overlap += glycosite.CalculateOverlaps(overlapType, moleculeType, recordOverlap, printOverlap);
    }
    return overlap;
}

std::vector<GlycosylationSite*> GlycoproteinBuilder::DetermineSitesWithOverlap(double tolerance, OverlapType overlapType)
{
    std::vector<GlycosylationSite*> sites_to_return;
    double overlap = 0.0;
    for (std::vector<GlycosylationSite>::iterator current_glycosite = glycosites_.begin(); current_glycosite != glycosites_.end(); ++current_glycosite)
    {
        overlap = current_glycosite->CalculateOverlaps(overlapType);
        if ( overlap > tolerance)
        {
            sites_to_return.push_back(&(*current_glycosite));
        }
    }
    return sites_to_return;
}

std::vector<GlycosylationSite*> GlycoproteinBuilder::GetSitesWithOverlap(double tolerance)
{
	//std::stringstream logss;
    std::vector<GlycosylationSite*> sites_to_return;
    double overlap = 0.0;
//  logss << "      Site        |  Total | Protein | Glycan " << std::endl;
    for (std::vector<GlycosylationSite>::iterator current_glycosite = glycosites_.begin(); current_glycosite != glycosites_.end(); ++current_glycosite)
    {
        overlap = current_glycosite->GetOverlap();
        if ( overlap > tolerance)
        {
//            logss << "Site " << current_glycosite->GetResidue()->GetId() << " is over tolerance with " << overlap << "\n";
            sites_to_return.push_back(&(*current_glycosite));
        }
    }
	//gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    return sites_to_return;
}

void GlycoproteinBuilder::DeleteSitesIterativelyWithAtomicOverlapAboveTolerance(std::vector<GlycosylationSite> &glycosites, double tolerance)
{
	std::stringstream logss;
    logss << "Atomic overlap before deleting sites is " << this->CalculateOverlaps(ATOMIC) << "\n";
    bool continue_deleting = true;
    // While overlap for any site is > tolerance delete site with highest overlap then re-calculate overlaps as there may be glycan-glycan overlap.
    while (continue_deleting)
    {
    	GlycosylationSite *worst_site = glycosites.data(); // Pointer to the first glycosite. Remember an erase/remove "advances"
    	for (std::vector<GlycosylationSite>::iterator current_glycosite = glycosites.begin(); current_glycosite != glycosites.end(); ++current_glycosite)
        {
        	if  ( current_glycosite->CalculateOverlaps(ATOMIC) > worst_site->CalculateOverlaps(ATOMIC))
            {
                worst_site = &(*current_glycosite); // The C is strong with this one.
            }
        }
        double worst_site_overlap = worst_site->CalculateOverlaps(ATOMIC);
        if ( worst_site_overlap > tolerance)
        {
            continue_deleting = true;
            //worst_site->Remove(this->GetGlycoproteinAssemblyPtr());
            worst_site->Rename_Protein_Residue_From_GLYCAM_To_Standard();
            logss << "Site " << worst_site->GetResidueId() << ": " << worst_site_overlap << " :" << "Removed\n";
            glycosites.erase(std::remove(glycosites.begin(), glycosites.end(), *worst_site), glycosites.end());
            this->SetOtherGlycosites(); // This needs to be updated after you delete each site, so the others no longer have a pointer to it. Great design Oliver well done you nailed it high five.
            this->Set_Other_Glycan_Beads(glycosites); // After "erasing", the actual atoms still exist and pointers to them are valid. Need to reset what beads are part of "other".
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

void GlycoproteinBuilder::UpdateAtomsThatMoveInLinkages()
{
    for (auto &glycosite : this->GetGlycosites())
    {
        glycosite.UpdateAtomsThatMoveInLinkages();
    }
    return;
}

//void GlycoproteinBuilder::StashCoordinates()
//{
//    for(auto &glycosite : this->GetGlycosites())
//    {
//        glycosite.StashCoordinates();
//    }
//    return;
//}

//void GlycoproteinBuilder::SetStashedCoordinatesWithLowestOverlap()
//{
//    for(auto &glycosite : this->GetGlycosites())
//    {
//        glycosite.SetStashedCoordinates();
//    }
//}

void GlycoproteinBuilder::Add_Beads(Assembly &glycoprotein, std::vector<GlycosylationSite> &glycosites)
{	// This function creates atoms with large radii called "beads" for the protein and glycans
    // The beads are added directly in glycoprotein and put in the AtomNode connection network
    // Within each glycosite an atomvector (vector of pointers to each bead) is added.
    // Each glycosite will have atomvectors of protein beads, glycan beads, and other glycan beads.
    // Other glycan beads are beads from glycans attached to other glycosites.
	std::vector<Atom*> proteinBeads = beads::Add_Beads_To_Protein(glycoprotein);

	for (auto &glycosite : glycosites)
	{ // Go through all glycosite glycans, add bead to each residue, attach it to one other atom in residue.
		glycosite.AddBeads(proteinBeads);
	} // Now find beads from other glycans and add them to list of other_glycan_beads for each glycosite
	this->Set_Other_Glycan_Beads(glycosites);
	return;
}

void GlycoproteinBuilder::Set_Other_Glycan_Beads(std::vector<GlycosylationSite> &glycosites)
{
    for (auto &glycosite1 : glycosites)
    {
    	std::vector<Atom*> other_glycan_beads;
        for (auto &glycosite2 : glycosites)
        {
            if(glycosite1 != glycosite2) // Check if same site
            {
                std::vector<Atom*> temp = glycosite2.GetSelfGlycanBeads();
                other_glycan_beads.insert(std::end(other_glycan_beads), std::begin(temp), std::end(temp));
                std::cout << "Adding beads of glycosite " << glycosite2.GetResidueId() << " to " << glycosite1.GetResidueId() << std::endl;
            }
        }
        glycosite1.SetOtherGlycanBeads(other_glycan_beads);
    }
    return;
}
