#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp" // calculateCoordinateFromInternalCoords
#include "includes/CentralDataStructure/Selections/atomSelections.hpp" //cdsSelections
#include "includes/CentralDataStructure/Editors/superimposition.hpp"
#include "includes/CentralDataStructure/Overlaps/beadResidues.hpp"
#include "includes/CodeUtils/templatedSelections.hpp"
#include "includes/CodeUtils/logging.hpp"
#include <bits/std_abs.h>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>
#include <iomanip> // For setting precision and formating in std::cout
#include <algorithm> //  std::erase, std::remove
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
GlycosylationSite::GlycosylationSite(Residue* residue, std::vector<Residue*> otherProteinResidues, std::string glycanInputString, std::string prepFileLocation) : residue_(residue),  glycan_(Carbohydrate(glycanInputString, prepFileLocation)), otherProteinResidues_(otherProteinResidues)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Attaching glycan!");
    this->AttachGlycan();
    std::cout << "Done attach glycan" << std::endl;

}
//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
void GlycosylationSite::AddBeads(std::vector<Atom*> proteinBeads)
{		// Add bead to each residue, attach it to one other atom in residue.
    Atom *cb_atom = this->GetResidue()->FindAtom("CB");
    if (cb_atom == nullptr)
    {
        throw "No CB atom in attachment residue. What da hell?";
    }
    double glycanLength = (cds::CalculateMaxDistanceBetweenCoordinates(cdsSelections::getCoordinates(this->GetAttachedGlycan()->getAtoms())) + 5); // added 5 to account for CB-C1(glycan) distance
    std::cout << "She's a smidge over: " << glycanLength << std::endl;
    this->SetProteinBeads(cdsSelections::FindAtomsWithinDistance(cb_atom, proteinBeads, glycanLength ));
    std::vector<Atom*> self_glycan_beads = beads::Add_Beads_To_Glycan(this->GetAttachedGlycan()->getResidues());
    std::cout << "I gots the self beads: " << self_glycan_beads.size() << std::endl;
    this->SetSelfGlycanBeads(&self_glycan_beads);
    std::cout << "Peace" << std::endl;
    return;
}

double GlycosylationSite::GetWeightedOverlap(double glycan_weight, double protein_weight)
{
    return ( (glycan_overlap_ * glycan_weight) + (protein_overlap_ * protein_weight) );
}

// Just to reproduce old code. Switch to distance based by residue and remove this crap later.
std::vector<Atom*> GlycosylationSite::GetProteinAtoms()
{
    std::vector<Atom*> foundAtoms;
    for (auto & residue : this->GetOtherProteinResidues())
    {
        for (auto &atom : residue->getAtoms())
        {
            foundAtoms.push_back(atom);
        }
    }
    return foundAtoms;
}
//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void GlycosylationSite::AttachGlycan()
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Start of AttachGlycan. Residue ID is: " + this->GetResidue()->getId());
	this->Prepare_Glycans_For_Superimposition_To_Particular_Residue(this->GetResidue()->getName());
    gmml::log(__LINE__, __FILE__, gmml::INF, "Superimpose prep done");
	this->Superimpose_Glycan_To_Glycosite(this->GetResidue());
    gmml::log(__LINE__, __FILE__, gmml::INF, "SuperimposedGlycanToGlycosite");
	this->Rename_Protein_Residue_To_GLYCAM_Nomenclature();
    gmml::log(__LINE__, __FILE__, gmml::INF, "Setting internal bond count to check if more form later");
	this->SetInternalBondCount(cdsSelections::CountInternalHeavyAtomBonds(this->GetAttachedGlycan()->getAtoms()));
    gmml::log(__LINE__, __FILE__, gmml::INF, "Attach glycan done");
}
// This function prepares the glycan molecule in the glycan_ assembly for superimpostion onto an amino acid in the protein
// It does this by "growing" the atoms of the amino acid side chain (e.g. Asn, Thr or Ser) out from the glycan reducing terminal
// Another function will use these additional atoms to superimpose the glycan onto residue
void GlycosylationSite::Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name)
{
    //Dear future self, the order that you add the atoms to the residue matters for superimposition ie N, CA, CB , not CB, CA, N.
    Residue* reducing_Residue = this->GetAttachedGlycan()->GetReducingResidue();
    // Want: residue.FindAtomByTag("anomeric-carbon"); The below is risky as it uses atoms names, i.e. would break for Sialic acid.
    // ToDo Ok so I reckon the below is just assuming alpha or beta depending on the concext. Need to fix a lot, but need to reproduce functionality after refactor first.
//    Atom* anomericAtom = cdsSelections::guessAnomericAtom(reducing_Residue);
//   This won't work as sometimes want alpha, sometimes beta. i.e. a CreateCoordinateForCenterAwayFromNeighbors function
    // This needs to be abstracted so it works for C2 reducing residues:
    Coordinate *coordC5 = reducing_Residue->FindAtom("C5")->getCoordinate();
    Coordinate *coordO5 = reducing_Residue->FindAtom("O5")->getCoordinate();
    Coordinate *coordC1 = reducing_Residue->FindAtom("C1")->getCoordinate();
    Atom* anomericAtom = reducing_Residue->FindAtom("C1"); // For adding bond.
    // Delete aglycon atoms from glycan.
    Residue * aglycon = this->GetAttachedGlycan()->GetAglycone();
    for (auto &atom : aglycon->getAtoms())
    {
       aglycon->deleteAtom(atom);
    }
    // Ok so going to set it so that the new "superimposition residue" is the old aglycon residue
    // This avoids having to delete the algycon residue object from assembly and adding the super residue to assembly.
    Residue* superimposition_residue = aglycon; // "renaming" so the below reads better.
    superimposition_residue->setName("SUP");
    // I put both the regular name and the O/N-linked glycam name here, as I'm not sure when it will be renamed.
    if ( (amino_acid_name == "ASN") || (amino_acid_name == "NLN") )
    {
        Atom *atomND2 = superimposition_residue->addAtom(std::make_unique<Atom>("ND2", (cds::calculateCoordinateFromInternalCoords(*coordC5, *coordO5, *coordC1, 109.3, 180, 1.53))));
        Atom *atomCG = superimposition_residue->addAtom(std::make_unique<Atom>("CG", (cds::calculateCoordinateFromInternalCoords(*coordO5, *coordC1, *atomND2->getCoordinate(), 109.3, 261, 1.325))));
        superimposition_residue->addAtom(std::make_unique<Atom>("OD1", (cds::calculateCoordinateFromInternalCoords(*coordC1, *atomND2->getCoordinate(), *atomCG->getCoordinate(), 126, 0, 1.22))));
        anomericAtom->addBond(atomND2); // This is so findAnomericAtom works later, needs a foreign residue neighbor.
    }
    else if ( (amino_acid_name.compare("THR")==0) || (amino_acid_name.compare("SER")==0) || (amino_acid_name.compare("OLT")==0) || (amino_acid_name.compare("OLS")==0) )
    {
        Atom *atomOG1 = superimposition_residue->addAtom(std::make_unique<Atom>("OG", (cds::calculateCoordinateFromInternalCoords(*coordC5, *coordO5, *coordC1, 112, 68, 1.46))));
        Atom *atomCB = superimposition_residue->addAtom(std::make_unique<Atom>("CB", (cds::calculateCoordinateFromInternalCoords(*coordO5, *coordC1, *atomOG1->getCoordinate(), 109.3, 75, 1.53))));
        superimposition_residue->addAtom(std::make_unique<Atom>("CA", (cds::calculateCoordinateFromInternalCoords(*coordC1, *atomOG1->getCoordinate(), *atomCB->getCoordinate(), 109.3, 125, 1.53))));
        if ( (amino_acid_name.compare("THR")==0) || (amino_acid_name.compare("OLT")==0) )
        {
            atomOG1->setName("OG1"); // It's OG in Ser.
        }
        anomericAtom->addBond(atomOG1); // This is so findAnomericAtom works later, needs a foreign residue neighbor.
    }
    else if ( (amino_acid_name.compare("TYR")==0) || (amino_acid_name.compare("OLY")==0) )
    {
        Atom *atomOH = superimposition_residue->addAtom(std::make_unique<Atom>("OH", (cds::calculateCoordinateFromInternalCoords(*coordC5, *coordO5, *coordC1, 112, 68, 1.46))));
        Atom *atomCZ = superimposition_residue->addAtom(std::make_unique<Atom>("CZ", (cds::calculateCoordinateFromInternalCoords(*coordO5, *coordC1, *atomOH->getCoordinate(), 117, 60, 1.35))));
        superimposition_residue->addAtom(std::make_unique<Atom>("CE1", (cds::calculateCoordinateFromInternalCoords(*coordC1, *atomOH->getCoordinate(), *atomCZ->getCoordinate(), 120, 180, 1.37))));
        anomericAtom->addBond(atomOH); // This is so findAnomericAtom works later, needs a foreign residue neighbor.
    }
    else
    {
        std::string message = "Problem creating glycosylation site. The amino acid requested: " + amino_acid_name + " has name that isn't supported. Currently you can glycosylate ASN, THR, SER or TYR. Email us to request others. Ideally include examples of 3D structures we can use as a template.";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    return;
}

void GlycosylationSite::Superimpose_Glycan_To_Glycosite(Residue *glycosite_residue)
{
    // Get the 3 target atoms from protein residue.
    std::vector<Coordinate*> targetCoords;
   // superimposition_atoms_ points to three atoms that were added to the glycan. Based on their names e.g. CG, ND2, we will superimpose them onto
   // the correspoinding "target" atoms in the protein residue (glycosite_residue).
    for (auto &superimposition_atom : this->GetAttachedGlycan()->GetAglycone()->getAtoms())
    {
        std::cout << "Superimposition aglycone atom is named " << superimposition_atom->getName() << "\n";
        for(auto &protein_atom : glycosite_residue->getAtoms())
        {
            if (protein_atom->getName() == superimposition_atom->getName())
            {
                std::cout << "Adding " << protein_atom->getName() << " to superimposition atoms\n";
                targetCoords.push_back(protein_atom->getCoordinate());
            }
        }
    }
    std::cout << "\n\n\n";
    std::vector<Coordinate*> aglyconeCoords = cdsSelections::getCoordinates(this->GetAttachedGlycan()->GetAglycone()->getAtoms());
    std::cout << "\n\n\n";
    std::vector<Coordinate*> glycanCoords = cdsSelections::getCoordinates(this->GetAttachedGlycan()->getAtoms());
    std::cout << "\n\n\n";
    std::cout << "Number of moving coords: " << glycanCoords.size() << "vs" << this->GetAttachedGlycan()->getAtoms().size() << "\n";
    std::cout << "Number of aglycone coords: " << aglyconeCoords.size() << " vs " << this->GetAttachedGlycan()->GetAglycone()->getAtoms().size() << "\n";
    std::cout << "Number of target coords: " << targetCoords.size() << "\n";
    for(auto &coord : aglyconeCoords)
    {
        std::cout << coord->ToString() << "\n";
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Superimposing via the aglycone.");
    cds::Superimpose(aglyconeCoords, targetCoords, glycanCoords);
    //Connect the glycan and protein atoms to each other.
    gmml::log(__LINE__, __FILE__, gmml::INF, "Connecting the aglycone to the protein");
    Atom *protein_connection_atom = this->GetConnectingProteinAtom(glycosite_residue->getName());
    protein_connection_atom->addBond(this->GetAttachedGlycan()->GetAnomericAtom());
    gmml::log(__LINE__, __FILE__, gmml::INF, "Deleting the aglycone");
    this->GetAttachedGlycan()->deleteResidue(this->GetAttachedGlycan()->GetAglycone());
    gmml::log(__LINE__, __FILE__, gmml::INF, "Completed superimposition to " + glycosite_residue->getId());
    return;
}

void GlycosylationSite::Rename_Protein_Residue_To_GLYCAM_Nomenclature()
{
    std::string amino_acid_name = this->GetResidue()->getName();
    if (amino_acid_name == "ASN") {this->GetResidue()->setName("NLN");}
    if (amino_acid_name == "SER") {this->GetResidue()->setName("OLS");}
    if (amino_acid_name == "THR") {this->GetResidue()->setName("OLT");}
    if (amino_acid_name == "TYR") {this->GetResidue()->setName("OLY");}
}

void GlycosylationSite::Rename_Protein_Residue_From_GLYCAM_To_Standard()
{
    std::string amino_acid_name = this->GetResidue()->getName();
    if (amino_acid_name == "NLN") {this->GetResidue()->setName("ASN");}
    if (amino_acid_name == "OLS") {this->GetResidue()->setName("SER");}
    if (amino_acid_name == "OLT") {this->GetResidue()->setName("THR");}
    if (amino_acid_name == "OLY") {this->GetResidue()->setName("TYR");}
}

double GlycosylationSite::CalculateOverlaps(OverlapType overlapType, MoleculeType moleculeType, bool recordOverlap, bool printOverlap)
{
    double overlap = 0.0;
    if (moleculeType == ALL) // glycan and protein overlap are calculated and stored sepatately. They get combined only when reporting total overlaps. Combined value is never stored.
    {
        double proteinOverlap = this->CalculateOverlaps(overlapType, PROTEIN, recordOverlap, printOverlap);
        double glycanOverlap = this->CalculateOverlaps(overlapType, GLYCAN, recordOverlap, printOverlap);
        overlap = proteinOverlap + glycanOverlap;
        if (printOverlap)
        {
        	this->PrintOverlaps();
        }
    }
    else
    {
        if(overlapType == ATOMIC)
        {
            overlap = this->CalculateAtomicOverlaps(moleculeType, printOverlap);
        }
        else if (overlapType == BEAD)
        {
            std::cout << "Nah here" << std::endl;
            overlap = this->CalculateBeadOverlaps(moleculeType);
        }
        if (recordOverlap)
        {
            this->SetOverlap(moleculeType, overlap);
        }
    }
    return overlap;
}

double GlycosylationSite::CalculateAtomicOverlaps(MoleculeType moleculeType, bool print)
{
    double overlap = 0.0;
    if(moleculeType == ALL)
    {
        overlap = ( this->CalculateAtomicOverlaps(PROTEIN) + this->CalculateAtomicOverlaps(GLYCAN) );
    }
    if(moleculeType == PROTEIN)
    {
        overlap = cds::CalculateAtomicOverlaps(this->GetProteinAtoms(), this->GetAttachedGlycan()->getAtoms(), print);
    }
    if(moleculeType == GLYCAN)
    {
        for(auto &other_glycosite : other_glycosites_)
        {
        	double current_overlap = cds::CalculateAtomicOverlaps(other_glycosite->GetAttachedGlycan()->getAtoms(), this->GetAttachedGlycan()->getAtoms(), print);
            overlap += current_overlap;
            overlap += cds::CalculateAtomicOverlaps(other_glycosite->GetAttachedGlycan()->getAtoms(), this->GetAttachedGlycan()->getAtoms(), print);
        }
    }
    return overlap;
}

void GlycosylationSite::PrintOverlaps()
{
    std::stringstream logss;
    logss << std::fixed; // Formating ouput
    logss << std::setprecision(2); // Formating ouput
    logss << std::setw(17) << this->GetResidue()->getId() << " | "
        << std::setw(6) << this->GetOverlap() << " |  "
        << std::setw(6) << this->GetProteinOverlap() << " | "
        << std::setw(6) << this->GetGlycanOverlap() << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());    
}

double GlycosylationSite::CalculateBeadOverlaps(MoleculeType moleculeType)
{
    double overlap = 0.0;
    if(moleculeType == ALL)
    {
        overlap = ( this->CalculateBeadOverlaps(PROTEIN) + this->CalculateBeadOverlaps(GLYCAN) );
    }
    if(moleculeType == PROTEIN)
    {
        overlap = this->CalculateBeadOverlaps(self_glycan_beads_, protein_beads_);
    }
    if(moleculeType == GLYCAN)
    {
        overlap = this->CalculateBeadOverlaps(self_glycan_beads_, other_glycan_beads_);
    }
    return overlap;
}

double GlycosylationSite::CalculateBeadOverlaps(std::vector<Atom*> &atomsA, std::vector<Atom*> &atomsB)
{
    double radius = 3.0; //Using same radius for all beads.
    double overlap = 0.0, current_overlap = 0.0;
    if (atomsA.empty() || atomsB.empty() )
    {
        throw std::runtime_error("Something very wrong in CalculateBeadOverlaps");
    }
    for(std::vector<Atom*>::iterator it1 = atomsA.begin(); it1 != atomsA.end(); ++it1)
    {
        Atom *atomA = *it1;
        for(std::vector<Atom*>::iterator it2 = atomsB.begin(); it2 != atomsB.end(); ++it2)
        {
            Atom *atomB = *it2;
            if ( (atomA->getCoordinate()->GetX() - atomB->getCoordinate()->GetX()) < (radius * 2) ) // This is faster than calulating distance, and rules out tons of atom pairs.
            {
                if (atomA->getCoordinate()->withinDistance(atomB->getCoordinate(), (radius + radius)) )
                {
                    current_overlap = cds::CalculateAtomicOverlaps(atomA, atomB, radius, radius); // This calls the version with radius values
                    overlap += current_overlap;
                    //std::cout << atomA->GetResidue()->getId() << " overlaping with " << atomB->GetResidue()->getId() << ": " << current_overlap << "\n";
                }
            }
        }
    }
    return (overlap / constants::CARBON_SURFACE_AREA); //Normalise to area of a buried carbon
}

void GlycosylationSite::Wiggle(OverlapType overlapType, bool firstLinkageOnly, double tolerance, int interval)
{ // I want to find the lowest overlap as close to each bonds default as possible. So code is a bit more complicated.
	if (firstLinkageOnly)
	{
		this->WiggleOneLinkage(all_residue_linkages_.front(), overlapType, tolerance, interval);
	}
	else
	{
		for(auto &linkage : all_residue_linkages_)
		{
			this->WiggleOneLinkage(linkage, overlapType, tolerance, interval);
		}
	}
    return;
}
//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void GlycosylationSite::SetOverlap(MoleculeType moleculeType, double overlap)
{
    switch (moleculeType)
    {
    case PROTEIN:
        this->SetProteinOverlap(overlap);
        break;
    case GLYCAN:
        this->SetGlycanOverlap(overlap);
        break;
    default:
        std::cerr << "ERROR, No type specified in GlycosylationSite::SetOverlap. Required.\n";
        break;
    }
}

// Should just pass a vec of beads that doesn't have the residue in them. The erase is expensive.
void GlycosylationSite::SetProteinBeads(std::vector<Atom*> beads)
{
    //Remove beads from attachment point residue. Don't want to count overlaps between glycan and residue it is attached to.
    std::vector<std::string> beadAtomNames = {"mfat", "sfat"};
    std::vector<Atom*> myProteinBeadAtoms = codeUtils::getElementsWithNames(this->GetProteinAtoms(), beadAtomNames);
    std::cout << "I HAVE THIS MANY BEASDS::::::::::::::::::::;" << myProteinBeadAtoms.size() << std::flush << std::endl;
    protein_beads_ = codeUtils::findElementsNotInVector(beads, myProteinBeadAtoms);
    std::cout << "which is now HAVE THIS MANY BEASDS::::::::::::::::::::;" << protein_beads_.size() << std::flush << std::endl;
}

void GlycosylationSite::SetDefaultDihedralAnglesUsingMetadata()
{
    for(auto &linkage : all_residue_linkages_)
    {
        linkage.SetDefaultShapeUsingMetadata();
    }
    if( ! this->NoNewInternalCloseContacts() )
    {
        this->ResetDihedralAngles();
    }
    return;
}

void GlycosylationSite::SetRandomDihedralAnglesUsingMetadata()
{
    for(auto &linkage : all_residue_linkages_)
    {
        linkage.SetRandomShapeUsingMetadata();
    }
    if( ! this->NoNewInternalCloseContacts() )
    {
        this->ResetDihedralAngles();
    }
    return;
}

void GlycosylationSite::SetRandomDihedralAnglesUsingMetadataForNthLinkage(long unsigned int linkage_number)
{
    if(linkage_number < all_residue_linkages_.size())
    {
        all_residue_linkages_.at(linkage_number).SetRandomShapeUsingMetadata();
    }
    if( ! this->NoNewInternalCloseContacts() )
    {
        this->ResetDihedralAngles();
    }
    return;
}

void GlycosylationSite::ResetDihedralAngles()
{
    for(auto &linkage : all_residue_linkages_)
    {
        linkage.SetShapeToPrevious();
    }
    return;
}

void GlycosylationSite::UpdateAtomsThatMoveInLinkages()
{
    for(auto &linkage : all_residue_linkages_)
    {
        linkage.DetermineAtomsThatMove();
    }
    return;
}

//void GlycosylationSite::StashCoordinates()
//{
//    std::cout << "Stashing coordinates for glycosite " << this->GetResidueNumber() << std::endl;
//    std::vector<Atom*> atoms = this->GetAttachedGlycan()->GetAllAtomsOfAssembly();
//    std::vector<Atom*> sidechain_atoms = this->GetResidue()->GetAtoms();
//    atoms.insert(atoms.end(), sidechain_atoms.begin(), sidechain_atoms.end() );
//    for(auto &atom : atoms)
//    {
//        atom->AddCoordinate(new GeometryTopology::Coordinate(atom->GetCoordinate())); // push back the currect coordinates onto the end of the coordinates vector.
//    }
//    return;
//}
//
//void GlycosylationSite::SetStashedCoordinates() // When a lower overlap is found, the coords are pushed back. Last one pushed back is the lowest overlap.
//{
//    std::cout << "Setting stashed coordinates for site " << this->GetResidueNumber() << "\n";
//    std::vector<Atom*> atoms = this->GetAttachedGlycan()->GetAllAtomsOfAssembly();
//    std::vector<Atom*> sidechain_atoms = this->GetResidue()->GetAtoms();
//    atoms.insert(atoms.end(), sidechain_atoms.begin(), sidechain_atoms.end() );
//    for(auto &atom : atoms)
//    {
//        GeometryTopology::Coordinate *first_coords = atom->GetCoordinate();
//        GeometryTopology::Coordinate *last_coords = atom->GetCoordinates().back();
//        *first_coords = *last_coords;
//    }
//    return;
//}

//std::vector<GlycosylationSite> GlycosylationSite::GetXClosestSitesWithinOverlapDistanceY(std::vector<GlycosylationSite> &glycosites, int maxNumberOfSitesToConsider)
//{
//    typedef std::pair<double, GlycosylationSite> siteDistancePair;
//    typedef std::vector<siteDistancePair> siteDistancePairVector;
//    siteDistancePairVector sitesWithinOverlapRange;
//    std::vector<GlycosylationSite> sitesToReturn;
//    double maxLengthOfThisGlycan = (selection::GetMaxDistanceBetweenAtoms(this->GetAttachedGlycan()->GetAllAtomsOfAssembly()) + 5);
//    Atom *this_cb_atom = this->GetResidue()->GetAtom("CB");
//    for(auto &glycosite : glycosites)
//    {
//        if (this->GetResidue()->getId().compare(glycosite.GetResidue()->getId())==0) // if not this site. Oly you should overload the = operator?
//        {
//            continue; // I think this skips to next item in for loop.
//        }
//        double maxLengthOfThatGlycan = (selection::GetMaxDistanceBetweenAtoms(glycosite.GetAttachedGlycan()->GetAllAtomsOfAssembly()) + 5);
//        Atom *that_cb_atom = glycosite.GetResidue()->GetAtom("CB");
//        double distanceBetweenCBAtoms = this_cb_atom->GetDistanceToAtom(that_cb_atom);
//        std::cout << distanceBetweenCBAtoms;
//        if ( ( distanceBetweenCBAtoms ) <= (maxLengthOfThisGlycan + maxLengthOfThatGlycan) )
//        {
//            sitesWithinOverlapRange.emplace_back(distanceBetweenCBAtoms, glycosite);
//        }
//    }
//    // Sort list by distance. Remember that you might have more or less sites than the maxNumberOfSitesToConsider
//    std::sort(sitesWithinOverlapRange.begin(), sitesWithinOverlapRange.end());
//    // Return up to maxNumberOfSitesToConsider
//    for(auto &pairInstance : sitesWithinOverlapRange)
//    {
//        if (sitesToReturn.size() < maxNumberOfSitesToConsider)
//        {
//            sitesToReturn.push_back(pairInstance.second);
//        }
//    }
//    return sitesToReturn;
//}

void GlycosylationSite::SetOtherGlycosites(std::vector<GlycosylationSite> &glycosites)
{
	other_glycosites_.clear();
    for(auto &glycosite : glycosites)
    {
        if(this->GetResidue()->getId() != glycosite.GetResidue()->getId())
        {
            other_glycosites_.push_back(&glycosite);
        }
    }
    return;
}
//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void GlycosylationSite::Print(std::string type)
{
    std::stringstream logss;
    if (type.compare("All")==0)
    {
        logss << "Residue ID: " << this->GetResidue()->getId() << ", overlap: " << this->GetOverlap();
        logss << ", Gly--Pro " << this->GetGlycanOverlap() << "--" << this->GetProteinOverlap() << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    }
}
//////////////////////////////////////////////////////////
//                   PRIVATE FUNCTIONS                 //
//////////////////////////////////////////////////////////
Atom* GlycosylationSite::GetConnectingProteinAtom(std::string residue_name)
{
    if(residue_name == "NLN" || residue_name == "ASN")
    {
        return this->GetResidue()->FindAtom("ND2");
    }
    else if(residue_name == "OLT" || residue_name == "THR")
    {
        return this->GetResidue()->FindAtom("OG1");
    }
    else if(residue_name == "OLS" || residue_name == "SER")
    {
        return this->GetResidue()->FindAtom("OG");
    }
    else if(residue_name == "OLY" || residue_name == "TYR")
    {
        return this->GetResidue()->FindAtom("OH");
    }
    else
    {
        std::string message = "Problem in GetConnectingProteinAtom. The amino acid requested: " + residue_name + " has name that isn't supported. Currently you can glycosylate ASN, THR, SER or TYR. Email us to request others. Ideally include examples of 3D structures we can use as a template.";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
}

// OG re-reading. I'm pretty sure this belongs in Residue_linkage. linkage.wiggle(output_pdb_id, tolerance, interval
// OG re-re-reading. It's because you need to calculate overlaps.
void GlycosylationSite::WiggleOneLinkage(ResidueLinkage &linkage, OverlapType overlapType, double tolerance, int interval)
{
    double current_overlap = this->CalculateOverlaps(overlapType);
    double lowest_overlap = current_overlap;
    // Reverse as convention is Glc1-4Gal and I want to wiggle in opposite direction i.e. from first rotatable bond in Asn outwards
    std::vector<RotatableDihedral> reversed_rotatable_bond_vector = linkage.GetRotatableDihedrals();
    std::reverse(reversed_rotatable_bond_vector.begin(), reversed_rotatable_bond_vector.end());
    for(auto &rotatable_dihedral : reversed_rotatable_bond_vector)
    {
        double best_dihedral_angle = rotatable_dihedral.CalculateDihedralAngle();
       // std::cout << "Starting new rotatable dihedral with best angle as " << best_dihedral_angle << "\n";
        gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata_entries = rotatable_dihedral.GetMetadata();
        for(auto &metadata : metadata_entries)
        {
            double lower_bound = (metadata.default_angle_value_ - metadata.lower_deviation_);
            double upper_bound = (metadata.default_angle_value_ + metadata.upper_deviation_);
            double current_dihedral = lower_bound;
            while(current_dihedral <= upper_bound )
            {
                rotatable_dihedral.SetDihedralAngle(current_dihedral);
                current_overlap = this->CalculateOverlaps(overlapType);
          //      std::cout << this->GetResidueNumber() << ": current dihedral : overlap " << current_dihedral << " : " << current_overlap << ". Best dihedral : overlap: " << best_dihedral_angle << " : "<< lowest_overlap << "\n";
                if (lowest_overlap >= (current_overlap + 0.01)) // 0.01 otherwise rounding errors
                {
                	if (this->NoNewInternalCloseContacts())
                	{
                		lowest_overlap = current_overlap;
                		best_dihedral_angle = current_dihedral;
//                		std::stringstream ss;
//                		ss << "betterOverlap_" << lowest_overlap << "_";
                		std::cout << "Wiggler: site " << this->GetResidue()->getNumber() << " has overlap: " << lowest_overlap << "\n";
                		//gmml::WritePDBFile(*(this->GetResidue()->GetAssembly()), "", "best_");
            //    		std::cout << "Best angle is now " << best_dihedral_angle << "\n";
                	}
                }
                // Perfer angles closer to default.
                else if ( (lowest_overlap == current_overlap) &&
                          (abs(metadata.default_angle_value_ - best_dihedral_angle ) > abs(metadata.default_angle_value_ - current_dihedral)))
                {
                	if (this->NoNewInternalCloseContacts())
                	{
                		best_dihedral_angle = current_dihedral;
                	}
                }
                current_dihedral += interval; // increment
            }
        }
        //std::cout << "Setting best angle as " << best_dihedral_angle << "\n";
        rotatable_dihedral.SetDihedralAngle(best_dihedral_angle);
        if(lowest_overlap <= tolerance) return;
    }
    return; // Note possibility of earlier return above
}

bool GlycosylationSite::NoNewInternalCloseContacts()
{
    unsigned long int newCount = cdsSelections::CountInternalHeavyAtomBonds(this->GetAttachedGlycan()->getAtoms());
	if (newCount > this->GetInternalBondCount())
	{
		return false;
	}
	return true;
}
