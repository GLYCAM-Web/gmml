#include <bits/std_abs.h>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include "includes/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/InternalPrograms/functionsForGMML.hpp"
#include "includes/InternalPrograms/beadResidues.hpp"
#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/GeometryTopology/geometrytopology.hpp"
#include "includes/GeometryTopology/ResidueLinkages/residue_linkage.hpp"
#include "includes/GeometryTopology/ResidueLinkages/rotatable_dihedral.hpp"
#include "includes/InternalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/MolecularModeling/assembly.hpp"
#include "includes/MolecularModeling/atom.hpp"
#include "includes/MolecularModeling/atomnode.hpp"
#include "includes/MolecularModeling/overlaps.hpp"
#include "includes/MolecularModeling/residue.hpp"
#include "includes/MolecularModeling/superimposition.hpp"
#include "includes/MolecularModeling/Selections/selections.hpp"
#include "includes/CodeUtils/logging.hpp"

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
GlycosylationSite::GlycosylationSite(Assembly* glycoprotein, std::string residueNumber, std::string glycanInputType, std::string glycanInput, std::string prepFileLocation)
{
    this->SetGlycanName(glycanInput);
    this->SetResidueNumber(residueNumber);
    this->SetGlycanOverlap(0.0);
    this->SetProteinOverlap(0.0);
    gmml::log(__LINE__, __FILE__, gmml::INF, "Finding protein residue: " + residueNumber);
    this->FindSetProteinResidue(residueNumber, glycoprotein->GetResidues()); // FindProteinResidue should be in Assembly? This does both finding and setting.
    gmml::log(__LINE__, __FILE__, gmml::INF, "Setting assembly for residue.");
    this->GetResidue()->SetAssembly(glycoprotein); // gawd-dern you gmml.
    gmml::log(__LINE__, __FILE__, gmml::INF, "Determining overlap atoms in protein");
    this->DetermineOverlapProteinAtoms(glycoprotein, this->GetResidue());
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building glycan!");
    this->BuildGlycan(glycanInputType, glycanInput, prepFileLocation);
    gmml::log(__LINE__, __FILE__, gmml::INF, "Attaching glycan!");
    this->AttachGlycan(glycoprotein);
    gmml::log(__LINE__, __FILE__, gmml::INF, "Added " + this->GetGlycanName() + " to " + this->GetResidueNumber());
}
//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
void GlycosylationSite::FindSetProteinResidue(std::string residueName, ResidueVector residues)
{
    std::stringstream logss;
	for (auto &residue : residues)
	{
		std::string formattedResidueNumber = "_" + residueName + "_";
		if( residue->GetId().compare(3, formattedResidueNumber.size(), formattedResidueNumber) == 0)
		{
			logss << "glycosite id:" << residue->GetId() << std::endl;
            gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
			this->SetResidue(residue);
			return;
		}
	}
    logss << "Error in glycosylation site class, did not find protein residue: " + residueName + " in 3D structure";
    gmml::log(__LINE__, __FILE__, gmml::ERR, logss.str());
    throw logss.str();
	return;
}

void GlycosylationSite::AddBeads(AtomVector proteinBeads)
{		// Add bead to each residue, attach it to one other atom in residue.
		Atom *cb_atom = this->GetResidue()->GetAtom("CB");
		double distance = (selection::GetMaxDistanceBetweenAtoms(this->GetAttachedGlycan()->GetAllAtomsOfAssembly()) + 5); // added 5 to account for CB-C1(glycan) distance
		AtomVector close_protein_beads = selection::AtomsWithinDistanceOf(cb_atom, distance, proteinBeads);
		this->SetProteinBeads(&close_protein_beads);
		AtomVector self_glycan_beads = beads::Add_Beads_To_Glycan(this->GetAttachedGlycan()->GetResidues());
		this->SetSelfGlycanBeads(&self_glycan_beads);
}

void GlycosylationSite::Remove(Assembly *glycoproteinAssembly)
{
	this->Rename_Protein_Residue_From_GLYCAM_To_Standard();
	for(auto &residue : this->GetAttachedGlycan()->GetResidues())
	{
		glycoproteinAssembly->RemoveResidue(residue);
	}
}

void GlycosylationSite::BuildGlycan(std::string glycanInputType, std::string glycanInput, std::string prepFileLocation)
{
	if (glycanInputType == "Sequence")
	{
	    gmml::log(__LINE__, __FILE__, gmml::INF, "Generating assembly from sequence with " + glycanInput);
	    gmml::log(__LINE__, __FILE__, gmml::INF, "Prepfile location is: " + prepFileLocation);
		CondensedSequence::carbohydrateBuilder carbBuilder(glycanInput, prepFileLocation);
		carbBuilder.SetDefaultShapeUsingMetadata();
		carbBuilder.ResolveOverlaps();
		this->SetGlycan(carbBuilder.GetAssembly());
	}
	else if (glycanInputType == "Library")
	{
	    gmml::log(__LINE__, __FILE__, gmml::INF, "Generating assembly from PDB with " + glycanInput);
		MolecularModeling::Assembly glycan(glycanInput, gmml::InputFileType::PDB);
		glycan.BuildStructureByDistance();
		this->SetGlycan(glycan);
	}
	else
	{
		throw "Error in GlycosylationSite class: the glycanInputType must be either \"Library\" or \"Sequence\"";
	}
	return;
}

double GlycosylationSite::GetWeightedOverlap(double glycan_weight, double protein_weight)
{
    return ( (glycan_overlap_ * glycan_weight) + (protein_overlap_ * protein_weight) );
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
// Only need glycoprotein so can merge assemblies and set bonding for the connecting atom
// Bond by distance wouldn't work as may have overlaps after superimposition.
void GlycosylationSite::AttachGlycan(Assembly* glycoprotein)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Start of AttachGlycan");
	this->Prepare_Glycans_For_Superimposition_To_Particular_Residue(residue_->GetName());
    gmml::log(__LINE__, __FILE__, gmml::INF, "Superimpose prep done");
	this->Superimpose_Glycan_To_Glycosite(residue_);
    gmml::log(__LINE__, __FILE__, gmml::INF, "Residue ID is: " + this->GetResidue()->GetId());
	this->RenumberGlycanToMatch(*glycoprotein);
    gmml::log(__LINE__, __FILE__, gmml::INF, "SuperimposedGlycanToGlycosite");
	this->Rename_Protein_Residue_To_GLYCAM_Nomenclature();
	gmml::WritePDBFile(glycan_,"./" , this->GetResidue()->GetId() + "_glycan", false);
	glycoprotein->MergeAssembly(&glycan_); // Add glycan to glycoprotein assembly, allows SetDihedral later. May not be necessary anymore with new Rotatable Dihedral class.
	gmml::log(__LINE__, __FILE__, gmml::INF, "Merge done");
	all_residue_linkages_.emplace_back(glycan_.GetResidues().at(0), residue_);
    gmml::log(__LINE__, __FILE__, gmml::INF, "Figuring out residue linkages");
	this->FigureOutResidueLinkagesInGlycan(glycan_.GetResidues().at(0), glycan_.GetResidues().at(0), &all_residue_linkages_);
    gmml::log(__LINE__, __FILE__, gmml::INF, "Setting internal bond count to check if more form later");
	this->SetInternalBondCount(gmml::CountInternalBonds(glycan_));
    gmml::log(__LINE__, __FILE__, gmml::INF, "Attach glycan done");
}

/*
 * This function prepares the glycan molecule in the glycan_ assembly for superimpostion onto an amino acid in the protein
 * It does this by "growing" the atoms of the amino acid side chain (e.g. Asn, Thr or Ser) out from the glycan reducing terminal
 * Another function will use these additional atoms to superimpose the glycan onto residue_
 * This function assumes that the glycan_ assembly for this glycosylation site has already been set
*/
void GlycosylationSite::Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name)
{
    //Dear future self, the order that you add the atoms to the residue matters for superimposition ie N, CA, CB , not CB, CA, N.
    // Want: assembly.FindResidueByTag("reducing-residue");
    Residue* reducing_Residue = glycan_.GetAllResiduesOfAssembly().at(1); // I assume I assumed something stupid here.
    // Want: residue.FindAtomByTag("anomeric-carbon"); The below is risky as it uses atoms names, i.e. would break for Sialic acid.
    Atom *atomC5 = reducing_Residue->GetAtom("C5");
    Atom *atomO5 = reducing_Residue->GetAtom("O5");
    Atom *atomC1 = reducing_Residue->GetAtom("C1");
    // Delete aglycon atoms from glycan.
    Residue * aglycon = glycan_.GetAllResiduesOfAssembly().at(0); // Oh jeez these assumptions are really building up.
    AtomVector aglycon_Atoms = aglycon->GetAtoms();
    for(AtomVector::iterator it = aglycon_Atoms.begin(); it != aglycon_Atoms.end(); ++it)
    {
       Atom* atom = *it;
       aglycon->RemoveAtom(atom); // Note only removes from residue. Atoms still exist and can be found through AtomNodes. They aren't written out in a PDB file though.
    }
    // Ok so going to set it so that the new "superimposition residue" is the old aglycon residue i.e. .at(0)
    // This avoids having to delete the algycon residue object from assembly and adding the super residue to assembly.
    // Deleting the residue is actually hard as the atoms still exist and are referenced from other places.
    Residue* superimposition_residue = aglycon; // "renaming" so the below reads better.
    superimposition_residue->SetName("SUP");
    superimposition_residue->SetId("SUP_?_1_?_?_1");
    // I put both the regular name and the O/N-linked glycam name here, as I'm not sure when it will be renamed.
    if ( (amino_acid_name.compare("ASN")==0) || (amino_acid_name.compare("NLN")==0) )
    {
        Atom *atomND2 = new Atom(superimposition_residue, "ND2", (GeometryTopology::get_cartesian_point_from_internal_coords(atomC5, atomO5, atomC1, 109.3, 180, 1.53)));
        Atom *atomCG = new Atom(superimposition_residue, "CG", (GeometryTopology::get_cartesian_point_from_internal_coords(atomO5, atomC1, atomND2, 109.3, 261, 1.325)));
        Atom *atomOD1 = new Atom(superimposition_residue, "OD1", (GeometryTopology::get_cartesian_point_from_internal_coords(atomC1, atomND2, atomCG, 126, 0, 1.22)));
        superimposition_residue->AddAtom(atomCG);
        superimposition_residue->AddAtom(atomOD1);
        superimposition_residue->AddAtom(atomND2);
        superimposition_atoms_ = superimposition_residue->GetAtoms();
    }
    else if ( (amino_acid_name.compare("THR")==0) || (amino_acid_name.compare("SER")==0) || (amino_acid_name.compare("OLT")==0) || (amino_acid_name.compare("OLS")==0) )
    {
        Atom *atomOG1 = new Atom(superimposition_residue, "OG", (GeometryTopology::get_cartesian_point_from_internal_coords(atomC5, atomO5, atomC1, 112, 68, 1.46)));
        Atom *atomCB = new Atom(superimposition_residue, "CB", (GeometryTopology::get_cartesian_point_from_internal_coords(atomO5, atomC1, atomOG1, 109.3, 75, 1.53)));
        Atom *atomCA = new Atom(superimposition_residue, "CA", (GeometryTopology::get_cartesian_point_from_internal_coords(atomC1, atomOG1, atomCB, 109.3, 125, 1.53)));
        superimposition_residue->AddAtom(atomCA);
        superimposition_residue->AddAtom(atomCB);
        superimposition_residue->AddAtom(atomOG1);
        superimposition_atoms_ = superimposition_residue->GetAtoms();
        if ( (amino_acid_name.compare("THR")==0) || (amino_acid_name.compare("OLT")==0) )
        {
            atomOG1->SetName("OG1"); // It's OG in Ser.
        }
    }
    else if ( (amino_acid_name.compare("TYR")==0) || (amino_acid_name.compare("OLY")==0) )
    {
        Atom *atomOH = new Atom(superimposition_residue, "OH", (GeometryTopology::get_cartesian_point_from_internal_coords(atomC5, atomO5, atomC1, 112, 68, 1.46)));
        Atom *atomCZ = new Atom(superimposition_residue, "CZ", (GeometryTopology::get_cartesian_point_from_internal_coords(atomO5, atomC1, atomOH, 117, 60, 1.35)));
        Atom *atomCE1 = new Atom(superimposition_residue, "CE1", (GeometryTopology::get_cartesian_point_from_internal_coords(atomC1, atomOH, atomCZ, 120, 180, 1.37)));
        superimposition_residue->AddAtom(atomCE1);
        superimposition_residue->AddAtom(atomCZ);
        superimposition_residue->AddAtom(atomOH);
        superimposition_atoms_ = superimposition_residue->GetAtoms();
    }
    else
    {
        // OK I DON'T KNOW HOW TO HANDLE EXCEPTIONS. This will never happen though I promise.
        std::cerr << "Problem in glycosylationsite::Prepare_Glycans_For_Superimposition_To_Particular_Residue(). Expect Segfault soon." << std::endl;
    }
    return;
}

void GlycosylationSite::Superimpose_Glycan_To_Glycosite(Residue *glycosite_residue)
{
    // Get the 3 target atoms from protein residue.
    AtomVector target_atoms;
   // superimposition_atoms_ points to three atoms that were added to the glycan. Based on their names e.g. CG, ND2, we will superimpose them onto
   // the correspoinding "target" atoms in the protein residue (glycosite_residue).
    for (auto &superimposition_atom : superimposition_atoms_)
    {
        for(auto &protein_atom : glycosite_residue->GetAtoms())
        {
            if (protein_atom->GetName() == superimposition_atom->GetName())
            {
                target_atoms.push_back(protein_atom);
            }
        }
    }
    AtomVector glycan_atoms = glycan_.GetAllAtomsOfAssembly();
    gmml::Superimpose(superimposition_atoms_, target_atoms, glycan_atoms);
    Residue* reducing_Residue = glycan_.GetResidues().at(1); // I assume I assumed something stupid here.
    AtomVector reducing_Atoms = reducing_Residue->GetAtoms();
    Atom* atomC1;
    for(AtomVector::iterator it = reducing_Atoms.begin(); it != reducing_Atoms.end(); it++)
    {
        Atom* atom = *it;
        if(atom->GetName().compare("C1")==0)
        {
            atomC1 = atom;
        }
    }
    //Connect the glycan and protein atoms to each other.
    Atom *protein_connection_atom = this->GetConnectingProteinAtom(glycosite_residue->GetName());
    protein_connection_atom->GetNode()->AddNodeNeighbor(atomC1);
    atomC1->GetNode()->AddNodeNeighbor(protein_connection_atom);
    //Delete the atoms used to superimpose the glycan onto the protein. Remove the residue.
    Residue *superimposition_residue = glycan_.GetAllResiduesOfAssembly().at(0);
    glycan_.RemoveResidue(superimposition_residue);
}

void GlycosylationSite::RenumberGlycanToMatch(Assembly &glycoprotein)
{
    std::string chainID = this->GetResidue()->GetChainID();
    gmml::log(__LINE__, __FILE__, gmml::INF, "Setting chainID to " + chainID + "\nNote a ? is fine and means no chain ID here.");
    int highestResidueNumber = selection::FindHighestResidueNumber(glycoprotein, chainID);
    for (auto &glycanResidue : this->GetAttachedGlycan()->GetResidues())
    {
        glycanResidue->SetChainID(chainID);
        glycanResidue->SetResidueNumber(std::to_string(++highestResidueNumber));
    }
    return;
}

void GlycosylationSite::Rename_Protein_Residue_To_GLYCAM_Nomenclature()
{
    std::string amino_acid_name = this->GetResidue()->GetName();
    if (amino_acid_name.compare("ASN")==0) {this->GetResidue()->SetName("NLN");}
    if (amino_acid_name.compare("SER")==0) {this->GetResidue()->SetName("OLS");}
    if (amino_acid_name.compare("THR")==0) {this->GetResidue()->SetName("OLT");}
    if (amino_acid_name.compare("TYR")==0) {this->GetResidue()->SetName("OLY");}
}

void GlycosylationSite::Rename_Protein_Residue_From_GLYCAM_To_Standard()
{
    std::string amino_acid_name = this->GetResidue()->GetName();
    if (amino_acid_name.compare("NLN")==0) {this->GetResidue()->SetName("ASN");}
    if (amino_acid_name.compare("OLS")==0) {this->GetResidue()->SetName("SER");}
    if (amino_acid_name.compare("OLT")==0) {this->GetResidue()->SetName("THR");}
    if (amino_acid_name.compare("OLY")==0) {this->GetResidue()->SetName("TYR");}
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
        overlap = gmml::CalculateAtomicOverlaps(this->GetProteinAtoms(), this->GetAttachedGlycan()->GetAllAtomsOfAssembly(), print);
    }
    if(moleculeType == GLYCAN)
    {
        for(auto &other_glycosite : other_glycosites_)
        {
        	double current_overlap = gmml::CalculateAtomicOverlaps(other_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly(), this->GetAttachedGlycan()->GetAllAtomsOfAssembly(), print);
            overlap += current_overlap;
            overlap += gmml::CalculateAtomicOverlaps(other_glycosite->GetAttachedGlycan()->GetAllAtomsOfAssembly(), this->GetAttachedGlycan()->GetAllAtomsOfAssembly(), print);
        }
    }
    return overlap;
}

void GlycosylationSite::PrintOverlaps()
{
    std::stringstream logss;
    logss << std::fixed; // Formating ouput
    logss << std::setprecision(2); // Formating ouput
    logss << std::setw(17) << this->GetResidue()->GetId() << " | "
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

double GlycosylationSite::CalculateBeadOverlaps(AtomVector &atomsA, AtomVector &atomsB)
{
    double radius = 3.0; //Using same radius for all beads.
    double distance = 0.0, overlap = 0.0, current_overlap = 0.0;
    for(AtomVector::iterator it1 = atomsA.begin(); it1 != atomsA.end(); ++it1)
    {
        Atom *atomA = *it1;
        for(AtomVector::iterator it2 = atomsB.begin(); it2 != atomsB.end(); ++it2)
        {
            Atom *atomB = *it2;
            if ( (atomA->GetCoordinates().at(0)->GetX() - atomB->GetCoordinates().at(0)->GetX()) < (radius * 2) ) // This is faster than calulating distance, and rules out tons of atom pairs.
            {
                distance = atomA->GetDistanceToAtom(atomB);
                if ( ( distance < (radius + radius) ) && ( distance > 0.0 ) ) //Close enough to overlap, but not the same atom
                {
                    current_overlap = gmml::CalculateAtomicOverlaps(atomA, atomB, radius, radius); // This calls the version with radius values
                    overlap += current_overlap;
                    //std::cout << atomA->GetResidue()->GetId() << " overlaping with " << atomB->GetResidue()->GetId() << ": " << current_overlap << "\n";
                }
            }
        }
    }
    return (overlap / gmml::CARBON_SURFACE_AREA); //Normalise to area of a buried carbon
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

void GlycosylationSite::SetProteinBeads(AtomVector *beads)
{
    protein_beads_ = *beads; // This creates a copy
    //Remove beads from attachment point residue. Don't want to count overlaps between glycan and residue it is attached to.
    for(AtomVector::iterator it1 = protein_beads_.begin(); it1 != protein_beads_.end(); /* Not incrementing here as erasing increments*/)
    {
        Atom *atom = *it1;
        if (atom->GetResidue() == this->GetResidue()) // this->GetResidue returns the glycosites protein residue.
        {
            protein_beads_.erase(std::remove(protein_beads_.begin(), protein_beads_.end(), *it1), protein_beads_.end());
        }
        else
        {
            ++it1; // "erase" increments it1, so if no erase happens this does an increment instead.
        }
    }
}

void GlycosylationSite::SetDefaultDihedralAnglesUsingMetadata()
{
    for(auto &linkage : all_residue_linkages_)
    {
        linkage.SetDefaultShapeUsingMetadata();
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

void GlycosylationSite::SetRandomDihedralAnglesUsingMetadataForNthLinkage(int linkage_number)
{
    if(linkage_number < all_residue_linkages_.size())
    {
        all_residue_linkages_.at(linkage_number).SetRandomShapeUsingMetadata();
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
//    AtomVector atoms = this->GetAttachedGlycan()->GetAllAtomsOfAssembly();
//    AtomVector sidechain_atoms = this->GetResidue()->GetAtoms();
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
//    AtomVector atoms = this->GetAttachedGlycan()->GetAllAtomsOfAssembly();
//    AtomVector sidechain_atoms = this->GetResidue()->GetAtoms();
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
//        if (this->GetResidue()->GetId().compare(glycosite.GetResidue()->GetId())==0) // if not this site. Oly you should overload the = operator?
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
        if(this->GetResidue()->GetId() != glycosite.GetResidue()->GetId())
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
        logss << "Residue ID: " << this->GetResidue()->GetId() << ", overlap: " << this->GetOverlap();
        logss << ", Gly--Pro " << this->GetGlycanOverlap() << "--" << this->GetProteinOverlap() << std::endl;
        gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    }
}
//////////////////////////////////////////////////////////
//                   PRIVATE FUNCTIONS                 //
//////////////////////////////////////////////////////////
// These functions do not fit here. Also, as ResidueNodes are not set, they are overly complex (nested recursive functions anyone?).
// They traverse the residues. I want to ignore the algycon for now.
void GlycosylationSite::FigureOutResidueLinkagesInGlycan(Residue *from_this_residue1, Residue *to_this_residue2, ResidueLinkageVector *residue_linkages)
{
    Atom *start_atom = to_this_residue2->GetAtoms().at(0); // Need to start somewhere.
    ResidueVector neighbors;
    RecursivelyGetAllNeighboringResidues(start_atom, &neighbors);
    ResidueVector glycan_residues = glycan_.GetResidues();
    for(auto &neighbor : neighbors)
    {
        if( (neighbor->GetIndex() != from_this_residue1->GetIndex()) && (std::find(glycan_residues.begin(), glycan_residues.end(), neighbor) != glycan_residues.end())  )
        {
            residue_linkages->emplace_back(neighbor, to_this_residue2);
        }
    }
    for(auto &neighbor : neighbors)
    {
        if( (neighbor->GetIndex() != from_this_residue1->GetIndex()) && (std::find(glycan_residues.begin(), glycan_residues.end(), neighbor) != glycan_residues.end())  )
        {
            this->FigureOutResidueLinkagesInGlycan(to_this_residue2, neighbor, residue_linkages);
        }
    }
    return;
}

void GlycosylationSite::RecursivelyGetAllNeighboringResidues(Atom* current_atom, ResidueVector* neighbors)
{
    current_atom->SetDescription("VisitedByRecursivelyGetAllNeighboringResidues");
    for(auto &neighboring_atom : current_atom->GetNode()->GetNodeNeighbors())
    {
        unsigned long long neighbor_index = neighboring_atom->GetResidue()->GetIndex();
        if(neighbor_index != current_atom->GetResidue()->GetIndex())
        {
            neighbors->push_back(neighboring_atom->GetResidue());
//            std::cout << "Foreign neighbor of " << current_atom->GetId() << " is " << neighboring_atom->GetId() << "\n";
//            std::cout << "size is now " << neighbors->size() << "\n";
        }
        else if (neighboring_atom->GetDescription().compare("VisitedByRecursivelyGetAllNeighboringResidues")!=0)
        { // If in current residue, and not visited already:
            this->RecursivelyGetAllNeighboringResidues(neighboring_atom, neighbors);
        }
    }
}

Atom* GlycosylationSite::GetConnectingProteinAtom(std::string residue_name)
{
    if(residue_name == "NLN" || residue_name == "ASN")
    {
        return residue_->GetAtom("ND2");
    }
    else if(residue_name == "OLT" || residue_name == "THR")
    {
        return residue_->GetAtom("OG1");
    }
    else if(residue_name == "OLS" || residue_name == "SER")
    {
        return residue_->GetAtom("OG");
    }
    else if(residue_name == "OLY" || residue_name == "TYR")
    {
        return residue_->GetAtom("OH");
    }
    else
    {
        std::cerr << "In GlycosylationSite::GetConectinProteinAtom you passed in a residue (" << residue_name << ") that isn't ON THE LIST. Exiting early: " << std::endl;
        exit(1);
    }
}

// OG re-reading. I'm pretty sure this belongs in Residue_linkage. linkage.wiggle(output_pdb_id, tolerance, interval
// OG re-re-reading. It's because you need to calculate overlaps.
void GlycosylationSite::WiggleOneLinkage(Residue_linkage &linkage, OverlapType overlapType, double tolerance, int interval)
{
    double current_overlap = this->CalculateOverlaps(overlapType);
    double lowest_overlap = current_overlap;
    // Reverse as convention is Glc1-4Gal and I want to wiggle in opposite direction i.e. from first rotatable bond in Asn outwards
    std::vector<Rotatable_dihedral> reversed_rotatable_bond_vector = linkage.GetRotatableDihedrals();
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
                		std::cout << "Wiggler: site " << this->GetResidueNumber() << " has overlap: " << lowest_overlap << "\n";
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
	int newCount = gmml::CountInternalBonds(*this->GetAttachedGlycan());
	if (newCount > this->GetInternalBondCount())
	{
		//std::cerr << "Internal glycan clash detected, rejecting shape.\n";
		return false;
	}
	return true;
}

void GlycosylationSite::DetermineOverlapProteinAtoms(Assembly* glycoprotein, Residue* attachmentResidue)
{ // Finds protein atoms not in the attachment residue
    AtomVector allProteinAtoms = glycoprotein->GetAllAtomsOfAssemblyWithinProteinResidues();
    AtomVector attachmentResidueAtoms = attachmentResidue->GetAtoms();
    this->SetProteinAtoms(selection::GetAtomsin_a_Notin_b_AtomVectors(allProteinAtoms, attachmentResidueAtoms));
}
