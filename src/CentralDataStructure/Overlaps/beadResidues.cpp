#include "includes/CentralDataStructure/Overlaps/beadResidues.hpp"
#include "includes/CodeUtils/templatedSelections.hpp"
#include "includes/CodeUtils/biology.hpp"

std::vector<Atom*> beads::Add_Beads_To_Protein(Assembly &assembly)
{
	// sfat stands for sidechain fat atom. mfat stands for mainchain fat atom. fat atoms are now called beads, but the name stayed.
	// Different names allow for easy selection in VMD visualisation problem. They serve no other purpose.
	// The beads should completely envelope all atoms in the protein. I first did all CA atoms in the mainchain, then visualized in VMD
	// and manually selected sidechain atoms that would cover the rest of the protein.
	// The code below adds to CA atoms and then selected other atoms based on their name and sometimes also the residue name.
	std::vector<Atom*> protein_beads;
	std::vector<Residue*> proteinResidues = codeUtils::getElementsWithNames(assembly.getResidues(), biology::proteinResidueNames);
	for (auto &residue : proteinResidues)
	{
		for (auto &atom : residue->getAtoms())
		{
			if (atom->getName().compare("CA")==0) // Main chain (mfat) CA atoms
			{
				//std::cout << "Adding bead to protein " << residue->GetId() << std::endl;
				residue->addAtom(std::make_unique<Atom>("mfat", *(atom->getCoordinate())));
				protein_beads.push_back(residue->getAtoms().back());
			}
			else if ( (atom->getName().compare("NZ")==0) || // Sidechain (sfat) atoms I've manually selected
					(atom->getName().compare("CZ")==0) ||
					(atom->getName().compare("NE2")==0) ||
					(atom->getName().compare("OD1")==0) ||
					(atom->getName().compare("SD")==0)  ||
					( (atom->getName().compare("CE2")==0) && residue->getName().compare("TRP")==0 ) ||
					( (atom->getName().compare("CD1")==0) && ( residue->getName().compare("LEU")==0 || residue->getName().compare("ILE")==0 ) ) ||
					( (atom->getName().compare("CD")==0) && residue->getName().compare("GLU")==0 )
			)
			{ // sfats should move when a chi1, chi2 is moved, so make sure they are connected to something for the SetDihedral function to move them.
				residue->addAtom(std::make_unique<Atom>("sfat", *(atom->getCoordinate())));
                protein_beads.push_back(residue->getAtoms().back());
			}
		}
	}
	return protein_beads;
}

std::vector<Atom*> beads::Add_Beads_To_Glycan(std::vector<Residue*> glycan_residues)
{
	std::vector<Atom*> glycan_beads;
	for (auto &residue : glycan_residues)
	{
		if ( residue->getName().at(1) != 'S' ) // don't add a bead to a sialic acid in this section (see "else" below). Middle character of resname is always S for sialic acid.
		{
			residue->addAtom(std::make_unique<Atom>("gfat", *(residue->getGeometricCenter())));
			glycan_beads.push_back(residue->getAtoms().back());
			//Bond bead_atom to any other atom in residue so when glycan is moved, bead_atom moves too.
			Atom *any_atom = residue->getAtoms().front();
			any_atom->addBond(residue->getAtoms().back());
			for (auto &atom : residue->getAtoms())
			{
				if ( (atom->getName().compare("C2N") == 0) || (atom->getName().compare("C6") == 0) )
				{
					residue->addAtom(std::make_unique<Atom>("gfat", *(atom->getCoordinate())));
					glycan_beads.push_back(residue->getAtoms().back());
					any_atom->addBond(residue->getAtoms().back());
				}
			}
		}
		else // if it is sialic acid
		{
			for (auto &atom : residue->getAtoms())
			{
				if ( (atom->getName().compare("C2") == 0) || (atom->getName().compare("N5") == 0) || (atom->getName().compare("C8") == 0) )
				{
				    residue->addAtom(std::make_unique<Atom>("gfat", *(atom->getCoordinate())));
				    glycan_beads.push_back(residue->getAtoms().back());
				    Atom *any_atom = residue->getAtoms().front();
				    any_atom->addBond(residue->getAtoms().back());
				}
			}
		}
	}
	return glycan_beads;
}

void beads::Remove_Beads(Assembly &ass)
{
	// Removes all bead atoms from the assembly.
	// Based on having "fat" in the name.
	// When finished, assembly should look the same as before Add_Beads was called.
	std::vector<Residue*> all_residues = ass.getResidues();
	for (auto &residue : all_residues)
	{
		for (auto &atom : residue->getAtoms())
		{
			if (atom->getName().find("fat") == 1)
			{
				residue->deleteAtom(atom);
			}
		}
	}
	return;
}
