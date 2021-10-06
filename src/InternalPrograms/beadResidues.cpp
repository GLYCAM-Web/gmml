#include "includes/InternalPrograms/beadResidues.hpp"
//#include "includes/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/MolecularModeling/atomnode.hpp"

AtomVector beads::Add_Beads_To_Protein(MolecularModeling::Assembly &assembly)
{
    // sfat stands for sidechain fat atom. mfat stands for mainchain fat atom. fat atoms are now called beads, but the name stayed.
    // Different names allow for easy selection in VMD visualisation problem. They serve no other purpose.
    // The beads should completely envelope all atoms in the protein. I first did all CA atoms in the mainchain, then visualized in VMD
    // and manually selected sidechain atoms that would cover the rest of the protein.
    // The code below adds to CA atoms and then selected other atoms based on their name and sometimes also the residue name.
    AtomVector protein_beads;
    for (auto &residue : assembly.GetAllProteinResiduesOfAssembly())
    {
        for (auto &atom : residue->GetAtoms())
        {
            if (atom->GetName().compare("CA")==0) // Main chain (mfat) CA atoms
            {
                //std::cout << "Adding bead to protein " << residue->GetId() << std::endl;
                Atom* bead_atom = new Atom(residue, "mfat", atom->GetCoordinates());
                residue->AddAtom(bead_atom);
                protein_beads.push_back(bead_atom);
            }
            else if ( (atom->GetName().compare("NZ")==0) || // Sidechain (sfat) atoms I've manually selected
                      (atom->GetName().compare("CZ")==0) ||
                      (atom->GetName().compare("NE2")==0) ||
                      (atom->GetName().compare("OD1")==0) ||
                      (atom->GetName().compare("SD")==0)  ||
                      ( (atom->GetName().compare("CE2")==0) && residue->GetName().compare("TRP")==0 ) ||
                      ( (atom->GetName().compare("CD1")==0) && ( residue->GetName().compare("LEU")==0 || residue->GetName().compare("ILE")==0 ) ) ||
                      ( (atom->GetName().compare("CD")==0) && residue->GetName().compare("GLU")==0 )
                    )
            { // sfats should move when a chi1, chi2 is moved, so make sure they are connected to something for the SetDihedral function to move them.
                Atom* bead_atom = new Atom(residue, "sfat", atom->GetCoordinates());
                residue->AddAtom(bead_atom);
                protein_beads.push_back(bead_atom);
            }
        }
    }
    return protein_beads;
}

AtomVector beads::Add_Beads_To_Glycan(ResidueVector glycan_residues)
{
    AtomVector glycan_beads;
    for (auto &residue : glycan_residues)
    {
        if ( residue->GetName().at(1) != 'S' ) // don't add a bead to a sialic acid in this section (see "else" below). Middle character of resname is always S for sialic acid.
        {
            Atom* bead_atom = new Atom(residue, "gfat", residue->GetGeometricCenter());
            residue->AddAtom(bead_atom);
            glycan_beads.push_back(bead_atom);
            //Bond bead_atom to any other atom in residue so when glycan is moved, bead_atom moves too.
            Atom *any_atom = residue->GetAtoms().at(0); // 0 is arbitrary, any atom would do.
            any_atom->GetNode()->AddNodeNeighbor(bead_atom);
            AtomVector temp = {any_atom};
            MolecularModeling::AtomNode *node = new MolecularModeling::AtomNode(); // DELETE IS FOR LOSERS.
            bead_atom->SetNode(node);
            bead_atom->GetNode()->SetNodeNeighbors(temp);
            for (auto &atom : residue->GetAtoms())
            {
                if ( (atom->GetName().compare("C2N") == 0) || (atom->GetName().compare("C6") == 0) )
                {
                    bead_atom = new Atom(residue, "gfat", atom->GetCoordinates().at(0));
                    residue->AddAtom(bead_atom);
                    glycan_beads.push_back(bead_atom);
                    any_atom->GetNode()->AddNodeNeighbor(bead_atom);
                    MolecularModeling::AtomNode *node1 = new MolecularModeling::AtomNode(); // DELETE IS FOR LOSERS.
                    bead_atom->SetNode(node1);
                    bead_atom->GetNode()->SetNodeNeighbors(temp);
                }
            }
        }
        else // if it is sialic acid
        {
            for (auto &atom : residue->GetAtoms())
            {
                if ( (atom->GetName().compare("C2") == 0) || (atom->GetName().compare("N5") == 0) || (atom->GetName().compare("C8") == 0) )
                {
                    Atom* bead_atom = new Atom(residue, "gfat", atom->GetCoordinates().at(0));
                    residue->AddAtom(bead_atom);
                    glycan_beads.push_back(bead_atom);
                    Atom *any_atom = residue->GetAtoms().at(0);
                    any_atom->GetNode()->AddNodeNeighbor(bead_atom);
                    AtomVector temp = {any_atom};
                    MolecularModeling::AtomNode *node = new MolecularModeling::AtomNode(); // DELETE IS FOR LOSERS.
                    bead_atom->SetNode(node);
                    bead_atom->GetNode()->SetNodeNeighbors(temp);
                }
            }
        }
    }
    return glycan_beads;
}

void beads::Remove_Beads(MolecularModeling::Assembly &ass)
{
    // Removes all bead atoms from the assembly.
    // Based on having "fat" in the name.
    // When finished, assembly should look the same as before Add_Beads was called.
    ResidueVector all_residues = ass.GetAllResiduesOfAssembly();
    for (auto &residue : all_residues)
    {
        for (auto &atom : residue->GetAtoms())
        {
            if (atom->GetName().find("fat")==1)
            {
                residue->RemoveAtom(atom);
//                std::cout << "RemovedBead: " << atom->GetId() << "\n";
            }
        }
    }
    return;
}





