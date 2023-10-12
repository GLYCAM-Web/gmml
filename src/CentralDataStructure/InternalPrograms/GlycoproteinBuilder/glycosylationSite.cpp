#include <bits/std_abs.h>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>
#include <iomanip> // For setting precision and formating in output
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/CondensedSequence/carbohydrate.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/CentralDataStructure/Measurements/measurements.hpp" // calculateCoordinateFromInternalCoords
#include "includes/CentralDataStructure/Selections/atomSelections.hpp" //cdsSelections
#include "includes/CentralDataStructure/Editors/superimposition.hpp"
// #include "includes/CentralDataStructure/Overlaps/beadResidues.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CentralDataStructure/Selections/shaperSelections.hpp" // For the ClearAtomLabels sillyness.
#include "includes/MolecularMetadata/glycoprotein.hpp"
#include "includes/CentralDataStructure/Selections/templatedSelections.hpp"

// #include <algorithm> //  std::erase, std::remove
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
GlycosylationSite::GlycosylationSite(Residue* residue, Carbohydrate* carbohydrate,
                                     std::vector<Residue*> otherProteinResidues, unsigned int glycanStartResidueNumber)
    : residue_(residue), glycan_(carbohydrate), otherProteinResidues_(otherProteinResidues)
{
    gmml::log(__LINE__, __FILE__, gmml::INF, "Attaching glycan!");
    this->AttachGlycan(glycanStartResidueNumber);
    cdsSelections::ClearAtomLabels(carbohydrate->GetReducingResidue()); // jfc
    cdsSelections::ClearAtomLabels(this->GetResidue());
    // this->SetInternalBondCount(cdsSelections::CountInternalHeavyAtomBonds(this->GetGlycan()->getAtoms()));
    gmml::log(__LINE__, __FILE__, gmml::INF, "Done attach glycan to " + this->GetResidueId());
    gmml::log(__LINE__, __FILE__, gmml::INF, "Here are the linkages:");
    for (auto& linkage : carbohydrate->GetGlycosidicLinkages())
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, linkage.GetName());
        linkage.DetermineResiduesForOverlapCheck(); // Now that the protein is attached.
    }
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void GlycosylationSite::AttachGlycan(unsigned int glycanResidueStartNumber)
{
    gmml::log(__LINE__, __FILE__, gmml::INF,
              "Start of AttachGlycan. Residue ID is: " + this->GetResidue()->getStringId());
    this->Prepare_Glycans_For_Superimposition_To_Particular_Residue(this->GetResidue()->getName());
    gmml::log(__LINE__, __FILE__, gmml::INF, "Superimpose prep done");
    this->Superimpose_Glycan_To_Glycosite(this->GetResidue());
    gmml::log(__LINE__, __FILE__, gmml::INF, "SuperimposedGlycanToGlycosite");
    gmml::log(__LINE__, __FILE__, gmml::INF, "Setting internal bond count to check if more form later");
    this->RenumberGlycanToMatch(glycanResidueStartNumber);
    gmml::log(__LINE__, __FILE__, gmml::INF, "Attach glycan done");
}

void GlycosylationSite::RenumberGlycanToMatch(unsigned int startNumber)
{
    for (auto& glycanResidue : this->GetGlycan()->getResidues())
    {
        glycanResidue->setNumber(++startNumber);
    }
    return;
}

// This function prepares the glycan molecule in the glycan_ assembly for superimpostion onto an amino acid in the
// protein It does this by "growing" the atoms of the amino acid side chain (e.g. Asn, Thr or Ser) out from the glycan
// reducing terminal Another function will use these additional atoms to superimpose the glycan onto residue
void GlycosylationSite::Prepare_Glycans_For_Superimposition_To_Particular_Residue(std::string amino_acid_name)
{
    // Dear future self, the order that you add the atoms to the residue matters for superimposition ie N, CA, CB , not
    // CB, CA, N.
    Residue* reducing_Residue = this->GetGlycan()->GetReducingResidue();
    // Want: residue.FindAtomByTag("anomeric-carbon"); The below is risky as it uses atoms names, i.e. would break for
    // Sialic acid. ToDo Ok so I reckon the below is just assuming alpha or beta depending on the concext. Need to fix a
    // lot, but need to reproduce functionality after refactor first.
    //    Atom* anomericAtom = cdsSelections::guessAnomericAtom(reducing_Residue);
    //   This won't work as sometimes want alpha, sometimes beta. i.e. a CreateCoordinateForCenterAwayFromNeighbors
    //   function
    // This needs to be abstracted so it works for C2 reducing residues:
    Coordinate* coordC5       = reducing_Residue->FindAtom("C5")->getCoordinate();
    Coordinate* coordO5       = reducing_Residue->FindAtom("O5")->getCoordinate();
    Coordinate* coordC1       = reducing_Residue->FindAtom("C1")->getCoordinate();
    Atom* anomericAtom        = reducing_Residue->FindAtom("C1"); // For adding bond.
    // Delete aglycon atoms from glycan.
    Residue* aglycon          = this->GetGlycan()->GetAglycone();
    for (auto& atom : aglycon->getAtoms())
    {
        aglycon->deleteAtom(atom);
    }
    // Ok so going to set it so that the new "superimposition residue" is the old aglycon residue
    // This avoids having to delete the algycon residue object from assembly and adding the super residue to assembly.
    Residue* superimposition_residue = aglycon; // "renaming" so the below reads better.
    superimposition_residue->setName("SUP");
    // I put both the regular name and the O/N-linked glycam name here, as I'm not sure when it will be renamed.
    if ((amino_acid_name == "ASN") || (amino_acid_name == "NLN"))
    {
        Atom* atomND2 = superimposition_residue->addAtom(std::make_unique<Atom>(
            "ND2", (cds::calculateCoordinateFromInternalCoords(*coordC5, *coordO5, *coordC1, 109.3, 180, 1.53))));
        Atom* atomCG  = superimposition_residue->addAtom(
            std::make_unique<Atom>("CG", (cds::calculateCoordinateFromInternalCoords(
                                             *coordO5, *coordC1, *atomND2->getCoordinate(), 109.3, 261, 1.325))));
        superimposition_residue->addAtom(std::make_unique<Atom>(
            "OD1", (cds::calculateCoordinateFromInternalCoords(*coordC1, *atomND2->getCoordinate(),
                                                               *atomCG->getCoordinate(), 126, 0, 1.22))));
        anomericAtom->addBond(atomND2); // This is so findAnomericAtom works later, needs a foreign residue neighbor.
    }
    else if ((amino_acid_name.compare("THR") == 0) || (amino_acid_name.compare("SER") == 0) ||
             (amino_acid_name.compare("OLT") == 0) || (amino_acid_name.compare("OLS") == 0))
    {
        Atom* atomOG1 = superimposition_residue->addAtom(std::make_unique<Atom>(
            "OG", (cds::calculateCoordinateFromInternalCoords(*coordC5, *coordO5, *coordC1, 112, 68, 1.46))));
        Atom* atomCB  = superimposition_residue->addAtom(
            std::make_unique<Atom>("CB", (cds::calculateCoordinateFromInternalCoords(
                                             *coordO5, *coordC1, *atomOG1->getCoordinate(), 109.3, 75, 1.53))));
        superimposition_residue->addAtom(std::make_unique<Atom>(
            "CA", (cds::calculateCoordinateFromInternalCoords(*coordC1, *atomOG1->getCoordinate(),
                                                              *atomCB->getCoordinate(), 109.3, 125, 1.53))));
        if ((amino_acid_name.compare("THR") == 0) || (amino_acid_name.compare("OLT") == 0))
        {
            atomOG1->setName("OG1"); // It's OG in Ser.
        }
        anomericAtom->addBond(atomOG1); // This is so findAnomericAtom works later, needs a foreign residue neighbor.
    }
    else if ((amino_acid_name.compare("TYR") == 0) || (amino_acid_name.compare("OLY") == 0))
    {
        Atom* atomOH = superimposition_residue->addAtom(std::make_unique<Atom>(
            "OH", (cds::calculateCoordinateFromInternalCoords(*coordC5, *coordO5, *coordC1, 112, 68, 1.46))));
        Atom* atomCZ = superimposition_residue->addAtom(std::make_unique<Atom>(
            "CZ",
            (cds::calculateCoordinateFromInternalCoords(*coordO5, *coordC1, *atomOH->getCoordinate(), 117, 60, 1.35))));
        superimposition_residue->addAtom(std::make_unique<Atom>(
            "CE1", (cds::calculateCoordinateFromInternalCoords(*coordC1, *atomOH->getCoordinate(),
                                                               *atomCZ->getCoordinate(), 120, 180, 1.37))));
        anomericAtom->addBond(atomOH); // This is so findAnomericAtom works later, needs a foreign residue neighbor.
    }
    else
    {
        std::string message =
            "Problem creating glycosylation site. The amino acid requested: " + amino_acid_name +
            " has name that isn't supported. Currently you can glycosylate ASN, THR, SER or TYR. Email us to request " +
            "others. Ideally include examples of 3D structures we can use as a template.";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    return;
}

void GlycosylationSite::Superimpose_Glycan_To_Glycosite(Residue* glycosite_residue)
{
    // Get the 3 target atoms from protein residue.
    std::vector<Coordinate*> targetCoords;
    // superimposition_atoms_ points to three atoms that were added to the glycan. Based on their names e.g. CG, ND2, we
    // will superimpose them onto the correspoinding "target" atoms in the protein residue (glycosite_residue).
    for (auto& superimposition_atom : this->GetGlycan()->GetAglycone()->getAtoms())
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Aglycone atom is named: " + superimposition_atom->getName());
        for (auto& protein_atom : glycosite_residue->getAtoms())
        {
            if (protein_atom->getName() == superimposition_atom->getName())
            {
                gmml::log(__LINE__, __FILE__, gmml::INF,
                          "Adding " + protein_atom->getName() + " to superimposition atoms\n");
                targetCoords.push_back(protein_atom->getCoordinate());
            }
        }
    }
    std::vector<Coordinate*> aglyconeCoords =
        cdsSelections::getCoordinates(this->GetGlycan()->GetAglycone()->getAtoms());
    std::vector<Coordinate*> glycanCoords = cdsSelections::getCoordinates(this->GetGlycan()->getAtoms());
    //        std::cout << "Number of moving coords: " << glycanCoords.size() << "vs" <<
    //        this->GetGlycan()->getAtoms().size() << "\n"; std::cout << "Number of aglycone coords: " <<
    //        aglyconeCoords.size() << " vs " << this->GetGlycan()->GetAglycone()->getAtoms().size() << "\n"; std::cout
    //        << "Number of target coords: " << targetCoords.size() << "\n";
    // Sanity checks:
    if (aglyconeCoords.size() < 3)
    {
        throw std::runtime_error("The aglycone does not contain enough atoms to perform the requested "
                                 "superimposition to " +
                                 this->GetResidueId() + ".\nCheck your input structure!\n");
    }
    if (targetCoords.size() < 3)
    {
        throw std::runtime_error("Did not find the correctly named atoms in target residue to perform the requested "
                                 "superimposition to " +
                                 this->GetResidueId() + ".\nCheck your input structure!\n");
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Superimposing via the aglycone.");
    cds::Superimpose(aglyconeCoords, targetCoords, glycanCoords);
    // Connect the glycan and protein atoms to each other.
    Atom* protein_connection_atom = this->GetConnectingProteinAtom(glycosite_residue->getName());
    gmml::log(__LINE__, __FILE__, gmml::INF,
              "Connecting the reducing atom of " + this->GetGlycan()->GetReducingResidue()->getStringId() +
                  " to the protein:" + this->GetGlycan()->GetAnomericAtom()->getId() + " bonded to " +
                  protein_connection_atom->getId());
    protein_connection_atom->addBond(this->GetGlycan()->GetAnomericAtom()); // Atom connectivity
    this->Rename_Protein_Residue_To_GLYCAM_Nomenclature();                  // e.g. ASN to NLN
    gmml::log(__LINE__, __FILE__, gmml::INF, "Replacing the aglycone");
    // this->GetGlycan()->deleteResidue(this->GetGlycan()->GetAglycone());
    this->GetGlycan()->replaceAglycone(this->GetResidue());
    gmml::log(__LINE__, __FILE__, gmml::INF, "Completed superimposition to " + glycosite_residue->getStringId());
    return;
}

void GlycosylationSite::Rename_Protein_Residue_To_GLYCAM_Nomenclature()
{
    this->GetResidue()->setName(glycoproteinMetadata::ConvertGlycosylatedResidueName(this->GetResidue()->getName()));
}

void GlycosylationSite::Rename_Protein_Residue_From_GLYCAM_To_Standard()
{
    this->GetResidue()->setName(glycoproteinMetadata::ConvertGlycosylatedResidueName(this->GetResidue()->getName()));
}

unsigned int GlycosylationSite::CountOverlaps(MoleculeType moleculeType)
{
    unsigned int overlap = 0;
    if (moleculeType == ALL)
    {
        unsigned int proteinOverlap = this->CountOverlaps(MoleculeType::PROTEIN);
        unsigned int glycanOverlap  = this->CountOverlaps(MoleculeType::GLYCAN);
        return (proteinOverlap + glycanOverlap);
    }
    if (moleculeType == PROTEIN)
    {
        // This should not be necessary, either return a ref for the get, or accept a value in the function.
        std::vector<Residue*> proteinResidues = this->GetOtherProteinResidues();
        std::vector<Residue*> glycanResidues  = this->GetGlycan()->getResidues();
        unsigned int numberOfOverlaps         = this->CountOverlaps(proteinResidues, glycanResidues);
        return numberOfOverlaps;
    }
    if (moleculeType == GLYCAN)
    {
        std::vector<Residue*> glycanResidues = this->GetGlycan()->getResidues();
        for (auto& other_glycosite : other_glycosites_)
        {
            std::vector<Residue*> otherGlycanResidues = other_glycosite->GetGlycan()->getResidues();
            overlap                                   += this->CountOverlaps(otherGlycanResidues, glycanResidues);
        }
    }
    return overlap;
}

unsigned int GlycosylationSite::CountOverlaps(const std::vector<Residue*> residuesA,
                                              const std::vector<Residue*> residuesB)
{
    return cds::CountOverlappingResidues(residuesA, residuesB);
}

void GlycosylationSite::PrintOverlaps()
{
    std::stringstream logss;
    logss << std::fixed << std::setprecision(2) << std::setw(17) << this->GetResidue()->getStringId() << " | "
          << std::setw(6) << this->CountOverlaps() << std::endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
}

void GlycosylationSite::Wiggle(bool firstLinkageOnly, int interval)
{ // I want to find the lowest overlap as close to each bonds default as possible. So code is a bit more complicated.
    this->WiggleOneLinkage(this->GetProteinGlycanLinkage(), interval);
    if (!firstLinkageOnly)
    {
        for (auto& linkage : this->GetGlycan()->GetGlycosidicLinkages())
        {
            this->WiggleOneLinkage(linkage, interval);
        }
    }
    return;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
// void GlycosylationSite::SetDefaultDihedralAnglesUsingMetadata()
//{
//     for(auto &linkage : this->GetGlycan()->GetGlycosidicLinkages())
//     {
//         linkage.SetDefaultShapeUsingMetadata();
//     }
//     if( ! this->NoNewInternalCloseContacts() )
//     {
//         this->ResetDihedralAngles();
//     }
//     return;
// }

void GlycosylationSite::SetRandomDihedralAnglesUsingMetadata()
{
    if (this->GetGlycan()->GetGlycosidicLinkages().empty())
    {
        std::string message =
            "Asked to set a dihedral for glycosite with no glycosidic linkages; " + this->GetResidueId();
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    for (auto& linkage : this->GetGlycan()->GetGlycosidicLinkages())
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Setting random shape for linkage: " + linkage.GetName());
        linkage.SetRandomShapeUsingMetadata();
    }
    //    if (!this->NoNewInternalCloseContacts())
    //    {
    //        gmml::log(__LINE__, __FILE__, gmml::INF, "Nope that caused an internal contact");
    //        this->ResetDihedralAngles();
    //    }
    return;
}

void GlycosylationSite::ResetDihedralAngles()
{
    for (auto& linkage : this->GetGlycan()->GetGlycosidicLinkages())
    {
        linkage.SetShapeToPrevious();
    }
    return;
}

void GlycosylationSite::SetOtherGlycosites(std::vector<GlycosylationSite>& glycosites)
{
    other_glycosites_.clear();
    for (auto& glycosite : glycosites)
    {
        if (this != &glycosite)
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
    if (type.compare("All") == 0)
    {
        logss << "Residue ID: " << this->GetResidue()->getStringId() << ", overlap: " << this->CountOverlaps();
        gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    }
}

//////////////////////////////////////////////////////////
//                   PRIVATE FUNCTIONS                 //
//////////////////////////////////////////////////////////
Atom* GlycosylationSite::GetConnectingProteinAtom(const std::string residue_name) const
{
    std::string connectionAtomName = glycoproteinMetadata::GetGlycositeConnectionAtomName(residue_name);
    if (connectionAtomName == "")
    {
        std::string message =
            "Problem in GetConnectingProteinAtom. The amino acid requested: " + residue_name +
            " has name that isn't supported. Currently you can glycosylate ASN, THR, SER or TYR. Email us to request "
            "others. Ideally include examples of 3D structures we can use as a template.";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    return this->GetResidue()->FindAtom(connectionAtomName);
}

// OG re-reading. I'm pretty sure this belongs in Residue_linkage. linkage.wiggle(output_pdb_id, tolerance, interval
// OG re-re-reading. It's because you need to calculate overlaps.
void GlycosylationSite::WiggleOneLinkage(ResidueLinkage& linkage, int interval)
{
    int current_overlap                                           = this->CountOverlaps();
    // int lowest_overlap                                            = current_overlap;
    //  Reverse as convention is Glc1-4Gal and I want to wiggle in opposite direction i.e. from first rotatable bond in
    //  Asn outwards
    std::vector<RotatableDihedral> reversed_rotatable_bond_vector = linkage.GetRotatableDihedralsRef();
    std::reverse(reversed_rotatable_bond_vector.begin(), reversed_rotatable_bond_vector.end());
    for (auto& rotatable_dihedral : reversed_rotatable_bond_vector)
    {
        //        gmml::log(__LINE__, __FILE__, gmml::INF,
        //                  "About to wiggle all the rotamers of " + rotatable_dihedral.GetName() +
        //                      " in this linkage: " + linkage.GetName());
        rotatable_dihedral.WiggleUsingAllRotamers(linkage.GetMovingResidues(), linkage.GetFixedResidues(), interval);
    }
    return; // Note possibility of earlier return above
}
