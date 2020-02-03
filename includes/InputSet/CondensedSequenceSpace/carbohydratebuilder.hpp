#ifndef CARBOHYDRATEBUILDER_H
#define CARBOHYDRATEBUILDER_H
#include "condensedsequence.hpp"
#include "../../MolecularModeling/assembly.hpp"


/*
 * The current but probably very naive plan is that I use my awesome Residue Linkage class and not touch the condensedsequence.cc monster at all.
 * It's a freaking mess, but if I do this correctly all I need from that function is a 3D structure with the bonding set ok.
 * Good luck if there's a bug in there. The only problem I foresee (ha!) is that I won't know the sequence name of each linkage, e.g.:
 * DGalpb1-4DGlcpa- is linkage index 4. I'll need to pass that info back up to Dan in the JSON object, so my plan is to write a function
 * here that will figure that out from the 3D structure... Yeah I really don't want to touch the condensedsequence class. Have you read it?
 * */

namespace CondensedSequenceSpace
{
class carbohydrateBuilder
{
public:

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTORS                   //
    //////////////////////////////////////////////////////////

  //  carbohydrateBuilder();
    carbohydrateBuilder(std::string condensedSequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH", std::string prepFilePath = "../dat/prep/GLYCAM_06j-1.prep");

    //////////////////////////////////////////////////////////
    //                       ACCESSORS                      //
    //////////////////////////////////////////////////////////

    CondensedSequence GetCondensedSequence();
    std::string GetOfficialSequenceString();
    std::string GetInputSequenceString();
    MolecularModeling::Assembly* GetAssembly();
    ResidueLinkageVector* GetGlycosidicLinkages();
    bool GetSequenceIsValid();

    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////



    //////////////////////////////////////////////////////////
    //                        FUNCTIONS                     //
    //////////////////////////////////////////////////////////

    void GenerateSingle3DStructure(std::string outputFileNaming = "default");
    std::string GenerateUserOptionsJSON();
    void ReadUserSelectionsJSON(std::string jsonInput);
    void GenerateRotamers(std::string jsonSelection = "");

private:
    void Write3DStructureFile(std::string type = "PDB", std::string filename = "output.pdb");
    void SetInputSequenceString(std::string sequence);
    void SetOfficialSequenceString(std::string sequence);
    void SetSequenceIsValid(bool isValid);
    void SetDefaultShapeUsingMetadata();
    void ResolveOverlaps();
    void FigureOutResidueLinkagesInGlycan(MolecularModeling::Residue *from_this_residue1, MolecularModeling::Residue *to_this_residue2, ResidueLinkageVector *residue_linkages);
    void InitializeClass(std::string condensedSequence, std::string prepFilePath);

    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    MolecularModeling::Assembly assembly_;
    CondensedSequence condensedSequence_;
    std::string officialSequenceString_;
    std::string inputSequenceString_;
    //GlycanMetadataContainer metadataInformation_; // need this class
    ResidueLinkageVector glycosidicLinkages_;
    bool sequenceIsValid_;

};
}
#endif // CARBOHYDRATEBUILDER_H



// Proposed usage examples:

//carbohydrateBuilder builder(std::string sequenceString, std::string prepFilePath, enum jobTypeOption);

// Scrap

/* A big problem I have now is that condensed sequence class is separate. It makes a 3D structure of a carbohydrate in an assembly
 * and I don't know what the linkage "name" is for my linkages. i.e. Galb1-4Glc corresponds to linkage 1. Hmm. I think I need to
 * add a name attribute to residue_linkage and bring the functionality for building the 3D structure out of the old condensed sequence class.
 * Also need to bring in the set angles functionality. Then I can delete all the old stuff and cheer.
 * The other thing I need to sort out is rotamer generation (easy) and selection (?).
 * Then naming of PDB files and tracking/reporting which is which rotamer permutant.

*/
