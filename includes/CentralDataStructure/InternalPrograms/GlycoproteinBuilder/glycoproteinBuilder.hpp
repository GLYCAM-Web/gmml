#ifndef GMML_INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINBUILDER_HPP
#define GMML_INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINBUILDER_HPP
#include "includes/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/Abstract/absBuilder.hpp"
#include <string>
#include <dirent.h>
#include <sys/stat.h>
class GlycoproteinBuilder : public Abstract::absBuilder
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
	GlycoproteinBuilder(std::string inputFile);
    GlycoproteinBuilder(std::string inputFile = "Default", std::string workingDirectory = "Default");
    GlycoproteinBuilder(GlycoproteinBuilderInputs inputStruct);
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    inline MolecularModeling::Assembly* GetGlycoproteinAssemblyPtr()	{return &glycoprotein_;}
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    void ResolveOverlaps();
    void WriteOutputFiles();
private:
    //////////////////////////////////////////////////////////
    //                  PRIVATE ACCESSORS                   //
    //////////////////////////////////////////////////////////
    inline std::vector<GlycosylationSite>& GetGlycosites() 			{return glycosites_;}
    inline MolecularModeling::Assembly& GetGlycoproteinAssembly()	{return glycoprotein_;}
    inline std::string GetWorkingDirectory() 						{return workingDirectory_;}
    inline std::string GetPrepFileLocation() 						{return prepFileLocation_;}
    inline std::string GetProteinPDBFileName() 						{return proteinPDBFileName_;}
    inline std::vector<GlycosylationSite> GetGlycosylationSites() 	{return glycosites_;}
    inline double GetOverlapTolerance() 							{return overlapTolerance_;}
    inline int GetNumberOfOutputStructures() 						{return numberOfOutputStructures_;}
    inline int GetPersistCycles()									{return persistCycles_;}
    inline bool GetIsDeterministic()								{return isDeterministic_;}
    double GetGlobalOverlap();
    //////////////////////////////////////////////////////////
    //                  PRIVATE MUTATORS                    //
    //////////////////////////////////////////////////////////
    void SetWorkingDirectory(std::string workingDirectory);
    void SetPrepFileLocation(std::string prepFileLocation);
    inline void SetPersistCycles(int i) 							{persistCycles_ = i;}
    inline void SetOverlapTolerance(double d)						{overlapTolerance_ = d;}
    inline void SetProteinPDBFileName(std::string s)				{proteinPDBFileName_ = s;}
    inline void SetNumberOfOutputStructures(int i) 					{numberOfOutputStructures_ = i;}
    inline void SetIsDeterministic(std::string isDeterministic);
    //////////////////////////////////////////////////////////
    //                  PRIVATE FUNCTIONS                   //
    //////////////////////////////////////////////////////////
    // Class instantiation
    void InitializeGlycoproteinBuilder(GlycoproteinBuilderInputs inputStruct);
    void ConvertInputStructEntries(GlycoproteinBuilderInputs inputStruct);
    void CreateGlycosites(std::vector<GlycositeInput> glycositesInputVector);
    // Overlap Resolution
    void Wiggle(OverlapType overlapType = ATOMIC, int persistCycles = 100, bool firstLinkageOnly = false, int interval = 5);
    void RandomDescent(OverlapType overlapType, int persistCycles, bool monte_carlo);
    void SetDefaultDihedralAnglesUsingMetadata();
    void SetRandomDihedralAnglesUsingMetadata();
    bool DumbRandomWalk(int maxCycles = 10);
    // I/O
    void PrintDihedralAnglesAndOverlapOfGlycosites();
    void CalculateAndPrintOverlaps();
    // Overlap Calculation
    double CalculateOverlaps(OverlapType overlapType = BEAD, MoleculeType moleculeType = ALL, bool recordOverlap = true, bool printOverlap = false);
    std::vector<GlycosylationSite*> DetermineSitesWithOverlap(double tolerance, OverlapType overlapType = BEAD);
    std::vector<GlycosylationSite*> GetSitesWithOverlap(double tolerance);
    void DeleteSitesIterativelyWithAtomicOverlapAboveTolerance(std::vector<GlycosylationSite> &glycosites, double tolerance);
    void UpdateAtomsThatMoveInLinkages();
    void SetOtherGlycosites();
    void SetChargesAndAtomTypes(); // This will be silly for now.
    // Selection
//    ResidueLinkageVector GetAllFirstAnd1_6Linkages();
//    ResidueLinkageVector GetAllFirstAnd2_XLinkages();
//    ResidueLinkageVector SplitLinkagesIntoPermutants(ResidueLinkageVector inputLinkages);
    // Stashing 3D structures in Assembly coordinate vector
//    void StashCoordinates();
//    void SetStashedCoordinatesWithLowestOverlap();
    // Beads.
    void Add_Beads(MolecularModeling::Assembly &glycoprotein, std::vector<GlycosylationSite> &glycosites);
    void Set_Other_Glycan_Beads(std::vector<GlycosylationSite> &glycosites);
 //   void Remove_Beads(MolecularModeling::Assembly &glycoprotein);
 //   AtomVector Add_Beads_To_Glycan(ResidueVector glycan_residues);
 //   AtomVector Add_Beads_To_Protein(MolecularModeling::Assembly &assembly);
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    std::vector<GlycosylationSite> glycosites_; 	// Info about each glycosylation site. See the class.
    MolecularModeling::Assembly glycoprotein_; 		// Generated by this code.
    double overlapTolerance_; 						// What amount of overlap to tolerate.
    std::string proteinPDBFileName_; 				// The protein pdb file to attach the glycans to.
    std::string prepFileLocation_; 					// Prep file for carbohydrate builder. Default is to figure out install directory + ../dat/prep/GLYCAM_06j-1_GAGS.prep.
    std::string workingDirectory_; 					// Where the outputs will go, and where to read the input file from. Default is current folder.
    int numberOfOutputStructures_;  				// Number of output 3D structures. Default 1.
    int persistCycles_;								// Algos continue for persistCycles, reset count if there is an overlap improvement.
    bool isDeterministic_;							// "true" means produce the same output for a given input each time. Good for testing.
};
#endif // GMML_INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINBUILDER_HPP
