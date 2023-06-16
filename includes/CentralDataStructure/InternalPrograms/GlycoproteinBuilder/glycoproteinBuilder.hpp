#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINBUILDER_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINBUILDER_HPP
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/glycosylationSite.hpp"
#include "includes/CentralDataStructure/InternalPrograms/GlycoproteinBuilder/gpInputStructs.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbPreprocessorInputs.hpp"
#include "includes/Abstract/absBuilder.hpp"
#include <string>
//#include <dirent.h>
//#include <sys/stat.h>

using cds::Assembly;
class GlycoproteinBuilder : public Abstract::absBuilder
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    GlycoproteinBuilder(glycoprotein::GlycoproteinBuilderInputs inputStruct, pdb::PreprocessorOptions preprocessingOptions = pdb::PreprocessorOptions());
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    inline Assembly* getGlycoprotein()	{return &glycoprotein_;}
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    void ResolveOverlaps();
    void WritePdbFile(const std::string prefix = "glycoprotein", const bool writeConectSection = true);
    void WriteOffFile(const std::string prefix = "glycoprotein");
    //void WriteOutputFiles(std::string prefix = "Glycoprotein_All_Resolved");
    void PrintDihedralAnglesAndOverlapOfGlycosites();
private:
    //////////////////////////////////////////////////////////
    //                  PRIVATE ACCESSORS                   //
    //////////////////////////////////////////////////////////
    inline std::vector<GlycosylationSite>& GetGlycosites() 			{return glycosites_;}
    inline std::vector<GlycosylationSite>& GetGlycosylationSites() 	{return glycosites_;}
    inline int GetOverlapTolerance() const			                {return overlapTolerance_;}
    inline int GetNumberOfOutputStructures() const 					{return numberOfOutputStructures_;}
    inline int GetPersistCycles() const						        {return persistCycles_;}
    inline bool GetIsDeterministic() const							{return isDeterministic_;}
    //////////////////////////////////////////////////////////
    //                  PRIVATE MUTATORS                    //
    //////////////////////////////////////////////////////////
    void SetWorkingDirectory(const std::string workingDirectory);
    void SetPrepFileLocation(const std::string prepFileLocation);
    inline void SetPersistCycles(const int i) 							{persistCycles_ = i;}
    inline void SetOverlapTolerance(const int i)						{overlapTolerance_ = i;}
    inline void SetProteinPDBFileName(const std::string s)				{proteinPDBFileName_ = s;}
    inline void SetNumberOfOutputStructures(const int i) 				{numberOfOutputStructures_ = i;}
    inline void SetIsDeterministic(const bool b)                        {isDeterministic_ = b;}
    //////////////////////////////////////////////////////////
    //                  PRIVATE FUNCTIONS                   //
    //////////////////////////////////////////////////////////
    // Class instantiation
    void CreateGlycosites(std::vector<glycoprotein::GlycositeInput> glycositesInputVector, const std::string prepFileLocation);
    Residue* SelectResidueFromInput(const std::string userSelection);
    // Overlap Resolution
    void Wiggle(Resolution resolutionLevel = RESIDUE, int persistCycles = 100, bool firstLinkageOnly = false, int interval = 5);
    void RandomDescent(Resolution resolutionLevel, int persistCycles, bool monte_carlo);
    void SetRandomDihedralAnglesUsingMetadata();
    bool DumbRandomWalk(int maxCycles = 10);
    // I/O
    // Overlap Calculation
    int CalculateOverlaps(Resolution resolutionLevel = RESIDUE, MoleculeType moleculeType = ALL);
    std::vector<GlycosylationSite*> DetermineSitesWithOverlap(Resolution resolutionLevel = RESIDUE, MoleculeType moleculeType = ALL);
    void DeleteSitesIterativelyWithAtomicOverlapAboveTolerance(std::vector<GlycosylationSite> &glycosites, int tolerance);
    void UpdateAtomsThatMoveInLinkages();
    void SetOtherGlycosites();
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    std::vector<GlycosylationSite> glycosites_; 	// Info about each glycosylation site. See the class.
    Assembly glycoprotein_; 		                // Generated by this code.
    int overlapTolerance_ = 1; 						    // What amount of overlap to tolerate.
    std::string proteinPDBFileName_; 				// The protein pdb file to attach the glycans to.
    int numberOfOutputStructures_;  				// Number of output 3D structures. Default 1.
    int persistCycles_;								// Algos continue for persistCycles, reset count if there is an overlap improvement.
    bool isDeterministic_;							// "true" means produce the same output for a given input each time. Good for testing.
};
#endif // GMML_INCLUDES_INTERNALPROGRAMS_GLYCOPROTEINBUILDER_GLYCOPROTEINBUILDER_HPP
