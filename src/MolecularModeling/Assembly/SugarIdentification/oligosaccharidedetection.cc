#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <stack>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <sys/stat.h>

#include "../../../../includes/MolecularModeling/assembly.hpp"
#include "../../../../includes/MolecularModeling/residue.hpp"
#include "../../../../includes/MolecularModeling/atom.hpp"
#include "../../../../includes/MolecularModeling/atomnode.hpp"
#include "../../../../includes/Glycan/monosaccharide.hpp"
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyfile.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologybond.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologybondtype.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyangletype.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyatompair.hpp"
#include "../../../../includes/InputSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbtitlesection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbconnectsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinksection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../../includes/utils.hpp"
#include "../../../../includes/common.hpp"
#include "../../../../includes/GeometryTopology/grid.hpp"
#include "../../../../includes/GeometryTopology/cell.hpp"
#include "../../../../includes/MolecularModeling/ring_shape_detection.hpp"
#include "../../../../includes/Glycan/monosaccharide.hpp"
#include "../../../../includes/Glycan/oligosaccharide.hpp"
#include "../../../../includes/Glycan/glycosidiclinkage.hpp"

#include <unistd.h>
#include <errno.h>
#include <string.h>

using MolecularModeling::Assembly;
using Glycan::Oligosaccharide;
using Glycan::Monosaccharide;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

// As of June 2018, DetectShape is called once in gmml, during the ExtractSugars function. Yao commented it out a month ago
// as it was causing segfaults and also removing all sidechain atoms from the assembly. Davis thinks this will be a problem
// when we go to build the ontology again, as it needs Detect Shape. I (Oliver) intend to comment out the below, and set up
// gmml to use the BFMP code I ported into gmml. See ring_shape_detection.hpp in MolecularModeling folder.
//void Assembly::DetectShape(MolecularModeling::AtomVector cycle, Glycan::Monosaccharide* mono)
//void Assembly::GetRingShapeBFMP_deprecatedExternalProgramCall(MolecularModeling::AtomVector cycle, Glycan::Monosaccharide* mono)
//{
//    ///Creating a new assembly only from the ring atoms for external detect shape program
//    Assembly* detect_shape_assembly = new Assembly();
//    detect_shape_assembly->AddResidue(cycle.at(0)->GetResidue());
//
//    Residue* detect_shape_residue = detect_shape_assembly->GetResidues().at(0);
//    for(unsigned int i = 0; i < cycle.size(); i++)
//    {
//        std::string name = cycle.at(i)->GetName();
//        std::string id = cycle.at(i)->GetId();
//
//        ///Preparing atom name and id for detect shape program. It doesn't work with atoms containing special characters: C', C* etc
//        //        replace( id.begin(), id.end(), '?', 'n'); // replace all '?' with 'n'
//        gmml::FindReplaceString(id, "\'", "");
//        gmml::FindReplaceString(id, ",", "");
//        gmml::FindReplaceString(name, "*", "");
//        replace( id.begin(), id.end(), '*', 's'); // replace all '*' with ''
//
//        //        replace( name.begin(), name.end(), '?', 'n'); // replace all '?' with 'n'
//        gmml::FindReplaceString(name, "\'", "");
//        gmml::FindReplaceString(name, ",", "");
//        gmml::FindReplaceString(name, "*", "");
//
//
//        cycle.at(i)->SetName(name);
//        cycle.at(i)->SetId(id);
//    }
//    detect_shape_residue->SetAtoms(cycle); //Commented out, suspected bug.
//    //detect_shape_residue->SetAtoms(detect_shape_assembly-> GetAllAtomsOfAssembly());
//
//		// First need to get environment variable GEMSHOME
//		std::string gemshome( getenv( "GEMSHOME" ) );
//		// Concatenate GEMSHOME to name of all files used for BFMP detect shape program.
//		std::string detect_shape = gemshome + "/apps/BFMP/detect_shape";
//		std::string temp_gmml_PDB = gemshome + "/temp_gmml_pdb.pdb";
//		std::string temp_detect_shape_PDB = gemshome + "/temp_detect_shape_pdb.pdb";
//		std::string temp_config = gemshome + "/temp_config";
//		std::string canonicals = gemshome + "/apps/BFMP/canonicals.txt";
//		// This file doesn't get concatenated to GEMSHOME because BFMP program generates it in the PWD
//		std::string ring_conformations = "ring_conformations.txt";
//
//    ///Write a new PDB file from the new assembly
//    PdbFileSpace::PdbFile* pdb = detect_shape_assembly->BuildPdbFileStructureFromAssembly();
//    pdb->Write( temp_gmml_PDB );
//
//    ///Converting the written PDB file to format readable by detect_shape program
//    std::string line = "";
//    std::ifstream gmml_pdb ( temp_gmml_PDB.c_str() );
//    std::ofstream detect_shape_pdb ( temp_detect_shape_PDB.c_str() );
//    int n = 0;
//    if (gmml_pdb.is_open())
//    {
//        while (!gmml_pdb.eof()) {
//            getline(gmml_pdb, line);
//            if(line.find("HETATM") != std::string::npos)
//            {
//                detect_shape_pdb << line << "\n";
//            }
//            n++;
//        }
//        gmml_pdb.close();
//        detect_shape_pdb.close();
//    }
//    else std::cout << "Unable to open " << temp_detect_shape_PDB << " file" << "\n";
//
//    ///Writing a configuration file for the second argument of the detect_sugar program
//    std::ofstream detect_shape_configuration ( temp_config.c_str() );
//    detect_shape_configuration << "Atom" << "\n";
//    for(unsigned int i = 0; i < cycle.size(); i++)
//    {
//        detect_shape_configuration << cycle.at(i)->GetName() << "\n";
//    }
//    detect_shape_configuration << "Residue" << "\n";
//    detect_shape_configuration << "1" << "\n";
//    detect_shape_configuration << "Path" << "\n";
//    detect_shape_configuration << canonicals << "\n";
//    detect_shape_configuration.close();
//
//    ///Calling detect_shape program
//    // The better way to do this would be to pass the bfmp code the six cartesian
//    // coords it needs and call it as a part of the library.  not sure how much
//    // work it would take to do this, but it would be better.
//		std::string command = detect_shape + " " + temp_detect_shape_PDB + " " + temp_config + " > /dev/null";
//		system( command.c_str() );
//
//    ///Adding the BFMP ring conformation information gained from the detect_sugar program to the monosaccharide
//    std::ifstream shape_detection_result ( ring_conformations.c_str() );
//    line = "";
//    if (shape_detection_result.is_open())
//    {
//        getline (shape_detection_result,line);
//        getline (shape_detection_result,line);
//        std::vector<std::string> line_tokens = gmml::Split(line, "\t");
//        if(line_tokens.at(1).compare("-") == 0)
//            mono->bfmp_ring_conformation_ = line_tokens.at(2);
//        else
//            mono->bfmp_ring_conformation_ = line_tokens.at(1);
//        shape_detection_result.close();
//    }
//    else std::cout << "Unable to open " << ring_conformations << " file from detect shape program" << "\n";
//
//    ///Deleting temporary files
//    remove( temp_detect_shape_PDB.c_str() );
//    remove( temp_gmml_PDB.c_str() );
//    remove( temp_config.c_str() );
//    remove( ring_conformations.c_str() );
//}

//This is a wrapper of the original ExtractSugars function. I want to make the vector of monosaccharides external, so I can do some manipulations on them.
//Other users call this wrapper.
std::vector< Glycan::Oligosaccharide* > Assembly::ExtractSugars( std::vector< std::string > amino_lib_files, bool glyprobity_report, bool populate_ontology, bool individualOntologies, std::string CCD_Path )
{

    std::vector<Glycan::Monosaccharide*> monos = std::vector<Glycan::Monosaccharide*>();

    OligosaccharideVector oligos = Assembly::ExtractSugars( amino_lib_files, monos, glyprobity_report,  populate_ontology, individualOntologies, CCD_Path );


    return oligos; // Oliver thinks Yao hosed ontology::analysis. This might fix it?
}

//A function for PDB integration of gmml where they can feed a stringstream of atom cards and get a list of oligosaccharide names
std::vector<MolecularModeling::Assembly::gmml_api_output> Assembly::PDBExtractSugars(std::vector< std::string > amino_lib_files, std::string CCD_Path)
{
  int local_debug = -1;
  std::vector<Glycan::Monosaccharide*> monos = std::vector<Glycan::Monosaccharide*>();
  OligosaccharideVector oligos = Assembly::ExtractSugars( amino_lib_files, monos, false, false, false, CCD_Path );
  std::vector<gmml_api_output> api_output;
  int i = 1;
  oligos.shrink_to_fit();
  for( OligosaccharideVector::iterator it = oligos.begin(); it != oligos.end(); it++ )
  {
    Glycan::Oligosaccharide* thisOligo = *it;
    MolecularModeling::Assembly::gmml_api_output thisOutput;
    thisOutput.linear_descriptors.push_back(std::make_pair("oligoName",  thisOligo->oligosaccharide_name_));
    thisOutput.linear_descriptors.push_back(std::make_pair("oligoSequenceName",  thisOligo->oligosaccharide_name_));
    thisOutput.linear_descriptors.push_back(std::make_pair("oligoIUPACname",  thisOligo->IUPAC_name_));
    thisOutput.linear_descriptors.push_back(std::make_pair("authorOligo",  thisOligo->author_IUPAC_name_));
    thisOligo->mono_nodes_.shrink_to_fit();
    for(std::vector<Glycan::Monosaccharide*>::reverse_iterator rit = thisOligo->mono_nodes_.rbegin(); rit != thisOligo->mono_nodes_.rend();rit++)
    {
      Glycan::Monosaccharide* thisMono = *rit;
      if(local_debug > 0)
      {
        std::string monoID =  thisMono->cycle_atoms_[0]->GetResidue()->GetId();
        std::cout << monoID << "\n";
      }
      //Mono Index
      thisOutput.indices.push_back(std::make_pair(std::to_string(thisMono->oligosaccharide_index_), thisMono->cycle_atoms_[0]->GetResidue()->GetId()));
      //Mono connectivity
      thisMono->mono_neighbors_.shrink_to_fit();
      for (std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> >::iterator it2 = thisMono->mono_neighbors_.begin(); it2 != thisMono->mono_neighbors_.end(); it2++)
      {
          Glycan::Monosaccharide* thisNeighbor = (*it2).second;
          std::string neighborMonoID = thisNeighbor->cycle_atoms_[0]->GetResidue()->GetId();
          Glycan::GlycosidicLinkage* thisLinkage = (*it2).first;
          int NeighborNum = thisNeighbor->oligosaccharide_index_;
          //residue ID, atom name, residue ID2, atom name 2.
          std::vector<std::string> residue_links_vector;
          residue_links_vector.push_back(thisMono->cycle_atoms_[0]->GetResidue()->GetId());
          std::string thisID;
          std::size_t endAtomName;
          if((thisLinkage->reducing_mono_ == thisMono) && (thisLinkage->glycosidic_oxygen_ != NULL))
          {
            thisID = thisLinkage->glycosidic_oxygen_->GetId();
            endAtomName = thisID.find("_");
            thisID = thisID.substr(0, endAtomName);
            residue_links_vector.push_back(thisID);
          }
          else if((thisLinkage->non_reducing_mono_ == thisMono) && thisLinkage->non_reducing_mono_carbon_ != NULL)
          {
            thisID = thisLinkage->non_reducing_mono_carbon_->GetId();
            endAtomName = thisID.find("_");
            thisID = thisID.substr(0, endAtomName);
            residue_links_vector.push_back(thisID);
          }
          else if((thisLinkage->non_reducing_mono_2_ == thisMono)  && thisLinkage->non_reducing_mono_2_carbon_ != NULL)
          {
            thisID = thisLinkage->non_reducing_mono_2_carbon_->GetId();
            endAtomName = thisID.find("_");
            thisID = thisID.substr(0, endAtomName);
            residue_links_vector.push_back(thisID);
          }
          residue_links_vector.push_back(thisNeighbor->cycle_atoms_[0]->GetResidue()->GetId());
          if((thisLinkage->reducing_mono_ == thisNeighbor) && (thisLinkage->glycosidic_oxygen_ != NULL))
          {
            thisID = thisLinkage->glycosidic_oxygen_->GetId();
            endAtomName = thisID.find("_");
            thisID = thisID.substr(0, endAtomName);
            residue_links_vector.push_back(thisID);
          }
          else if((thisLinkage->non_reducing_mono_ == thisNeighbor) && thisLinkage->non_reducing_mono_carbon_ != NULL)
          {
            thisID = thisLinkage->non_reducing_mono_carbon_->GetId();
            endAtomName = thisID.find("_");
            thisID = thisID.substr(0, endAtomName);
            residue_links_vector.push_back(thisID);
          }
          else if((thisLinkage->non_reducing_mono_2_ == thisNeighbor)  && thisLinkage->non_reducing_mono_2_carbon_ != NULL)
          {
            thisID = thisLinkage->non_reducing_mono_2_carbon_->GetId();
            endAtomName = thisID.find("_");
            thisID = thisID.substr(0, endAtomName);
            residue_links_vector.push_back(thisID);
          }
            thisOutput.residue_links.push_back(residue_links_vector);
      }
      //Errors at the mono level
      for(std::vector<Glycan::Note*>::iterator it3 = thisMono->mono_notes_.begin(); it3 != thisMono->mono_notes_.end(); it3++)
      {
        Glycan::Note* thisNote = (*it3);
        std::string thisNoteString = thisNote->type_ + ": " + thisNote->description_;
        thisOutput.error_warning_messages.push_back(thisNoteString);
      }
    }
    //Errors at the Oligo level
    for(std::vector<Glycan::Note*>::iterator it = thisOligo->oligo_notes_.begin(); it != thisOligo->oligo_notes_.end(); it++)
    {
      Glycan::Note* thisNote = (*it);
      std::string thisNoteString = thisNote->ConvertGlycanNoteType2String(thisNote->type_) + ": " + thisNote->description_;
      thisOutput.error_warning_messages.push_back(thisNoteString);
    }
    api_output.push_back(thisOutput);
    i++;
  }
  return api_output;
}

std::vector< Glycan::Oligosaccharide* > Assembly::ExtractSugars( std::vector< std::string > amino_lib_files, std::vector <Glycan::Monosaccharide*>& monos, bool glyprobity_report, bool populate_ontology, bool individualOntologies, std::string CCD_Path)
{
  int local_debug = -1;
  gmml::ResidueNameMap dataset_residue_names = GetAllResidueNamesFromMultipleLibFilesMap( amino_lib_files );
  if(local_debug > 0)
  {
    std::cout << "\n" << "Extracting Sugars\n";
  }
  ///CYCLE DETECTION
  // DetectCyclesByExhaustiveRingPerception gets stuck in infinite loops for some files, DFS doesn't.
  // CycleMap cycles = DetectCyclesByExhaustiveRingPerception();
  CycleMap cycles = DetectCyclesByDFS();
  ///PRINTING ALL DETECTED CYCLES
  if(local_debug > 0)
  {
    std::cout << "\n" << "All detected cycles\n";
  }
  for( CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++ )
  {
    std::string cycle_atoms_str = ( *it ).first;
    if(local_debug > 0)
    {
      std::cout << cycle_atoms_str << "\n";
    }
  }

  ///FILTERING OUT FUSED CYCLES. aka Cycles that are sharing an edge
  RemoveFusedCycles( cycles );

  ///FILTERING OUT OXYGENLESS CYCLES
  FilterAllCarbonCycles( cycles );
  if(local_debug > 0)
  {
    std::cout << "\n" << "Cycles after discarding rings that are all-carbon" << "\n";
    for( CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++ )
    {
      std::string cycle_atoms_str = ( *it ).first;
      std::cout << cycle_atoms_str << "\n";
    }
  }
  //FILTERING OUT CYCLES WITH DOUBLE BONDS and add them to new cyclemap
  CycleMap double_bond_cycles = FilterCyclesWithDoubleBonds( cycles );


  //TODO figure out what to do with double bonded cycles

  if(local_debug > 0)
  {
    std::cout << "\n" << "Cycles after discarding rings containing C-C double bonds" << "\n";
  }

  ///ANOMERIC CARBON DETECTION and SORTING
  std::vector< std::string > anomeric_carbons_status = std::vector< std::string >();
  std::vector< Glycan::Note* > anomeric_notes = std::vector< Glycan::Note* >();
  CycleMap sorted_cycles = CycleMap();
  MolecularModeling::AtomVector AllAnomericCarbons = MolecularModeling::AtomVector();
  for( CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++ )
  {
    MolecularModeling::AtomVector cycle_atoms = ( *it ).second;
    for (std::vector<MolecularModeling::Atom*>::iterator it = cycle_atoms.begin(); it != cycle_atoms.end(); it++ )
    {
      (*it)-> SetIsCycle(true);
    }
  }
  for( CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++ )
  {
    std::string cycle_atoms_str = ( *it ).first;
    MolecularModeling::AtomVector cycle_atoms = ( *it ).second;
    Glycan::Monosaccharide* mono = new Glycan::Monosaccharide(&cycle_atoms_str, cycle_atoms, this, CCD_Path);
    mono->assembly_ = this;
    monos.push_back(mono);
    if(local_debug > 0)
    {
      std::cout << cycle_atoms_str << "\n";
      ///e.g. C1_3810_NAG_A_1521_?_?_1-O5_3821_NAG_A_1521_?_?_1-
      ///C5_3814_NAG_A_1521_?_?_1-C4_3813_NAG_A_1521_?_?_1-C3_3812_NAG_A_1521_?_?_1-C2_3811_NAG_A_1521_?_?_1
    }
  }
  std::map <unsigned long long, Glycan::Monosaccharide* > ordered_monos_map =  std::map <unsigned long long, Glycan::Monosaccharide*>();
  for (std::vector <Glycan::Monosaccharide*>::iterator it = monos.begin(); it != monos.end(); it++)
  {
    Glycan::Monosaccharide* this_mono = *it;
    ordered_monos_map[this_mono->cycle_atoms_.at(0)->GetIndex()] = this_mono;
  }
  std::vector<Glycan::Monosaccharide*> ordered_monos;
  for(std::map <unsigned long long, Glycan::Monosaccharide* >::iterator it = ordered_monos_map.begin(); it != ordered_monos_map.end(); it++)
  {
    Glycan::Monosaccharide* this_mono = it->second;
    ordered_monos.push_back(this_mono);
    // std::cout << this_mono->cycle_atoms_str_ << "\n";
  }

  // ///CREATING MONOSACCHARIDE STRUCTURE. Ring atoms, side atoms, chemical code (Glycode), modifications/derivatives, names
  if(local_debug > 0)
  {
    std::cout << "\n" << "Detailed information of sorted cycles after discarding fused or oxygenless rings: " << "\n";
  }
  int mono_id_ = 0;
  for (std::vector<Glycan::Monosaccharide*>::iterator it = ordered_monos.begin(); it != ordered_monos.end(); it++)
  {
    Glycan::Monosaccharide* mono = *it;
    if(local_debug > 0)
    {
      std::cout << "Ring atoms: " << mono->cycle_atoms_str_ << "\n";

      ///PRINTING ASSIGNED SIDE ATOMS
      std::cout << "Side group atoms: " << "\n";
      for( std::vector< MolecularModeling::AtomVector >::iterator it1 = mono->side_atoms_.begin(); it1 != mono->side_atoms_.end(); it1++ )
      {
        MolecularModeling::AtomVector sides = ( *it1 );
        if( it1 == mono->side_atoms_.begin() ) //side atoms of anomeric carbon
        {
          if( sides.at( 0 ) != NULL && sides.at( 1 ) != NULL )
          {
            std::cout << "[1] -> " << sides.at( 0 )->GetId() << ", " << sides.at( 1 )->GetId() << "\n";
          }
          else if( sides.at( 1 ) != NULL )
          {
            std::cout << "[1] -> " << sides.at( 1 )->GetId() << "\n";
          }
          else if( sides.at( 0 ) != NULL )
          {
            std::cout << "[1] -> " << sides.at( 0 )->GetId() << "\n";
          }
        }
        else if( it1 == mono->side_atoms_.end() - 1 ) //side atoms of last carbon of the
        {
          std::cout << "[" << mono->cycle_atoms_.size() - 1 << "] -> ";
          if( sides.at( 0 ) != NULL )
          {
            std::cout << sides.at( 0 )->GetId() << "\n";
          }
        }
        else if( sides.at( 1 ) != NULL )
        {
          int cycle_atom_index = distance( mono->side_atoms_.begin(), it1 );
          std::cout << "[" << cycle_atom_index + 1 << "] -> " << sides.at( 1 )->GetId() << "\n";
        }
      }

      ///PRINTING ANOMERIC STATUS
      std::cout << mono->anomeric_status_ << mono->cycle_atoms_.at( 0 )->GetId() << "\n";
      std::cout << "\n" << "Ring Stereo chemistry chemical code:"  << "\n";
      std::cout << mono->chemical_code_->toString() << "\n";
      mono->chemical_code_->Print( std::cout );
      std::cout << "\n";
      // Actually, it only works for hexoses, so this should read "!=6" rather than ">5"
      if( mono->cycle_atoms_.size() > 5 )
      {
        if( mono->bfmp_ring_conformation_.compare( "" ) != 0 )
        {
          std::cout << "BFMP ring conformation: " << mono->bfmp_ring_conformation_ << "\n" << "\n"; ///Part of Glyprobity report
        }
      }

      // ///IF NO DERIVATIVES THEN NO NEED TO ITERATE THROUGH WHATS NOT THERE
      // /// Also, this now prints Ring Position instead of Position. e.g: -1, a( 1 ), 2, 3, 4, 5, +1, +2, +3
      if( !mono->derivatives_map_.empty() )
      {
        for( std::vector<std::pair< std::string, std::string> >::iterator it1 = mono->derivatives_map_.begin(); it1 != mono->derivatives_map_.end(); it1++ )
        {
          std::string key = ( *it1 ).first;
          std::string value = ( *it1 ).second;
          std::string trimmedValue;
          if(value.substr(0,2) == "xC")
          {
            trimmedValue = value.substr(2,value.size()-1);
          }
          else
          {
            trimmedValue = value;
          }
          if(value != "")
            std::cout << "Carbon at Ring Position " << key << " is attached to " << trimmedValue << "\n";
        }
      }
    }

    ///GLYPROBITY REPORT (GEOMETRY OUTLIERS)
    if( glyprobity_report ) {
      CalculateGlyprobityGeometryOutliers( mono );
    }

    ///PRINTING NAMES OF MONOSACCHARIDE
    if(local_debug > 0)
    {
      std::cout << "Stereochemistry name: " << mono->sugar_name_.monosaccharide_stereochemistry_name_ << "\n";
      std::cout << "Stereochemistry short name: " << mono->sugar_name_.monosaccharide_stereochemistry_short_name_ << "\n";
      std::cout << "Complete name: " << mono->sugar_name_.monosaccharide_name_ << "\n";
      std::cout << "Short name: " << mono->sugar_name_.monosaccharide_short_name_ << "\n";
    }

    //Check author's naming vs what's detected
    if(CCD_Path != " ")
    {
      if(mono->sugar_name_.pdb_code_.find(mono->residue_name_) == std::string::npos)
      {
        if(this->source_file_.find(mono->residue_name_) ==  std::string::npos) //For CCD Lookup, if the residue name matches the file name we are already looking it up.  Prevents infinite loop for unidentified CCD sugars
        {
          if(local_debug > 0)
          {
            std::cout << this->source_file_ << ": " << mono->residue_name_ << "\n";
          }
          GetAuthorNaming(amino_lib_files, mono, CCD_Path);
          mono->createAuthorSNFGname();

          if(mono->sugar_name_.monosaccharide_name_ != mono->author_sugar_name_.monosaccharide_name_)
          {
            Glycan::Note* mismatch_note = new Glycan::Note();
            mismatch_note->type_ = Glycan::ERROR;
            mismatch_note->category_ = Glycan::DER_MOD;
            std::stringstream note;
            note << "Residue name, " << mono->cycle_atoms_[0]->GetResidue()->GetName() << " (" << mono->cycle_atoms_[0]->GetResidue()->GetId() << "), in input PDB file for " << mono->sugar_name_.monosaccharide_short_name_ << " does not match GlyFinder residue code: " << mono->sugar_name_.pdb_code_ << ". "/* << mono->sugar_name_.monosaccharide_name_ << " vs. " << mono->author_sugar_name_.monosaccharide_name_*/;
            mismatch_note->description_ = note.str();
            // this->AddNote(mismatch_note);
            mono->mono_notes_.push_back(mismatch_note);
          }
        }
      }
      else
      {
        if(this->source_file_.find(mono->residue_name_) ==  std::string::npos) //For CCD Lookup, if the residue name matches the file name we are already looking it up.  Prevents infinite loop for unidentified CCD sugars
        {
          if(local_debug > 0)
          {
            std::cout << this->source_file_ << ": " << mono->residue_name_ << "\n";
          }
          GetAuthorNaming(amino_lib_files, mono, CCD_Path);
          mono->createAuthorSNFGname();
        }
        // mono->author_sugar_name_ = mono->sugar_name_;
        // mono->createAuthorSNFGname();
      }
      //replace KDO/N with Kdo/n for the PDB (CCD is always included for the PDB runs)
      if (mono->sugar_name_.monosaccharide_short_name_.find("KDO") != std::string::npos)
      {
        mono->sugar_name_.monosaccharide_short_name_.replace(1, 3, "Kdo");
      }
      else if (mono->sugar_name_.monosaccharide_short_name_.find("KDN") != std::string::npos)
      {
        mono->sugar_name_.monosaccharide_short_name_.replace(1, 3, "Kdn");
      }
    }
    if(local_debug > 0)
    {
      std::cout << "SNFG Name: " << mono->SNFG_name_ << "\n";
      std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------\n";
    }
    mono_id_++;
    mono->mono_id_ = mono_id_;
    *it = mono;
  }

  //Checking if any sugar named residue is not detected
  MolecularModeling::ResidueVector all_residues = GetAllResiduesOfAssembly();
  for(MolecularModeling::ResidueVector::iterator it = all_residues.begin(); it != all_residues.end(); it++)
  {
    MolecularModeling::Residue* this_residue = *it;
    std::string residue_name = this_residue->GetName();
    Glycan::SugarName residue_sugar_name = gmml::ResidueSugarNameLookup(residue_name);
    Glycan::SugarName residue_complex_sugar_name = gmml::ResidueComplexSugarNameLookup(residue_name);
    if((residue_sugar_name.pdb_code_ != "") || (residue_complex_sugar_name.pdb_code_ != "")) //This is a sugar residue
    {
      if (this_residue->GetIsSugar() == false)
      {
        if(this_residue->GetId()[12] == '?')//alternate atomic locations caused new residues that weren't assigned as sugars.  A ? at 12 in the ID
                                            //means it isnt an atom with alternate coordinates
        {
          Glycan::Note* undetected_note = new Glycan::Note();
          undetected_note->type_ = Glycan::ERROR;
          undetected_note->category_ = Glycan::RESIDUE_NAME;
          std::stringstream note;
          note << "Residue " << this_residue->GetName() << " (" << this_residue->GetId();
          note << "), in PDB file indicates a sugar residue, but no sugar was detected.";
          undetected_note->description_ = note.str();
          this->AddNote(undetected_note);
        }
      }
    }
  }

  ///CREATING TREE-LIKE STRUCTURE OF OLIGOSACCHARIDE
  if(local_debug > 0)
  {
    std::cout << "\n" << "Oligosaccharides:" << "\n";
  }
  int number_of_covalent_links = 0;
  int number_of_probable_non_covalent_complexes = 0;
  // std::vector< Glycan::Oligosaccharide* > oligosaccharides = ExtractOligosaccharides( ordered_monos, dataset_residue_names, number_of_covalent_links, number_of_probable_non_covalent_complexes );
  createOligosaccharideGraphs(ordered_monos, dataset_residue_names, number_of_covalent_links, number_of_probable_non_covalent_complexes);
  // Glycan::Oligosaccharide* testOligo = new Glycan::Oligosaccharide(ordered_monos, dataset_residue_names, this);
  // Glycan::Oligosaccharide* testOligo = new Glycan::Oligosaccharide();
  std::vector<Glycan::Oligosaccharide*> testOligos = createOligosaccharides(ordered_monos);
  testOligos = ExtractOligosaccharides( ordered_monos, dataset_residue_names, number_of_covalent_links, number_of_probable_non_covalent_complexes );
  if(local_debug > 0)
  {
    std::cout << std::to_string(testOligos.size()) << "\n";
  }
  // testOligo.createOligosaccharideGraphs(ordered_monos);

  ///BUILDING OLIGOSACCHARIDE SEQUENCE
  int number_of_oligosaccharides = 0;
  int number_of_monosaccharides = 0;

  for(std::vector< Glycan::Oligosaccharide* >::iterator it = testOligos.begin(); it != testOligos.end(); it++)
  {
    Glycan::Oligosaccharide* thisOligo = *it;
    if(local_debug > 0)
    {
      std::cout << "Oligo IUPAC Name:       " << thisOligo->IUPAC_name_ << "\n";
      std::cout << "Oligo author IUPAC Name:" << thisOligo->author_IUPAC_name_ << "\n";
      std::cout << "Oligo Name:             ";
    }
    thisOligo->Print( std::cout );// This for some reason does a ton of stuff instead of printing....
  }

  ///PRINTING NOTES AND ISSUES FOUND WITH THE INPUT FILE IF THERE ARE ANY NOTES
  std::vector< Glycan::Note* > notes = this->GetNotes();
  if(local_debug > 0)
  {
    if( !notes.empty() ) {
      std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << "\n";
      std::cout << "\n" << "NOTES/ISSUES:" << "\n";
      for( std::vector< Glycan::Note* >::iterator note_it = notes.begin(); note_it != notes.end(); note_it++ ) {
        Glycan::Note* note = ( *note_it );
        std::cout << "\n" << "Category: " << note->ConvertGlycanNoteCat2String( note->category_ ) << "\n";
        std::cout << "Type: " << note->ConvertGlycanNoteType2String( note->type_ ) << "\n";
        std::cout << "Description: " << note->description_ << "\n";
      }
      std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << "\n";
    }

    ///PRINTING STATISTICAL REPORT OF GLYPROBITY
    if( glyprobity_report ) {
      std::cout << "\n" << "GLYPROBITY REPORT" << "\n";
      std::cout << "<-------Topology------>" << "\n";
      std::cout << "Monosaccharide detected: " << monos.size() << "\n";
      std::cout << "Residue Distribution " << "\n";
      std::cout << " Monosaccharides: " << number_of_monosaccharides << "\n";
      std::cout << " Oligosaccharides: " << number_of_oligosaccharides << "\n";
      std::cout << "Carbohydrate Context " << "\n";
      std::cout << " Covalently linked to protein: " << number_of_covalent_links << "\n";
      std::cout << " Non-covalent complex: " << number_of_probable_non_covalent_complexes << "\n";
      std::cout << "<--------------------->" << "\n";
    }
  }
  ///POPULATING GMMO ONTOLOGY
  if( populate_ontology )
  {
    if( testOligos.size() > 0 )
    {
      std::ofstream out_file;
      if(individualOntologies)
      {
        std::string gmmo = this->GetSourceFile();
        gmmo = gmmo.substr(0,gmmo.size()-4) + ".ttl";
        
        // gmmo.insert(gmmo.size()-8, gmmo.substr(gmmo.size()-7, 2));
        // gmmo.insert(gmmo.size()-8, "/");
        // std::string gmmoDirectory = gmmo.substr(0, gmmo.size()-8);
        std::string ontologyDirectory = "Ontologies";
        std::string gmmoDirectory = ontologyDirectory + "/" + gmmo.substr(gmmo.size()-7, 2);
        gmmo = ontologyDirectory + "/" + gmmo.substr(gmmo.size()-7, 2) + "/" + gmmo;        
        mkdir(ontologyDirectory.c_str(),  S_IRWXU | S_IRWXG | S_IRWXO);
        mkdir(gmmoDirectory.c_str(),  S_IRWXU | S_IRWXG | S_IRWXO);
        out_file.open( gmmo.c_str(), std::fstream::out | std::fstream::trunc);
        out_file << Ontology::TTL_FILE_PREFIX << "\n";
      }
      else
      {
        std::string gmmo = "gmmo.ttl";
        out_file.open( gmmo.c_str(), std::fstream::app );
        std::ifstream in( gmmo );///Checking if the file is empty
        size_t out_file_size = 0;
        in.seekg( 0,std::ios_base::end );
        out_file_size = in.tellg();
        in.close();
        if( out_file_size == 0 ) ///If the file is empty add the prefixes first
        {
          out_file << Ontology::TTL_FILE_PREFIX << "\n";
        }
      }
      // this->PopulateOntology( out_file, oligosaccharides );
      this->PopulateOntology( out_file, testOligos );
      out_file.close();
    }
  }

  return testOligos;
}

bool Assembly::MatchDisaccharide(std::queue<Glycan::Oligosaccharide*> oligo_queue, double &phi_angle, double &psi_angle, std::string first_mono, char mono1_carbon_index, std::string second_mono, char mono2_carbon_index)
{
    Glycan::Oligosaccharide* oligo = oligo_queue.front();
    oligo_queue.pop();
    bool found_disaccharide = false;

    Glycan::Oligosaccharide* corresponding_second_oligo = NULL;
    ///second mono of the input disaccharide comes first in the tree structure of oligo (parent of first mono)
    if(oligo->root_->sugar_name_.monosaccharide_short_name_.compare(second_mono) == 0 )///found a mono with the same name as the right side mono of the disaccharide
        corresponding_second_oligo = oligo;

    OligosaccharideVector child_oligos = oligo->child_oligos_;
    for(OligosaccharideVector::iterator it1 = child_oligos.begin(); it1 != child_oligos.end(); it1++)
    {
        oligo_queue.push(*it1);

        if(corresponding_second_oligo != NULL) ///if current mono matches the disaccharide, look for a linked mono that macthes the other mono in the disaccharide
        {
            Glycan::Oligosaccharide* corresponding_first_oligo = (*it1);
            oligo_queue.push(corresponding_first_oligo);

            if(corresponding_first_oligo->root_->sugar_name_.monosaccharide_short_name_.compare(first_mono) == 0)///found a mono with the same name as the left side mono of the disaccharide
            {
                std::vector<std::string> child_links = corresponding_second_oligo->child_oligos_linkages_;///links from right mono of the disaccharide to child monos

                std::vector<std::string> mono2_cycle_atom_tokens = gmml::Split(corresponding_second_oligo->root_->cycle_atoms_str_, "-");
                for(std::vector<std::string>::iterator it2 = child_links.begin(); it2 != child_links.end(); it2++)
                {
                    //                int index = distance(child_links.begin(), it2);
                    std::string link = (*it2);
                    std::vector<std::string> link_tokens = gmml::Split(link, "-");
                    int parent_c_index = 0;
                    int child_c_index = 0;

                    if(corresponding_second_oligo->root_->side_atoms_.at(0).at(0) != NULL)
                        parent_c_index++;
                    for(unsigned int i = 0; i < mono2_cycle_atom_tokens.size(); i++)
                    {
                        parent_c_index++;
                        if(mono2_cycle_atom_tokens.at(i).compare(link_tokens.at(0)) == 0)
                            break;
                    }
                    std::vector<std::string> mono1_cycle_atom_tokens = gmml::Split(corresponding_first_oligo->root_->cycle_atoms_str_, "-");
                    if(corresponding_first_oligo->root_->side_atoms_.at(0).at(0) != NULL)
                        child_c_index++;
                    for(unsigned int i = 0; i < mono1_cycle_atom_tokens.size(); i++)
                    {
                        child_c_index++;
                        if(mono1_cycle_atom_tokens.at(i).compare(link_tokens.at(2)) == 0)
                            break;
                    }

                    if(mono2_carbon_index == parent_c_index + '0' && mono1_carbon_index == child_c_index + '0')///indeces matched the indeces of the given disaccharide
                    {
                        found_disaccharide = true;
                        MolecularModeling::AtomVector mono1_ring_atoms = corresponding_first_oligo->root_->cycle_atoms_;

                        ///Preparing atoms for phi and psi angle
                        MolecularModeling::Atom* phi_atom1 = mono1_ring_atoms.at(mono1_ring_atoms.size() - 1); ///O5 or N5
                        MolecularModeling::Atom* phi_atom2 = new MolecularModeling::Atom(); ///C1
                        phi_atom2 = corresponding_first_oligo->root_->cycle_atoms_.at(0); ///anomeric carbon
                        MolecularModeling::Atom* phi_atom3 = NULL;
                        MolecularModeling::AtomVector phi_atom2_neighbors = phi_atom2->GetNode()->GetNodeNeighbors();
                        for(MolecularModeling::AtomVector::iterator it3 = phi_atom2_neighbors.begin(); it3 != phi_atom2_neighbors.end(); it3++)
                        {
                            MolecularModeling::Atom* atom2_neighbor= (*it3);
                            ///If the neighbor id is the same as the linkage intermediate atom id
                            if(atom2_neighbor->GetId().compare(link_tokens.at(1)) == 0)
                            {
                                phi_atom3 = atom2_neighbor; ///Ox
                                break;
                            }
                        }
                        if(phi_atom3 != NULL)
                        {
                            MolecularModeling::Atom* phi_atom4 = NULL;
                            MolecularModeling::AtomVector phi_atom3_neighbors = phi_atom3->GetNode()->GetNodeNeighbors();
                            for(MolecularModeling::AtomVector::iterator it4 = phi_atom3_neighbors.begin(); it4 != phi_atom3_neighbors.end(); it4++)
                            {
                                MolecularModeling::Atom* atom3_neighbor= (*it4);
                                ///If the neighbor id is the same as the linkage carbon at second mono side (first linkage atom from second mono side)
                                if(atom3_neighbor->GetId().compare(link_tokens.at(0)) == 0)
                                {
                                    phi_atom4 = atom3_neighbor; ///Cx
                                    phi_angle = CalculateTorsionAngleByAtoms(phi_atom1, phi_atom2, phi_atom3, phi_atom4); /// ϕ (O5′-C1′-Ox-Cx)
                                    gmml::log(__LINE__, __FILE__, gmml::INF, std::to_string(phi_angle));
                                    MolecularModeling::AtomVector phi_atom4_neighbors = phi_atom4->GetNode()->GetNodeNeighbors();

                                    int atom4_index = gmml::ConvertString<int>(gmml::Split(phi_atom4->GetName(), "C*,\'").at(0));
                                    for(MolecularModeling::AtomVector::iterator it5 = phi_atom4_neighbors.begin(); it5 != phi_atom4_neighbors.end(); it5++)
                                    {
                                        MolecularModeling::Atom* atom4_neighbor = (*it5);
                                        std::string neighbor_name = atom4_neighbor->GetName();
                                        if(neighbor_name.find("C") != std::string::npos)
                                        {
                                            int neighbor_index = gmml::ConvertString<int>(gmml::Split(neighbor_name, "C*,\'").at(0));
                                            if(atom4_index > neighbor_index) ///Cx-1
                                            {
                                                psi_angle = CalculateTorsionAngleByAtoms(phi_atom2, phi_atom3, phi_atom4, atom4_neighbor); /// ψ (C1′-Ox-Cx-Cx−1)
                                                break;
                                            }
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                    if(found_disaccharide)
                        break;
                }
                if(found_disaccharide)
                    break;
            }
        }
    }
    if(!found_disaccharide && !oligo_queue.empty())
        return MatchDisaccharide(oligo_queue, phi_angle, psi_angle, first_mono, mono1_carbon_index, second_mono, mono2_carbon_index);
    else
        return found_disaccharide;
}

// void Assembly::ExtractRingAtomsInformation()
// {
//     ///CYCLE DETECTION
//     CycleMap cycles = DetectCyclesByExhaustiveRingPerception();
//     ///FILTERING OUT FUSED CYCLES
//     RemoveFusedCycles(cycles);
//     ///FILTERING OUT OXYGENLESS CYCLES
//     FilterAllCarbonCycles(cycles);
//     CycleMap sorted_cycles = CycleMap();
//     std::vector<std::string> anomeric_carbons_status = std::vector<std::string>();
//     ///ANOMERIC DETECTION and SORTING
//     for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
//     {
//         std::string cycle_atoms_str = (*it).first;
//         MolecularModeling::AtomVector cycle_atoms = (*it).second;
//         Glycan::Note* anomeric_note = new Glycan::Note();
//         MolecularModeling::Atom* anomeric = Glycan::Monosaccharide::FindAnomericCarbon(anomeric_note, anomeric_carbons_status, cycle_atoms, cycle_atoms_str);
//
//         if(anomeric != NULL)
//         {
//             MolecularModeling::AtomVector sorted_cycle_atoms = MolecularModeling::AtomVector();
//             std::stringstream sorted_cycle_stream;
//             sorted_cycle_atoms = SortCycle( cycle_atoms, anomeric, sorted_cycle_stream );
//             sorted_cycles[sorted_cycle_stream.str()] = sorted_cycle_atoms;
//         }
//     }
//     cycles = sorted_cycles;
//     std::cout << "\n" << "Detailed information of sorted cycles after discarding fused or oxygenless rings: " << "\n";
//     gmml::log(__LINE__, __FILE__,  gmml::INF, "Detailed information of sorted cycles after discarding fused or oxygenless rings: ");
//
//     for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
//     {
//         Glycan::Monosaccharide* mono = new Glycan::Monosaccharide();
//         int status_index = distance(cycles.begin(), it);
//         mono->anomeric_status_ = anomeric_carbons_status.at(status_index);
//         std::string cycle_atoms_str = (*it).first;
//         MolecularModeling::AtomVector cycle = (*it).second;
//
//         std::stringstream ring_atoms;
//         ring_atoms << "Ring atoms: " << cycle_atoms_str;
//         std::cout << ring_atoms.str() << "\n";
//         gmml::log(__LINE__, __FILE__,  gmml::INF, ring_atoms.str());
//
//         std::stringstream ring_atom_names;
//         ring_atom_names << "Ring atom names: " ;
//         for(MolecularModeling::AtomVector::iterator it1 = cycle.begin(); it1 != cycle.end(); it1++)
//         {
//             MolecularModeling::Atom* atom = (*it1);
//             std::vector<std::string> atom_tokens = gmml::Split(atom->GetId(), "_");
//             if(it1 == cycle.end() - 1)
//                 ring_atom_names << atom_tokens.at(0);
//             else
//                 ring_atom_names << atom_tokens.at(0) << ",";
//         }
//         std::cout << ring_atom_names.str() << "\n";
//         gmml::log(__LINE__, __FILE__,  gmml::INF, ring_atom_names.str());
//     }
// }

void Assembly::ReturnCycleAtoms(std::string src_id, MolecularModeling::Atom *current_atom, AtomIdAtomMap &atom_parent_map, MolecularModeling::AtomVector &cycle, std::stringstream &cycle_stream)
{
    cycle.push_back(current_atom);
    cycle_stream << current_atom->GetId() << "-";
    MolecularModeling::Atom* parent = atom_parent_map[current_atom->GetId()];
    if(src_id.compare(parent->GetId()) == 0)
    {
        cycle.push_back(parent);
        cycle_stream << parent->GetId();
        return;
    }
    ReturnCycleAtoms(src_id, parent, atom_parent_map, cycle, cycle_stream);
}
Assembly::CycleMap Assembly::FilterCyclesWithDoubleBonds(CycleMap &cycles)
{
  int local_debug = -1;
  if (local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__,  gmml::INF, "In Filter Double bond function");
  }
  std::map<std::string, bool> to_be_deleted_cycles = std::map<std::string, bool>();
    bool all_single_bonds = true;
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        std::string cycle_str = (*it).first;
        MolecularModeling::AtomVector cycle_atoms = (*it).second;
        all_single_bonds = true;
        std::stringstream debugStr;
        debugStr << cycle_atoms.size();
        if (local_debug > 0)
        {
          gmml::log(__LINE__, __FILE__,  gmml::INF, cycle_str);
          gmml::log(__LINE__, __FILE__,  gmml::INF, debugStr.str());
        }
        for(int it1 = 0; it1 != (int) cycle_atoms.size() -1; it1++)
        {

            MolecularModeling::Atom* atom1 = cycle_atoms[it1];
            if (local_debug > 0)
            {
              gmml::log(__LINE__, __FILE__,  gmml::INF, atom1->GetId());
            }
            MolecularModeling::Atom* atom2;
            bool doubleBond = false;
            if(it1 == 0)
            {
              atom2 = cycle_atoms[cycle_atoms.size()-1];
            }
            else
            {
              atom2 = cycle_atoms[it1+1];
            }
            doubleBond = guessIfC_CDoubleBond(atom1, atom2);
            if(doubleBond == true)
            {
                all_single_bonds = false;
                break;
            }
        }
        if(all_single_bonds == false)
            to_be_deleted_cycles[cycle_str] = true;
    }
    CycleMap all_single_bond_filtered_cycles = CycleMap();
    CycleMap all_double_bond_cycles = CycleMap();
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        std::string cycle_str = (*it).first;
        MolecularModeling::AtomVector cycle_atoms = (*it).second;
        if(to_be_deleted_cycles.find(cycle_str) == to_be_deleted_cycles.end())
        {
          all_single_bond_filtered_cycles[cycle_str] = cycle_atoms;
        }

        else
        {
          all_double_bond_cycles[cycle_str] = cycle_atoms;
        }
    }
    cycles.clear();
    cycles = all_single_bond_filtered_cycles;
    return all_double_bond_cycles;
}


void Assembly::FilterAllCarbonCycles(CycleMap &cycles)
{
    std::map<std::string, bool> to_be_deleted_cycles = std::map<std::string, bool>();
    bool all_carbons = true;
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        std::string cycle_str = (*it).first;
        MolecularModeling::AtomVector cycle_atoms = (*it).second;
        all_carbons = true;
        for(MolecularModeling::AtomVector::iterator it1 = cycle_atoms.begin(); it1 != cycle_atoms.end(); it1++)
        {
            MolecularModeling::Atom* atom = (*it1);
            if(atom->GetElementSymbol() != "C")
            {
                all_carbons = false;
                break;
            }
        }
        if(all_carbons)
            to_be_deleted_cycles[cycle_str] = true;
    }
    CycleMap all_carbons_filtered_cycles = CycleMap();
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        std::string cycle_str = (*it).first;
        MolecularModeling::AtomVector cycle_atoms = (*it).second;
        if(to_be_deleted_cycles.find(cycle_str) == to_be_deleted_cycles.end())
            all_carbons_filtered_cycles[cycle_str] = cycle_atoms;
    }
    cycles.clear();
    cycles = all_carbons_filtered_cycles;
}

void Assembly::RemoveFusedCycles(CycleMap &cycles)
{
    std::map<std::string, bool> to_be_deleted_cycles = std::map<std::string, bool>();
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        std::string cycle_i_str = (*it).first; ///cycle i will be compared with all the other cycles
        MolecularModeling::AtomVector cycle_i_atoms = (*it).second;
        for(CycleMap::iterator it1 = cycles.begin(); it1 != cycles.end(); it1++)
        {
            if(it != it1)
            {
                std::string cycle_j_str = (*it1).first; ///cycle j to be compared with cycle i
                for(unsigned int i = 0; i < cycle_i_atoms.size(); i++)
                {
                    std::stringstream mutual_edge;
                    std::stringstream mutual_edge_reverse;
                    MolecularModeling::Atom* a1 = new MolecularModeling::Atom();
                    MolecularModeling::Atom* a2 = new MolecularModeling::Atom();
                    if(i == cycle_i_atoms.size() - 1)
                    {
                        a1 = cycle_i_atoms.at(i);
                        a2 = cycle_i_atoms.at(0);
                    }
                    else
                    {
                        a1 = cycle_i_atoms.at(i);
                        a2 = cycle_i_atoms.at(i + 1);
                    }
                    mutual_edge << a1->GetId() << "-" << a2->GetId();
                    mutual_edge_reverse << a2->GetId() << "-" << a1->GetId();
                    if(cycle_j_str.find(mutual_edge.str()) != std::string::npos || cycle_j_str.find(mutual_edge_reverse.str()) != std::string::npos)///mutual edge found
                    {
                        to_be_deleted_cycles[cycle_i_str] = true;
                        to_be_deleted_cycles[cycle_j_str] = true;
                        break;
                    }
                }
            }
        }
    }
    CycleMap fused_filtered_cycles = CycleMap();
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        std::string cycle_str = (*it).first;
        MolecularModeling::AtomVector cycle_atoms = (*it).second;
        if(to_be_deleted_cycles.find(cycle_str) == to_be_deleted_cycles.end())
            fused_filtered_cycles[cycle_str] = cycle_atoms;
    }
    cycles.clear();
    cycles = fused_filtered_cycles;
}

// MolecularModeling::Atom* Assembly::FindAnomericCarbon( Glycan::Note* anomeric_note, std::vector< std::string > & anomeric_carbons_status, MolecularModeling::AtomVector cycle, std::string cycle_atoms_str ) {
//   MolecularModeling::Atom* anomeric_carbon = new MolecularModeling::Atom();
//   for( MolecularModeling::AtomVector::iterator it = cycle.begin(); it != cycle.end(); it++ ) {
//     MolecularModeling::Atom* cycle_atom = ( *it );
//
//     if( ( cycle_atom->GetName().substr( 0, 1 ).compare( "O" ) == 0 ) ) {///find oxygen in ring
//       //                && isdigit(gmml::ConvertString<char>(cycle_atom->GetName().substr(1,1)))))
//       MolecularModeling::AtomNode* node = cycle_atom->GetNode();
//       MolecularModeling::AtomVector neighbors = node->GetNodeNeighbors();
//
//       ///Check the first neighbor of oxygen
//       MolecularModeling::Atom* o_neighbor1 = neighbors.at( 0 );
//       MolecularModeling::AtomNode* o_neighbor1_node = o_neighbor1->GetNode();
//       MolecularModeling::AtomVector o_neighbor1_neighbors = o_neighbor1_node->GetNodeNeighbors();
//
//       for( MolecularModeling::AtomVector::iterator it1 = o_neighbor1_neighbors.begin(); it1 != o_neighbor1_neighbors.end(); it1++ ) {///check if neighbor1 of oxygen has another oxygen or nitrogen neighbor
//         MolecularModeling::Atom* o_neighbor1_neighbor = ( *it1 );
//
//         if( cycle_atoms_str.find( o_neighbor1_neighbor->GetId() ) == std::string::npos ///if the neighbor is not one of the cycle atoms
//                 && ( o_neighbor1_neighbor->GetName().substr( 0, 1 ).compare( "O" ) == 0 || o_neighbor1_neighbor->GetName().substr( 0, 1 ).compare( "N" ) == 0 ) ) { ///if first element is "O" or "N"
//           //                        && isdigit(gmml::ConvertString<char>(neighbor1_neighbor->GetName().substr(1,1))))///if second element is a digit
//           anomeric_carbon = o_neighbor1;
//           anomeric_carbons_status.push_back( "Anomeric carbon: "  );
//           anomeric_note->description_ = "";
//
//           return anomeric_carbon;
//         }
//       }
//
//       ///Check the second neighbor of oxygen
//       MolecularModeling::Atom* o_neighbor2 = neighbors.at( 1 );
//       MolecularModeling::AtomNode* o_neighbor2_node = o_neighbor2->GetNode();
//       MolecularModeling::AtomVector o_neighbor2_neighbors = o_neighbor2_node->GetNodeNeighbors();
//
//       for( MolecularModeling::AtomVector::iterator it2 = o_neighbor2_neighbors.begin(); it2 != o_neighbor2_neighbors.end(); it2++ ) {///check if neighbor2 of oxygen has another oxygen or nitrogen neighbor
//         MolecularModeling::Atom* o_neighbor2_neighbor = ( *it2 );
//
//         if( cycle_atoms_str.find( o_neighbor2_neighbor->GetId() ) == std::string::npos
//                 && ( o_neighbor2_neighbor->GetName().substr( 0, 1 ).compare( "O" ) == 0 || o_neighbor2_neighbor->GetName().substr( 0, 1 ).compare( "N" ) == 0 ) ) {
//           //                        && isdigit(gmml::ConvertString<char>(neighbor2_neighbor->GetName().substr(1,1))))
//           anomeric_carbon = o_neighbor2;
//           anomeric_carbons_status.push_back( "Anomeric carbon: "  );
//           anomeric_note->description_ = "";
//
//           return anomeric_carbon;
//         }
//       }
//
//       ///Seems redundant to have this multiple times, if it is going to say the same thing and be generated for all cases after the first logic check.
//       ///Specially if we ever want to change this Note to say something different. Like I am doing now. :)
//       std::stringstream ss;
//       ss << "Could not find glycosidic oxygen or nitrogen within " << gmml::dCutOff << " Angstroms of a ring carbon bonded to the ring oxygen(" << node->GetAtom()->GetName() << ").";
//       anomeric_note->description_ = ss.str();
//
//       ///Check the order of the carbons based on their names to locate the anomeric
//       std::stringstream ss1;
//       for( unsigned int i = 0; i < o_neighbor1->GetName().size(); i++ ) {
//         if( isdigit( o_neighbor1->GetName().at( i )) != 0 )
//           ss1 << o_neighbor1->GetName().at( i );
//       }
//       std::stringstream ss2;
//       for( unsigned int i = 0; i < o_neighbor2->GetName().size(); i++ ) {
//         if( isdigit( o_neighbor2->GetName().at( i )) != 0 )
//           ss2 << o_neighbor2->GetName().at(i);
//       }
//       if( gmml::ConvertString< int >( ss1.str() ) < gmml::ConvertString< int >( ss2.str() ) ) {
//         anomeric_note->type_ = Glycan::WARNING;
//         anomeric_note->category_ = Glycan::ANOMERIC;
//         anomeric_carbons_status.push_back( "Anomeric carbon assigned based on its atom name index (" + ss1.str() + "), Anomeric Carbon is: " + anomeric_carbon->GetName() );
//
//         return o_neighbor1;
//       }
//       if( gmml::ConvertString< int >( ss2.str() ) < gmml::ConvertString< int >( ss1.str() ) ) {
//         anomeric_note->type_ = Glycan::WARNING;
//         anomeric_note->category_ = Glycan::ANOMERIC;
//         anomeric_carbons_status.push_back( "Anomeric carbon assigned based on its atom name index (" + ss1.str() + "), Anomeric Carbon is: " + anomeric_carbon->GetName() );
//
//         return o_neighbor2;
//       }
//
//       ///Check non-ring neighbors of oxygen neighbors (the one without non-ring carbon is anomeric)
//       bool neighbor2_is_anomeric = false;
//       for( MolecularModeling::AtomVector::iterator it1 = o_neighbor1_neighbors.begin(); it1 != o_neighbor1_neighbors.end(); it1++ ) {///check if neighbor1 of oxygen has non-ring carbon neighbor
//         MolecularModeling::Atom* neighbor1_neighbor = ( *it1 );
//         if( cycle_atoms_str.find( neighbor1_neighbor->GetId()) == std::string::npos ///if the neighbor is not one of the cycle atoms
//                 && ( neighbor1_neighbor->GetName().substr( 0, 1 ).compare( "C" ) == 0 ) ) {///if first element is "C"
//           neighbor2_is_anomeric = true;
//           break;
//         }
//       }
//       bool neighbor1_is_anomeric = false;
//       for( MolecularModeling::AtomVector::iterator it1 = o_neighbor2_neighbors.begin(); it1 != o_neighbor2_neighbors.end(); it1++ ) {///check if neighbor1 of oxygen has non-ring carbon neighbor
//         MolecularModeling::Atom* neighbor2_neighbor = ( *it1 );
//         if( cycle_atoms_str.find( neighbor2_neighbor->GetId()) == std::string::npos ///if the neighbor is not one of the cycle atoms
//                 && ( neighbor2_neighbor->GetName().substr( 0, 1 ).compare( "C" ) == 0 ) ) {///if first element is "C"
//           neighbor1_is_anomeric = true;
//         }
//       }
//       if( !neighbor1_is_anomeric ) {
//         anomeric_note->type_ = Glycan::WARNING;
//         anomeric_note->category_ = Glycan::ANOMERIC;
//         anomeric_carbons_status.push_back( "Anomeric carbon assigned to the ring carbon neighboring the ring oxygen but is not attached to an exocyclic carbon (that is, assuming the monosaccharide is an aldose), Anomeric Carbon is: " + anomeric_carbon->GetName() );
//         return o_neighbor2;
//       }
//       else if( !neighbor2_is_anomeric ) {
//         anomeric_note->type_ = Glycan::WARNING;
//         anomeric_note->category_ = Glycan::ANOMERIC;
//         anomeric_carbons_status.push_back( "Anomeric carbon assigned to the ring carbon neighboring the ring oxygen but is not attached to an exocyclic carbon (that is, assuming the monosaccharide is an aldose), Anomeric Carbon is: " + anomeric_carbon->GetName() );
//         return o_neighbor1;
//       }
//       // @TODO August 8, 2017 Davis/Lachele - Once we get to a point where we read in mmcif definitions, we need to use it to determine the Anomeric Carbon.
//       anomeric_note->type_ = Glycan::WARNING;
//       anomeric_note->category_ = Glycan::ANOMERIC;
//       anomeric_carbons_status.push_back( "Anomeric carbon assigned to the first ring carbon (" + o_neighbor1->GetName() + "), Anomeric Carbon is: " + anomeric_carbon->GetName() );
//       return o_neighbor1;
//     }
//   }
//   return NULL;
// }

MolecularModeling::AtomVector Assembly::SortCycle(MolecularModeling::AtomVector cycle, MolecularModeling::Atom *anomeric_atom, std::stringstream& sorted_cycle_stream)
{
    MolecularModeling::AtomVector sorted_cycle = MolecularModeling::AtomVector();
    for(MolecularModeling::AtomVector::iterator it = cycle.begin(); it != cycle.end(); it++)
    {
        MolecularModeling::Atom* atom = (*it);
        unsigned int index = distance(cycle.begin(), it);
        if(atom->GetId().compare(anomeric_atom->GetId()) == 0)
        {
            if(index == cycle.size() - 1)///anomeric atom is at the end of the cycle
            {
                sorted_cycle.push_back(anomeric_atom);///anomeric atom as the first atom of the cycle
                sorted_cycle_stream << anomeric_atom->GetId() << "-";
                anomeric_atom->SetIsRing(true);
                MolecularModeling::Atom* a0 = cycle.at(0);
                if(a0->GetName().substr(0,1).compare("O") == 0)///a0 is oxygen so the std::vector is in reverse order
                {
                    for(MolecularModeling::AtomVector::iterator it1 = it - 1; it1 != cycle.begin(); it1--)///atoms before the anomeric atom in reverse order
                    {
                        MolecularModeling::Atom* a = (*it1);
                        sorted_cycle.push_back(a);
                        sorted_cycle_stream << a->GetId() << "-";
                        a->SetIsRing(true);
                    }
                    sorted_cycle.push_back((*cycle.begin()));
                    sorted_cycle_stream << (*cycle.begin())->GetId();
                    (*cycle.begin())->SetIsRing(true);
                }
                else
                {
                    for(MolecularModeling::AtomVector::iterator it1 = cycle.begin(); it1 != it; it1++)///atoms before the anomeric atom from beginning of std::vector
                    {
                        MolecularModeling::Atom* a = (*it1);
                        sorted_cycle.push_back(a);
                        a->SetIsRing(true);
                        if(it1 == it - 1)
                          sorted_cycle_stream << a->GetId();
                        else
                          sorted_cycle_stream << a->GetId() << "-";
                    }
                }
            }
            else///anomeric is not at the end of the cycle
            {
                MolecularModeling::Atom* next_atom = cycle.at(index + 1);
                if(next_atom->GetName().substr(0,1).compare("O") == 0)///next atom is oxygen so the std::vector is in reverse order
                {
                    for(MolecularModeling::AtomVector::iterator it1 = it; it1 != cycle.begin(); it1--) ///atoms befor anomeric atom down to beginning of the std::vector
                    {
                        MolecularModeling::Atom* a_before = (*it1);
                        sorted_cycle.push_back(a_before);
                        sorted_cycle_stream << a_before->GetId() << "-";
                        a_before->SetIsRing(true);
                    }
                    sorted_cycle.push_back((*cycle.begin()));
                    sorted_cycle_stream << (*cycle.begin())->GetId() << "-";
                    (*cycle.begin())->SetIsRing(true);
                    for(MolecularModeling::AtomVector::iterator it2 = cycle.end() - 1; it2 != it; it2--)///atoms from end of the std::vector down to anomeric atom
                    {
                        MolecularModeling::Atom* atom_after = (*it2);
                        sorted_cycle.push_back(atom_after);
                        atom_after->SetIsRing(true);
                        if(it2 == it + 1)
                          sorted_cycle_stream << atom_after->GetId();
                        else
                          sorted_cycle_stream << atom_after->GetId() << "-";
                    }
                }
                else///oxygen is before the anomeric atom so the std::vector is in normal order
                {
                    for(MolecularModeling::AtomVector::iterator it1 = it; it1 != cycle.end(); it1++) ///atoms after anomeric atom to the end of the std::vector
                    {
                        MolecularModeling::Atom* atom_after = (*it1);
                        sorted_cycle.push_back(atom_after);
                        atom_after->SetIsRing(true);
                        sorted_cycle_stream << atom_after->GetId() << "-";
                    }
                    for(MolecularModeling::AtomVector::iterator it2 = cycle.begin(); it2 != it; it2++)///atoms befor the anomeric atom from beginning of std::vector
                    {
                        MolecularModeling::Atom* atom_before = (*it2);
                        sorted_cycle.push_back(atom_before);
                        atom_before->SetIsRing(true);
                        if(it2 == it - 1)
                          sorted_cycle_stream << atom_before->GetId();
                        else
                          sorted_cycle_stream << atom_before->GetId() << "-";
                    }
                }
            }
        }
    }
    return sorted_cycle;
}

// std::vector<std::string> Assembly::GetSideGroupOrientations(Glycan::Monosaccharide* mono, std::string cycle_atoms_str)
// {
//   //9/14/18 Removed side atom initialization in this function and either moved it to detectSideGroups() or used Yao's InitiateDetectionOfCompleteSideGroupAtoms()
//   // Dave
//   int local_debug = -1;
//     std::vector<std::string> orientations = std::vector<std::string>();
//     std::vector<AtomVector> side_atoms = std::vector<AtomVector>();
//     MolecularModeling::AtomVector default_atom_vector = MolecularModeling::AtomVector(3, NULL);
//
//     for(MolecularModeling::AtomVector::iterator it = mono->cycle_atoms_.begin(); it != mono->cycle_atoms_.end() - 1; it++) ///iterate on cycle atoms except the oxygen in the ring
//     {
//       orientations.push_back("N");
//       side_atoms.push_back(default_atom_vector);
//       unsigned int index = distance(mono->cycle_atoms_.begin(), it);
//       MolecularModeling::Atom* prev_atom = new MolecularModeling::Atom();
//       MolecularModeling::Atom* current_atom = (*it);
//       MolecularModeling::Atom* next_atom = new MolecularModeling::Atom();
//       if(index == 0)///if the current atom is the anomeric atom
//           prev_atom = mono->cycle_atoms_.at(mono->cycle_atoms_.size() - 1); ///previous atom is the oxygen(last atom of the sorted cycle)
//       else
//           prev_atom = mono->cycle_atoms_.at(index - 1);
//       next_atom = mono->cycle_atoms_.at(index + 1);
//
//       ///Calculating the plane based on the two ring neighbors of the current atom
//       GeometryTopology::Coordinate prev_atom_coord = GeometryTopology::Coordinate(*prev_atom->GetCoordinates().at(model_index_));
//       GeometryTopology::Coordinate current_atom_coord = GeometryTopology::Coordinate(*current_atom->GetCoordinates().at(model_index_));
//       GeometryTopology::Coordinate next_atom_coord = GeometryTopology::Coordinate(*next_atom->GetCoordinates().at(model_index_));
//       prev_atom_coord.operator -(current_atom_coord) ;
//       next_atom_coord.operator -(current_atom_coord) ;
//       GeometryTopology::Plane plane = GeometryTopology::Plane();
//       plane.SetV1(prev_atom_coord);
//       plane.SetV2(next_atom_coord);
//       GeometryTopology::Coordinate normal_v = plane.GetUnitNormalVector();
//
//       ///Calculating the orientation of the side atoms
//       MolecularModeling::AtomNode* node = current_atom->GetNode();
//       MolecularModeling::AtomVector neighbors = node->GetNodeNeighbors();
//       int not_h_neighbors = 0;
//       for(MolecularModeling::AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
//       {
//         MolecularModeling::Atom* neighbor = (*it1);
//         std::string neighbor_id = neighbor->GetId();
//         if(cycle_atoms_str.find(neighbor_id) == std::string::npos) ///if not one of the cycle atoms
//         {
//           if(neighbor->GetName().at(0) != 'H') ///deoxy check
//               not_h_neighbors++;
//         }
//       }
//       if(not_h_neighbors != 0)
//       {
//         for(MolecularModeling::AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
//         {
//           MolecularModeling::Atom* neighbor = (*it1);
//           std::string neighbor_id = neighbor->GetId();
//           if(cycle_atoms_str.find(neighbor_id) == std::string::npos) ///if not one of the cycle atoms
//           {
//             // if(neighbor->GetName().at(0) != 'H') ///deoxy check
//             //     not_h_neighbors++;
//             GeometryTopology::Coordinate side_atom_coord = GeometryTopology::Coordinate(*neighbor->GetCoordinates().at(model_index_));
//             side_atom_coord.operator -(current_atom_coord);
//             side_atom_coord.Normalize();
//             double theta = acos(normal_v.DotProduct(side_atom_coord));
//
//             if(index == 0 && neighbor_id.at(0) == 'C')///if anomeric atom has a non-ring carbon neighbor
//             {
//               if(orientations.at(index).compare("N") == 0) ///if the position of non-ring oxygen or nitrogen hasn't been set yet
//               {
//                 if(theta > (gmml::PI_RADIAN/2))
//                 {
//                   orientations.at(index) = "-1D";
//                   side_atoms.at(index).at(0) = neighbor;
//                   neighbor -> SetIsSideChain(true);
//                 }
//                 else
//                 {
//                   orientations.at(index) = "-1U";
//                   side_atoms.at(index).at(0) = neighbor;
//                   neighbor -> SetIsSideChain(true);
//                 }
//                 continue;
//               }
//               else
//               { ///position of non-ring oxygen or nitrogen + the non-ring carbon
//                 std::stringstream ss;
//                 if(theta > (gmml::PI_RADIAN/2))
//                 {
//                   ss << orientations.at(index) << "-1D";
//                   orientations.at(index) = ss.str();
//                   side_atoms.at(index).at(0) = neighbor;
//                   neighbor -> SetIsSideChain(true);
//                 }
//                 else
//                 {
//                   ss << orientations.at(index) << "-1U";
//                   orientations.at(index) = ss.str();
//                   side_atoms.at(index).at(0) = neighbor;
//                   neighbor -> SetIsSideChain(true);
//                 }
//                 break;
//               }
//             }
//             else if(neighbor_id.at(0) == 'O' || neighbor_id.at(0) == 'N')///if neighbor is a non-ring oxygen or nitrogen
//             {
//               if(index == 0)///current atom is anomeric
//               {
//                 if(orientations.at(index).compare("N") == 0) ///if the position of non-ring carbon neighbor (if exist) hasn't been set yet
//                 {
//                   if(theta > (gmml::PI_RADIAN/2))
//                   {
//                     orientations.at(index) = "D";
//                     side_atoms.at(index).at(1) = neighbor;
//                     neighbor -> SetIsSideChain(true);
//                   }
//                   else
//                   {
//                     orientations.at(index) = "U";
//                     side_atoms.at(index).at(1) = neighbor;
//                     neighbor -> SetIsSideChain(true);
//                   }
//                   continue;
//                 }
//                 else
//                 { ///position of non-ring oxygen or nitrogen + the non-ring carbon
//                   std::stringstream ss;
//                   if(theta > (gmml::PI_RADIAN/2))
//                   {
//                     ss << "D" << orientations.at(index);
//                     orientations.at(index) = ss.str();
//                     side_atoms.at(index).at(1) = neighbor;
//                     neighbor -> SetIsSideChain(true);
//                   }
//                   else
//                   {
//                     ss << "U" << orientations.at(index);
//                     orientations.at(index) = ss.str();
//                     side_atoms.at(index).at(1) = neighbor;
//                     neighbor -> SetIsSideChain(true);
//                   }
//                   break;
//                 }
//               }
//               else
//               {
//                 if(theta > (gmml::PI_RADIAN/2))
//                 {
//                   orientations.at(index) = "D";
//                   side_atoms.at(index).at(1) = neighbor;
//                   neighbor -> SetIsSideChain(true);
//                 }
//                 else
//                 {
//                   orientations.at(index) = "U";
//                   side_atoms.at(index).at(1) = neighbor;
//                   neighbor -> SetIsSideChain(true);
//                 }
//                 break;
//               }
//             }
//             else if(index == mono->cycle_atoms_.size() - 2 && neighbor_id.at(0) == 'C')///if the last ring carbon has a non-ring carbon neighbor
//             {
//               ///Check if neighbor of neighbor is oxygen or nitrogen
//               MolecularModeling::AtomNode* neighbor_node = neighbor->GetNode();
//               MolecularModeling::AtomVector neighbors_of_neighbor = neighbor_node->GetNodeNeighbors();
//               int o_neighbors = 0;
//               int not_h_neighbors = 0;
//               for(MolecularModeling::AtomVector::iterator it2 = neighbors_of_neighbor.begin(); it2 != neighbors_of_neighbor.end(); it2++)
//               {
//                 MolecularModeling::Atom* neighbor_of_neighbor = (*it2);
//                 std::string neighbor_of_neighbor_id = neighbor_of_neighbor->GetId();
//                 if((cycle_atoms_str.find(neighbor_of_neighbor_id) == std::string::npos) &&
//                   ((neighbor_of_neighbor_id.at(0) == 'O') || (neighbor_of_neighbor_id.at(0) == 'N'))) ///if neighbor of neighbor is a non-ring oxygen or nitrogen
//                 {
//                   o_neighbors++;
//                 }
//                 if(cycle_atoms_str.find(neighbor_of_neighbor_id) == std::string::npos &&
//                   (neighbor_of_neighbor_id.at(0) != 'H'))///if neighbor of neighbor is any non-ring atom other than hydrogen
//                 {
//                   not_h_neighbors++;
//                 }
//               }
//               if (o_neighbors >= 1)
//               {
//                 if(theta > (gmml::PI_RADIAN/2))
//                 {
//                   orientations.at(index) = "D";
//                   side_atoms.at(index).at(0) = neighbor;
//                   neighbor -> SetIsSideChain(true);
//                 }
//                 else
//                 {
//                   orientations.at(index) = "U";
//                   side_atoms.at(index).at(0) = neighbor;
//                   neighbor -> SetIsSideChain(true);
//                 }
//               }
//               else if(not_h_neighbors == 0)///Type Deoxy
//               {
//                 if(theta > (gmml::PI_RADIAN/2))
//                 {
//                   orientations.at(index) = "Dd";
//                   side_atoms.at(index).at(0) = neighbor;
//                   neighbor -> SetIsSideChain(true);
//                 }
//                 else
//                 {
//                   orientations.at(index) = "Ud";
//                   side_atoms.at(index).at(0) = neighbor;
//                   neighbor -> SetIsSideChain(true);
//                 }
//                 break;
//               }
//               if(orientations.at(index).compare("N") != 0)
//                 break;
//             }
//           }
//         }
//       }
//     }
//     mono->side_atoms_ = side_atoms;
//     // if ((local_debug > 0) && (side_atoms.length() > 0))
//     // {
//     //   for (std::vector<AtomVector>::iterator sideAtomVectorIT = side_atoms.begin(); sideAtomVectorIT != side_atoms.end(); sideAtomVectorIT ++)
//     //   {
//     //     MolecularModeling::AtomVector sideAtomVector = (*sideAtomVectorIT);
//     //     for (MolecularModeling::AtomVector::iterator sideAtomIT = sideAtomVector.begin(); sideAtomIT != sideAtomVector.end(); sideAtomIT++)
//     //     {
//     //       MolecularModeling::Atom* sideAtom = (*sideAtomIT);
//     //       std::stringstream debugStr;
//     //       debugStr << "Side atom: " << sideAtom->GetName();
//     //       gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
//     //     }
//     //   }
//     // }
//
//     return orientations;
// }
//
// void Assembly::InitiateDetectionOfCompleteSideGroupAtoms (std::vector<Glycan::Monosaccharide*> monos)
// {
//     for (std::vector<Glycan::Monosaccharide*>::iterator it1= monos.begin(); it1 != monos.end(); it1++){
//         Glycan::Monosaccharide* mono = *it1;
//         MolecularModeling::AtomVector& cycle_atoms = mono->cycle_atoms_;
//
//         for (std::vector<AtomVector>::iterator it2 = mono->side_atoms_.begin(); it2 != mono->side_atoms_.end(); it2++){
//             MolecularModeling::AtomVector& SideAtomArm = *it2;
//             MolecularModeling::AtomVector all_plus_one_side_atoms = MolecularModeling::AtomVector();
//             for (MolecularModeling::AtomVector::iterator it3 = SideAtomArm.begin(); it3 != SideAtomArm.end(); it3++){
//                 if ((*it3) != NULL){
//                     all_plus_one_side_atoms.push_back(*it3);
//                 }
//             }
//
//             for (MolecularModeling::AtomVector::iterator it4= all_plus_one_side_atoms.begin(); it4 != all_plus_one_side_atoms.end(); it4++ ){
//                 Atom* plus_one_side_atom = *it4;
//                 MolecularModeling::AtomVector visited_atoms = MolecularModeling::AtomVector();
//                 if ( this->CheckIfPlusOneSideAtomBelongsToCurrentMonosaccharide(SideAtomArm, cycle_atoms, plus_one_side_atom) ){
//                       this->SetCompleteSideGroupAtoms (SideAtomArm, plus_one_side_atom, cycle_atoms ,visited_atoms);
//                 }
//
//             }
//         }
//     }//for
// }//InitiateDetectionOfCompleteSideGroupAtoms
//
// bool Assembly::CheckIfPlusOneSideAtomBelongsToCurrentMonosaccharide (MolecularModeling::AtomVector& SideAtomArm, MolecularModeling::AtomVector& cycle_atoms, Atom* working_atom )
// {
//     bool belongs_to_current_monosaccharide = true;
//         MolecularModeling::AtomVector working_node_neighbors = working_atom -> GetNode() -> GetNodeNeighbors();
//         MolecularModeling::AtomVector attached_anomeric_carbons = MolecularModeling::AtomVector();
//         for (MolecularModeling::AtomVector::iterator it = working_node_neighbors.begin(); it != working_node_neighbors.end(); it++){
//             if (*it != NULL){
// 	        Atom* working_node_neighbor = *it;
// 	        if (working_node_neighbor->GetIsAnomericCarbon()){
// 	            attached_anomeric_carbons.push_back(working_node_neighbor);
// 	        }
//             }
//         }//for
//
//         if (attached_anomeric_carbons.size() == 2){
//         //TODO:figure out which branch is shorter
//         }
//
//         if (attached_anomeric_carbons.size() == 1){
//             if (std::find(cycle_atoms.begin(),cycle_atoms.end(),attached_anomeric_carbons.at(0)) != cycle_atoms.end()){ //if the anomeric carbon is from the current monosaccharide
// 	        bool sidechain_is_terminal = true;
// 	        MolecularModeling::AtomVector cycle_and_visisted_atoms = cycle_atoms;
// 	        CheckIfSideChainIsTerminal(working_atom,cycle_and_visisted_atoms,sidechain_is_terminal);
//
// 	        if (!sidechain_is_terminal ){
// 	            SideAtomArm.erase(std::find(SideAtomArm.begin(), SideAtomArm.end(), working_atom));	//Unless it's terminal, this oxygen belongs to the other monosaccharide.
// 		    belongs_to_current_monosaccharide = false;
// 	        }
//             }
//         }
//     //If number of attached anomeric carbon is zero, then this atom belongs to current side chain, and it is already. No need to do anything!
//
//     return belongs_to_current_monosaccharide;
// }
//
//
// void Assembly::SetCompleteSideGroupAtoms(MolecularModeling::AtomVector& SideAtomArm, Atom* working_atom, MolecularModeling::AtomVector& cycle_atoms, MolecularModeling::AtomVector& visited_atoms)
// {
//     MolecularModeling::AtomVector working_node_neighbors = working_atom -> GetNode() -> GetNodeNeighbors();
//     for (MolecularModeling::AtomVector::iterator it = working_node_neighbors.begin(); it != working_node_neighbors.end(); it++){  //Detect atoms that should be added to the current side chain arm.
// 	Atom* working_node_neighbor = *it;
//         if (working_node_neighbor != NULL){
// 	    if ( !working_node_neighbor->GetIsCycle() &&
// 		 std::find(visited_atoms.begin(), visited_atoms.end(), working_node_neighbor) == visited_atoms.end() &&
// 		 !working_node_neighbor->GetResidue()->CheckIfProtein() ){	//If not part of another side chain && not part of cycle && not visited && is not protein (NLN,LNK etc)
// 		AtomVector working_node_neighbor_neighbors = MolecularModeling::AtomVector();
//     		std::string working_neighbor_name = working_node_neighbor-> GetName();
//
// 	        if ( working_neighbor_name.substr(0,1) == "O" || working_neighbor_name.substr(0,1) == "N" || working_neighbor_name.substr(0,1) == "S") { //if this neighbor is an O,N or S
// 		    MolecularModeling::AtomVector attached_anomeric_carbons = MolecularModeling::AtomVector();
// 		    for (MolecularModeling::AtomVector::iterator it2 = working_node_neighbor_neighbors.begin(); it2 != working_node_neighbor_neighbors.end(); it2++){
// 		        if (*it2 != NULL ){
// 			    if ( (*it2)-> GetIsAnomericCarbon() ){
// 		                attached_anomeric_carbons.push_back(*it2); //Check how many anomeric carbons this working neighbor is attacheed to, should be either 0,1,2.
// 			    }
// 		        }
// 	            }
//
//
// 	            if (attached_anomeric_carbons.size() == 1){  // if attached to only one anomeric carbon
// 		        if (attached_anomeric_carbons.front() != working_atom){  // if the working atom IS NOT the anomeric carbon the neighbor is attached to,the current side chain gets this oxygen.
// 			    SideAtomArm.push_back(working_node_neighbor);
// 		    	    working_node_neighbor-> SetIsSideChain(true);
// 		        }
// 			// if the workin atom IS INDEED the anomeric carbon the neighbor is attached to, the current side chain doesn't get this oxygen, unless this is a terminal hydroxyl.
//
// 			/*else {
// 			    MolecularModeling::AtomVector cycle_and_visited_atoms = MolecularModeling::AtomVector();
// 			    cycle_and_visited_atoms.push_back(attached_anomeric_carbons.front() );
// 			    bool side_chain_is_terminal = true;
// 			    CheckIfSideChainIsTerminal (working_node_neighbor, cycle_and_visited_atoms, side_chain_is_terminal);
// 			    if (side_chain_is_terminal){
// 			        SideAtomArm.push_back(working_node_neighbor);
// 		    	        working_node_neighbor-> SetIsSideChain(true);
// 			    }
//
// 			}*/
//
// 		    }
//
// 		    if (attached_anomeric_carbons.size() == 0){
// 		        SideAtomArm.push_back(working_node_neighbor);
// 			working_node_neighbor-> SetIsSideChain(true);
// 		    }
//
// 		} // if is oxygen,nitrogen or sulfur
//
// 		else {  // if this neighbor is not oxygen,nitrogen or sulfur, the current side chain should get this atom
// 		    SideAtomArm.push_back(working_node_neighbor);
// 		    working_node_neighbor-> SetIsSideChain(true);
// 		}
// 	    }//if not cycle
// 	}
//     }//for
//     visited_atoms.push_back(working_atom);
//
//     for (std::vector<Atom*>::iterator it3 = working_node_neighbors.begin(); it3 != working_node_neighbors.end(); it3++){
// 	if (*it3 != NULL){
// 	    Atom* working_node_neighbor = *it3;
// 	    //if a node neighbor is not part of a ring,not visited, and not part on protein, change working atom to this neighbor, and start a new recursion call.
// 	    if ( working_node_neighbor != working_atom && !working_node_neighbor->GetIsCycle() &&
//                  std::find(visited_atoms.begin(), visited_atoms.end(), working_node_neighbor) == visited_atoms.end() &&
// 		!working_node_neighbor->GetResidue()->CheckIfProtein()) {
//                 working_atom = working_node_neighbor;
// 	        SetCompleteSideGroupAtoms(SideAtomArm, working_atom, cycle_atoms ,visited_atoms);
// 	    }//if
// 	}//if
//     }//for
//
//     return;
// }//SetCompleteSideGroupAtoms
//
// void Assembly::CheckIfSideChainIsTerminal(Atom* starting_atom, std::vector<Atom*> & cycle_and_visited_atoms, bool & is_terminal)
// {
//     cycle_and_visited_atoms.push_back(starting_atom);
//     std::vector<Atom*> starting_atom_node_neighbors = starting_atom->GetNode()->GetNodeNeighbors();
//     for (std::vector<Atom*>::iterator it = starting_atom_node_neighbors.begin(); it != starting_atom_node_neighbors.end(); it++){
// 	if (*it != NULL){
// 	    Atom* current_atom_node_neighbor = *it;
// 	    if (std::find(cycle_and_visited_atoms.begin(), cycle_and_visited_atoms.end(), current_atom_node_neighbor) == cycle_and_visited_atoms.end() ){
// 	        if (current_atom_node_neighbor->GetIsCycle() || current_atom_node_neighbor->GetResidue()->CheckIfProtein() ){
// 		    is_terminal = false;
// 	        }
// 	    }
// 	}
//     }//for
//
//     for (std::vector<Atom*>::iterator it = starting_atom_node_neighbors.begin(); it != starting_atom_node_neighbors.end(); it++){
// 	if ( *it != NULL){
// 	    if(std::find(cycle_and_visited_atoms.begin(), cycle_and_visited_atoms.end(), *it) == cycle_and_visited_atoms.end() &&
// 	       !(*it)->GetResidue()->CheckIfProtein()  && is_terminal ){
// 	        starting_atom = *it;
// 	        CheckIfSideChainIsTerminal(starting_atom,cycle_and_visited_atoms,is_terminal);
// 	    }//if
// 	}
//     }//for
//
//     return;
// } //CheckIfSideChainIsTerminal

// Glycan::ChemicalCode* Assembly::BuildChemicalCode(std::vector<std::string> orientations)
// {
//     Glycan::ChemicalCode* code = new Glycan::ChemicalCode();
//     if(orientations.size() == 5 )
//         code->base_ = "P";
//     else if(orientations.size() == 4 )
//         code->base_ = "F";
//     else
//         code->base_ = "?";
//
//     ///Side atom(s) of anomeric
//     ///Has only non-ring oxygen neighbor
//     if(orientations.at(0).compare("U") == 0)
//         code->right_up_.push_back("a");
//     else if(orientations.at(0).compare("D") == 0)
//         code->right_down_.push_back("a");
//
//     ///Has non-ring oxygen and carbon neighbors
//     else if(orientations.at(0).compare("U-1U") == 0)
//     {
//         code->right_up_.push_back("a");
//         code->right_up_.push_back("-1");
//     }
//     else if(orientations.at(0).compare("D-1U") == 0)
//     {
//         code->right_down_.push_back("a");
//         code->right_up_.push_back("-1");
//     }
//     else if(orientations.at(0).compare("U-1D") == 0)
//     {
//         code->right_up_.push_back("a");
//         code->right_down_.push_back("-1");
//     }
//     else if(orientations.at(0).compare("D-1D") == 0)
//     {
//         code->right_down_.push_back("a");
//         code->right_down_.push_back("-1");
//     }
//
//     ///Has only non-ring carbon neighbor
//     else if(orientations.at(0).compare("-1U") == 0)
//         code->right_up_.push_back("-1");
//     else if(orientations.at(0).compare("-1D") == 0)
//         code->right_down_.push_back("-1");
//
//     ///Side atom of other carbons of the ring
//     for(std::vector<std::string>::iterator it = orientations.begin() + 1; it != orientations.end() - 1; it++)
//     {
//         std::string orientation = (*it);
//         int index = distance(orientations.begin(), it);
//         if(orientation.compare("U") == 0)
//             code->left_up_.push_back(gmml::ConvertT(index + 1));
//         else if(orientation.compare("D") == 0)
//             code->left_down_.push_back(gmml::ConvertT(index + 1));
//         else if(orientation.compare("N") == 0)
//         {
//             std::stringstream ss;
//             ss << gmml::ConvertT(index + 1) << "d";
//             code->left_middle_.push_back(ss.str() );
//         }
//     }
//
//     ///Side atom(s) of last carbon
//     if(orientations.at(orientations.size() - 1).compare("U") == 0)
//         code->right_up_.push_back("+1");
//     else if(orientations.at(orientations.size() - 1).compare("D") == 0)
//         code->right_down_.push_back("+1");
//     ///Type Deoxy
//     else if(orientations.at(orientations.size() - 1).compare("Ud") == 0)
//         code->right_up_.push_back("+1d");
//     else if(orientations.at(orientations.size() - 1).compare("Dd") == 0)
//         code->right_down_.push_back("+1d");
//
//     return code;
// }

// MolecularModeling::AtomVector Assembly::ExtractAdditionalSideAtoms(Glycan::Monosaccharide *mono)
// {
//     MolecularModeling::AtomVector plus_sides = MolecularModeling::AtomVector();
//     if(mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(0) != NULL)///if there exist a +1 carbon atom. in side_atoms_ structure (std::vector<AtomVector>) the first index of the last element is dedicated to +1 atom
//     {
//         plus_sides.push_back(mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(0));
//         MolecularModeling::AtomVector plus_one_atom_neighbors = mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(0)->GetNode()->GetNodeNeighbors();
//         for(MolecularModeling::AtomVector::iterator it1 = plus_one_atom_neighbors.begin(); it1 != plus_one_atom_neighbors.end(); it1++)
//         {
//             if((*it1)->GetName().at(0) == 'C' && mono->cycle_atoms_str_.find((*it1)->GetId()) == std::string::npos)///+2 carbon atom found
//             {
//                 MolecularModeling::Atom* plus_two = (*it1);
//                 plus_sides.push_back(plus_two);
//                 mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(1) = plus_two;///in side_atoms_ structure (std::vector<AtomVector>) the second index of the last element is dedicated to +2 atom
//
//                 MolecularModeling::AtomVector plus_two_atom_neighbors = plus_two->GetNode()->GetNodeNeighbors();
//                 for(MolecularModeling::AtomVector::iterator it2 = plus_two_atom_neighbors.begin(); it2 != plus_two_atom_neighbors.end(); it2++)
//                 {
//                     MolecularModeling::Atom* plus_three = (*it2);
//                     if(plus_three->GetName().at(0) == 'C' && plus_three->GetId().compare(mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(0)->GetId()) != 0)///+3 carbon atom found
//                     {
//                         plus_sides.push_back(plus_three);
//                         mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(2) = plus_three;///in side_atoms_ structure (std::vector<AtomVector>) the third index of the last element is dedicated to +3 atom
//                         break;
//                     }
//                 }
//                 break;
//             }
//         }
//     }
//     return plus_sides;
// }

// void Assembly::ExtractDerivatives(Glycan::Monosaccharide * mono)
// {
//     int local_debug = -1;
//     if (local_debug > 0)
//     {
//       std::stringstream debugStr;
//       debugStr << "On monosaccharide: " << mono->mono_id_;
//       gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
//     }
//     for(MolecularModeling::AtomVector::iterator it = mono->cycle_atoms_.begin(); it != mono->cycle_atoms_.end() - 1; it++) ///iterate on cycle atoms except the oxygen of the ring
//     {
//         int index = distance(mono->cycle_atoms_.begin(), it);
//         MolecularModeling::Atom* target = (*it);
//         std::string key = "";
//         std::string value = "";
//         MolecularModeling::AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
//         if (local_debug > 0)
//         {
//           std::stringstream debugStr;
//           debugStr << "On atom: " << target->GetName();
//           gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
//         }
//         for(MolecularModeling::AtomVector::iterator it1 = t_neighbors.begin(); it1 != t_neighbors.end(); it1++)
//         {
//
//             MolecularModeling::AtomVector pattern_atoms = MolecularModeling::AtomVector();
//             MolecularModeling::Atom* t_neighbor = (*it1);
//             if(local_debug > 0)
//             {
//               std::stringstream debugStr;
//               debugStr << "On neighbor: " << t_neighbor->GetName();
//               gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
//             }
//             if(t_neighbor->GetName().at(0) == 'N' && mono->cycle_atoms_str_.find(t_neighbor->GetId()) == std::string::npos)///check formulas with nitrogen
//             {
//                 if((value = CheckxC_N(target, mono->cycle_atoms_str_/*, pattern_atoms*/)).compare("") != 0)///xCH-N
//                     break;
//                 if((value = CheckxC_NxO_CO_C(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH3
//                     break;
//                 if((value = CheckxC_NxO_CO_CO(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH2OH
//                     break;
//                 if((value = CheckxC_NxO_SO3(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-SO3
//                     break;
//                 if((value = CheckxC_NxO_PO3(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-PO3
//                     break;
//                 if((value = CheckxC_NxO_C(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-CH3
//                     break;
//             }
//             if(t_neighbor->GetName().at(0) == 'O' && mono->cycle_atoms_str_.find(t_neighbor->GetId()) == std::string::npos)///check formulas with oxygen
//             {
//                 if((value = CheckxC_NxO_CO_C(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH3
//                     break;
//                 if((value = CheckxC_NxO_CO_CO(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH2OH
//                     break;
//                 if((value = CheckxC_NxO_SO3(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-SO3
//                     break;
//                 if((value = CheckxC_NxO_PO3(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-PO3
//                     break;
//                 if((value = CheckxC_NxO_C(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-CH3
//                     break;
//                 if((value = CheckxCOO(target, mono->cycle_atoms_str_/*, pattern_atoms*/)).compare("") != 0)///xC-(O,O) and xC-(O,OH)
//                     break;
//             }
//         }
//         if(value.compare("") != 0)///if any pattern matched add it to the index-derivative map
//         {
//             if(index == 0)
//                 key = "a";
//             else
//                 key = gmml::ConvertT(index + 1);
//             mono->derivatives_map_[key] = value;
//         }
//     }
//     for(std::vector<AtomVector>::iterator it = mono->side_atoms_.begin(); it != mono->side_atoms_.end(); it++) ///iterate on side atoms
//     {
//         int index = distance(mono->side_atoms_.begin(), it);
//         int side_branch_last_carbon_index = 0;
//         MolecularModeling::AtomVector sides = (*it);
//         MolecularModeling::Atom* target = NULL;
//         std::string key = "";
//         std::string value = "";
//         // if(it == mono->side_atoms_.begin())///side atoms of anomeric carbon
//         // {
//         //     if(sides.at(0) != NULL)
//         //         target = sides.at(0);///first index of each side is for carbon atoms in the std::vector<AtomVector> structure
//         // }
//         // //WHERE ARE ALL OF THE OTHER SIDE ATOMS?!?!?!?!?!?!?!?!?!??!
//         // if(it == mono->side_atoms_.end() - 1)//side atoms of last carbon of the ring
//         // {
//         //     if(sides.at(0) != NULL)
//         //     {
//         //         for(side_branch_last_carbon_index = sides.size() - 1; sides.at(side_branch_last_carbon_index) == NULL; side_branch_last_carbon_index-- )
//         //         {
//         //           target = sides.at(side_branch_last_carbon_index);
//         //         }
//         //
//         //     }
//         // }
//         // else
//         // {
//           if(sides.at(0) != NULL)
//               target = sides.at(0);
//         // }
//         if ((local_debug > 0) && (target != NULL))
//         {
//           std::stringstream debugStr;
//           debugStr << "On side atom: " << target->GetName();
//           gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
//         }
//         if(target != NULL)
//         {
//             MolecularModeling::AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
//             for(MolecularModeling::AtomVector::iterator it1 = t_neighbors.begin(); it1 != t_neighbors.end(); it1++)
//             {
//                 MolecularModeling::AtomVector pattern_atoms = MolecularModeling::AtomVector();
//                 MolecularModeling::Atom* t_neighbor = (*it1);
//                 if (local_debug > 0)
//                 {
//                   std::stringstream debugStr;
//                   debugStr << "On side atom neighbor: " << t_neighbor->GetName();
//                   gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
//                 }
//                 if(t_neighbor->GetName().at(0) == 'N' && mono->cycle_atoms_str_.find(t_neighbor->GetId()) == std::string::npos)///check formulas with nitrogen
//                 {
//                     if((value = CheckxC_N(target, mono->cycle_atoms_str_/*, pattern_atoms*/)).compare("") != 0)///xCH-N
//                         break;
//                     if((value = CheckxC_NxO_CO_C(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH3
//                         break;
//                     if((value = CheckxC_NxO_CO_CO(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH2OH
//                         break;
//                     if((value = CheckxC_NxO_SO3(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-SO3
//                         break;
//                     if((value = CheckxC_NxO_PO3(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-PO3
//                         break;
//                     if((value = CheckxC_NxO_C(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-CH3
//                         break;
//                 }
//                 if(t_neighbor->GetName().at(0) == 'O' && mono->cycle_atoms_str_.find(t_neighbor->GetId()) == std::string::npos)///check formulas with oxygen
//                 {
//                     if((value = CheckxC_NxO_CO_C(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH3
//                         break;
//                     if((value = CheckxC_NxO_CO_CO(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH2OH
//                         break;
//                     if((value = CheckxC_NxO_SO3(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-SO3
//                         break;
//                     if((value = CheckxC_NxO_PO3(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-PO3
//                         break;
//                     if((value = CheckxC_NxO_C(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-CH3
//                         break;
//                     if((value = CheckxCOO(target, mono->cycle_atoms_str_/*, pattern_atoms*/)).compare("") != 0)///xC-(O,O) and xC-(O,OH)
//                         break;
//                 }
//                 gmml::log(__LINE__, __FILE__, gmml::INF, "Done checking formulas");
//             }
//             if(value.compare("") != 0)///if any pattern matched add it to the index-derivative map
//             {
//               gmml::log(__LINE__, __FILE__, gmml::INF, "Makes it here");
//                 if(index == 0)
//                     key = "-1";
//                 else
//                 {
//                     switch (side_branch_last_carbon_index)
//                     {
//                         case 0:
//                             key = "+1";
//                             break;
//                         case 1:
//                             key = "+2";
//                             break;
//                         case 2:
//                             key = "+3";
//                             break;
//                     }
//
//                 }
//                 gmml::log(__LINE__, __FILE__, gmml::INF, "And then here");
//                 mono->derivatives_map_[key] = value;
//             }
//         }
//     }
// }

std::string Assembly::CheckxC_N(MolecularModeling::Atom* target, std::string cycle_atoms_str/*, MolecularModeling::AtomVector& pattern_atoms*/)
{
    int local_debug = -1;
    std::stringstream pattern;
    pattern << "xC";
    MolecularModeling::Atom* N = NULL;
    MolecularModeling::AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(MolecularModeling::AtomVector::iterator it1 = t_neighbors.begin(); it1 != t_neighbors.end(); it1++)
    {
        MolecularModeling::Atom* t_neighbor = (*it1);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == std::string::npos)
        {
            if(t_neighbor->GetName().at(0) != 'N')
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == 'N' && N != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == 'N' && N == NULL)
                N = t_neighbor;
        }
    }
    if(N != NULL)
    {
        pattern << "-N";
        MolecularModeling::AtomVector n_neighbors = N->GetNode()->GetNodeNeighbors();
        for(MolecularModeling::AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            MolecularModeling::Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0)
                pattern << n_neighbor->GetName().at(0);
        }
    }
    if (local_debug > 0)
    {
      std::stringstream debugStr;
      debugStr << "CheckxC_N:" << pattern.str() << "\n";
      gmml::log(__LINE__, __FILE__,  gmml::INF, debugStr.str());
    }
    if(pattern.str().compare("xCH-NHH") == 0 || pattern.str().compare("xC-N") == 0 || pattern.str().compare("xCHH-NHH") == 0 || pattern.str().compare("xCH-N") == 0  ||
            pattern.str().compare("xCHH-N") == 0 || pattern.str().compare("xC-NHH") == 0)
        return "xCH-N";
    else
        return "";
}

std::string Assembly::CheckxC_NxO_CO_C(MolecularModeling::Atom *target, std::string cycle_atoms_str, char NxO, MolecularModeling::AtomVector& pattern_atoms)
{
    int local_debug = -1;
    std::stringstream pattern;
    pattern << "xC";
    MolecularModeling::Atom* N_or_O = NULL;
    MolecularModeling::AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(MolecularModeling::AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        MolecularModeling::Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == std::string::npos)
        {
            if(t_neighbor->GetName().at(0) != NxO)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
                N_or_O = t_neighbor;
        }
    }
    if(N_or_O != NULL)
    {
        MolecularModeling::Atom* C = NULL;
        pattern << "-" << NxO;
        MolecularModeling::AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
        for(MolecularModeling::AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            MolecularModeling::Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'C')
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C != NULL)
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C == NULL)
                C = n_neighbor;
        }
        if(C != NULL)
        {
            pattern_atoms.push_back(C);
            MolecularModeling::Atom* CC = NULL;
            MolecularModeling::Atom* CO = NULL;
            pattern << "-C";
            MolecularModeling::AtomVector c_neighbors = C->GetNode()->GetNodeNeighbors();
            for(MolecularModeling::AtomVector::iterator it = c_neighbors.begin(); it != c_neighbors.end(); it++)
            {
                MolecularModeling::Atom* c_neighbor = (*it);
                if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 && c_neighbor->GetName().at(0) != 'C' && c_neighbor->GetName().at(0) != 'O')
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'C' && CC != NULL)
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'O' && CO != NULL)
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'O' && CO == NULL)
                    CO = c_neighbor;
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'C' && CC == NULL)
                    CC = c_neighbor;
            }
            if(CO != NULL)
            {
                pattern_atoms.push_back(CO);
                pattern << "O";
                MolecularModeling::AtomVector co_neighbors = CO->GetNode()->GetNodeNeighbors();
                for(MolecularModeling::AtomVector::iterator it = co_neighbors.begin(); it != co_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* co_neighbor = (*it);
                    if(co_neighbor->GetId().compare(C->GetId()) != 0)
                        pattern << co_neighbor->GetName().at(0);
                }
            }
            if(CC != NULL)
            {
                pattern_atoms.push_back(CC);
                pattern << "-C";
                MolecularModeling::AtomVector cc_neighbors = CC->GetNode()->GetNodeNeighbors();
                for(MolecularModeling::AtomVector::iterator it = cc_neighbors.begin(); it != cc_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* cc_neighbor = (*it);
                    if(cc_neighbor->GetId().compare(C->GetId()) != 0)
                        pattern << cc_neighbor->GetName().at(0);
                }
            }
        }
    }
    if (local_debug > 0)
    {
      std::stringstream debugStr;
      debugStr << "CheckxC_NxO_CO_C:" << pattern.str() << "\n";
      gmml::log(__LINE__, __FILE__,  gmml::INF, debugStr.str());
    }
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-CO-CHHH") == 0 || pattern.str().compare("xC-N-CO-C") == 0 || pattern.str().compare("xCHH-NH-CO-CHHH") == 0 || pattern.str().compare("xC-NH-CO-CHHH") == 0 ||
                pattern.str().compare("xC-N-CO-CHHH") == 0 || pattern.str().compare("xC-NH-CO-C") == 0 || pattern.str().compare("xCHH-N-CO-CHHH") == 0 || pattern.str().compare("xCHH-NH-CO-C") == 0 ||
                pattern.str().compare("xCHH-N-CO-C") == 0 || pattern.str().compare("xCH-N-CO-CHHH") == 0 || pattern.str().compare("xCH-NH-CO-C") == 0 || pattern.str().compare("xCH-N-CO-C") == 0 )
            return "xC-N-C=OCH3";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-OH-CO-CHHH") == 0 || pattern.str().compare("xC-O-CO-C") == 0 || pattern.str().compare("xCHH-OH-CO-CHHH") == 0 || pattern.str().compare("xC-OH-CO-CHHH") == 0 ||
                pattern.str().compare("xC-O-CO-CHHH") == 0 || pattern.str().compare("xC-OH-CO-C") == 0 || pattern.str().compare("xCHH-O-CO-CHHH") == 0 || pattern.str().compare("xCHH-OH-CO-C") == 0 ||
                pattern.str().compare("xCHH-O-CO-C") == 0 || pattern.str().compare("xCH-O-CO-CHHH") == 0 || pattern.str().compare("xCH-OH-CO-C") == 0 || pattern.str().compare("xCH-O-CO-C") == 0)
            return "xC-O-C=OCH3";
        else
            return "";
    }
    else
        return "";
}

std::string Assembly::CheckxC_NxO_CO_CO(MolecularModeling::Atom *target, std::string cycle_atoms_str, char NxO, MolecularModeling::AtomVector& pattern_atoms)
{
    int local_debug = -1;
    std::stringstream pattern;
    pattern << "xC";
    MolecularModeling::Atom* N_or_O = NULL;
    MolecularModeling::AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(MolecularModeling::AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        MolecularModeling::Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == std::string::npos)
        {
            if(t_neighbor->GetName().at(0) != NxO)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
                N_or_O = t_neighbor;
        }
    }
    if(N_or_O != NULL)
    {
        MolecularModeling::Atom* C = NULL;
        pattern << "-" << NxO;
        MolecularModeling::AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
        for(MolecularModeling::AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            MolecularModeling::Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'C')
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C != NULL)
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C == NULL)
                C = n_neighbor;
        }
        if(C != NULL)
        {
            pattern_atoms.push_back(C);
            MolecularModeling::Atom* CC = NULL;
            MolecularModeling::Atom* CO = NULL;
            pattern << "-C";
            MolecularModeling::AtomVector c_neighbors = C->GetNode()->GetNodeNeighbors();
            for(MolecularModeling::AtomVector::iterator it = c_neighbors.begin(); it != c_neighbors.end(); it++)
            {
                MolecularModeling::Atom* c_neighbor = (*it);
                if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 && c_neighbor->GetName().at(0) != 'C' && c_neighbor->GetName().at(0) != 'O')
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'C' && CC != NULL)
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'O' && CO != NULL)
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'O' && CO == NULL)
                    CO = c_neighbor;
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'C' && CC == NULL)
                    CC = c_neighbor;
            }
            if(CO != NULL)
            {
                pattern_atoms.push_back(CO);
                pattern << "O";
                MolecularModeling::AtomVector co_neighbors = CO->GetNode()->GetNodeNeighbors();
                for(MolecularModeling::AtomVector::iterator it = co_neighbors.begin(); it != co_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* co_neighbor = (*it);
                    if(co_neighbor->GetId().compare(C->GetId()) != 0)
                        pattern << co_neighbor->GetName().at(0);
                }
            }
            if(CC != NULL)
            {
                pattern_atoms.push_back(CC);
                MolecularModeling::Atom* O = NULL;
                pattern << "-C";
                MolecularModeling::AtomVector cc_neighbors = CC->GetNode()->GetNodeNeighbors();
                for(MolecularModeling::AtomVector::iterator it = cc_neighbors.begin(); it != cc_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* cc_neighbor = (*it);
                    if(cc_neighbor->GetId().compare(C->GetId()) != 0 && cc_neighbor->GetName().at(0) != 'O')
                        pattern << cc_neighbor->GetName().at(0);
                    else if(cc_neighbor->GetId().compare(C->GetId()) != 0 && cc_neighbor->GetName().at(0) == 'O' && O != NULL)
                        pattern << cc_neighbor->GetName().at(0);
                    else if(cc_neighbor->GetId().compare(C->GetId()) != 0 && cc_neighbor->GetName().at(0) == 'O' && O == NULL)
                        O = cc_neighbor;
                }
                if(O != NULL)
                {
                    pattern_atoms.push_back(O);
                    pattern << "-O";
                    MolecularModeling::AtomVector o_neighbors = O->GetNode()->GetNodeNeighbors();
                    for(MolecularModeling::AtomVector::iterator it = o_neighbors.begin(); it != o_neighbors.end(); it++)
                    {
                        MolecularModeling::Atom* o_neighbor = (*it);
                        if(o_neighbor->GetId().compare(CC->GetId()) != 0)
                            pattern << o_neighbor->GetName().at(0);
                    }
                }
            }
        }
    }
    if (local_debug > 0)
    {
      std::stringstream debugStr;
      debugStr << "CheckxC_NxO_CO_CO:" << pattern.str() << "\n";
      gmml::log(__LINE__, __FILE__,  gmml::INF, debugStr.str());
    }
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-CO-CHH-OH") == 0 || pattern.str().compare("xCH-N-CO-CHH-OH") == 0 || pattern.str().compare("xCH-NH-CO-C-OH") == 0 ||
                pattern.str().compare("xCH-NH-CO-CHH-O") == 0 || pattern.str().compare("xCH-N-CO-C-OH") == 0 || pattern.str().compare("xCH-N-CO-CHH-O") == 0 ||
                pattern.str().compare("xCH-NH-CO-C-OH") == 0 || pattern.str().compare("xCH-N-CO-C-O") == 0 || pattern.str().compare("xC-N-CO-C-O") == 0 ||
                pattern.str().compare("xC-NH-CO-C-O") == 0 || pattern.str().compare("xC-NH-CO-C-OH") == 0 || pattern.str().compare("xC-NH-CO-CHH-O") == 0 ||
                pattern.str().compare("xC-N-CO-CHH-O") == 0 || pattern.str().compare("xC-N-CO-C-OH") == 0 || pattern.str().compare("xC-NH-CO-CHH-OH") == 0 ||
                pattern.str().compare("xC-N-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-NH-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-N-CO-C-O") == 0 ||
                pattern.str().compare("xCHH-N-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-N-CO-CHH-O") == 0 || pattern.str().compare("xCHH-N-CO-C-OH") == 0 ||
                pattern.str().compare("xCHH-NH-CO-C-OH") == 0 || pattern.str().compare("xCHH-NH-CO-CHH-O") == 0 || pattern.str().compare("xCHH-NH-CO-C-O") == 0 )
            return "xC-N-C=OCH2OH";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-OH-CO-CHH-OH") == 0 || pattern.str().compare("xCH-O-CO-CHH-OH") == 0 || pattern.str().compare("xCH-OH-CO-C-OH") == 0 ||
                pattern.str().compare("xCH-OH-CO-CHH-O") == 0 || pattern.str().compare("xCH-O-CO-C-OH") == 0 || pattern.str().compare("xCH-O-CO-CHH-O") == 0 ||
                pattern.str().compare("xCH-OH-CO-C-OH") == 0 || pattern.str().compare("xCH-O-CO-C-O") == 0 || pattern.str().compare("xC-O-CO-C-O") == 0 ||
                pattern.str().compare("xC-OH-CO-C-O") == 0 || pattern.str().compare("xC-OH-CO-C-OH") == 0 || pattern.str().compare("xC-OH-CO-CHH-O") == 0 ||
                pattern.str().compare("xC-O-CO-CHH-O") == 0 || pattern.str().compare("xC-O-CO-C-OH") == 0 || pattern.str().compare("xC-OH-CO-CHH-OH") == 0 ||
                pattern.str().compare("xC-O-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-OH-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-O-CO-C-O") == 0 ||
                pattern.str().compare("xCHH-O-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-O-CO-CHH-O") == 0 || pattern.str().compare("xCHH-O-CO-C-OH") == 0 ||
                pattern.str().compare("xCHH-OH-CO-C-OH") == 0 || pattern.str().compare("xCHH-OH-CO-CHH-O") == 0 || pattern.str().compare("xCHH-OH-CO-C-O") == 0 )

            return "xC-O-C=OCH2OH";
        else
            return "";
    }
    else
        return "";
}

std::string Assembly::CheckxC_NxO_SO3(MolecularModeling::Atom *target, std::string cycle_atoms_str, char NxO, MolecularModeling::AtomVector& pattern_atoms)
{
    int local_debug = -1;
    std::stringstream pattern;
    pattern << "xC";
    MolecularModeling::Atom* N_or_O = NULL;
    MolecularModeling::AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(MolecularModeling::AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        MolecularModeling::Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == std::string::npos)
        {
            if(t_neighbor->GetName().at(0) != NxO)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
                N_or_O = t_neighbor;
        }
    }
    if(N_or_O != NULL)
    {
        MolecularModeling::Atom* S = NULL;
        pattern << "-" << NxO;
        MolecularModeling::AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
        for(MolecularModeling::AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            MolecularModeling::Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'S')
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'S' && S != NULL)
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'S' && S == NULL)
                S = n_neighbor;
        }
        if(S != NULL)
        {
            pattern_atoms.push_back(S);
            MolecularModeling::Atom* O1 = NULL;
            MolecularModeling::Atom* O2 = NULL;
            MolecularModeling::Atom* O3 = NULL;
            pattern << "-S";
            MolecularModeling::AtomVector s_neighbors = S->GetNode()->GetNodeNeighbors();
            for(MolecularModeling::AtomVector::iterator it = s_neighbors.begin(); it != s_neighbors.end(); it++)
            {
                MolecularModeling::Atom* s_neighbor = (*it);
                if(s_neighbor->GetId().compare(N_or_O->GetId()) != 0 && s_neighbor->GetName().at(0) != 'O')
                    pattern << s_neighbor->GetName().at(0);
                else if(s_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  s_neighbor->GetName().at(0) == 'O' && O1 == NULL)
                    O1 = s_neighbor;
                else if(s_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  s_neighbor->GetName().at(0) == 'O' && O2 == NULL)
                    O2 = s_neighbor;
                else if(s_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  s_neighbor->GetName().at(0) == 'O' && O3 == NULL)
                    O3 = s_neighbor;
            }
            if(O1 != NULL && O2 != NULL && O3 != NULL)
            {
                pattern_atoms.push_back(O1);
                pattern_atoms.push_back(O2);
                pattern_atoms.push_back(O3);
                pattern << "OOO";
                MolecularModeling::AtomVector o1_neighbors = O1->GetNode()->GetNodeNeighbors();
                for(MolecularModeling::AtomVector::iterator it = o1_neighbors.begin(); it != o1_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* o1_neighbor = (*it);
                    if(o1_neighbor->GetId().compare(S->GetId()) != 0)
                        pattern << o1_neighbor->GetName().at(0);
                }
                MolecularModeling::AtomVector o2_neighbors = O2->GetNode()->GetNodeNeighbors();
                for(MolecularModeling::AtomVector::iterator it = o2_neighbors.begin(); it != o2_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* o2_neighbor = (*it);
                    if(o2_neighbor->GetId().compare(S->GetId()) != 0)
                        pattern << o2_neighbor->GetName().at(0);
                }
                MolecularModeling::AtomVector o3_neighbors = O3->GetNode()->GetNodeNeighbors();
                for(MolecularModeling::AtomVector::iterator it = o3_neighbors.begin(); it != o3_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* o3_neighbor = (*it);
                    if(o3_neighbor->GetId().compare(S->GetId()) != 0)
                        pattern << o3_neighbor->GetName().at(0);
                }
            }
        }
    }
    if (local_debug > 0)
    {
      std::stringstream debugStr;
      debugStr << "CheckxC_NxO_SO3:" << pattern.str() << "\n";
      gmml::log(__LINE__, __FILE__,  gmml::INF, debugStr.str());
    }
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-SOOOH") == 0 || pattern.str().compare("xCH-N-SOOOH") == 0 || pattern.str().compare("xCH-NH-SOOO") == 0 || pattern.str().compare("xCH-N-SOOO") == 0 ||
                pattern.str().compare("xCHH-NH-SOOOH") == 0 || pattern.str().compare("xCHH-NH-SOOO") == 0 || pattern.str().compare("xCHH-N-SOOOH") == 0 || pattern.str().compare("xCHH-N-SOOO") == 0 ||
                pattern.str().compare("xC-N-SOOO") == 0 || pattern.str().compare("xC-NH-SOOOH") == 0 || pattern.str().compare("xC-N-SOOOH") == 0 || pattern.str().compare("xC-NH-SOOO") == 0 )
            return "xC-N-SO3";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-OH-SOOOH") == 0 || pattern.str().compare("xCH-O-SOOOH") == 0 || pattern.str().compare("xCH-OH-SOOO") == 0 || pattern.str().compare("xCH-O-SOOO") == 0 ||
                pattern.str().compare("xCHH-OH-SOOOH") == 0 || pattern.str().compare("xCHH-OH-SOOO") == 0 || pattern.str().compare("xCHH-O-SOOOH") == 0 || pattern.str().compare("xCHH-O-SOOO") == 0 ||
                pattern.str().compare("xC-O-SOOO") == 0 || pattern.str().compare("xC-OH-SOOOH") == 0 || pattern.str().compare("xC-O-SOOOH") == 0 || pattern.str().compare("xC-OH-SOOO") == 0 )
            return "xC-O-SO3";
        else
            return "";
    }
    else
        return "";
}

std::string Assembly::CheckxC_NxO_PO3(MolecularModeling::Atom *target, std::string cycle_atoms_str, char NxO, MolecularModeling::AtomVector& pattern_atoms)
{
    int local_debug = -1;
    std::stringstream pattern;
    pattern << "xC";
    MolecularModeling::Atom* N_or_O = NULL;
    MolecularModeling::AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(MolecularModeling::AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        MolecularModeling::Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == std::string::npos)
        {
            if(t_neighbor->GetName().at(0) != NxO)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
                N_or_O = t_neighbor;
        }
    }
    if(N_or_O != NULL)
    {
        MolecularModeling::Atom* P = NULL;
        pattern << "-" << NxO;
        MolecularModeling::AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
        for(MolecularModeling::AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            MolecularModeling::Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'P')
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'P' && P != NULL)
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'P' && P == NULL)
                P = n_neighbor;
        }
        if(P != NULL)
        {
            pattern_atoms.push_back(P);
            MolecularModeling::Atom* O1 = NULL;
            MolecularModeling::Atom* O2 = NULL;
            MolecularModeling::Atom* O3 = NULL;
            pattern << "-P";
            MolecularModeling::AtomVector p_neighbors = P->GetNode()->GetNodeNeighbors();
            for(MolecularModeling::AtomVector::iterator it = p_neighbors.begin(); it != p_neighbors.end(); it++)
            {
                MolecularModeling::Atom* p_neighbor = (*it);
                if(p_neighbor->GetId().compare(N_or_O->GetId()) != 0 && p_neighbor->GetName().at(0) != 'O')
                    pattern << p_neighbor->GetName().at(0);
                else if(p_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  p_neighbor->GetName().at(0) == 'O' && O1 == NULL)
                    O1 = p_neighbor;
                else if(p_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  p_neighbor->GetName().at(0) == 'O' && O2 == NULL)
                    O2 = p_neighbor;
                else if(p_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  p_neighbor->GetName().at(0) == 'O' && O3 == NULL)
                    O3 = p_neighbor;
            }
            if(O1 != NULL && O2 != NULL && O3 != NULL)
            {
                pattern_atoms.push_back(O1);
                pattern_atoms.push_back(O2);
                pattern_atoms.push_back(O3);
                pattern << "OOO";
                MolecularModeling::AtomVector o1_neighbors = O1->GetNode()->GetNodeNeighbors();
                for(MolecularModeling::AtomVector::iterator it = o1_neighbors.begin(); it != o1_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* o1_neighbor = (*it);
                    if(o1_neighbor->GetId().compare(P->GetId()) != 0)
                        pattern << o1_neighbor->GetName().at(0);
                }
                MolecularModeling::AtomVector o2_neighbors = O2->GetNode()->GetNodeNeighbors();
                for(MolecularModeling::AtomVector::iterator it = o2_neighbors.begin(); it != o2_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* o2_neighbor = (*it);
                    if(o2_neighbor->GetId().compare(P->GetId()) != 0)
                        pattern << o2_neighbor->GetName().at(0);
                }
                MolecularModeling::AtomVector o3_neighbors = O3->GetNode()->GetNodeNeighbors();
                for(MolecularModeling::AtomVector::iterator it = o3_neighbors.begin(); it != o3_neighbors.end(); it++)
                {
                    MolecularModeling::Atom* o3_neighbor = (*it);
                    if(o3_neighbor->GetId().compare(P->GetId()) != 0)
                        pattern << o3_neighbor->GetName().at(0);
                }
            }
        }
    }
    if (local_debug > 0)
    {
      std::stringstream debugStr;
      debugStr << "CheckxC_NxO_PO3:" << pattern.str() << "\n";
      gmml::log(__LINE__, __FILE__,  gmml::INF, debugStr.str());
    }
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-POOOH") == 0 || pattern.str().compare("xCH-N-POOOH") == 0 || pattern.str().compare("xCH-NH-POOO") == 0 || pattern.str().compare("xCH-N-POOO") == 0 ||
                pattern.str().compare("xCHH-NH-POOOH") == 0 || pattern.str().compare("xCHH-NH-POOO") == 0 || pattern.str().compare("xCHH-N-POOOH") == 0 || pattern.str().compare("xCHH-N-POOO") == 0 ||
                pattern.str().compare("xC-N-POOO") == 0 || pattern.str().compare("xC-NH-POOOH") == 0 || pattern.str().compare("xC-N-POOOH") == 0 || pattern.str().compare("xC-NH-POOO") == 0 )
            return "xC-N-PO3";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-OH-POOOH") == 0 || pattern.str().compare("xCH-O-POOOH") == 0 || pattern.str().compare("xCH-OH-POOO") == 0 || pattern.str().compare("xCH-O-POOO") == 0 ||
                pattern.str().compare("xCHH-OH-POOOH") == 0 || pattern.str().compare("xCHH-OH-POOO") == 0 || pattern.str().compare("xCHH-O-POOOH") == 0 || pattern.str().compare("xCHH-O-POOO") == 0 ||
                pattern.str().compare("xC-O-POOO") == 0 || pattern.str().compare("xC-OH-POOOH") == 0 || pattern.str().compare("xC-O-POOOH") == 0 || pattern.str().compare("xC-OH-POOO") == 0 )            return "xC-O-PO3";
        else
            return "";
    }
    else
        return "";
}

std::string Assembly::CheckxC_NxO_C(MolecularModeling::Atom *target, std::string cycle_atoms_str, char NxO, MolecularModeling::AtomVector& pattern_atoms)
{
    int local_debug = -1;
    std::stringstream pattern;
    pattern << "xC";
    MolecularModeling::Atom* N_or_O = NULL;
    MolecularModeling::AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(MolecularModeling::AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        MolecularModeling::Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == std::string::npos)
        {
            if(t_neighbor->GetName().at(0) != NxO)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
                N_or_O = t_neighbor;
        }
    }
    if(N_or_O != NULL)
    {
        MolecularModeling::Atom* C = NULL;
        pattern << "-" << NxO;
        MolecularModeling::AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
        for(MolecularModeling::AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            MolecularModeling::Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'C')
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C != NULL)
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C == NULL)
                C = n_neighbor;
        }
        if(C != NULL)
        {
            pattern_atoms.push_back(C);
            pattern << "-C";
            MolecularModeling::AtomVector c_neighbors = C->GetNode()->GetNodeNeighbors();
            for(MolecularModeling::AtomVector::iterator it = c_neighbors.begin(); it != c_neighbors.end(); it++)
            {
                MolecularModeling::Atom* c_neighbor = (*it);
                if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0)
                    pattern << c_neighbor->GetName().at(0);
            }
        }
    }
    if (local_debug > 0)
    {
      std::stringstream debugStr;
      debugStr << "CheckxC_NxO_C:" << pattern.str() << "\n";
      gmml::log(__LINE__, __FILE__,  gmml::INF, debugStr.str());
    }
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-N-CHHH") == 0 ||  pattern.str().compare("xCH-NH-CHHH") == 0 || pattern.str().compare("xCH-NH-C") == 0 || pattern.str().compare("xCH-N-C") == 0 ||
                pattern.str().compare("xCHH-N-CHHH") == 0 || pattern.str().compare("xCHH-NH-CHHH") == 0 || pattern.str().compare("xCHH-NH-C") == 0 || pattern.str().compare("xCHH-N-C") == 0 ||
                pattern.str().compare("xC-NH-CHHH") == 0 || pattern.str().compare("xC-NH-C") == 0 || pattern.str().compare("xC-N-CHHH") == 0 || pattern.str().compare("xC-N-C") == 0)
            return "xC-N-CH3";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-O-CHHH") == 0 ||  pattern.str().compare("xCH-OH-CHHH") == 0 || pattern.str().compare("xCH-OH-C") == 0 || pattern.str().compare("xCH-O-C") == 0 ||
                pattern.str().compare("xCHH-O-CHHH") == 0 || pattern.str().compare("xCHH-OH-CHHH") == 0 || pattern.str().compare("xCHH-OH-C") == 0 || pattern.str().compare("xCHH-O-C") == 0 ||
                pattern.str().compare("xC-OH-CHHH") == 0 || pattern.str().compare("xC-OH-C") == 0 || pattern.str().compare("xC-O-CHHH") == 0 || pattern.str().compare("xC-O-C") == 0)
            return "xC-O-CH3";
        else
            return "";
    }
    else
        return "";
}

std::string Assembly::CheckxCOO(MolecularModeling::Atom *target, std::string cycle_atoms_str/*, MolecularModeling::AtomVector& pattern_atoms*/)
{
  int local_debug = -1;
    std::stringstream pattern;
    pattern << "xC";
    MolecularModeling::Atom* O1 = NULL;
    MolecularModeling::Atom* O2 = NULL;
    MolecularModeling::AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(MolecularModeling::AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        MolecularModeling::Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == std::string::npos)
        {
            if(t_neighbor->GetName().at(0) != 'O')
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == 'O' && O1 == NULL)
                O1 = t_neighbor;
            else if(t_neighbor->GetName().at(0) == 'O' && O2 == NULL)
                O2 = t_neighbor;
        }
    }
    if(O1 != NULL && O2 != NULL)
    {
        pattern << "OO";
        MolecularModeling::AtomVector o1_neighbors = O1->GetNode()->GetNodeNeighbors();
        for(MolecularModeling::AtomVector::iterator it = o1_neighbors.begin(); it != o1_neighbors.end(); it++)
        {
            MolecularModeling::Atom* o1_neighbor = (*it);
            if(o1_neighbor->GetId().compare(target->GetId()) != 0)
                pattern << o1_neighbor->GetName().at(0);
        }
        MolecularModeling::AtomVector o2_neighbors = O2->GetNode()->GetNodeNeighbors();
        for(MolecularModeling::AtomVector::iterator it = o2_neighbors.begin(); it != o2_neighbors.end(); it++)
        {
            MolecularModeling::Atom* o2_neighbor = (*it);
            if(o2_neighbor->GetId().compare(target->GetId()) != 0)
                pattern << o2_neighbor->GetName().at(0);
        }
    }
    if (local_debug > 0)
    {
      std::stringstream debugStr;
      debugStr << "CheckxCOO:" << pattern.str() << "\n";
      gmml::log(__LINE__, __FILE__,  gmml::INF, debugStr.str());
    }
    if(pattern.str().compare("xCOO") == 0 || pattern.str().compare("xCHOO") == 0)
        return "xC-(O,O)";
    else if(pattern.str().compare("xCOOH") == 0 || pattern.str().compare("xCHOOH") == 0)
        return "xC-(O,OH)";
    else
        return "";
}

// std::string Assembly::CheckxC_NxO_CH3C_COO(MolecularModeling::Atom *target, std::string cycle_atoms_str, MolecularModeling::AtomVector& pattern_atoms)
// {//isopropionic acid
  // int local_debug = -1;
  // std::stringstream pattern;
  // pattern << "xC";
  // MolecularModeling::Atom* N_or_O = NULL;
  // MolecularModeling::AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
  // for(MolecularModeling::AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
  // {
  //   MolecularModeling::Atom* t_neighbor = (*it);
  //   if(cycle_atoms_str.find(t_neighbor->GetId()) == std::string::npos)
  //   {
  //     if(t_neighbor->GetName().at(0) != NxO)
  //       pattern << t_neighbor->GetName().at(0);
  //     else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
  //       pattern << t_neighbor->GetName().at(0);
  //     else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
  //       N_or_O = t_neighbor;
  //   }
  // }
  // if(N_or_O != NULL)
  // {
  //   MolecularModeling::Atom* C = NULL;
  //   pattern << "-" << NxO;
  //   MolecularModeling::AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
  //   for(MolecularModeling::AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
  //   {
  //     MolecularModeling::Atom* n_neighbor = (*it);
  //     if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'C')
  //       pattern << n_neighbor->GetName().at(0);
  //     else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C != NULL)
  //       pattern << n_neighbor->GetName().at(0);
  //     else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C == NULL)
  //       C = n_neighbor;
  //   }
  //   if(C != NULL)
  //   {
  //     pattern_atoms.push_back(C);
  //     pattern << "-C";
  //     MolecularModeling::AtomVector c_neighbors = C->GetNode()->GetNodeNeighbors();
  //     for(MolecularModeling::AtomVector::iterator it = c_neighbors.begin(); it != c_neighbors.end(); it++)
  //     {
  //       MolecularModeling::Atom* c_neighbor = (*it);
  //       if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0)
  //         pattern << c_neighbor->GetName().at(0);
  //     }
  //   }
  // }
  // if (local_debug > 0)
  // {
  //   std::stringstream debugStr;
  //   debugStr << "CheckxC_NxO_C:" << pattern.str() << "\n";
  //   gmml::log(__LINE__, __FILE__,  gmml::INF, debugStr.str());
  // }
  // if(NxO == 'N')
  // {
  //   if(pattern.str().compare("xCH-N-CHHH") == 0 ||
  //     pattern.str().compare("xCH-NH-CHHH") == 0 ||
  //     pattern.str().compare("xCH-NH-C") == 0 ||
  //     pattern.str().compare("xCH-N-C") == 0 ||
  //     pattern.str().compare("xCHH-N-CHHH") == 0 ||
  //     pattern.str().compare("xCHH-NH-CHHH") == 0 ||
  //     pattern.str().compare("xCHH-NH-C") == 0 ||
  //     pattern.str().compare("xCHH-N-C") == 0 ||
  //     pattern.str().compare("xC-NH-CHHH") == 0 ||
  //     pattern.str().compare("xC-NH-C") == 0 ||
  //     pattern.str().compare("xC-N-CHHH") == 0 ||
  //     pattern.str().compare("xC-N-C") == 0)
  //   {
  //     return "xC-N-CH3";
  //   }
  //   else
  //     return "";
  // }
  // else if(NxO == 'O')
  // {
  //   if(pattern.str().compare("xCH-O-CHHH") == 0 ||
  //     pattern.str().compare("xCH-OH-CHHH") == 0 ||
  //     pattern.str().compare("xCH-OH-C") == 0 ||
  //     pattern.str().compare("xCH-O-C") == 0 ||
  //     pattern.str().compare("xCHH-O-CHHH") == 0 ||
  //     pattern.str().compare("xCHH-OH-CHHH") == 0 ||
  //     pattern.str().compare("xCHH-OH-C") == 0 ||
  //     pattern.str().compare("xCHH-O-C") == 0 ||
  //     pattern.str().compare("xC-OH-CHHH") == 0 ||
  //     pattern.str().compare("xC-OH-C") == 0 ||
  //     pattern.str().compare("xC-O-CHHH") == 0 ||
  //     pattern.str().compare("xC-O-C") == 0)
  //   {
  //     return "xC-O-CH3";
  //   }
  //   else
  //     return "";
  // }
  // else
  //   return "";
// }

// void Assembly::GenerateCompleteSugarName(Glycan::Monosaccharide *mono)
// {
//     std::stringstream in_bracket;
//     std::stringstream head;
//     std::stringstream tail;
//     bool minus_one = false;
//     if(mono->derivatives_map_.find("-1") != mono->derivatives_map_.end())
//     {
//         minus_one = true;
//     }
//     for(std::map<std::string, std::string>::iterator it1 = mono->derivatives_map_.begin(); it1 != mono->derivatives_map_.end(); it1++)
//     {
//         std::string key = (*it1).first;
//         std::string value = (*it1).second;
//         std::string long_name_pattern = "";
//         std::string cond_name_pattern = "";
//         std::string long_name_pattern_at_minus_one = "";
//         std::string long_name_pattern_at_plus_one = "";
//         std::string pattern = "";
//         if(value.compare("xCH-N") == 0)
//         {
//             long_name_pattern = "-osamine";
//             cond_name_pattern = "N";
//             pattern = "CH-N";
//             AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
//         }
//         if(value.compare("xC-N-C=OCH3") == 0)
//         {
//             long_name_pattern = "N-acetyl-";
//             cond_name_pattern = "NAc";
//             pattern = "C-N-C=OCH3";
//             AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
//         }
//         if(value.compare("xC-N-C=OCH2OH") == 0)
//         {
//             long_name_pattern = "N-glycolyl-";
//             cond_name_pattern = "NGc";
//             pattern = "C-N-C=OCH2OH";
//             AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
//         }
//         if(value.compare("xC-N-SO3") == 0)
//         {
//             long_name_pattern = "N-sulfo-";
//             cond_name_pattern = "NS";
//             pattern = "C-N-SO3";
//             AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
//         }
//         if(value.compare("xC-N-PO3") == 0)
//         {
//             long_name_pattern = "N-phospho-";
//             cond_name_pattern = "NP";
//             pattern = "C-N-PO3";
//             AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
//         }
//         if(value.compare("xC-N-CH3") == 0)
//         {
//             long_name_pattern = "N-methyl-";
//             cond_name_pattern = "NMe";
//             pattern = "C-N-CH3";
//             AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
//         }
//         if(value.compare("xC-O-C=OCH3") == 0)
//         {
//             long_name_pattern = "-acetyl-";
//             cond_name_pattern = "Ac";
//             pattern = "C-O-C=OCH3";
//             AddDerivativeRuleInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
//         }
//         if(value.compare("xC-O-C=OCH2OH") == 0)
//         {
//             long_name_pattern = "-glycolyl-";
//             cond_name_pattern = "Gc";
//             pattern = "C-O-C=OCH2OH";
//             AddDerivativeRuleInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
//         }
//         if(value.compare("xC-O-SO3") == 0)
//         {
//             long_name_pattern = "-sulfo-";
//             cond_name_pattern = "S";
//             pattern = "C-O-SO3";
//             AddDerivativeRuleInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
//         }
//         if(value.compare("xC-O-PO3") == 0)
//         {
//             long_name_pattern = "-phospho-";
//             cond_name_pattern = "P";
//             pattern = "C-O-PO3";
//             AddDerivativeRuleInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
//         }
//         if(value.compare("xC-O-CH3") == 0)
//         {
//             long_name_pattern = "-methyl-";
//             cond_name_pattern = "Me";
//             pattern = "C-O-CH3";
//             AddDerivativeRuleInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
//         }
//         if(value.compare("xC-(O,OH)") == 0)
//         {
//             long_name_pattern_at_minus_one = "-ulosonic acid";
//             long_name_pattern_at_plus_one = "-uronic acid";
//             cond_name_pattern = "AH";
//             pattern = "C-(O,OH)";
//             AddModificationRuleTwoInfo(key, pattern, mono, long_name_pattern_at_minus_one, long_name_pattern_at_plus_one, cond_name_pattern, tail, minus_one, in_bracket);
//         }
//         if(value.compare("xC-(O,O)") == 0)
//         {
//             long_name_pattern_at_minus_one = "-ulosonate";
//             long_name_pattern_at_plus_one = "-uronate";
//             cond_name_pattern = "A";
//             pattern = "C-(O,O)";
//             AddModificationRuleTwoInfo(key, pattern, mono, long_name_pattern_at_minus_one, long_name_pattern_at_plus_one, cond_name_pattern, tail, minus_one, in_bracket);
//         }
//     }
//     if(in_bracket.str().size() != 0)
//     {
//         std::stringstream short_name;
//         std::string sn;
//         if(mono->sugar_name_.monosaccharide_short_name_.compare("") != 0)
//             sn = mono->sugar_name_.monosaccharide_short_name_;
//         else
//             sn = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
//
//         ///moving a, b or x to after the bracket: short-name + [...] + a/b/x and removing ", " from the end of bracket stream
//         int condensed_name_size = sn.size();
//         std::string condensed_name = sn;
//         std::string new_name_part1 = "";
//         char new_name_part2 = ' ';
//         if(condensed_name_size > 0)
//         {
//           new_name_part1 = condensed_name.substr(0, (condensed_name_size - 1));///short_name
//           new_name_part2 = condensed_name.at(condensed_name_size - 1);///a/b/x
//         }
//         short_name << new_name_part1 << "[" << in_bracket.str().substr(0, in_bracket.str().size() - 1) << "]" << new_name_part2;
//
//         mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
//     }
//     else if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0 && mono->sugar_name_.monosaccharide_short_name_.compare("") == 0)
//     {
//         mono->sugar_name_.monosaccharide_short_name_ = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
//     }
//     std::stringstream long_name;
//     if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
//     {
//         long_name << head.str() << mono->sugar_name_.monosaccharide_stereochemistry_name_ << tail.str();
//         mono->sugar_name_.monosaccharide_name_ = long_name.str();
//     }
// }

void Assembly::AddModificationRuleOneInfo(std::string key, std::string pattern, Glycan::Monosaccharide* mono, std::string long_name_pattern,
                                          std::string cond_name_pattern, std::stringstream& head,std::stringstream& tail, bool minus_one,
                                          std::stringstream& in_bracket)
{
  int local_debug = -1;
  std::stringstream ss;
  ss << pattern;
  if(key.compare("a") == 0)
  {
    ss << " is at warning position: anomeric";
    gmml::log(__LINE__, __FILE__,  gmml::WAR, ss.str());
    Glycan::Note* der_mod_note = new Glycan::Note();
    der_mod_note->type_ = Glycan::WARNING;
    der_mod_note->category_ = Glycan::DER_MOD;
    std::stringstream note;
    note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
    der_mod_note->description_ = note.str();
    // this->AddNote(der_mod_note);
    mono->mono_notes_.push_back(der_mod_note);
    if(local_debug > 0)
    {
      std::cout << ss.str() << "\n";
    }
  }
  else if(key.compare("2") == 0 && mono->sugar_name_.ring_type_.compare("P") == 0 &&
          find(mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), "-1") == mono->chemical_code_->right_down_.end() &&
          find(mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), "-1") == mono->chemical_code_->right_up_.end() &&
          mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
  {
    if(long_name_pattern.compare("-osamine") == 0)
        tail << long_name_pattern;
    else
        head << long_name_pattern;
    std::stringstream short_name;
    if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
    {
      ///moving a, b or x to after the N expression: short-name + Condensed name pattern + a/b/x
      int stereo_condensed_name_size = mono->sugar_name_.monosaccharide_stereochemistry_short_name_.size();
      std::string stereo_condensed_name = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
      std::string new_name_part1 = stereo_condensed_name.substr(0, (stereo_condensed_name_size - 1));///short_name
      char new_name_part2 = stereo_condensed_name.at(stereo_condensed_name_size - 1);///a/b/x
      short_name << new_name_part1 << cond_name_pattern << new_name_part2;
      mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
    }
  }
  else if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
  {
    if(!minus_one)
      ss << " is at error position: 4";
    else
      ss << " is at error position: 5";
    gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
    Glycan::Note* der_mod_note = new Glycan::Note();
    der_mod_note->type_ = Glycan::ERROR;
    der_mod_note->category_ = Glycan::DER_MOD;
    std::stringstream note;
    note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
    der_mod_note->description_ = note.str();
    // this->AddNote(der_mod_note);
    mono->mono_notes_.push_back(der_mod_note);
    if(local_debug > 0)
    {
      std::cout << ss.str() << "\n";
    }
  }
  else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
  {
    if(!minus_one)
      ss << " is at error position: 5";
    else
      ss << " is at error position: 6";
    gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
    Glycan::Note* der_mod_note = new Glycan::Note();
    der_mod_note->type_ = Glycan::ERROR;
    der_mod_note->category_ = Glycan::DER_MOD;
    std::stringstream note;
    note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
    der_mod_note->description_ = note.str();
    // this->AddNote(der_mod_note);
    mono->mono_notes_.push_back(der_mod_note);
    if(local_debug > 0)
    {
      std::cout << ss.str() << "\n";
    }
  }
  else
  {
    if(!minus_one)
    {
      if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
        in_bracket << mono->cycle_atoms_.size() - 1 + gmml::ConvertString<int>(key) << cond_name_pattern << ",";
      else
        in_bracket << key << cond_name_pattern << ",";
    }
    else
    {
      if(key.compare("-1") == 0)
        in_bracket << "1" << cond_name_pattern << ",";
      else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
        in_bracket << mono->cycle_atoms_.size() + gmml::ConvertString<int>(key) << cond_name_pattern << ",";
      else
        in_bracket << gmml::ConvertString<int>(key) + 1 << cond_name_pattern << ",";
    }
  }
}

void Assembly::AddUnknownDerivativeRuleInfo(std::string key, std::string pattern, Glycan::Monosaccharide *mono, std::string long_name_pattern,
                                  std::string cond_name_pattern, std::stringstream &head, bool minus_one, std::stringstream &in_bracket)
{
  std::stringstream ss;
  ss << pattern;
  if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
  {
    if(!minus_one)
    {
      if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
        head << mono->cycle_atoms_.size() - 1 + gmml::ConvertString<int>(key) << long_name_pattern;
      else if(key.compare("a") != 0)
        head << gmml::ConvertString<int>(key) << long_name_pattern;
    }
    else
    {
      if(key.compare("-1") == 0)
        head << "1" << long_name_pattern;
      else if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
        head << mono->cycle_atoms_.size() + gmml::ConvertString<int>(key) << long_name_pattern;
      else if(key.compare("a") != 0)
        head << gmml::ConvertString<int>(key) + 1 << long_name_pattern;
    }
  }
  if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
  {
    if(!minus_one)
    {
      if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
      {
        in_bracket << mono->cycle_atoms_.size() - 1 + gmml::ConvertString<int>(key) << cond_name_pattern << ",";
        // gmml::log(__LINE__, __FILE__, gmml::INF, in_bracket.str());
      }
      else if(key.compare("a") != 0)
      {
        in_bracket << gmml::ConvertString<int>(key) << cond_name_pattern << ",";
        // gmml::log(__LINE__, __FILE__, gmml::INF, in_bracket.str());
      }
    }
    else
    {
      if(key.compare("-1") == 0)
      {
        in_bracket << "1" << cond_name_pattern << ",";
        // gmml::log(__LINE__, __FILE__, gmml::INF, in_bracket.str());
      }
      else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
      {
        in_bracket << mono->cycle_atoms_.size() + gmml::ConvertString<int>(key) << cond_name_pattern << ",";
        // gmml::log(__LINE__, __FILE__, gmml::INF, in_bracket.str());
      }
      else if(key.compare("a") != 0)
      {
        in_bracket << gmml::ConvertString<int>(key) + 1<< cond_name_pattern << ",";
        // gmml::log(__LINE__, __FILE__, gmml::INF, in_bracket.str());
      }
    }
  }
}


void Assembly::AddDerivativeRuleInfo(std::string key, std::string pattern, Glycan::Monosaccharide *mono, std::string long_name_pattern, std::string cond_name_pattern, std::stringstream &head,
                                     bool minus_one, std::stringstream &in_bracket)
{
  int local_debug = -1;
    std::stringstream ss;
    ss << pattern;
    // std::cout << key << ": " << pattern << "\n";
    if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
    {
        if(!minus_one)
        {
            if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                head << mono->cycle_atoms_.size() - 1 + std::stoi(key.substr(1,1)) << long_name_pattern;
           else if(key.compare("a") == 0)
               head << "1" << long_name_pattern;
            else if(key.compare("a") != 0)
                head << gmml::ConvertString<int>(key) << long_name_pattern;
        }
        else
        {
            if(key.compare("-1") == 0)
                head << "1" << long_name_pattern;
            else if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                head << mono->cycle_atoms_.size()  +  std::stoi(key.substr(1,1)) << long_name_pattern;
           else if(key.compare("a") == 0)
               head << "2" << long_name_pattern;
            else if(key.compare("a") != 0)
                head << gmml::ConvertString<int>(key) + 1 << long_name_pattern;
        }
    }
    if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
    {
        if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
        {
            if(!minus_one)
                ss << " is at error position: 4";
            else
                ss << " is at error position: 5";
            gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
            Glycan::Note* der_mod_note = new Glycan::Note();
            der_mod_note->type_ = Glycan::ERROR;
            der_mod_note->category_ = Glycan::DER_MOD;
            std::stringstream note;
            note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
            der_mod_note->description_ = note.str();
            // this->AddNote(der_mod_note);
            mono->mono_notes_.push_back(der_mod_note);
            std::cout << ss.str() << "\n";
        }
        else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
        {
            if(!minus_one)
                ss << " is at error position: 5";
            else
                ss << " is at error position: 6";
            gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
            Glycan::Note* der_mod_note = new Glycan::Note();
            der_mod_note->type_ = Glycan::ERROR;
            der_mod_note->category_ = Glycan::DER_MOD;
            std::stringstream note;
            note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
            der_mod_note->description_ = note.str();
            // this->AddNote(der_mod_note);
            mono->mono_notes_.push_back(der_mod_note);
            std::cout << ss.str() << "\n";
        }
        else
        {
            if(!minus_one)
            {
                if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                {
                  in_bracket << mono->cycle_atoms_.size() - 1 + gmml::ConvertString<int>(key) << cond_name_pattern << ",";
                  gmml::log(__LINE__, __FILE__, gmml::INF, in_bracket.str());
                }
               else if(key.compare("a") == 0)
                   in_bracket << "1" << cond_name_pattern << ",";
                else if(key.compare("a") != 0)
                {
                  in_bracket << gmml::ConvertString<int>(key) << cond_name_pattern << ",";
                  gmml::log(__LINE__, __FILE__, gmml::INF, in_bracket.str());
                }
            }
            else
            {
                if(key.compare("-1") == 0)
                {
                  in_bracket << "1" << cond_name_pattern << ",";
                  gmml::log(__LINE__, __FILE__, gmml::INF, in_bracket.str());
                }
                else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                {
                  in_bracket << mono->cycle_atoms_.size() + gmml::ConvertString<int>(key) << cond_name_pattern << ",";
                  gmml::log(__LINE__, __FILE__, gmml::INF, in_bracket.str());
                }
               else if(key.compare("a") == 0)
                   in_bracket << "2" << cond_name_pattern << ",";
                else if(key.compare("a") != 0)
                {
                  in_bracket << gmml::ConvertString<int>(key) + 1 << cond_name_pattern << ",";
                  gmml::log(__LINE__, __FILE__, gmml::INF, in_bracket.str());
                }
            }
        }
    }
}
void Assembly::AddModificationRuleTwoInfo(std::string key, std::string pattern, Glycan::Monosaccharide *mono, std::string long_name_pattern_at_minus_one, std::string long_name_pattern_at_plus_one,
                                          std::string cond_name_pattern, std::stringstream &tail, bool minus_one, std::stringstream &in_bracket)
{
    Glycan::Note* der_mod_note = new Glycan::Note();
    std::stringstream ss;
    ss << pattern;
    if((key.compare("-1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0) && mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
    {
        if(key.compare("-1") == 0)
        {
            ss << " is at warning position: " << gmml::ConvertString<int>(key) + 1;
            der_mod_note->type_ = Glycan::WARNING;
        }
        if(!minus_one)
        {
            ss << " is at error position: " << key;
            der_mod_note->type_ = Glycan::ERROR;
        }
        else
        {
            ss << " is at error position: " << gmml::ConvertString<int>(key) + 1;
            der_mod_note->type_ = Glycan::ERROR;
        }
        std::cout << ss.str() << "\n";
        gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
        tail << long_name_pattern_at_minus_one;

        if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
        {
            if(!minus_one)
            {
                if( key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() - 1 + gmml::ConvertString<int>(key) << cond_name_pattern << ",";
            }
            else
            {
                if(key.compare("-1") == 0)
                    in_bracket << "1" << cond_name_pattern << ",";
                else if( key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() + gmml::ConvertString<int>(key) << cond_name_pattern << ",";
            }
        }
    }
    else if(key.compare("+1") == 0 && mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
    {
        tail << long_name_pattern_at_plus_one;
        std::stringstream short_name;
        if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
        {
            ///moving a, b or x to after the AH expression: short-name + AH + a/b/x
            int stereo_condensed_name_size = mono->sugar_name_.monosaccharide_stereochemistry_short_name_.size();
            std::string stereo_condensed_name = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
            std::string new_name_part1 = stereo_condensed_name.substr(0, (stereo_condensed_name_size - 1));///short_name
            char new_name_part2 = stereo_condensed_name.at(stereo_condensed_name_size - 1);///a/b/x
            short_name << new_name_part1 << cond_name_pattern << new_name_part2;

            mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
        }
    }
    else
    {
        if(!minus_one)
        {
            if(key.compare("a") != 0)
                ss << " is at warning position: " << key;
            else
                ss << " is at warning position: 1";
        }
        else
        {
            if(key.compare("a") != 0)
                ss << " is at warning position: " << gmml::ConvertString<int>(key) + 1;
            else
                ss << " is at warning position: 2";
        }
        der_mod_note->type_ = Glycan::WARNING;
    }
    der_mod_note->category_ = Glycan::DER_MOD;
    std::stringstream note;
    note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
    der_mod_note->description_ = note.str();
    // this->AddNote(der_mod_note);
    mono->mono_notes_.push_back(der_mod_note);
    // std::cout << ss.str() << "\n";
    gmml::log(__LINE__, __FILE__,  gmml::WAR, ss.str());
}

// void Assembly::UpdateComplexSugarChemicalCode( Glycan::Monosaccharide * mono ) {
//     for( std::map< std::string, std::string >::iterator it1 = mono->derivatives_map_.begin(); it1 != mono->derivatives_map_.end(); it1++ ) {
//       std::string key = ( *it1 ).first;
//       std::string value = ( *it1 ).second;
//       if( value.compare( "" ) != 0 ) {
//         std::string code = "";
//         if( value.compare( "xCH-N" ) == 0 ) {
//           code = key + "N";
//         } else if( value.compare( "xC-N-C=OCH3" ) == 0 ) {
//           code = key + "NAc";
//         } else if( value.compare( "xC-N-C=OCH2OH" ) == 0 ) {
//           code = key + "NGc";
//         } else if( value.compare( "xC-N-SO3" ) == 0 ) {
//           code = key + "NS";
//         } else if( value.compare( "xC-N-PO3" ) == 0 ) {
//           code = key + "NP";
//         } else if( value.compare( "xC-N-CH3" ) == 0 ) {
//           code = key + "NMe";
//         } else if( value.compare( "xC-O-C=OCH3" ) == 0 ) {
//           code = key + "Ac";
//         } else if( value.compare( "xC-O-C=OCH3OH" ) == 0 ) {
//           code = key + "Gc";
//         } else if( value.compare( "xC-O-SO3" ) == 0 ) {
//           code = key + "S";
//         } else if( value.compare( "xC-O-PO3" ) == 0 ) {
//           code = key + "P";
//         } else if( value.compare( "xC-O-CH3" ) == 0 ) {
//           code = key + "Me";
//         } else if( value.compare( "xC-(O,OH)" ) == 0 ) {
//           code = key + "AH";
//         } else if( value.compare( "xC-(O,O)" ) == 0 ) {
//           code = key + "A";
//         }
//         std::vector< std::string >::iterator index_it;
//         if( ( key.compare( "a" ) == 0 ) || ( key.compare( "-1" ) == 0 ) || ( key.compare( "+1" ) == 0 ) || ( key.compare( "+2" )== 0 ) || ( key.compare( "+3" ) == 0 ) ) {
//           if( ( index_it = find( mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), key ) ) != mono->chemical_code_->right_up_.end() ) {
//               ( *index_it ) = code;
//           } else if( ( index_it = find( mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), key ) ) != mono->chemical_code_->right_down_.end() ) {
//               ( *index_it ) = code;
//           }
//         } else if( ( key.compare( "2" ) == 0 ) || ( key.compare( "3" ) == 0 ) || ( key.compare( "4" ) == 0 ) ) {
//           if( ( index_it = find( mono->chemical_code_->left_up_.begin(), mono->chemical_code_->left_up_.end(), key ) ) != mono->chemical_code_->left_up_.end() ) {
//               ( *index_it ) = code;
//           } else if( ( index_it = find( mono->chemical_code_->left_down_.begin(), mono->chemical_code_->left_down_.end(), key ) ) != mono->chemical_code_->left_down_.end() ) {
//               ( *index_it ) = code;
//           }
//         }
//       }
//     }
// }
//
// void Assembly::UpdatePdbCode( Glycan::Monosaccharide * mono )
// {
//   std::string code = mono->chemical_code_->toString();
//   for( int i = 0; i < gmml::COMPLEXSUGARNAMELOOKUPSIZE; i++ )
//   {
//     if( code.compare( gmml::COMPLEXSUGARNAMELOOKUP[ i ].chemical_code_string_ ) == 0 )
//     {
//       mono->sugar_name_.pdb_code_ = gmml::COMPLEXSUGARNAMELOOKUP[ i ].pdb_code_;
//       return;
//     }
//   }
//   for( int i = 0; i < gmml::SUGARNAMELOOKUPSIZE; i++ )
//   {
//     if( code.compare( gmml::SUGARNAMELOOKUP[ i ].chemical_code_string_ ) == 0 )
//     {
//       mono->sugar_name_.pdb_code_ = gmml::SUGARNAMELOOKUP[ i ].pdb_code_;
//       return;
//     }
//   }
// }
void Assembly::createOligosaccharideGraphs(std::vector<Glycan::Monosaccharide*> detected_monos,
                                                            gmml::ResidueNameMap dataset_residue_names,
                                                            int& number_of_covalent_links,
                                                            int& number_of_probable_non_covalent_complexes)
{
  std::string terminal_residue_name = "";
  gmml::ResidueNameMap common_terminal_residues = gmml::InitializeCommonTerminalResidueMap();
  std::map<Glycan::Monosaccharide*, std::vector<Glycan::Monosaccharide*> > monos_table = std::map<Glycan::Monosaccharide*, std::vector<Glycan::Monosaccharide*> >();
  std::map<Glycan::Monosaccharide*, std::vector<std::string> > monos_table_linkages = std::map<Glycan::Monosaccharide*, std::vector<std::string> >();
  // gmml::log(__LINE__, __FILE__,  gmml::INF, " Start iterating on list ..." );
  ///Iterating on list of monos to check if there is a connection to another mono in the list
  for(std::vector<Glycan::Monosaccharide*>::iterator it = detected_monos.begin(); it != detected_monos.end(); it++)
  {
    Glycan::Monosaccharide* mono1 = (*it);
    monos_table[mono1] = std::vector<Glycan::Monosaccharide*>();
    monos_table_linkages[mono1] = std::vector<std::string>();
    for(std::vector<std::vector<MolecularModeling::Atom*> >::iterator it1 = mono1->side_atoms_.begin(); it1 != mono1->side_atoms_.end(); it1++) ///iterate on side atoms
    {
      int index = distance(mono1->side_atoms_.begin(), it1);
      std::vector<MolecularModeling::Atom*>  sides = (*it1);
      std::map<MolecularModeling::Atom*, MolecularModeling::Atom*> target_parent_map = std::map<MolecularModeling::Atom*, MolecularModeling::Atom*>(); /// A map of target atom to it's parent atom. Target atom is a non ring oxygen or nitrogen
      if(it1 == mono1->side_atoms_.begin())///side atoms of anomeric
      {
        if(sides.at(1) != NULL)
          target_parent_map[sides.at(1)] = mono1->cycle_atoms_.at(0);
      }
      else if(it1 == mono1->side_atoms_.end() - 1) ///side atoms of last carbon of the ring
      {
        for(std::vector<MolecularModeling::Atom*> ::iterator last_c_side_it = sides.begin(); last_c_side_it != sides.end(); last_c_side_it++)
        {
          MolecularModeling::Atom* side_of_last_carbon = (*last_c_side_it);
          if(side_of_last_carbon != NULL)
          {
            std::vector<MolecularModeling::Atom*>  last_c_side_neighbors = side_of_last_carbon->GetNode()->GetNodeNeighbors();
            for(std::vector<MolecularModeling::Atom*> ::iterator it2 = last_c_side_neighbors.begin(); it2 != last_c_side_neighbors.end(); it2++)
            {
              if((*it2)->GetId().at(0) == 'O' || (*it2)->GetId().at(0) == 'N')
              {
                target_parent_map[(*it2)] = side_of_last_carbon;
                break;
              }
            }
          }
        }
      }
      else
      {
        if(sides.at(1) != NULL)
          target_parent_map[sides.at(1)] = mono1->cycle_atoms_.at(index);///index 1 of each side is for non-carbon side atoms in the std::vector<std::vector<MolecularModeling::Atom*> > structure
      }
      ///Examine neighbors of each target atom to check if they can be found in other monos side/ring atoms
      for(std::map<MolecularModeling::Atom*, MolecularModeling::Atom*>::iterator map_it = target_parent_map.begin(); map_it != target_parent_map.end(); map_it++)
      {
        bool found_in_other_mono = false;
        MolecularModeling::Atom* target = (*map_it).first;
        MolecularModeling::Atom* target_parent = (*map_it).second;
        std::vector<MolecularModeling::Atom*>  t_neighbors = target->GetNode()->GetNodeNeighbors();
        for(std::vector<MolecularModeling::Atom*> ::iterator it2 = t_neighbors.begin(); it2 != t_neighbors.end(); it2++)
        {
          MolecularModeling::Atom* t_neighbor = (*it2);
          if(t_neighbor->GetId().compare(target_parent->GetId()) != 0)///making sure neighbor is not the parent of target atom
          {
            for(std::vector<Glycan::Monosaccharide*>::iterator it3 = detected_monos.begin(); it3 != detected_monos.end(); it3++)
            {
              if(it3 != it)///Cheking monos other than the current mono
              {
                Glycan::Monosaccharide* mono2 = (*it3);
                std::vector<MolecularModeling::Atom*>  mono2_sides = mono2->side_atoms_.at(mono2->side_atoms_.size() - 1); ///side of last ring carbon
                bool found_in_side = false;
                for(std::vector<MolecularModeling::Atom*> ::iterator mono2_last_c_side_it = mono2_sides.begin(); mono2_last_c_side_it != mono2_sides.end(); mono2_last_c_side_it++)
                {
                  MolecularModeling::Atom* mono2_last_c_side = (*mono2_last_c_side_it);
                  if(mono2_last_c_side != NULL)
                  {
                    if(t_neighbor->GetId().compare(mono2_last_c_side->GetId()) == 0) ///target atom has been attached to another cycle's side atom
                      found_in_side = true;
                  }
                }
                if(found_in_side || mono2->cycle_atoms_str_.find(t_neighbor->GetId()) != std::string::npos) //if target's neighbor found in another mono's side or ring atoms
                {
                  found_in_other_mono = true;
                  monos_table[mono1].push_back(mono2);
                  std::string mono1_carbon = target_parent->GetId();
                  std::string mono1_name = "";
                  std::string mono2_carbon = t_neighbor->GetId();
                  std::string mono2_name = "";
                  if(mono1->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                    mono1_name = mono1->sugar_name_.monosaccharide_short_name_;
                  else
                    mono1_name = mono1->sugar_name_.monosaccharide_stereochemistry_short_name_;
                  if(mono2->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                    mono2_name = mono2->sugar_name_.monosaccharide_short_name_;
                  else
                    mono2_name = mono2->sugar_name_.monosaccharide_stereochemistry_short_name_;
                  std::stringstream linkage;
                  linkage << mono1_carbon << "-" << target->GetId() << "-" << mono2_carbon;
                  monos_table_linkages[mono1].push_back(linkage.str());
                  // std::string mono1_carbon_name(1, mono1_carbon[1]);
                  // std::string mono2_carbon_name(1, mono2_carbon[1]);
                  Glycan::GlycosidicLinkage* thisLinkage = new Glycan::GlycosidicLinkage(mono1, mono2, mono1_carbon, mono2_carbon);
                  mono1->mono_neighbors_.push_back(std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*>(thisLinkage, mono2));
                  // std::stringstream ss;
                  // ss << mono1->cycle_atoms_[0]->GetResidue()->GetId() << " is connected to " << mono2->cycle_atoms_[0]->GetResidue()->GetId() << " via " << mono1_carbon_name << "-" << mono2_carbon_name;

                  // gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
                  break;
                }
              }
            }
          }
        }
        if(found_in_other_mono)
          break;
      }
    }
    // for(std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> >::iterator monoNeighbor = mono1->mono_neighbors_.begin(); monoNeighbor != mono1->mono_neighbors_.end(); monoNeighbor++)
    // {
    //   std::stringstream ss;
    //   Glycan::Monosaccharide* thisNeighbor = (monoNeighbor->second);
    //   Glycan::GlycosidicLinkage* thisLinkage = (monoNeighbor->first);
    //   ss << mono1->cycle_atoms_[0]->GetResidue()->GetId() << " is connected to " << thisNeighbor->cycle_atoms_[0]->GetResidue()->GetId() << " via " << thisLinkage->linkage_type_;
    //   gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
    // }

  }
  // gmml::log(__LINE__, __FILE__,  gmml::INF, " Done iterating list" );


  std::vector<int> visited_monos = std::vector<int>();
  std::vector<Glycan::Oligosaccharide*> oligosaccharides = std::vector<Glycan::Oligosaccharide*>();

  std::vector<std::string> checked_linkages = std::vector<std::string>();
  // gmml::log(__LINE__, __FILE__,  gmml::INF, " Start for loop ..." );
  for(std::map<Glycan::Monosaccharide*, std::vector<Glycan::Monosaccharide*> >::iterator it = monos_table.begin(); it != monos_table.end(); it++)
  {
    Glycan::Monosaccharide* key = (*it).first;
    std::vector<Glycan::Monosaccharide*> values = (*it).second;

    std::vector<std::string> visited_linkages = std::vector<std::string>();
    if(find(visited_monos.begin(), visited_monos.end(), key->mono_id_) == visited_monos.end())///if the mono is not visited
    {
      bool isRoot = false;
      std::stringstream anomeric_linkage;
      anomeric_linkage << key->cycle_atoms_.at(0)->GetId() << "-";
      std::vector<std::string> mono_linkages = monos_table_linkages[key];
      std::vector<MolecularModeling::Atom*>  terminal_atoms = std::vector<MolecularModeling::Atom*> ();
      if(values.size() == 0) ///mono is not attached to any other mono
      {
        MolecularModeling::Atom* anomeric_o = NULL;
        if(key->side_atoms_.at(0).at(1) != NULL)
          anomeric_o = key->side_atoms_.at(0).at(1);
        if(anomeric_o != NULL)
        {
          if(dataset_residue_names.find(anomeric_o->GetResidue()->GetName()) != dataset_residue_names.end() ||
              common_terminal_residues.find(anomeric_o->GetResidue()->GetName()) != common_terminal_residues.end())///check if there is any terminal
          {
            terminal_residue_name = anomeric_o->GetResidue()->GetName();
          }
          else
          {
            terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
          }
        }
        isRoot = true;
      }
      else if (values.size() == 1 && ///mono is attached to one other mono
             (monos_table_linkages[values.at(0)].size() == 1))///the other mono is only attached to this mono
      {
        ///CHECKING LINKAGE ISSUES, e.g. C1-O3-C4 is an issue
        CheckLinkageNote(key, values.at(0), mono_linkages.at(0), checked_linkages);
        std::stringstream other_mono_anomeric_linkage_as_right_side;
        other_mono_anomeric_linkage_as_right_side << "-" << values.at(0)->cycle_atoms_.at(0)->GetId();///atom id on the right side of the linkage c-o-c

        MolecularModeling::Atom* anomeric_o = NULL;
        MolecularModeling::Atom* o_neighbor_1 = NULL;
        MolecularModeling::Atom* o_neighbor_2 = NULL;
        std::vector<MolecularModeling::Atom*>  o_neighbors = std::vector<MolecularModeling::Atom*> ();
        if(key->side_atoms_.at(0).at(1) != NULL)///Getting the information of anomeric oxygen's neighbors is needed for choosing the root
        {
          anomeric_o = key->side_atoms_.at(0).at(1);
          o_neighbors = anomeric_o->GetNode()->GetNodeNeighbors();
          if(o_neighbors.size() > 1)
          {
            o_neighbor_1 = o_neighbors.at(0);
            o_neighbor_2 = o_neighbors.at(1);
          }
        }
        if(anomeric_o != NULL)
        {
          ///RULE1: anomeric to anomeric linkage
          if(((mono_linkages.at(0)).find(anomeric_linkage.str()) != std::string::npos) && ///this mono is attached to other mono through anomeric
              (mono_linkages.at(0).find(other_mono_anomeric_linkage_as_right_side.str()) != std::string::npos))///the other mono is only attached to this mono through anomeric
            isRoot = true;
          ///RULE2: Directed graph
          else if(((mono_linkages.at(0)).find(anomeric_linkage.str()) == std::string::npos) && ///this mono is not attached to other mono through anomeric
                  (mono_linkages.at(0).find(other_mono_anomeric_linkage_as_right_side.str()) != std::string::npos)) ///the other mono is attached to this mono through anomeric
          {
            isRoot = true;
            terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
          }
          ///RULE3: Terminal
          else if(o_neighbors.size() == 1) ///anomeric oxygen is not attached to anything else except the carbon of the ring
          {
            isRoot = true;
            terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
          }
          else if(o_neighbors.size() == 2 && (((o_neighbor_1->GetDescription().find("Het;") != std::string::npos) && (o_neighbor_2->GetDescription().find("Het;") == std::string::npos)) ||
                                              ((o_neighbor_2->GetDescription().find("Het;") != std::string::npos) && (o_neighbor_1->GetDescription().find("Het;") == std::string::npos))) )
          {
            ///anomeric oxygen is attached to protein
            isRoot = true;
            terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
            number_of_covalent_links++;
            if(terminal_residue_name.compare("NLN") != 0 && terminal_residue_name.compare("OLS") != 0 && terminal_residue_name.compare("OLT") != 0)
            {
              std::stringstream ss;
              ss << "Root anomeric atom is attached to a non-standard " << terminal_residue_name << " protein residue!";
              gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
              std::cout << ss.str() << "\n";
              terminal_residue_name = "";
            }
          }
          else if(dataset_residue_names.find(anomeric_o->GetResidue()->GetName()) != dataset_residue_names.end() ||
                  common_terminal_residues.find(anomeric_o->GetResidue()->GetName()) != common_terminal_residues.end())///mono is attached to a terminal through anomeric oxygen
          {
            terminal_residue_name = anomeric_o->GetResidue()->GetName();
            isRoot = true;
          }
          else if((terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms)).compare("") != 0)
            isRoot = true;
        }
        ///RULE2: Directed graph
        else if((mono_linkages.at(0).find(other_mono_anomeric_linkage_as_right_side.str()) != std::string::npos)) ///this mono doesn't have anomeric oxygen and the other mono is attached to this mono through anomeric
        {
          isRoot = true;
          terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
          number_of_probable_non_covalent_complexes++;
        }
      }
      else
      {
        MolecularModeling::Atom* anomeric_o = NULL;
        MolecularModeling::Atom* o_neighbor_1 = NULL;
        MolecularModeling::Atom* o_neighbor_2 = NULL;
        std::vector<MolecularModeling::Atom*>  o_neighbors = std::vector<MolecularModeling::Atom*> ();
        if(key->side_atoms_.at(0).at(1) != NULL)///Getting the information of anomeric oxygen's neighbors is needed for choosing the root
        {
          anomeric_o = key->side_atoms_.at(0).at(1);
          o_neighbors = anomeric_o->GetNode()->GetNodeNeighbors();
          if(o_neighbors.size() > 1)
          {
            o_neighbor_1 = o_neighbors.at(0);
            o_neighbor_2 = o_neighbors.at(1);
          }
        }
        //testing
        if(anomeric_o != NULL)
        {
          ///RULE1: anomeric to anomeric linkage
          for(unsigned int i = 0; i < values.size(); i++)
          {
            CheckLinkageNote(key, values.at(i), mono_linkages.at(i), checked_linkages);
            std::stringstream other_mono_anomeric_linkage_as_right_side;
            other_mono_anomeric_linkage_as_right_side << "-" << values.at(i)->cycle_atoms_.at(0)->GetId();///atom id on the right side of the linkage c-o-c
            if(((mono_linkages.at(i)).find(anomeric_linkage.str()) != std::string::npos) && ///this mono is attached to another mono through anomeric
                (mono_linkages.at(i).find(other_mono_anomeric_linkage_as_right_side.str()) != std::string::npos))///the other mono is attached to this mono through anomeric
            {
              isRoot = true;
              break;
            }
          }
          if(!isRoot) ///RULE2: Directed graph
          {
            for(unsigned int i = 0; i < values.size(); i++)
            {
              CheckLinkageNote(key, values.at(i), mono_linkages.at(i), checked_linkages);
              std::stringstream other_mono_anomeric_linkage_as_right_side;
              other_mono_anomeric_linkage_as_right_side << "-" << values.at(i)->cycle_atoms_.at(0)->GetId();///atom id on the right side of the linkage c-o-c
              if(((mono_linkages.at(i)).find(anomeric_linkage.str()) != std::string::npos)) ///this mono is attached to other mono through anomeric
              {
                isRoot = false;
                break;
              }
              else if((mono_linkages.at(i).find(other_mono_anomeric_linkage_as_right_side.str()) != std::string::npos)) ///the other mono is attached to this mono through anomeric
              {
                isRoot = true;
                terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
              }
            }
          }
          else if(!isRoot)///RULE3: Terminal
          {
            if(o_neighbors.size() == 1) ///anomeric oxygen is not attached to anything else, except the carbon of the ring
            {
              isRoot = true;
              terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
            }
            else if(o_neighbors.size() == 2 && (((o_neighbor_1->GetDescription().find("Het;") != std::string::npos) && (o_neighbor_2->GetDescription().find("Het;") == std::string::npos)) ||
                                                ((o_neighbor_2->GetDescription().find("Het;") != std::string::npos) && (o_neighbor_1->GetDescription().find("Het;") == std::string::npos))) )
            {
              ///anomeric oxygen is attached to protein
              isRoot = true;
              terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
              number_of_covalent_links++;
              if(terminal_residue_name.compare("NLN") != 0 && terminal_residue_name.compare("OLS") != 0 && terminal_residue_name.compare("OLT") != 0)
              {
                std::stringstream ss;
                ss << "Root anomeric atom is attached to a non-standard " << terminal_residue_name << " protein residue!";
                gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
                std::cout << ss.str() << "\n";
                terminal_residue_name = "";
              }
            }
            else if(dataset_residue_names.find(anomeric_o->GetResidue()->GetName()) != dataset_residue_names.end() ||
                    common_terminal_residues.find(anomeric_o->GetResidue()->GetName()) != common_terminal_residues.end())///mono is attached to a terminal through anomeric oxygen
            {
              terminal_residue_name = anomeric_o->GetResidue()->GetName();
              isRoot = true;
            }
            else if((terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms)).compare("") != 0)
              isRoot = true;
          }
        }
        //this mono doesn't have anomeric oxygen
        else ///RULE2: Directed graph
        {
          for(unsigned int i = 0; i < values.size(); i++)
          {
            CheckLinkageNote(key, values.at(i), mono_linkages.at(i), checked_linkages);
            std::vector<std::string> other_mono_linkage = monos_table_linkages[values.at(i)];
            std::stringstream other_mono_anomeric_linkage;
            other_mono_anomeric_linkage << values.at(i)->cycle_atoms_.at(0)->GetId() << "-";///atom id on the left side of the linkage c-o-c
            if((other_mono_linkage.at(0).find(other_mono_anomeric_linkage.str()) != std::string::npos)) ///the other mono is attached to this mono through anomeric
            {
              isRoot = true;
              terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
              number_of_probable_non_covalent_complexes++;
              break;
            }
          }
        }
      }
      if(isRoot)
      {
        Glycan::Oligosaccharide* oligo = new Glycan::Oligosaccharide(this);
        CalculateOligosaccharideBFactor(oligo,  oligo->mono_nodes_);
        BuildOligosaccharideTreeStructure(key, values, oligo, visited_monos, monos_table, monos_table_linkages, visited_linkages);
        oligo->terminal_ = terminal_residue_name;
        oligosaccharides.push_back(oligo);
      }
    }
  }
  // gmml::log(__LINE__, __FILE__,  gmml::INF, " End for loop ..." );
  // gmml::log(__LINE__, __FILE__,  gmml::INF, " Another for loop ..." );
  for(std::map<Glycan::Monosaccharide*, std::vector<Glycan::Monosaccharide*> >::iterator it = monos_table.begin(); it != monos_table.end(); it++)
  {
    Glycan::Monosaccharide* key = (*it).first;
    std::vector<Glycan::Monosaccharide*> values = (*it).second;
    if(values.size() > 1)
    {
      std::vector<std::string> visited_linkages = std::vector<std::string>();
      if(find(visited_monos.begin(), visited_monos.end(), key->mono_id_) == visited_monos.end())///if the mono is not visited
      {
        std::vector<std::string> mono_linkages = monos_table_linkages[key];
        std::stringstream anomeric_linkage;
        anomeric_linkage << key->cycle_atoms_.at(0)->GetId() << "-";
        for(std::vector<std::string>::iterator it1 = mono_linkages.begin(); it1 != mono_linkages.end(); it1++)
        {
          if((*it1).find(anomeric_linkage.str()) != std::string::npos)///mono is attached to another mono through anomeric
          {
            Glycan::Oligosaccharide* oligo = new Glycan::Oligosaccharide(this);
            CalculateOligosaccharideBFactor(oligo,  oligo->mono_nodes_);
            BuildOligosaccharideTreeStructure(key, values, oligo, visited_monos, monos_table, monos_table_linkages, visited_linkages);
            oligosaccharides.push_back(oligo);
            break;
          }
        }
      }
    }
  }
    // gmml::log(__LINE__, __FILE__,  gmml::INF, "Done with that too ..." );
    // return oligosaccharides;
}

std::vector<Glycan::Oligosaccharide*> Assembly::createOligosaccharides(std::vector<Glycan::Monosaccharide*> detected_monos)
{
  // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
  std::vector<Glycan::Oligosaccharide*> detected_oligos;
  for(std::vector<Glycan::Monosaccharide*>::iterator it = detected_monos.begin(); it != detected_monos.end(); it++)
  {

    Glycan::Monosaccharide* this_mono = *it;
    // gmml::log(__LINE__, __FILE__,  gmml::INF, this_mono->sugar_name_.monosaccharide_short_name_);
    if(this_mono->mono_neighbors_.empty())
    {
      this_mono->is_root_ = true;
      this_mono->is_visited_ = false;
    }
    else
    {
      this_mono->is_root_ = true;
      for(std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> >::iterator monoNeighbor = this_mono->mono_neighbors_.begin(); monoNeighbor != this_mono->mono_neighbors_.end(); monoNeighbor++)
      {
        Glycan::GlycosidicLinkage* thisLinkage = (*monoNeighbor).first;
        Glycan::Monosaccharide* thisNeighbor= (*monoNeighbor).second;
        std::stringstream ss;
        if(thisLinkage->reducing_mono_ != NULL)
        {
          ss << this_mono->cycle_atoms_[0]->GetResidue()->GetId() << " is being compared to " << thisLinkage->non_reducing_mono_->cycle_atoms_[0]->GetResidue()->GetId() << " and the linkage type is " << thisLinkage->inverse_linkage_type_;
          // gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
          if(this_mono->cycle_atoms_[0]->GetResidue()->GetId() == thisLinkage->non_reducing_mono_->cycle_atoms_[0]->GetResidue()->GetId())
          {//if this mono has a mono neighbor at the anomeric carbon, it can't be the root (attached to terminal)
            // gmml::log(__LINE__, __FILE__,  gmml::INF, "This mono is the non reducing mono");
            this_mono->is_root_ = false;
          }
          else
          {
            // gmml::log(__LINE__, __FILE__,  gmml::INF, "This mono is the reducing mono");
          }
        }
        else
        {
          ss << this_mono->cycle_atoms_[0]->GetResidue()->GetId() << " and " << thisNeighbor->cycle_atoms_[0]->GetResidue()->GetId() << " have a linkage of " << thisLinkage->linkage_type_;
          // gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
          // gmml::log(__LINE__, __FILE__,  gmml::INF, "This linkage is anomeric-anomeric");
          //TODO handle anomeric-anomeric; they both think they are the root and you get stuck in an infinite loop
          if((thisNeighbor->mono_neighbors_.size() == 1) && (this_mono->mono_neighbors_.size() == 1))//just a disaccharide
          {
            if(thisNeighbor->is_root_)
            {
              this_mono->is_root_ = false;
            }
            if(thisLinkage->linkage_type_ == "1-2")
            {
              // gmml::log(__LINE__, __FILE__,  gmml::INF, "This is a 1-2 anomeric linkage, so the other mono is the root");
              this_mono->is_root_ = false;
            }
          }
          else if((thisNeighbor->mono_neighbors_.size() == 1) && (this_mono->mono_neighbors_.size() > 1))//other neighbor is the 'terminal'
          {
            this_mono->is_root_ = false;
          }
        }
      }
    }
    if((this_mono->is_root_) && (!this_mono->is_visited_))
    {
      // gmml::log(__LINE__, __FILE__,  gmml::INF, "This mono is the root");
      Glycan::Oligosaccharide* this_Oligo = new Glycan::Oligosaccharide(this);
      // gmml::log(__LINE__, __FILE__, gmml::INF, this_mono->sugar_name_.monosaccharide_short_name_);
      this_Oligo->traverseGraph(this_mono, this_Oligo);
      this_Oligo->reindexRGroups(this_Oligo);
      this_Oligo->indexMonosaccharides();
      detected_oligos.push_back(this_Oligo);
      std::string iupac = "Oligo IUPAC Name: " + this_Oligo->IUPAC_name_;
      // gmml::log(__LINE__, __FILE__, gmml::INF, iupac);
      std::string oligoname =  "Oligo Name: " +this_Oligo->oligosaccharide_name_;
      // gmml::log(__LINE__, __FILE__, gmml::INF, oligoname);
    }
  }
  for(std::vector<Glycan::Monosaccharide*>::iterator it = detected_monos.begin(); it != detected_monos.end(); it++)
  {
    Glycan::Monosaccharide* this_mono = *it;
    if(this_mono->oligo_parent_ == NULL)
    {
      this_mono->is_root_ = true;
      Glycan::Oligosaccharide* this_Oligo = new Glycan::Oligosaccharide(this);
      //TODO the below function does not work correctly for cyclic oligosaccharides; write a new one or fix it
      this_Oligo->traverseGraph(this_mono, this_Oligo);
      this_Oligo->reindexRGroups(this_Oligo);
      this_Oligo->indexMonosaccharides();
      detected_oligos.push_back(this_Oligo);
    }
  }
  return detected_oligos;
}

std::vector<Glycan::Oligosaccharide*> Assembly::ExtractOligosaccharides( std::vector<Glycan::Monosaccharide*> monos,
                                                            gmml::ResidueNameMap dataset_residue_names,
                                                            int& number_of_covalent_links,
                                                            int& number_of_probable_non_covalent_complexes ) {
    std::string terminal_residue_name = "";
    gmml::ResidueNameMap common_terminal_residues = gmml::InitializeCommonTerminalResidueMap();
    std::map<Glycan::Monosaccharide*, std::vector<Glycan::Monosaccharide*> > monos_table = std::map<Glycan::Monosaccharide*, std::vector<Glycan::Monosaccharide*> >();
    std::map<Glycan::Monosaccharide*, std::vector<std::string> > monos_table_linkages = std::map<Glycan::Monosaccharide*, std::vector<std::string> >();
    // gmml::log(__LINE__, __FILE__,  gmml::INF, " Start iterating on list ..." );
    ///Iterating on list of monos to check if there is a connection to another mono in the list
    for(std::vector<Glycan::Monosaccharide*>::iterator it = monos.begin(); it != monos.end(); it++)
    {
        Glycan::Monosaccharide* mono1 = (*it);

        monos_table[mono1] = std::vector<Glycan::Monosaccharide*>();
        monos_table_linkages[mono1] = std::vector<std::string>();

        for(std::vector<AtomVector>::iterator it1 = mono1->side_atoms_.begin(); it1 != mono1->side_atoms_.end(); it1++) ///iterate on side atoms
        {
            int index = distance(mono1->side_atoms_.begin(), it1);
            MolecularModeling::AtomVector sides = (*it1);
            std::map<MolecularModeling::Atom*, MolecularModeling::Atom*> target_parent_map = std::map<MolecularModeling::Atom*, MolecularModeling::Atom*>(); /// A map of target atom to it's parent atom. Target atom is a non ring oxygen or nitrogen

            if(it1 == mono1->side_atoms_.begin())///side atoms of anomeric
            {
                if(sides.at(1) != NULL)
                    target_parent_map[sides.at(1)] = mono1->cycle_atoms_.at(0);
            }
            else if(it1 == mono1->side_atoms_.end() - 1) ///side atoms of last carbon of the ring
            {
                for(MolecularModeling::AtomVector::iterator last_c_side_it = sides.begin(); last_c_side_it != sides.end(); last_c_side_it++)
                {
                    MolecularModeling::Atom* side_of_last_carbon = (*last_c_side_it);
                    if(side_of_last_carbon != NULL)
                    {
                        MolecularModeling::AtomVector last_c_side_neighbors = side_of_last_carbon->GetNode()->GetNodeNeighbors();
                        for(MolecularModeling::AtomVector::iterator it2 = last_c_side_neighbors.begin(); it2 != last_c_side_neighbors.end(); it2++)
                        {
                            if((*it2)->GetId().at(0) == 'O' || (*it2)->GetId().at(0) == 'N')
                            {
                                target_parent_map[(*it2)] = side_of_last_carbon;
                                break;
                            }
                        }
                    }
                }
            }
            else
            {
                if(sides.at(1) != NULL)
                    target_parent_map[sides.at(1)] = mono1->cycle_atoms_.at(index);///index 1 of each side is for non-carbon side atoms in the std::vector<AtomVector> structure
            }
            ///Examine neighbors of each target atom to check if they can be found in other monos side/ring atoms
            for(std::map<MolecularModeling::Atom*, MolecularModeling::Atom*>::iterator map_it = target_parent_map.begin(); map_it != target_parent_map.end(); map_it++)
            {
                bool found_in_other_mono = false;
                MolecularModeling::Atom* target = (*map_it).first;
                MolecularModeling::Atom* target_parent = (*map_it).second;
                MolecularModeling::AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
                for(MolecularModeling::AtomVector::iterator it2 = t_neighbors.begin(); it2 != t_neighbors.end(); it2++)
                {
                    MolecularModeling::Atom* t_neighbor = (*it2);
                    if(t_neighbor->GetId().compare(target_parent->GetId()) != 0)///making sure neighbor is not the parent of target atom
                    {
                        for(std::vector<Glycan::Monosaccharide*>::iterator it3 = monos.begin(); it3 != monos.end(); it3++)
                        {
                            if(it3 != it)///Cheking monos other than the current mono
                            {
                                Glycan::Monosaccharide* mono2 = (*it3);
                                MolecularModeling::AtomVector mono2_sides = mono2->side_atoms_.at(mono2->side_atoms_.size() - 1); ///side of last ring carbon

                                bool found_in_side = false;
                                for(MolecularModeling::AtomVector::iterator mono2_last_c_side_it = mono2_sides.begin(); mono2_last_c_side_it != mono2_sides.end(); mono2_last_c_side_it++)
                                {
                                    MolecularModeling::Atom* mono2_last_c_side = (*mono2_last_c_side_it);
                                    if(mono2_last_c_side != NULL)
                                    {
                                        if(t_neighbor->GetId().compare(mono2_last_c_side->GetId()) == 0) ///target atom has been attached to another cycle's side atom
                                            found_in_side = true;
                                    }
                                }
                                if(found_in_side || mono2->cycle_atoms_str_.find(t_neighbor->GetId()) != std::string::npos) //if target's neighbor found in another mono's side or ring atoms
                                {
                                    found_in_other_mono = true;
                                    monos_table[mono1].push_back(mono2);

                                    std::string mono1_carbon = target_parent->GetId();
                                    std::string mono1_name = "";
                                    std::string mono2_carbon = t_neighbor->GetId();
                                    std::string mono2_name = "";
                                    if(mono1->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                                        mono1_name = mono1->sugar_name_.monosaccharide_short_name_;
                                    else
                                        mono1_name = mono1->sugar_name_.monosaccharide_stereochemistry_short_name_;

                                    if(mono2->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                                        mono2_name = mono2->sugar_name_.monosaccharide_short_name_;
                                    else
                                        mono2_name = mono2->sugar_name_.monosaccharide_stereochemistry_short_name_;

                                    std::stringstream linkage;
                                    linkage << mono1_carbon << "-" << target->GetId() << "-" << mono2_carbon;
                                    monos_table_linkages[mono1].push_back(linkage.str());
                                    break;
                                }
                            }
                        }
                    }
                }
                if(found_in_other_mono)
                    break;
            }
        }
    }
    // gmml::log(__LINE__, __FILE__,  gmml::INF, " Done iterating list" );
    std::vector<int> visited_monos = std::vector<int>();
    std::vector<Glycan::Oligosaccharide*> oligosaccharides = std::vector<Glycan::Oligosaccharide*>();

    std::vector<std::string> checked_linkages = std::vector<std::string>();
    // gmml::log(__LINE__, __FILE__,  gmml::INF, " Start for loop ..." );
    for(std::map<Glycan::Monosaccharide*, std::vector<Glycan::Monosaccharide*> >::iterator it = monos_table.begin(); it != monos_table.end(); it++)
    {
        Glycan::Monosaccharide* key = (*it).first;
        std::vector<Glycan::Monosaccharide*> values = (*it).second;

        std::vector<std::string> visited_linkages = std::vector<std::string>();
        if(find(visited_monos.begin(), visited_monos.end(), key->mono_id_) == visited_monos.end())///if the mono is not visited
        {
            bool isRoot = false;
            std::stringstream anomeric_linkage;
            anomeric_linkage << key->cycle_atoms_.at(0)->GetId() << "-";
            std::vector<std::string> mono_linkages = monos_table_linkages[key];
            MolecularModeling::AtomVector terminal_atoms = MolecularModeling::AtomVector();
            if(values.size() == 0) ///mono is not attached to any other mono
            {
                MolecularModeling::Atom* anomeric_o = NULL;
                if(key->side_atoms_.at(0).at(1) != NULL)
                    anomeric_o = key->side_atoms_.at(0).at(1);
                if(anomeric_o != NULL)
                {
                    if(dataset_residue_names.find(anomeric_o->GetResidue()->GetName()) != dataset_residue_names.end() ||
                            common_terminal_residues.find(anomeric_o->GetResidue()->GetName()) != common_terminal_residues.end())///check if there is any terminal
                    {
                        terminal_residue_name = anomeric_o->GetResidue()->GetName();
                    }
                    else
                    {
                        terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                    }
                }
                isRoot = true;
            }
            else if (values.size() == 1 && ///mono is attached to one other mono
                     (monos_table_linkages[values.at(0)].size() == 1))///the other mono is only attached to this mono
            {
                ///CHECKING LINKAGE ISSUES, e.g. C1-O3-C4 is an issue
                CheckLinkageNote(key, values.at(0), mono_linkages.at(0), checked_linkages);
                std::stringstream other_mono_anomeric_linkage_as_right_side;
                other_mono_anomeric_linkage_as_right_side << "-" << values.at(0)->cycle_atoms_.at(0)->GetId();///atom id on the right side of the linkage c-o-c

                MolecularModeling::Atom* anomeric_o = NULL;
                MolecularModeling::Atom* o_neighbor_1 = NULL;
                MolecularModeling::Atom* o_neighbor_2 = NULL;
                MolecularModeling::AtomVector o_neighbors = MolecularModeling::AtomVector();
                if(key->side_atoms_.at(0).at(1) != NULL)///Getting the information of anomeric oxygen's neighbors is needed for choosing the root
                {
                    anomeric_o = key->side_atoms_.at(0).at(1);
                    o_neighbors = anomeric_o->GetNode()->GetNodeNeighbors();
                    if(o_neighbors.size() > 1)
                    {
                        o_neighbor_1 = o_neighbors.at(0);
                        o_neighbor_2 = o_neighbors.at(1);
                    }
                }
                if(anomeric_o != NULL)
                {
                    ///RULE1: anomeric to anomeric linkage
                    if(((mono_linkages.at(0)).find(anomeric_linkage.str()) != std::string::npos) && ///this mono is attached to other mono through anomeric
                            (mono_linkages.at(0).find(other_mono_anomeric_linkage_as_right_side.str()) != std::string::npos))///the other mono is only attached to this mono through anomeric
                        isRoot = true;
                    ///RULE2: Directed graph
                    else if(((mono_linkages.at(0)).find(anomeric_linkage.str()) == std::string::npos) && ///this mono is not attached to other mono through anomeric
                            (mono_linkages.at(0).find(other_mono_anomeric_linkage_as_right_side.str()) != std::string::npos)) ///the other mono is attached to this mono through anomeric
                    {
                        isRoot = true;
                        terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                    }
                    ///RULE3: Terminal
                    else if(o_neighbors.size() == 1) ///anomeric oxygen is not attached to anything else except the carbon of the ring
                    {
                        isRoot = true;
                        terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                    }
                    else if(o_neighbors.size() == 2 && (((o_neighbor_1->GetDescription().find("Het;") != std::string::npos) && (o_neighbor_2->GetDescription().find("Het;") == std::string::npos)) ||
                                                        ((o_neighbor_2->GetDescription().find("Het;") != std::string::npos) && (o_neighbor_1->GetDescription().find("Het;") == std::string::npos))) )
                    {
                        ///anomeric oxygen is attached to protein
                        isRoot = true;
                        terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                        number_of_covalent_links++;
                        if(terminal_residue_name.compare("NLN") != 0 && terminal_residue_name.compare("OLS") != 0 && terminal_residue_name.compare("OLT") != 0)
                        {
                            std::stringstream ss;
                            ss << "Root anomeric atom is attached to a non-standard " << terminal_residue_name << " protein residue!";
                            gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
                            std::cout << ss.str() << "\n";
                            terminal_residue_name = "";
                        }
                    }
                    else if(dataset_residue_names.find(anomeric_o->GetResidue()->GetName()) != dataset_residue_names.end() ||
                            common_terminal_residues.find(anomeric_o->GetResidue()->GetName()) != common_terminal_residues.end())///mono is attached to a terminal through anomeric oxygen
                    {
                        terminal_residue_name = anomeric_o->GetResidue()->GetName();
                        isRoot = true;
                    }
                    else if((terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms)).compare("") != 0)
                        isRoot = true;
                }
                ///RULE2: Directed graph
                else if((mono_linkages.at(0).find(other_mono_anomeric_linkage_as_right_side.str()) != std::string::npos)) ///this mono doesn't have anomeric oxygen and the other mono is attached to this mono through anomeric
                {
                    isRoot = true;
                    terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                    number_of_probable_non_covalent_complexes++;
                }
            }
            else
            {
                MolecularModeling::Atom* anomeric_o = NULL;
                MolecularModeling::Atom* o_neighbor_1 = NULL;
                MolecularModeling::Atom* o_neighbor_2 = NULL;
                MolecularModeling::AtomVector o_neighbors = MolecularModeling::AtomVector();
                if(key->side_atoms_.at(0).at(1) != NULL)///Getting the information of anomeric oxygen's neighbors is needed for choosing the root
                {
                    anomeric_o = key->side_atoms_.at(0).at(1);
                    o_neighbors = anomeric_o->GetNode()->GetNodeNeighbors();
                    if(o_neighbors.size() > 1)
                    {
                        o_neighbor_1 = o_neighbors.at(0);
                        o_neighbor_2 = o_neighbors.at(1);
                    }
                }
		//testing
                if(anomeric_o != NULL)
                {
                    ///RULE1: anomeric to anomeric linkage
                    for(unsigned int i = 0; i < values.size(); i++)
                    {
                        CheckLinkageNote(key, values.at(i), mono_linkages.at(i), checked_linkages);
                        std::stringstream other_mono_anomeric_linkage_as_right_side;
                        other_mono_anomeric_linkage_as_right_side << "-" << values.at(i)->cycle_atoms_.at(0)->GetId();///atom id on the right side of the linkage c-o-c

                        if(((mono_linkages.at(i)).find(anomeric_linkage.str()) != std::string::npos) && ///this mono is attached to another mono through anomeric
                                (mono_linkages.at(i).find(other_mono_anomeric_linkage_as_right_side.str()) != std::string::npos))///the other mono is attached to this mono through anomeric
                        {
                            isRoot = true;
                            break;
                        }
                    }
                    if(!isRoot) ///RULE2: Directed graph
                    {
                        for(unsigned int i = 0; i < values.size(); i++)
                        {
                            CheckLinkageNote(key, values.at(i), mono_linkages.at(i), checked_linkages);
                            std::stringstream other_mono_anomeric_linkage_as_right_side;
                            other_mono_anomeric_linkage_as_right_side << "-" << values.at(i)->cycle_atoms_.at(0)->GetId();///atom id on the right side of the linkage c-o-c
                            if(((mono_linkages.at(i)).find(anomeric_linkage.str()) != std::string::npos)) ///this mono is attached to other mono through anomeric
                            {
                                isRoot = false;
                                break;
                            }
                            else if((mono_linkages.at(i).find(other_mono_anomeric_linkage_as_right_side.str()) != std::string::npos)) ///the other mono is attached to this mono through anomeric
                            {
                                isRoot = true;
                                terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                            }
                        }
                    }
                    else if(!isRoot)///RULE3: Terminal
                    {
                        if(o_neighbors.size() == 1) ///anomeric oxygen is not attached to anything else, except the carbon of the ring
                        {
                            isRoot = true;
                            terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                        }
                        else if(o_neighbors.size() == 2 && (((o_neighbor_1->GetDescription().find("Het;") != std::string::npos) && (o_neighbor_2->GetDescription().find("Het;") == std::string::npos)) ||
                                                            ((o_neighbor_2->GetDescription().find("Het;") != std::string::npos) && (o_neighbor_1->GetDescription().find("Het;") == std::string::npos))) )
                        {
                            ///anomeric oxygen is attached to protein
                            isRoot = true;
                            terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                            number_of_covalent_links++;
                            if(terminal_residue_name.compare("NLN") != 0 && terminal_residue_name.compare("OLS") != 0 && terminal_residue_name.compare("OLT") != 0)
                            {
                                std::stringstream ss;
                                ss << "Root anomeric atom is attached to a non-standard " << terminal_residue_name << " protein residue!";
                                gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
                                std::cout << ss.str() << "\n";
                                terminal_residue_name = "";
                            }
                        }
                        else if(dataset_residue_names.find(anomeric_o->GetResidue()->GetName()) != dataset_residue_names.end() ||
                                common_terminal_residues.find(anomeric_o->GetResidue()->GetName()) != common_terminal_residues.end())///mono is attached to a terminal through anomeric oxygen
                        {
                            terminal_residue_name = anomeric_o->GetResidue()->GetName();
                            isRoot = true;
                        }
                        else if((terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms)).compare("") != 0)
                            isRoot = true;
                    }
                }
                //this mono doesn't have anomeric oxygen
                else ///RULE2: Directed graph
                {
                    for(unsigned int i = 0; i < values.size(); i++)
                    {
                        CheckLinkageNote(key, values.at(i), mono_linkages.at(i), checked_linkages);
                        std::vector<std::string> other_mono_linkage = monos_table_linkages[values.at(i)];
                        std::stringstream other_mono_anomeric_linkage;
                        other_mono_anomeric_linkage << values.at(i)->cycle_atoms_.at(0)->GetId() << "-";///atom id on the left side of the linkage c-o-c
                        if((other_mono_linkage.at(0).find(other_mono_anomeric_linkage.str()) != std::string::npos)) ///the other mono is attached to this mono through anomeric
                        {
                            isRoot = true;
                            terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                            number_of_probable_non_covalent_complexes++;
                            break;
                        }
		    }
                }
            }
            if(key->is_root_)
            {
                Glycan::Oligosaccharide* oligo = key->oligo_parent_;
                CalculateOligosaccharideBFactor(oligo, oligo->mono_nodes_);
                BuildOligosaccharideTreeStructure(key, values, oligo, visited_monos, monos_table, monos_table_linkages, visited_linkages);
                oligo->terminal_ = terminal_residue_name;
                oligosaccharides.push_back(oligo);
            }
        }
    }
// gmml::log(__LINE__, __FILE__,  gmml::INF, " End for loop ..." );
// gmml::log(__LINE__, __FILE__,  gmml::INF, " Another for loop ..." );
    for(std::map<Glycan::Monosaccharide*, std::vector<Glycan::Monosaccharide*> >::iterator it = monos_table.begin(); it != monos_table.end(); it++)
    {
        Glycan::Monosaccharide* key = (*it).first;
        std::vector<Glycan::Monosaccharide*> values = (*it).second;
        if(values.size() > 1)
        {
            std::vector<std::string> visited_linkages = std::vector<std::string>();
            if(find(visited_monos.begin(), visited_monos.end(), key->mono_id_) == visited_monos.end())///if the mono is not visited
            {
                std::vector<std::string> mono_linkages = monos_table_linkages[key];
                std::stringstream anomeric_linkage;
                anomeric_linkage << key->cycle_atoms_.at(0)->GetId() << "-";

                for(std::vector<std::string>::iterator it1 = mono_linkages.begin(); it1 != mono_linkages.end(); it1++)
                {
                    if((*it1).find(anomeric_linkage.str()) != std::string::npos)///mono is attached to another mono through anomeric
                    {
                        Glycan::Oligosaccharide* oligo = new Glycan::Oligosaccharide(this);
                        CalculateOligosaccharideBFactor(oligo, oligo->mono_nodes_);
                        BuildOligosaccharideTreeStructure(key, values, oligo, visited_monos, monos_table, monos_table_linkages, visited_linkages);
                        oligosaccharides.push_back(oligo);
                        break;
                    }
                }
            }
        }
    }
    // gmml::log(__LINE__, __FILE__,  gmml::INF, "Done with that too ..." );
    return oligosaccharides;
}//End ExtractOligosaccharides

void Assembly::UpdateMonosaccharides2Residues(std::vector<Glycan::Monosaccharide*>& monos)
{
  /*An old PDB feature puts multiple monosaccharides in a sigle residue.This is incompatible with downstream GMML codes.
    Solution: Each monosaccharide becomes an residue, replacing the corresponding old residue. Complete side group atoms were determined in function SetCompleteSideGroupAtoms().*/
  std::vector<Residue*> OldResidue2BeErasedFromAssembly = std::vector<Residue*>();
  int mono_index = 0;
  for (std::vector<Glycan::Monosaccharide*>::iterator it= monos.begin(); it != monos.end(); it++){
      mono_index++;
      std::stringstream mono_index_stream;
      mono_index_stream << mono_index;
      std::string mono_index_str = mono_index_stream.str();

      Glycan::Monosaccharide* mono = *it;
      Residue * NewResidueForThisMonosaccharide = new Residue();
      Residue * OldResidueForThisMonosaccharide = mono->cycle_atoms_.at(0)->GetResidue();
      NewResidueForThisMonosaccharide-> SetAssembly(this);
      std::vector<Residue*> AllResiduesInAssembly= this -> GetResidues();
      this->InsertResidue(OldResidueForThisMonosaccharide ,NewResidueForThisMonosaccharide);
      NewResidueForThisMonosaccharide-> SetName(OldResidueForThisMonosaccharide->GetName() );
      NewResidueForThisMonosaccharide-> SetId(OldResidueForThisMonosaccharide->GetId() + "-mon-" + mono_index_str);

      if ( std::find(OldResidue2BeErasedFromAssembly.begin(),OldResidue2BeErasedFromAssembly.end(),OldResidueForThisMonosaccharide) == OldResidue2BeErasedFromAssembly.end() ){
          OldResidue2BeErasedFromAssembly.push_back(OldResidueForThisMonosaccharide);
      }


      std::vector<Atom*> AllAtomsInThisMonosaccharide = std::vector<Atom*>();
      for (std::vector<Atom*>::iterator it2= mono->cycle_atoms_.begin(); it2!= mono->cycle_atoms_.end(); it2++){ //Add cycle atoms to monosaccharide
	  MolecularModeling::Atom* cycle_atom = *it2;
          AllAtomsInThisMonosaccharide.push_back(cycle_atom);
	  MolecularModeling::AtomVector cycle_atom_neighbors = cycle_atom->GetNode()->GetNodeNeighbors();
	  for (MolecularModeling::AtomVector::iterator neighbor_it = cycle_atom_neighbors.begin(); neighbor_it != cycle_atom_neighbors.end(); neighbor_it++){
	      MolecularModeling::Atom* neighbor = *neighbor_it;
	      //If a neighbor of a ring atom is neither on a cycle or a sidechain, and is a hydrogen, then it is the ring hydrogen. Add it to new residue as well. 
	      if (!neighbor->GetIsCycle() && !neighbor->GetIsSideChain() && neighbor->GetElementSymbol() == "H"){ 
	          AllAtomsInThisMonosaccharide.push_back(neighbor);
	      }
	  }
      }

      for (std::vector<std::vector<Atom*> >::iterator it2= mono->side_atoms_.begin(); it2!= mono->side_atoms_.end(); it2++){ //Add side chain atoms to monosaccharide
	  if ( !(*it2).empty() ){
	      for (std::vector<Atom*>::iterator it3= (*it2).begin(); it3 != (*it2).end(); it3++){
		  if ( (*it3) != NULL){
                      AllAtomsInThisMonosaccharide.push_back(*it3);
		  }
	      }
	  }
      }

      NewResidueForThisMonosaccharide-> SetAtoms(AllAtomsInThisMonosaccharide);

      for (std::vector<Atom*>::iterator it2= AllAtomsInThisMonosaccharide.begin(); it2!= AllAtomsInThisMonosaccharide.end(); it2++){
          (*it2)-> SetResidue(NewResidueForThisMonosaccharide);
      }

  }//for


  for (std::vector<Residue*>::iterator it= OldResidue2BeErasedFromAssembly.begin(); it!= OldResidue2BeErasedFromAssembly.end(); it++){ ////remove old residues
      this->RemoveResidue(*it);
  }

}//UpdateMonosaccharides2Residues

std::string Assembly::CheckOMETerminal(MolecularModeling::Atom* target, MolecularModeling::AtomVector& terminal_atoms)
{
    terminal_atoms = MolecularModeling::AtomVector();
    MolecularModeling::AtomVector atoms_1 = MolecularModeling::AtomVector();
    MolecularModeling::AtomVector atoms_2 = MolecularModeling::AtomVector();
    std::stringstream pattern;
    pattern << "O";
    atoms_1.push_back(target);
    atoms_2.push_back(target);
    MolecularModeling::Atom* C1 = NULL;
    MolecularModeling::Atom* C2 = NULL;
    MolecularModeling::AtomVector o_neighbors = target->GetNode()->GetNodeNeighbors();
    for(MolecularModeling::AtomVector::iterator it = o_neighbors.begin(); it != o_neighbors.end(); it++)
    {
        MolecularModeling::Atom* o_neighbor = (*it);

        if(o_neighbor->GetName().at(0) == 'C' && C1 == NULL && C2 == NULL)
            C1 = o_neighbor;
        else if(o_neighbor->GetName().at(0) == 'C' && C1 != NULL && C2 == NULL)
            C2 = o_neighbor;
    }
    std::stringstream temp;
    if(C1 != NULL)
    {
        temp << pattern.str() << "-C";
        atoms_1.push_back(C1);
        MolecularModeling::AtomVector c1_neighbors = C1->GetNode()->GetNodeNeighbors();
        for(MolecularModeling::AtomVector::iterator it = c1_neighbors.begin(); it != c1_neighbors.end(); it++)
        {
            MolecularModeling::Atom* c1_neighbor = (*it);
            if(c1_neighbor->GetId().compare(target->GetId()) != 0)
            {
                temp << c1_neighbor->GetName().at(0);
                atoms_1.push_back(c1_neighbor);
            }
        }
    }
    if(temp.str().compare("O-C") == 0 || temp.str().compare("O-CHHH") == 0)
    {
        terminal_atoms = atoms_1;
        return "OME";
    }
    if(C2 != NULL)
    {
        pattern << "-C";
        atoms_2.push_back(C2);
        MolecularModeling::AtomVector c2_neighbors = C2->GetNode()->GetNodeNeighbors();
        for(MolecularModeling::AtomVector::iterator it = c2_neighbors.begin(); it != c2_neighbors.end(); it++)
        {
            MolecularModeling::Atom* c2_neighbor = (*it);
            if(c2_neighbor->GetId().compare(target->GetId()) != 0)
            {
                pattern << c2_neighbor->GetName().at(0);
                atoms_2.push_back(c2_neighbor);
            }
        }
        if(pattern.str().compare("O-C") == 0 || pattern.str().compare("O-CHHH") == 0)
        {
            terminal_atoms = atoms_2;
            return "OME";
        }
    }
    return "";
}

std::string Assembly::CheckROHTerminal(MolecularModeling::Atom* target, MolecularModeling::AtomVector& terminal_atoms)
{
    terminal_atoms = MolecularModeling::AtomVector();
    MolecularModeling::AtomVector o_neighbors = target->GetNode()->GetNodeNeighbors();
    if(o_neighbors.size() == 1)
    {
        terminal_atoms.push_back(target);
        return "ROH";
    }
    //else if (o_neighbors.size() > 1)
    else if (o_neighbors.size() > 1 && o_neighbors.size()== 2 )
    //ROH oxygen has up to two neighbors. It's safer to limit number of neighbors to two.
    {
        if((o_neighbors.at(0)->GetName().at(0) == 'H' && o_neighbors.at(1)->GetName().at(0) != 'H'))
        {
            terminal_atoms.push_back(target);
            terminal_atoms.push_back(o_neighbors.at(0));
            return "ROH";
        }
        else if(o_neighbors.at(1)->GetName().at(0) == 'H' && o_neighbors.at(0)->GetName().at(0) != 'H')
        {
            terminal_atoms.push_back(target);
            terminal_atoms.push_back(o_neighbors.at(1));
            return "ROH";
        }
    }
    return "";
}

std::string Assembly::CheckTBTTerminal(MolecularModeling::Atom *target, MolecularModeling::AtomVector& terminal_atoms)
{
    terminal_atoms = MolecularModeling::AtomVector();
    MolecularModeling::AtomVector atoms_1 = MolecularModeling::AtomVector();
    MolecularModeling::AtomVector atoms_2 = MolecularModeling::AtomVector();
    std::stringstream pattern;
    pattern << "O";
    atoms_1.push_back(target);
    atoms_2.push_back(target);
    MolecularModeling::Atom* C1 = NULL;
    MolecularModeling::Atom* C2 = NULL;
    MolecularModeling::AtomVector o_neighbors = target->GetNode()->GetNodeNeighbors();
    for(MolecularModeling::AtomVector::iterator it = o_neighbors.begin(); it != o_neighbors.end(); it++)
    {
        MolecularModeling::Atom* o_neighbor = (*it);
        if(o_neighbor->GetName().at(0) == 'C' && C1 == NULL && C2 == NULL)
            C1 = o_neighbor;
        else if(o_neighbor->GetName().at(0) == 'C' && C1 != NULL && C2 == NULL)
            C2 = o_neighbor;
    }
    std::stringstream temp;
    if(C1 != NULL)
    {
        MolecularModeling::Atom* C1C1 = NULL;
        MolecularModeling::Atom* C1C2 = NULL;
        MolecularModeling::Atom* C1C3 = NULL;
        temp << pattern.str() << "-" << "C";
        atoms_1.push_back(C1);
        MolecularModeling::AtomVector c1_neighbors = C1->GetNode()->GetNodeNeighbors();
        for(MolecularModeling::AtomVector::iterator it = c1_neighbors.begin(); it != c1_neighbors.end(); it++)
        {
            MolecularModeling::Atom* c1_neighbor = (*it);
            if(c1_neighbor->GetId().compare(target->GetId()) != 0 && c1_neighbor->GetName().at(0) == 'C' && C1C1 == NULL)
                C1C1 = c1_neighbor;
            else if(c1_neighbor->GetId().compare(target->GetId()) != 0 && c1_neighbor->GetName().at(0) == 'C' && C1C1 != NULL && C1C2 == NULL)
                C1C2 = c1_neighbor;
            else if(c1_neighbor->GetId().compare(target->GetId()) != 0 && c1_neighbor->GetName().at(0) == 'C' && C1C1 != NULL && C1C2 != NULL && C1C3 == NULL)
                C1C3 = c1_neighbor;
        }
        if(C1C1 != NULL && C1C2 != NULL && C1C3 != NULL)
        {
            temp << "C";
            atoms_1.push_back(C1C1);
            MolecularModeling::AtomVector c1c1_neighbors = C1C1->GetNode()->GetNodeNeighbors();
            for(MolecularModeling::AtomVector::iterator it = c1c1_neighbors.begin(); it != c1c1_neighbors.end(); it++)
            {
                MolecularModeling::Atom* c1c1_neighbor = (*it);
                if(c1c1_neighbor->GetId().compare(C1->GetId()) != 0)
                {
                    temp << c1c1_neighbor->GetName().at(0);
                    atoms_1.push_back(c1c1_neighbor);
                }
            }
            temp << "C";
            atoms_1.push_back(C1C2);
            MolecularModeling::AtomVector c1c2_neighbors = C1C2->GetNode()->GetNodeNeighbors();
            for(MolecularModeling::AtomVector::iterator it = c1c2_neighbors.begin(); it != c1c2_neighbors.end(); it++)
            {
                MolecularModeling::Atom* c1c2_neighbor = (*it);
                if(c1c2_neighbor->GetId().compare(C1->GetId()) != 0)
                {
                    temp << c1c2_neighbor->GetName().at(0);
                    atoms_1.push_back(c1c2_neighbor);
                }
            }
            temp << "C";
            atoms_1.push_back(C1C3);
            MolecularModeling::AtomVector c1c3_neighbors = C1C3->GetNode()->GetNodeNeighbors();
            for(MolecularModeling::AtomVector::iterator it = c1c3_neighbors.begin(); it != c1c3_neighbors.end(); it++)
            {
                MolecularModeling::Atom* c1c3_neighbor = (*it);
                if(c1c3_neighbor->GetId().compare(C1->GetId()) != 0)
                {
                    temp << c1c3_neighbor->GetName().at(0);
                    atoms_1.push_back(c1c3_neighbor);
                }
            }
        }
        if(temp.str().compare("O-CCCC") == 0 || temp.str().compare("O-CCHHHCHHHCHHH") == 0 )
        {
            terminal_atoms = atoms_1;
            return "TBT";
        }
    }
    if(C2 != NULL)
    {
        MolecularModeling::Atom* C2C1 = NULL;
        MolecularModeling::Atom* C2C2 = NULL;
        MolecularModeling::Atom* C2C3 = NULL;
        pattern << "-" << "C";
        atoms_2.push_back(C2);
        MolecularModeling::AtomVector c2_neighbors = C2->GetNode()->GetNodeNeighbors();
        for(MolecularModeling::AtomVector::iterator it = c2_neighbors.begin(); it != c2_neighbors.end(); it++)
        {
            MolecularModeling::Atom* c2_neighbor = (*it);
            if(c2_neighbor->GetId().compare(target->GetId()) != 0 && c2_neighbor->GetName().at(0) == 'C' && C2C1 == NULL)
                C2C1 = c2_neighbor;
            else if(c2_neighbor->GetId().compare(target->GetId()) != 0 && c2_neighbor->GetName().at(0) == 'C' && C2C1 != NULL && C2C2 == NULL)
                C2C2 = c2_neighbor;
            else if(c2_neighbor->GetId().compare(target->GetId()) != 0 && c2_neighbor->GetName().at(0) == 'C' && C2C1 != NULL && C2C2 != NULL && C2C3 == NULL)
                C2C3 = c2_neighbor;
        }
        if(C2C1 != NULL && C2C2 != NULL && C2C3 != NULL)
        {
            pattern << "C";
            atoms_2.push_back(C2C1);
            MolecularModeling::AtomVector c2c1_neighbors = C2C1->GetNode()->GetNodeNeighbors();
            for(MolecularModeling::AtomVector::iterator it = c2c1_neighbors.begin(); it != c2c1_neighbors.end(); it++)
            {
                MolecularModeling::Atom* c2c1_neighbor = (*it);
                if(c2c1_neighbor->GetId().compare(C2->GetId()) != 0)
                {
                    pattern << c2c1_neighbor->GetName().at(0);
                    atoms_2.push_back(c2c1_neighbor);
                }
            }
            pattern << "C";
            atoms_2.push_back(C2C2);
            MolecularModeling::AtomVector c2c2_neighbors = C2C2->GetNode()->GetNodeNeighbors();
            for(MolecularModeling::AtomVector::iterator it = c2c2_neighbors.begin(); it != c2c2_neighbors.end(); it++)
            {
                MolecularModeling::Atom* c2c2_neighbor = (*it);
                if(c2c2_neighbor->GetId().compare(C2->GetId()) != 0)
                {
                    pattern << c2c2_neighbor->GetName().at(0);
                    atoms_2.push_back(c2c2_neighbor);
                }
            }
            pattern << "C";
            atoms_2.push_back(C2C3);
            MolecularModeling::AtomVector c2c3_neighbors = C2C3->GetNode()->GetNodeNeighbors();
            for(MolecularModeling::AtomVector::iterator it = c2c3_neighbors.begin(); it != c2c3_neighbors.end(); it++)
            {
                MolecularModeling::Atom* c2c3_neighbor = (*it);
                if(c2c3_neighbor->GetId().compare(C2->GetId()) != 0)
                {
                    pattern << c2c3_neighbor->GetName().at(0);
                    atoms_2.push_back(c2c3_neighbor);
                }
            }
        }
        if(pattern.str().compare("O-CCCC") == 0 || pattern.str().compare("O-CCHHHCHHHCHHH") == 0 )
        {
            terminal_atoms = atoms_2;
            return "TBT";
        }
    }
    return "";
}

std::string Assembly::CheckTerminals(MolecularModeling::Atom* target, MolecularModeling::AtomVector& terminal_atoms)
{
    if(target !=NULL)
    {
      MolecularModeling::AtomVector o_neighbors = target->GetNode()->GetNodeNeighbors();
    	//I have encounter the situation where a NLN is recognized as ROH, because the anomeric nitrogen has hydrogen and satisfies the criteria for ROH.
    	// So, I added codes to check if terminal is protein.
    	bool non_protein_terminal = true;
    	if (target->GetResidue()->CheckIfProtein())
      {
    	    non_protein_terminal = false;
    	}

      if (non_protein_terminal)
      {
            if(CheckROHTerminal(target, terminal_atoms).compare("") != 0)
                return "ROH";
            else if(CheckOMETerminal(target, terminal_atoms).compare("") != 0)
                return "OME";
            else if(CheckTBTTerminal(target, terminal_atoms).compare("") != 0)
                return "TBT";
            else
		            return "Unknown";
	    }
      //else if(o_neighbors.size() == 2)
      else if(o_neighbors.size() >= 2) //Not just size =2 ,if terminal is NLN && input pdb file contains hydrogen.the sidechain connecting nitrogen contain 3 atoms
      {
          Atom* target_o_neighbor = NULL;
          if(o_neighbors.at(0)->GetDescription().find("Het;") != std::string::npos && o_neighbors.at(1)->GetDescription().find("Het;") == std::string::npos)	//if one is het and the other is not
              target_o_neighbor = o_neighbors.at(1);
          else if(o_neighbors.at(0)->GetDescription().find("Het;") == std::string::npos && o_neighbors.at(1)->GetDescription().find("Het;") != std::string::npos)
              target_o_neighbor = o_neighbors.at(0);

        	    //My code for assigning target_o_neighbor:
        	    /*for (unsigned int i=0; i< o_neighbors.size(); i++){
        		if (o_neighbors[i] -> GetResidue() -> CheckIfProtein()){
        		    //assuming normal structure, all neighbor atoms should belong to the same protein.
        		    target_o_neighbor = o_neighbors[i];
        		}
        	    }*/
        	    // Yao Xiao: my code ends.
          if(target_o_neighbor != NULL)
          {
              ResidueVector residues = this->GetAllResiduesOfAssembly();
              Residue* target_residue = NULL;
              for(ResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
              {
                  Residue* residue = *it;
                  if(residue->GetId().compare(target_o_neighbor->GetResidue()->GetId()) == 0)
                  {
                      target_residue = residue;
                      break;
                  }
              }
              if(target_residue != NULL)
              {
                terminal_atoms = target_residue->GetAtoms();
              }
              gmml::AminoacidGlycamMap aminoacid_glycam = gmml::AminoacidGlycamLookup(target_o_neighbor->GetResidue()->GetName());
              gmml::AminoacidGlycamMap glycam_aminoacid = gmml::GlycamAminoacidLookup(target_o_neighbor->GetResidue()->GetName());

              if(aminoacid_glycam.aminoacid_name_.compare("") != 0)
              {
                  return aminoacid_glycam.aminoacid_name_;	//Why return the amino acid name instead of Glycam name?
                                                            //For glyfinder; will work out solution if both are needed (Dave)
                  // return aminoacid_glycam.glycam_name_;
              }
              else if(glycam_aminoacid.glycam_name_.compare("") != 0)
              {
                  return glycam_aminoacid.aminoacid_name_;	////Why return the amino acid name instead of Glycam name?
                                                            //For glyfinder; will work out solution if both are needed (Dave)
                  // return glycam_aminoacid.glycam_name_;
              }
              else
              {
	                std::cout << "This return." << "\n";
                  return target_o_neighbor->GetResidue()->GetName();
              }
          }
          else
          {
              return "Unknown";
          }
      }
      else
      {
          return "Unknown";
      }
    }
    else
    {
        return "Unknown";
    }
}

void Assembly::CheckLinkageNote(Glycan::Monosaccharide* mono1, Glycan::Monosaccharide* mono2, std::string linkage, std::vector<std::string>& checked_linkages)
{
    if(find(checked_linkages.begin(), checked_linkages.end(), linkage) == checked_linkages.end())///If this linkage hasn't been checked before by calling the function on other side of the linkage
    {
        std::vector<std::string> linkage_tokens = gmml::Split(linkage, "-");

        std::stringstream reverse_linkage;
        reverse_linkage << linkage_tokens.at(2) << "-" << linkage_tokens.at(1) << "-" << linkage_tokens.at(0);
        checked_linkages.push_back(linkage);
        checked_linkages.push_back(reverse_linkage.str());

	// Check if vector has a size and do some error handling
	std::vector<std::string> left_c_index_vector = gmml::Split(gmml::Split(linkage_tokens.at(0), "_").at(0), "C*,\'");
	std::vector<std::string> right_c_index_vector = gmml::Split(gmml::Split(linkage_tokens.at(2), "_").at(0), "C*,\'");
	std::vector<std::string> glycosidic_o_index_vector = gmml::Split(gmml::Split(linkage_tokens.at(1), "_").at(0), "ON*,\'");

	if( !left_c_index_vector.empty() && !right_c_index_vector.empty() && !glycosidic_o_index_vector.empty() ) {
		// Check if string in vector.at(0) is int
		std::vector <std::vector <std::string> > C_index_ErrCheck = std::vector <std::vector <std::string> >();
		C_index_ErrCheck.push_back(left_c_index_vector);
		C_index_ErrCheck.push_back(right_c_index_vector);
		C_index_ErrCheck.push_back(glycosidic_o_index_vector);
		bool AllIndexesAreInt= true;
		std::stringstream ss;
		int index;

		for (std::vector <std::vector <std::string> >::iterator it= C_index_ErrCheck.begin(); it!= C_index_ErrCheck.end(); it++){
			ss << (*it).at(0);
			ss >> index;
			if (ss.fail()){
				AllIndexesAreInt= false;
				ss.clear();
			}
			ss.str("");
		}

		if( AllIndexesAreInt ) {
			int left_c_index = gmml::ConvertString<int>(left_c_index_vector.at(0));
			int right_c_index = gmml::ConvertString<int>(right_c_index_vector.at(0));
			int glycosidic_o_index = gmml::ConvertString<int>(glycosidic_o_index_vector.at(0));
			if((left_c_index <= 9) && (left_c_index >= 0)) {
				if(left_c_index != glycosidic_o_index && right_c_index != glycosidic_o_index)
				{
            				Glycan::Note* linkage_note = new Glycan::Note();
            				linkage_note->type_ = Glycan::ERROR;
            				linkage_note->category_ = Glycan::GLYCOSIDIC;
            				std::stringstream n;
            				n << mono1->sugar_name_.monosaccharide_short_name_ << ": Glycosidic oxygen/nitrogen index does not conform to carbon index in the linkage to "
              		  	  	  << mono2->sugar_name_.monosaccharide_short_name_ << ". " << gmml::Split(linkage_tokens.at(0), "_").at(0) << "-" << gmml::Split(linkage_tokens.at(1), "_").at(0)
              		  	  	  << "-" << gmml::Split(linkage_tokens.at(2), "_").at(0);
            				linkage_note->description_ = n.str();
            				// this->AddNote(linkage_note);
                    mono1->mono_notes_.push_back(linkage_note);
                    mono2->mono_notes_.push_back(linkage_note);
        			}
			}
		}
	}
    }
}

void Assembly::BuildOligosaccharideTreeStructure(Glycan::Monosaccharide *key, std::vector<Glycan::Monosaccharide*> values, Glycan::Oligosaccharide *oligo,
                                                 std::vector<int>& visited_monos, std::map<Glycan::Monosaccharide*, std::vector<Glycan::Monosaccharide*> > monos_table,
                                                 std::map<Glycan::Monosaccharide*, std::vector<std::string> > monos_table_linkages, std::vector<std::string>& visited_linkages)
{
    oligo->root_ = key;
    if(values.size() == 0)
    {
        oligo->child_oligos_ = std::vector<Glycan::Oligosaccharide*>();
        oligo->child_oligos_linkages_ = std::vector<std::string>();
        visited_monos.push_back(key->mono_id_);
        return;
    }
    else
    {
        for(std::vector<Glycan::Monosaccharide*>::iterator it = values.begin(); it != values.end(); it++)
        {
            Glycan::Monosaccharide* value_mono = (*it);
            if(find(visited_monos.begin(), visited_monos.end(), value_mono->mono_id_) == visited_monos.end())
            {
                int it_index = distance(values.begin(), it);
                std::vector<std::string> key_mono_linkages = monos_table_linkages[key];
                /*for(int i = 0; i < key_mono_linkages.size(); i++ )  {
                    std::cout << "HELP US!" << key_mono_linkages.at( i ) << "\n";
                }*/
                std::string link = key_mono_linkages.at(it_index);
                std::stringstream reverse_link;
                reverse_link << gmml::Split(link, "-").at(2) << "-" << gmml::Split(link, "-").at(1) << "-" << gmml::Split(link, "-").at(0);
                if(find(visited_linkages.begin(), visited_linkages.end(), link) == visited_linkages.end() &&
                        find(visited_linkages.begin(), visited_linkages.end(), reverse_link.str()) == visited_linkages.end())
                {
                    //std::cout << "key id " << key->mono_id_  << ", value id " << value_mono->mono_id_ << "\n";
                    Glycan::Oligosaccharide* child_oligo = new Glycan::Oligosaccharide(this);
                    CalculateOligosaccharideBFactor(child_oligo, values);
                    std::vector<Glycan::Monosaccharide*> value_mono_values = monos_table[value_mono];
                    visited_linkages.push_back(link);
                    //std::cout << "call " << value_mono->mono_id_ << "\n";
                    BuildOligosaccharideTreeStructure(value_mono, value_mono_values, child_oligo, visited_monos, monos_table, monos_table_linkages, visited_linkages);
                    oligo->child_oligos_.push_back(child_oligo);
                    oligo->child_oligos_linkages_.push_back(link);
                }
            }
        }
        visited_monos.push_back(key->mono_id_);
        return;
    }
}

void Assembly::CalculateOligosaccharideBFactor(Glycan::Oligosaccharide* oligo, std::vector<Glycan::Monosaccharide*> monos)
{
  float total_b_factor = 0;
  int num_monos = 0;
  for (std::vector<Glycan::Monosaccharide*>::iterator it = monos.begin(); it != monos.end(); it++)
  {
    Glycan::Monosaccharide* mono = (*it);
    float this_b_factor = mono->b_factor_;
    total_b_factor = total_b_factor + this_b_factor;
    num_monos++;
  }
  oligo->oligosaccharide_b_factor_ = total_b_factor/num_monos;
}

//I know there is a better way to do this and assign Phi/Psi/Omega values to the oligosaccharide struct, but for now it just needs to be done quickly
//See Figure 1 in Vina Carb paper by our group if confused about which atoms are chosen and why
//Dave 9/27/18
/*! \todo Fix unused parent_oligo variable

Here is the error:

src/MolecularModeling/Assembly/SugarIdentification/oligosaccharidedetection.cc: At global scope:
src/MolecularModeling/Assembly/SugarIdentification/oligosaccharidedetection.cc:4014:61: warning: unused parameter 'parent_oligo' [-Wunused-parameter]
 double Assembly::CalculatePhiAngle(Glycan::Oligosaccharide* parent_oligo, Glycan::Oligosaccharide* child_oligo, std::string parent_atom_id, std::string child_atom_id, std::string glycosidic_atom_id)
                                                             ^
*/
double Assembly::CalculatePhiAngle(Glycan::Oligosaccharide* child_oligo, std::string parent_atom_id, std::string child_atom_id, std::string glycosidic_atom_id)
{
  int local_debug = -1;
  if(local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__, gmml::INF, "CalculatingPhiAngle");
  }
  //(O5-C1-O-Cx') {Ring oxygen of child_oligo}-{child_atom_id}-{glycosidic_atom_id}-{parent_atom_id}

  MolecularModeling::Atom* O5 = NULL;
  MolecularModeling::Atom* C1 = NULL;
  MolecularModeling::Atom* glycosidicO = NULL;
  MolecularModeling::Atom* Cx = NULL;

  //Get Anomeric Carbon and ring oxygen
  for(std::vector<MolecularModeling::Atom*>::iterator it = child_oligo->root_->cycle_atoms_.begin();
      it != child_oligo->root_->cycle_atoms_.end(); ++it)
  {
    MolecularModeling::Atom* this_cycle_atom = (*it);
    std::string this_element = this_cycle_atom->GetElementSymbol();
    std::string this_atom_id = this_cycle_atom->GetId();

    if (this_atom_id.find("O") != std::string::npos)
    {
      O5 = this_cycle_atom;
    }
    if (this_atom_id == child_atom_id)
    {
      C1 = this_cycle_atom;
    }
  }
  if((O5 == NULL)||(C1 == NULL))
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "O5 or C1 is null in calculate Phi Angle");
    return -9999;
  }
  if(local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__, gmml::INF, "O5");
    gmml::log(__LINE__, __FILE__, gmml::INF, O5->GetId());
    gmml::log(__LINE__, __FILE__, gmml::INF, "C1");
    gmml::log(__LINE__, __FILE__, gmml::INF, C1->GetId());
  }
  //Find glycosidicO
  MolecularModeling::AtomVector C1_neighbors = C1->GetNode()->GetNodeNeighbors();
  for (MolecularModeling::AtomVector::iterator it = C1_neighbors.begin(); it != C1_neighbors.end(); ++it)
  {
    MolecularModeling::Atom* this_neighbor = (*it);
    std::string this_atom_id = this_neighbor->GetId();
    if (this_atom_id == glycosidic_atom_id)
    {
      glycosidicO = this_neighbor;
    }
  }
  if(glycosidicO == NULL)
  {
    gmml::log(__LINE__, __FILE__, gmml::INF, glycosidicO->GetId());
    gmml::log(__LINE__, __FILE__, gmml::ERR, "glycosidicO is null in calculate Phi Angle");
    return -9999;
  }

  //Get Cx
  MolecularModeling::AtomVector glycosidicO_neighbors = glycosidicO->GetNode()->GetNodeNeighbors();
  for (MolecularModeling::AtomVector::iterator it = glycosidicO_neighbors.begin(); it != glycosidicO_neighbors.end(); ++it)
  {
    MolecularModeling::Atom* this_neighbor = (*it);
    std::string this_atom_id = this_neighbor->GetId();
    if (this_atom_id == parent_atom_id)
    {
      Cx = this_neighbor;
    }
  }
  if(Cx == NULL)
  {
    gmml::log(__LINE__, __FILE__, gmml::INF, Cx->GetId());
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Cx is null in calculate Phi Angle");
    return -9999;
  }

  double Phi = CalculateTorsionAngleByAtoms(O5, C1, glycosidicO, Cx);
  if (local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__, gmml::INF, std::to_string(Phi));
  }
  return Phi;
}

double Assembly::CalculatePsiAngle(Glycan::Oligosaccharide* child_oligo, std::string parent_atom_id, std::string child_atom_id, std::string glycosidic_atom_id)
{
  //(C1-O-Cx'-C[x-1]') {child_atom_id}-{glycosidic_atom_id}-{parent_atom_id}-{parent_atom_id - 1}

  MolecularModeling::Atom* C1 = NULL;
  MolecularModeling::Atom* glycosidicO = NULL;
  MolecularModeling::Atom* Cx = NULL;
  MolecularModeling::Atom* Cx_1 = NULL;

  //Get Anomeric Carbon
  for(std::vector<MolecularModeling::Atom*>::iterator it = child_oligo->root_->cycle_atoms_.begin();
      it != child_oligo->root_->cycle_atoms_.end(); ++it)
  {
    MolecularModeling::Atom* this_cycle_atom = (*it);
    std::string this_atom_id = this_cycle_atom->GetId();
    if (this_atom_id == child_atom_id)
    {
      C1 = this_cycle_atom;
    }
  }
  if(C1 == NULL)
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "C1 is null in calculate Psi Angle");
    return -9999;
  }

  //Find glycosidicO
  MolecularModeling::AtomVector C1_neighbors = C1->GetNode()->GetNodeNeighbors();
  for (MolecularModeling::AtomVector::iterator it = C1_neighbors.begin(); it != C1_neighbors.end(); ++it)
  {
    MolecularModeling::Atom* this_neighbor = (*it);
    std::string this_atom_id = this_neighbor->GetId();
    if (this_atom_id == glycosidic_atom_id)
    {
      glycosidicO = this_neighbor;
    }
  }
  if(glycosidicO == NULL)
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "glycosidicO is null in calculate Psi Angle");
    return -9999;
  }

  //Get Cx
  MolecularModeling::AtomVector glycosidicO_neighbors = glycosidicO->GetNode()->GetNodeNeighbors();
  for (MolecularModeling::AtomVector::iterator it = glycosidicO_neighbors.begin(); it != glycosidicO_neighbors.end(); ++it)
  {
    MolecularModeling::Atom* this_neighbor = (*it);
    std::string this_atom_id = this_neighbor->GetId();
    if (this_atom_id == parent_atom_id)
    {
      Cx = this_neighbor;
    }
  }
  if(Cx == NULL)
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Cx is null in calculate Psi Angle");
    return -9999;
  }

  //Get Cx-1
  MolecularModeling::AtomVector Cx_neighbors = Cx->GetNode()->GetNodeNeighbors();
  for (MolecularModeling::AtomVector::iterator it = Cx_neighbors.begin(); it != Cx_neighbors.end(); ++it)
  {
    MolecularModeling::Atom* this_neighbor = (*it);
    std::string this_atom_id = this_neighbor->GetId();
    int Cx_atom_number = parent_atom_id.at(1) - '0';
    int this_atom_number = this_atom_id.at(1) - '0';
    if ((Cx_atom_number != 1) && ((Cx_atom_number - 1) == this_atom_number))
    {
      Cx_1 = this_neighbor;
    }
    else if((Cx_atom_number == 1) && ((Cx_atom_number + 1) == this_atom_number))
    {
      Cx_1 = this_neighbor;
    }
  }
  if(Cx_1 == NULL)
  {
    gmml::log(__LINE__, __FILE__, gmml::INF, Cx->GetId());
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Cx_1 is null in calculate Psi Angle");
    return -9999;
  }

  return CalculateTorsionAngleByAtoms(C1, glycosidicO, Cx, Cx_1);
}

double Assembly::CalculateOmegaAngle(Glycan::Oligosaccharide* parent_oligo, std::string parent_atom_id, std::string glycosidic_atom_id)
{
  int local_debug = -1;
  //(O-C6'-C5'-O5') {glycosidic_atom_id}-{parent_atom_id}-{Carbon 5 in parent oligo}-{Ring oxygen in parent_oligo}

  MolecularModeling::Atom* glycosidicO = NULL;
  MolecularModeling::Atom* C6prime = NULL;
  MolecularModeling::Atom* C5prime = NULL;
  MolecularModeling::Atom* O5prime = NULL;

  //Get C5' and O5' first
  if (local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Get C5' and O5' first");
  }
  for(std::vector<MolecularModeling::Atom*>::iterator it = parent_oligo->root_->cycle_atoms_.begin();
      it != parent_oligo->root_->cycle_atoms_.end(); ++it)
  {
    MolecularModeling::Atom* this_cycle_atom = (*it);
    std::string this_atom_id = this_cycle_atom->GetId();
    int this_atom_number = this_atom_id.at(1) - '0';
    std::string this_element = this_cycle_atom->GetElementSymbol();
    if (this_atom_number == 5)
    {
      if (this_atom_id.find("C") != std::string::npos)
      {
        C5prime = this_cycle_atom;
      }
    }
    if (this_atom_id.find("O") != std::string::npos)
    {
      O5prime = this_cycle_atom;
    }
  }
  if ((C5prime == NULL)||(O5prime == NULL))
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "O5prime or C5prime is null in calculate Omega Angle");
    return -9999;
  }

  //Get C6'
  if (local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Get C6'");
  }
  MolecularModeling::AtomVector C5prime_neighbors = C5prime->GetNode()->GetNodeNeighbors();
  for (MolecularModeling::AtomVector::iterator it = C5prime_neighbors.begin(); it != C5prime_neighbors.end(); ++it)
  {
    MolecularModeling::Atom* this_neighbor = (*it);
    std::string this_atom_id = this_neighbor->GetId();
    if ( this_atom_id == parent_atom_id)
    {
      C6prime = this_neighbor;
    }
  }
  if(C6prime == NULL)
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "C6prime is null in calculate Omega Angle");
    return -9999;
  }
  //Get glycosidic oxygen
  if (local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Get glycosidic oxygen");
  }
  MolecularModeling::AtomVector C6prime_neighbors = C6prime->GetNode()->GetNodeNeighbors();
  for (MolecularModeling::AtomVector::iterator it = C6prime_neighbors.begin(); it != C6prime_neighbors.end(); ++it)
  {
    MolecularModeling::Atom* this_neighbor = (*it);
    std::string this_atom_id = this_neighbor->GetId();
    if ( this_atom_id == glycosidic_atom_id)
    {
      glycosidicO = this_neighbor;
    }
  }
  if(glycosidicO == NULL)
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "glycosidicO is null in calculate Omega Angle");
    return -9999;
  }
  if (local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__,  gmml::INF, "About to calcuate torsion");
  }
  return CalculateTorsionAngleByAtoms(glycosidicO, C6prime, C5prime, O5prime);
}

void Assembly::GetAuthorNaming(std::vector< std::string > amino_lib_files, Glycan::Monosaccharide* mono, std::string CCD_Path )
{
  int local_debug = -1;
  std::string residueName = mono->cycle_atoms_.at(0)->GetResidue()->GetName();
  char CCDhash = residueName[0];
  std::stringstream CCDfilepath;
  CCDfilepath << CCD_Path << "/" << CCDhash << "/" << residueName << "/" << residueName << ".pdb";

  // Initialize an Assembly from the PDB file.

  //The next two lines are the fastest way to check if a file exists according to https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
  struct stat buffer;
  if( stat( CCDfilepath.str().c_str(), &buffer ) == 0 )
  {
    MolecularModeling::Assembly CCDassembly(CCDfilepath.str(), gmml::PDB);

    //TODO Remove H atoms from CCD PDB files; They are screwing up too much in the sugar naming (AH vs A in the chemical_code_)
    std::vector<MolecularModeling::Residue*> residues = CCDassembly.GetAllResiduesOfAssembly();
    for(std::vector<MolecularModeling::Residue*>::iterator it = residues.begin(); it != residues.end(); it++)
    {
      MolecularModeling::Residue* thisResidue = (*it);
      std::vector<MolecularModeling::Atom*> atoms = thisResidue->GetAtoms();
      for (std::vector<MolecularModeling::Atom*>::iterator it = atoms.begin(); it != atoms.end(); it++)
      {
        MolecularModeling::Atom* thisAtom = *it;
        if(thisAtom->GetElementSymbol() == "H")
        {
          thisResidue->RemoveAtom(thisAtom);
        }
      }
    }

    // Build by Distance
    CCDassembly.BuildStructureByDistance(10);
    // Find the Sugars.
    std::vector<Glycan::Oligosaccharide*> authorOligos = CCDassembly.ExtractSugars(amino_lib_files, false, false, false, CCD_Path);
    //The vector is really just a monosaccharide
    if(authorOligos.size() > 0)
    {
      Glycan::Oligosaccharide* authorOligo = authorOligos[0];
      if(authorOligo->mono_nodes_.size() > 0)
      {
        if(authorOligos[0]->mono_nodes_[0]->sugar_name_.pdb_code_ == "")
        {
          if(mono->sugar_name_.monosaccharide_name_ == authorOligos[0]->mono_nodes_[0]->sugar_name_.monosaccharide_name_)//should have the formula for R groups so will be an "exact" match
          {
            mono->sugar_name_.pdb_code_ = residueName;
          }
          authorOligos[0]->mono_nodes_[0]->sugar_name_.pdb_code_ = residueName;
        }
        mono->author_sugar_name_ = authorOligos[0]->mono_nodes_[0]->sugar_name_;
      }
      if(local_debug > 0)
      {
        gmml::log(__LINE__, __FILE__,  gmml::INF, "Author Name");
        gmml::log(__LINE__, __FILE__,  gmml::INF, mono->author_sugar_name_.monosaccharide_name_);
      }
    }
  }
}
