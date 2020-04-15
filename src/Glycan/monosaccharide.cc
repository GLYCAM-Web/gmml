#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iterator>

#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/Glycan/chemicalcode.hpp"
#include "../../includes/Glycan/sugarname.hpp"
#include "../../includes/Glycan/note.hpp"
#include "../../includes/Glycan/monosaccharide.hpp"
#include "../../includes/MolecularModeling/assembly.hpp"
#include "../../includes/MolecularModeling/ring_shape_detection.hpp"

using Glycan::Monosaccharide;
using MolecularModeling::Assembly;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Monosaccharide::Monosaccharide(const Monosaccharide &mono)
{

}

Monosaccharide::Monosaccharide(std::string* cycle_atoms_str, std::vector<MolecularModeling::Atom*>& cycle_atoms, MolecularModeling::Assembly* this_assembly, std::string CCD_Path)
{
  int local_debug = -1;
  residue_name_ = cycle_atoms[0]->GetResidue()->GetName();
  cycle_atoms[0]->GetResidue()->SetIsSugar(true);
  // std::cout << residue_name_ << "\n";
  is_visited_ = false;
  is_root_ = false;
  is_counted_ = false;
  anomeric_carbon_pointer_ = NULL;
  sugar_name_ = {"","","","","","","","","",""};
  //Set cycle atoms
  for (std::vector<MolecularModeling::Atom*>::iterator it = cycle_atoms.begin(); it != cycle_atoms.end(); it++ )
  {
    (*it)-> SetIsCycle(true);
  }

  //Detect Anomeric Carbon, assign Cycle atoms
  Glycan::Note* anomeric_note = new Glycan::Note();
  std::vector< std::string > anomeric_carbons_status = std::vector< std::string >();
  anomeric_carbon_pointer_ = FindAnomericCarbon( anomeric_note, anomeric_carbons_status, cycle_atoms, *cycle_atoms_str );
  Glycan::Note anomeric_carbon_debug = *anomeric_note;
  if(anomeric_carbon_pointer_!=NULL)
  {
    anomeric_carbon_pointer_->SetIsAnomericCarbon(true);
    std::stringstream sorted_cycle_stream;
    cycle_atoms_ = this_assembly->SortCycle( cycle_atoms, anomeric_carbon_pointer_, sorted_cycle_stream );
    cycle_atoms_str_ = sorted_cycle_stream.str();
  }
  else
  {
    cycle_atoms_ = cycle_atoms;
  }
  if ((local_debug > 0) && (cycle_atoms_.empty()))
  {
    gmml::log(__LINE__, __FILE__, gmml::INF, "This monosaccharide has no cycle atoms!!!");
    gmml::log(__LINE__, __FILE__, gmml::INF, cycle_atoms_str_);
    gmml::log(__LINE__, __FILE__, gmml::INF, std::to_string(cycle_atoms_.size()));
    gmml::log(__LINE__, __FILE__, gmml::INF, *cycle_atoms_str);
    gmml::log(__LINE__, __FILE__, gmml::INF, std::to_string(cycle_atoms.size()));

  }
  //Get BFMP
  if( cycle_atoms_.size() > 5 )
  {
    bfmp_ring_conformation_ = glylib::CalculateRingShapeBFMP(this);
  }

  //Assign side_atoms_
  std::vector< std::string > orientations = GetSideGroupOrientations(this_assembly);//This also sets the first side atom for each ring position

  //this function sometimes deletes the anomeric carbon's neighbor...
  // this->InitiateDetectionOfCompleteSideGroupAtoms ();
  std::vector<MolecularModeling::Atom*> plus_sides = ExtractAdditionalSideAtoms();


  ///DERIVATIVE/MODIFICATION PATTERN EXTRACTION
  this->ExtractDerivatives(this_assembly);

  //Assign Chemical Code, and SugarNames
  chemical_code_ = BuildChemicalCode(orientations);
  sugar_name_ = gmml::SugarStereoChemistryNameLookup(chemical_code_->toString());

  ///CALCULATING B FACTOR
  float total_b_factor = 0;
  int num_atoms =0;
  for (std::vector< std::vector<MolecularModeling::Atom*> >::iterator it1 = side_atoms_.begin(); it1 != side_atoms_.end(); it1++ )
  {
    std::vector<MolecularModeling::Atom*> sides = ( *it1 );
    for(std::vector<MolecularModeling::Atom*>::iterator it2 = sides.begin(); it2 != sides.end(); it2++ )
    {
      MolecularModeling::Atom* atom = (*it2);
      if(atom !=NULL)
      {
        float this_b_factor = atom->GetBFactor();
        total_b_factor = total_b_factor + this_b_factor;
        num_atoms++;
        // std::stringstream test;
        // test << "Atom_b_factor_mono" << this_b_factor;
        // gmml::log(__LINE__, __FILE__, gmml::INF, test.str());
      }
    }
  }
  for (std::vector<MolecularModeling::Atom*>::iterator it = cycle_atoms_.begin(); it != cycle_atoms_.end(); it++ )  {
    MolecularModeling::Atom* atom = (*it);
    float this_b_factor = atom->GetBFactor();
    total_b_factor = total_b_factor + this_b_factor;
    num_atoms++;
  }
  b_factor_ = total_b_factor/num_atoms;

  GenerateCompleteName(plus_sides, this, this_assembly);

  if(sugar_name_.chemical_code_string_ == "")
  {//No match in lookup table, check if its from Deoxy groups

    Glycan::ChemicalCode new_code = *chemical_code_;
    // std::cout << "Old code" << chemical_code_->toString() << "\n";
    std::vector<std::string> v(1, "");
    new_code.left_middle_ = v;
    new_code.right_middle_ = v;

    std::vector<std::string> deoxy_locations;
    for(unsigned int i = 0; i < new_code.right_up_.size(); i++)
    {
      if(((new_code.right_up_[i].substr(0,1) == "+") ||
         (new_code.right_up_[i].substr(0,1) == "-")) &&
         (new_code.right_up_[i].find("d") != std::string::npos))
      {
        // std::cout << new_code.right_up_[i] << "\n";
        new_code.right_up_[i] = new_code.right_up_[i].substr(0,2);
        // std::cout << new_code.right_up_[i] << "\n";
        derivatives_map_.push_back(std::make_pair(new_code.right_up_[i].substr(0,2), "xCHH"));
      }
    }
    for(unsigned int i = 0; i < new_code.right_down_.size(); i++)
    {
      if(((new_code.right_down_[i].substr(0,1) == "+") ||
         (new_code.right_down_[i].substr(0,1) == "-") )&&
         (new_code.right_down_[i].find("d") != std::string::npos))
      {
        // std::cout << new_code.right_down_[i] << "\n";
        new_code.right_down_[i] = new_code.right_down_[i].substr(0,2);
        // std::cout << new_code.right_down_[i] << "\n";
        derivatives_map_.push_back(std::make_pair(new_code.right_down_[i].substr(0,2), "xCHH"));
      }
    }
    if(!chemical_code_->left_middle_.empty())
    {
      for(unsigned int i = 0; i < chemical_code_->left_middle_.size(); i++)
      {
        // std::cout << chemical_code_->left_middle_[i] << "\n";
        deoxy_locations.push_back(chemical_code_->left_middle_[i].substr(0,1));
      }
    }
    if(!chemical_code_->right_middle_.empty())
    {
      for(unsigned int i = 0; i < chemical_code_->right_middle_.size(); i++)
      {
        // std::cout << chemical_code_->right_middle_[i] << "\n";
        deoxy_locations.push_back(chemical_code_->right_middle_[i].substr(0,1));
      }
    }
    std::sort(deoxy_locations.begin(), deoxy_locations.end());
    for(unsigned int i = 0; i < deoxy_locations.size(); i++)
    {
      // std::cout << "Deoxy at: " << deoxy_locations[i] << "\n";
      new_code.left_up_.push_back(deoxy_locations[i]);
      derivatives_map_.push_back(std::make_pair(deoxy_locations[i], "xCHH"));
    }
    std::sort(derivatives_map_.begin(), derivatives_map_.end());
    std::sort(new_code.left_up_.begin(), new_code.left_up_.end());
    // for(std::vector<std::string>::iterator it = chemical_code_->left_down_.begin(); it != chemical_code_->left_down_.end(); it++)
    // {
    //     std::cout << "left down: " << (*it) << "\n";
    // }
    // for(std::vector<std::string>::iterator it = chemical_code_->left_up_.begin(); it != chemical_code_->left_up_.end(); it++)
    // {
    //     std::cout << "left up: " << (*it) << "\n";
    // }
    // for(std::vector<std::string>::iterator it = chemical_code_->right_down_.begin(); it != chemical_code_->right_down_.end(); it++)
    // {
    //     std::cout << "right down: " << (*it) << "\n";
    // }
    // for(std::vector<std::string>::iterator it = chemical_code_->right_up_.begin(); it != chemical_code_->right_up_.end(); it++)
    // {
    //     std::cout << "right up: " << (*it) << "\n";
    // }
    // std::cout << "new code: " << new_code.toString() << "\n";


    Glycan::SugarName base_name = gmml::SugarStereoChemistryNameLookup( new_code.toString() );
    // std::cout << "new name: " << base_name.monosaccharide_stereochemistry_short_name_ << "\n";
    sugar_name_ = base_name;
    this->UpdateComplexSugarChemicalCode();
    Glycan::SugarName updated_name = gmml::ComplexSugarNameLookup(chemical_code_->toString());
    if((updated_name.monosaccharide_name_ != "") && (updated_name.monosaccharide_name_ != sugar_name_.monosaccharide_name_))
    {
      sugar_name_ = updated_name;
    }
    this->UpdatePdbCode();
    this->GenerateCompleteSugarName(this_assembly);

    // addDeoxyToName(base_name, chemical_code_, deoxy_locations);

    //this will create the correct name for the sugar, but the chemical code will not match, and will be incorrect, as it will have ^n at the deoxy locations.
    // sugar_name_ = base_name;

  }
  if((sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0) && (sugar_name_.monosaccharide_name_.compare("") == 0))
  {
    sugar_name_.monosaccharide_name_ = sugar_name_.monosaccharide_stereochemistry_name_;
  }

  //Check if Residue name matches; if not use the CCD to create author_sugar_name_
  std::string original_residue = cycle_atoms_.at(0)->GetResidue()->GetName();
  std::string original_residue_id = cycle_atoms_.at(0)->GetResidue()->GetId();
  CheckMonoNaming(original_residue, original_residue_id);

  // std::vector< std::string > idealOrientations = GetSideGroupOrientations(CCD_cycle_atoms);
  // author_chemical_code_ = BuildChemicalCode(idealOrientations);
  // author_sugar_name_ = gmml::SugarStereoChemistryNameLookup(chemical_code_->toString());

  //Add notes
  std::stringstream n;
  if( anomeric_note->description_.compare( "" ) != 0 )
  {
    n << sugar_name_.monosaccharide_short_name_ << ": " << anomeric_note->description_;
    anomeric_note->description_ = n.str();
    mono_notes_.push_back(anomeric_note );
  }

  ///ADDING NOTES/ISSUES OF RESIDUE NAMING

  std::vector<std::string> pdb_codes = gmml::Split(sugar_name_.pdb_code_, ",");
  if( pdb_codes.size() > 0 )
  {
    std::string pdb_code = "";
    bool found_code = false;
    for( std::vector< std::string >::iterator codes_it = pdb_codes.begin(); codes_it != pdb_codes.end(); codes_it++ )
    {
      pdb_code = ( *codes_it );
      if( pdb_code.compare( original_residue ) == 0 )
      {
        found_code = true;
        break;
      }
    }
    if( !found_code )
    {
      Glycan::Note* residue_naming_note = new Glycan::Note();
      residue_naming_note->category_ = Glycan::RESIDUE_NAME;
      std::stringstream res_ss;
      if( sugar_name_.pdb_code_.compare( "" ) == 0 )
      {
        residue_naming_note->type_ = Glycan::WARNING;
        res_ss << "PDB 3 letter code not found for " << sugar_name_.monosaccharide_short_name_;
      }
      else
      {
        residue_naming_note->type_ = Glycan::ERROR;
        res_ss << "Residue name, " << original_residue << " (" << original_residue_id << "), in input PDB file for " << sugar_name_.monosaccharide_short_name_ << " does not match GlyFinder residue code: " << sugar_name_.pdb_code_;
      }
      residue_naming_note->description_ = res_ss.str();
      mono_notes_.push_back(residue_naming_note);
      // mono_notes_.push_back( residue_naming_note );
      // GetAuthorNaming(amino_lib_files, mono, CCD_Path);
    }
  }
  // createAuthorSNFGname();
  createSNFGname();
  if(sugar_name_.pdb_code_ == "")
  {
    std::stringstream ss;
    ss << "Chemical code: " << chemical_code_->toString() << " not found for " << residue_name_;
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
  }
}
//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////



void Glycan::Monosaccharide::createSNFGname()
{//Added for the PDB to create their own SNFGs
  std::string monoSNFGName =  sugar_name_.monosaccharide_short_name_;
  if(monoSNFGName.size() > 0)
  {
    if(monoSNFGName[0] == sugar_name_.isomer_[0])
    {
      monoSNFGName.erase(0,1);
    }
    // if(((monoSNFGName[monoSNFGName.length()-1] == sugar_name_.configuration_[0]) && (monoSNFGName[monoSNFGName.length()-2] != 'N'))
    //   ||
    //   || (monoSNFGName[monoSNFGName.length()-1] == 'x'))
    // {
    //   monoSNFGName.erase(monoSNFGName.length()-1,1);
    // }
    if ((monoSNFGName[monoSNFGName.length()-1] == 'x') ||
       (monoSNFGName[monoSNFGName.length()-1] == 'a') ||
       (monoSNFGName[monoSNFGName.length()-1] == 'b'))
    {
       monoSNFGName.erase(monoSNFGName.length()-1,1);
    }
    // std::cout << "Ring type: " << sugar_name_.ring_type_[0] << "\n";
    //There is a bug in the ring type code; furanose sugar has ring type of P (CCD code GZL)
    // if(monoSNFGName[3] == tolower(sugar_name_.ring_type_[0]))
    // {
    //   monoSNFGName.erase(3,1);
    // }
    if((monoSNFGName[3] == 'p') || (monoSNFGName[3] == 'f'))
    {
      monoSNFGName.erase(3,1);
    }
    if((monoSNFGName[monoSNFGName.length()-1] == 'H') && (monoSNFGName.length() > 3))
    {
      monoSNFGName.erase(monoSNFGName.length()-1,1);
    }
  }
  SNFG_name_ = monoSNFGName;
}

void Glycan::Monosaccharide::createAuthorSNFGname()
{//Added for the PDB to create their own SNFGs
  std::string monoSNFGName =  author_sugar_name_.monosaccharide_short_name_;
  if(monoSNFGName.size() > 0)
  {
    if(monoSNFGName[0] == author_sugar_name_.isomer_[0])
    {
      monoSNFGName.erase(0,1);
    }
    if(((monoSNFGName[monoSNFGName.length()-1] == author_sugar_name_.configuration_[0]) && (monoSNFGName[monoSNFGName.length()-2] != 'N')) || (monoSNFGName[monoSNFGName.length()-1] == 'x'))
    {
      monoSNFGName.erase(monoSNFGName.length()-1,1);
    }
    if(monoSNFGName[3] == tolower(author_sugar_name_.ring_type_[0]))
    {
      monoSNFGName.erase(3,1);
    }
  }
  author_SNFG_name_ = monoSNFGName;
}

MolecularModeling::Atom* Glycan::Monosaccharide::FindAnomericCarbon( Glycan::Note* anomeric_note, std::vector< std::string > & anomeric_carbons_status, std::vector<MolecularModeling::Atom*> cycle, std::string cycle_atoms_str )
{
  MolecularModeling::Atom* anomeric_carbon = new MolecularModeling::Atom();
  for( std::vector<MolecularModeling::Atom*>::iterator it = cycle.begin(); it != cycle.end(); it++ ) {
    MolecularModeling::Atom* cycle_atom = ( *it );

    if( ( cycle_atom->GetName().substr( 0, 1 ).compare( "O" ) == 0 ) ) {///find oxygen in ring
      //                && isdigit(gmml::ConvertString<char>(cycle_atom->GetName().substr(1,1)))))
      MolecularModeling::AtomNode* node = cycle_atom->GetNode();
      std::vector<MolecularModeling::Atom*> neighbors = node->GetNodeNeighbors();

      ///Check the first neighbor of oxygen
      MolecularModeling::Atom* o_neighbor1 = neighbors.at( 0 );
      MolecularModeling::AtomNode* o_neighbor1_node = o_neighbor1->GetNode();
      std::vector<MolecularModeling::Atom*> o_neighbor1_neighbors = o_neighbor1_node->GetNodeNeighbors();

      for( std::vector<MolecularModeling::Atom*>::iterator it1 = o_neighbor1_neighbors.begin(); it1 != o_neighbor1_neighbors.end(); it1++ ) {///check if neighbor1 of oxygen has another oxygen or nitrogen neighbor
        MolecularModeling::Atom* o_neighbor1_neighbor = ( *it1 );

        if( cycle_atoms_str.find( o_neighbor1_neighbor->GetId() ) == std::string::npos ///if the neighbor is not one of the cycle atoms
                && ( o_neighbor1_neighbor->GetName().substr( 0, 1 ).compare( "O" ) == 0 || o_neighbor1_neighbor->GetName().substr( 0, 1 ).compare( "N" ) == 0 )
                && (o_neighbor1->GetElementSymbol() != "H")) { ///if first element is "O" or "N"
          //                        && isdigit(gmml::ConvertString<char>(neighbor1_neighbor->GetName().substr(1,1))))///if second element is a digit
          anomeric_carbon = o_neighbor1;
          anomeric_carbons_status.push_back( "Anomeric carbon: "  );
          anomeric_note->description_ = "";

          return anomeric_carbon;
        }
      }

      ///Check the second neighbor of oxygen
      MolecularModeling::Atom* o_neighbor2 = neighbors.at( 1 );
      MolecularModeling::AtomNode* o_neighbor2_node = o_neighbor2->GetNode();
      std::vector<MolecularModeling::Atom*> o_neighbor2_neighbors = o_neighbor2_node->GetNodeNeighbors();

      for( std::vector<MolecularModeling::Atom*>::iterator it2 = o_neighbor2_neighbors.begin(); it2 != o_neighbor2_neighbors.end(); it2++ ) {///check if neighbor2 of oxygen has another oxygen or nitrogen neighbor
        MolecularModeling::Atom* o_neighbor2_neighbor = ( *it2 );

        if( cycle_atoms_str.find( o_neighbor2_neighbor->GetId() ) == std::string::npos
                && ( o_neighbor2_neighbor->GetName().substr( 0, 1 ).compare( "O" ) == 0 || o_neighbor2_neighbor->GetName().substr( 0, 1 ).compare( "N" ) == 0 )
                && (o_neighbor2->GetElementSymbol() != "H")) {
          //                        && isdigit(gmml::ConvertString<char>(neighbor2_neighbor->GetName().substr(1,1))))
          anomeric_carbon = o_neighbor2;
          anomeric_carbons_status.push_back( "Anomeric carbon: "  );
          anomeric_note->description_ = "";

          return anomeric_carbon;
        }
      }

      ///Seems redundant to have this multiple times, if it is going to say the same thing and be generated for all cases after the first logic check.
      ///Specially if we ever want to change this Note to say something different. Like I am doing now. :)
      std::stringstream ss;
      ss << "Could not find glycosidic oxygen or nitrogen within " << gmml::maxCutOff << " Angstroms of a ring carbon bonded to the ring oxygen(" << node->GetAtom()->GetName() << ").";
      anomeric_note->description_ = ss.str();

      ///Check the order of the carbons based on their names to locate the anomeric
      std::stringstream ss1;
      for( unsigned int i = 0; i < o_neighbor1->GetName().size(); i++ ) {
        if( isdigit( o_neighbor1->GetName().at( i )) != 0 )
          ss1 << o_neighbor1->GetName().at( i );
      }
      std::stringstream ss2;
      for( unsigned int i = 0; i < o_neighbor2->GetName().size(); i++ ) {
        if( isdigit( o_neighbor2->GetName().at( i )) != 0 )
          ss2 << o_neighbor2->GetName().at(i);
      }
      if( gmml::ConvertString< int >( ss1.str() ) < gmml::ConvertString< int >( ss2.str() ) ) {
        anomeric_note->type_ = Glycan::WARNING;
        anomeric_note->category_ = Glycan::ANOMERIC;
        anomeric_carbons_status.push_back( "Anomeric carbon assigned based on its atom name index (" + ss1.str() + "), Anomeric Carbon is: " + anomeric_carbon->GetName() );
        if (o_neighbor1->GetElementSymbol() != "H")
          return o_neighbor1;
      }
      if( gmml::ConvertString< int >( ss2.str() ) < gmml::ConvertString< int >( ss1.str() ) ) {
        anomeric_note->type_ = Glycan::WARNING;
        anomeric_note->category_ = Glycan::ANOMERIC;
        anomeric_carbons_status.push_back( "Anomeric carbon assigned based on its atom name index (" + ss1.str() + "), Anomeric Carbon is: " + anomeric_carbon->GetName() );
        if (o_neighbor2->GetElementSymbol() != "H")
          return o_neighbor2;
      }

      ///Check non-ring neighbors of oxygen neighbors (the one without non-ring carbon is anomeric)
      bool neighbor2_is_anomeric = false;
      for( std::vector<MolecularModeling::Atom*>::iterator it1 = o_neighbor1_neighbors.begin(); it1 != o_neighbor1_neighbors.end(); it1++ ) {///check if neighbor1 of oxygen has non-ring carbon neighbor
        MolecularModeling::Atom* neighbor1_neighbor = ( *it1 );
        if( cycle_atoms_str.find( neighbor1_neighbor->GetId()) == std::string::npos ///if the neighbor is not one of the cycle atoms
                && ( neighbor1_neighbor->GetName().substr( 0, 1 ).compare( "C" ) == 0 ) ) {///if first element is "C"
          neighbor2_is_anomeric = true;
          break;
        }
      }
      bool neighbor1_is_anomeric = false;
      for( std::vector<MolecularModeling::Atom*>::iterator it1 = o_neighbor2_neighbors.begin(); it1 != o_neighbor2_neighbors.end(); it1++ ) {///check if neighbor1 of oxygen has non-ring carbon neighbor
        MolecularModeling::Atom* neighbor2_neighbor = ( *it1 );
        if( cycle_atoms_str.find( neighbor2_neighbor->GetId()) == std::string::npos ///if the neighbor is not one of the cycle atoms
                && ( neighbor2_neighbor->GetName().substr( 0, 1 ).compare( "C" ) == 0 ) ) {///if first element is "C"
          neighbor1_is_anomeric = true;
        }
      }
      if( !neighbor1_is_anomeric ) {
        anomeric_note->type_ = Glycan::WARNING;
        anomeric_note->category_ = Glycan::ANOMERIC;
        anomeric_carbons_status.push_back( "Anomeric carbon assigned to the ring carbon neighboring the ring oxygen but is not attached to an exocyclic carbon (that is, assuming the monosaccharide is an aldose), Anomeric Carbon is: " + anomeric_carbon->GetName() );
        if (o_neighbor2->GetElementSymbol() != "H")
          return o_neighbor2;
      }
      else if( !neighbor2_is_anomeric ) {
        anomeric_note->type_ = Glycan::WARNING;
        anomeric_note->category_ = Glycan::ANOMERIC;
        anomeric_carbons_status.push_back( "Anomeric carbon assigned to the ring carbon neighboring the ring oxygen but is not attached to an exocyclic carbon (that is, assuming the monosaccharide is an aldose), Anomeric Carbon is: " + anomeric_carbon->GetName() );
        if (o_neighbor2->GetElementSymbol() != "H")
          return o_neighbor1;
      }
      // @TODO August 8, 2017 Davis/Lachele - Once we get to a point where we read in mmcif definitions, we need to use it to determine the Anomeric Carbon.
      anomeric_note->type_ = Glycan::WARNING;
      anomeric_note->category_ = Glycan::ANOMERIC;
      anomeric_carbons_status.push_back( "Anomeric carbon assigned to the first ring carbon (" + o_neighbor1->GetName() + "), Anomeric Carbon is: " + anomeric_carbon->GetName() );
      if (o_neighbor2->GetElementSymbol() != "H")
        return o_neighbor1;
    }
  }
  return NULL;
}

std::vector<std::string> Glycan::Monosaccharide::GetSideGroupOrientations(MolecularModeling::Assembly* this_assembly)
{
  //9/14/18 Removed side atom initialization in this function and either moved it to Yao's InitiateDetectionOfCompleteSideGroupAtoms()
  // Dave
  int local_debug = -1;
  std::vector<std::string> orientations = std::vector<std::string>();
  if(!cycle_atoms_.empty())
  {
    int model_index_ = this_assembly->GetModelIndex();
    std::vector<std::vector<MolecularModeling::Atom*>> side_atoms = std::vector<std::vector<MolecularModeling::Atom*>>();
    std::vector<MolecularModeling::Atom*> default_atom_vector = std::vector<MolecularModeling::Atom*>(3, NULL);
    for(std::vector<MolecularModeling::Atom*>::iterator it = cycle_atoms_.begin(); it != cycle_atoms_.end() - 1; it++) ///iterate on cycle atoms except the oxygen in the ring
    {
      orientations.push_back("N");
      side_atoms.push_back(default_atom_vector);
      unsigned int index = distance(cycle_atoms_.begin(), it);
      MolecularModeling::Atom* prev_atom;
      MolecularModeling::Atom* current_atom = (*it);
      MolecularModeling::Atom* next_atom;
      if(index == 0)///if the current atom is the anomeric atom
      {
        prev_atom = cycle_atoms_.at(cycle_atoms_.size() - 1); ///previous atom is the oxygen(last atom of the sorted cycle)
      }
      else
      {
        prev_atom = cycle_atoms_.at(index - 1);
      }
      next_atom = cycle_atoms_.at(index + 1);

      ///Calculating the plane based on the two ring neighbors of the current atom
      GeometryTopology::Coordinate prev_atom_coord = GeometryTopology::Coordinate(*prev_atom->GetCoordinates().at(model_index_));
      GeometryTopology::Coordinate current_atom_coord = GeometryTopology::Coordinate(*current_atom->GetCoordinates().at(model_index_));
      GeometryTopology::Coordinate next_atom_coord = GeometryTopology::Coordinate(*next_atom->GetCoordinates().at(model_index_));
      prev_atom_coord.operator -(current_atom_coord) ;
      next_atom_coord.operator -(current_atom_coord) ;
      GeometryTopology::Plane plane = GeometryTopology::Plane();
      plane.SetV1(prev_atom_coord);
      plane.SetV2(next_atom_coord);
      GeometryTopology::Coordinate normal_v = plane.GetUnitNormalVector();

      ///Calculating the orientation of the side atoms
      MolecularModeling::AtomNode* node = current_atom->GetNode();
      std::vector<MolecularModeling::Atom*> neighbors = node->GetNodeNeighbors();
      int not_h_neighbors = 0;
      for(std::vector<MolecularModeling::Atom*>::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
      {
        MolecularModeling::Atom* neighbor = (*it1);
        std::string neighbor_id = neighbor->GetId();
        if(cycle_atoms_str_.find(neighbor_id) == std::string::npos) ///if not one of the cycle atoms
        {
          if(neighbor->GetName().at(0) != 'H') ///deoxy check
          {
            not_h_neighbors++;
          }
        }
      }
      if(not_h_neighbors != 0)
      {
        for(std::vector<MolecularModeling::Atom*>::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
        {
          MolecularModeling::Atom* neighbor = (*it1);
          std::string neighbor_id = neighbor->GetId();
          if(cycle_atoms_str_.find(neighbor_id) == std::string::npos) ///if not one of the cycle atoms
          {
            // if(neighbor->GetName().at(0) != 'H') ///deoxy check
            //     not_h_neighbors++;
            GeometryTopology::Coordinate side_atom_coord = GeometryTopology::Coordinate(*neighbor->GetCoordinates().at(model_index_));
            side_atom_coord.operator -(current_atom_coord);
            side_atom_coord.Normalize();
            double theta = acos(normal_v.DotProduct(side_atom_coord));

            if(index == 0 && neighbor->GetElementSymbol() == "C")///if anomeric atom has a non-ring carbon neighbor
            {
              if(orientations.at(index).compare("N") == 0) ///if the position of non-ring oxygen or nitrogen hasn't been set yet
              {
                side_atoms.at(index).at(0) = neighbor;
                neighbor -> SetIsSideChain(true);
                neighbor->SetIsExocyclicCarbon(true);
                if(theta > (gmml::PI_RADIAN/2))
                {
                  orientations.at(index) = "-1D";
                }
                else
                {
                  orientations.at(index) = "-1U";
                }
                continue;
              }
              else
              { ///position of non-ring oxygen or nitrogen + the non-ring carbon
                std::stringstream ss;
                if(theta > (gmml::PI_RADIAN/2))
                {
                  ss << orientations.at(index) << "-1D";
                  orientations.at(index) = ss.str();
                  side_atoms.at(index).at(0) = neighbor;
                  neighbor -> SetIsSideChain(true);
                }
                else
                {
                  ss << orientations.at(index) << "-1U";
                  orientations.at(index) = ss.str();
                  side_atoms.at(index).at(0) = neighbor;
                  neighbor -> SetIsSideChain(true);
                }
                break;
              }
            }
            //else if(neighbor_id.at(0) == 'O' || neighbor_id.at(0) == 'N')///if neighbor is a non-ring oxygen or nitrogen
	    //Yao: After talking to Dave. Replace search based on atom id to element symbol
            else if(neighbor->GetElementSymbol() == "O" || neighbor->GetElementSymbol() == "N")///if neighbor is a non-ring oxygen or nitrogen
            {
              if(index == 0)///current atom is anomeric
              {
                if(orientations.at(index).compare("N") == 0) ///if the position of non-ring carbon neighbor (if exist) hasn't been set yet
                {
                  if(theta > (gmml::PI_RADIAN/2))
                  {
                    orientations.at(index) = "D";
                    side_atoms.at(index).at(1) = neighbor;
                    neighbor -> SetIsSideChain(true);
                  }
                  else
                  {
                    orientations.at(index) = "U";
                    side_atoms.at(index).at(1) = neighbor;
                    neighbor -> SetIsSideChain(true);
                  }
                  continue;
                }
                else
                { ///position of non-ring oxygen or nitrogen + the non-ring carbon
                  std::stringstream ss;
                  if(theta > (gmml::PI_RADIAN/2))
                  {
                    ss << "D" << orientations.at(index);
                    orientations.at(index) = ss.str();
                    side_atoms.at(index).at(1) = neighbor;
                    neighbor -> SetIsSideChain(true);
                  }
                  else
                  {
                    ss << "U" << orientations.at(index);
                    orientations.at(index) = ss.str();
                    side_atoms.at(index).at(1) = neighbor;
                    neighbor -> SetIsSideChain(true);
                  }
                  break;
                }
              }
              else
              {
                if(theta > (gmml::PI_RADIAN/2))
                {
                  orientations.at(index) = "D";
                  side_atoms.at(index).at(1) = neighbor;
                  neighbor -> SetIsSideChain(true);
                }
                else
                {
                  orientations.at(index) = "U";
                  side_atoms.at(index).at(1) = neighbor;
                  neighbor -> SetIsSideChain(true);
                }
                break;
              }
            }
            else if(index == cycle_atoms_.size() - 2 && neighbor->GetElementSymbol() == "C")///if the last ring carbon has a non-ring carbon neighbor
            {

              //Add +n exocyclic flags TODO for Dave



              ///Check if neighbor of neighbor is oxygen or nitrogen
              MolecularModeling::AtomNode* neighbor_node = neighbor->GetNode();
              std::vector<MolecularModeling::Atom*> neighbors_of_neighbor = neighbor_node->GetNodeNeighbors();
              int o_neighbors = 0;
              int not_h_neighbors = 0;
              for(std::vector<MolecularModeling::Atom*>::iterator it2 = neighbors_of_neighbor.begin(); it2 != neighbors_of_neighbor.end(); it2++)
              {
                MolecularModeling::Atom* neighbor_of_neighbor = (*it2);
                std::string neighbor_of_neighbor_id = neighbor_of_neighbor->GetId();
                if((cycle_atoms_str_.find(neighbor_of_neighbor_id) == std::string::npos) &&
                  //((neighbor->GetElementSymbol() == "O") || (neighbor->GetElementSymbol() == "N"))) ///if neighbor of neighbor is a non-ring oxygen or nitrogen  Yao: should be "neighbor_of_neighbor"
                  ((neighbor_of_neighbor->GetElementSymbol() == "O") || (neighbor_of_neighbor->GetElementSymbol() == "N"))) ///if neighbor of neighbor is a non-ring oxygen or nitrogen  Yao: Correction
                {
                  o_neighbors++;
                }
                if(cycle_atoms_str_.find(neighbor_of_neighbor_id) == std::string::npos &&
                  //(neighbor->GetElementSymbol() != "H"))///if neighbor of neighbor is any non-ring atom other than hydrogen  Yao: should be "neighbor_of_neighbor"
                  (neighbor_of_neighbor->GetElementSymbol() != "H"))///if neighbor of neighbor is any non-ring atom other than hydrogen
                {
                  not_h_neighbors++;
                }
              }
              if (o_neighbors >= 1)
              {
                if(theta > (gmml::PI_RADIAN/2))
                {
                  orientations.at(index) = "D";
                  side_atoms.at(index).at(0) = neighbor;
                  neighbor -> SetIsSideChain(true);
                }
                else
                {
                  orientations.at(index) = "U";
                  side_atoms.at(index).at(0) = neighbor;
                  neighbor -> SetIsSideChain(true);
                }
              }
              else if(not_h_neighbors == 0)///Type Deoxy
              {
                if(theta > (gmml::PI_RADIAN/2))
                {
                  orientations.at(index) = "Dd";
                  side_atoms.at(index).at(0) = neighbor;
                  neighbor -> SetIsSideChain(true);
                }
                else
                {
                  orientations.at(index) = "Ud";
                  side_atoms.at(index).at(0) = neighbor;
                  neighbor -> SetIsSideChain(true);
                }
                break;
              }
              if(orientations.at(index).compare("N") != 0)
                break;
            }
          }
        }
      }
    }
    side_atoms_ = side_atoms;


    // if ((local_debug > 0) && (side_atoms.length() > 0))
    // {
    //   for (std::vector<std::vector<MolecularModeling::Atom*>>::iterator sidestd::vector<MolecularModeling::Atom*>IT = side_atoms.begin(); sidestd::vector<MolecularModeling::Atom*>IT != side_atoms.end(); sidestd::vector<MolecularModeling::Atom*>IT ++)
    //   {
    //     std::vector<MolecularModeling::Atom*> sidestd::vector<MolecularModeling::Atom*> = (*sidestd::vector<MolecularModeling::Atom*>IT);
    //     for (std::vector<MolecularModeling::Atom*>::iterator sideAtomIT = sidestd::vector<MolecularModeling::Atom*>.begin(); sideAtomIT != sidestd::vector<MolecularModeling::Atom*>.end(); sideAtomIT++)
    //     {
    //       MolecularModeling::Atom* sideAtom = (*sideAtomIT);
    //       std::stringstream debugStr;
    //       debugStr << "Side atom: " << sideAtom->GetName();
    //       gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
    //     }
    //   }
    // }
  }
  return orientations;
}

void Glycan::Monosaccharide::InitiateDetectionOfCompleteSideGroupAtoms ()
{
    // for (std::vector<Glycan::Monosaccharide*>::iterator it1= monos.begin(); it1 != monos.end(); it1++){
    //     Glycan::Monosaccharide* mono = *it1;
    //     std::vector<MolecularModeling::Atom*>& cycle_atoms = cycle_atoms_;

  for (std::vector<std::vector<MolecularModeling::Atom*> >::iterator it = side_atoms_.begin(); it != side_atoms_.end(); it++)
  {
    std::vector<MolecularModeling::Atom*>& SideAtomArm = *it;
    std::vector<MolecularModeling::Atom*> all_plus_one_side_atoms = std::vector<MolecularModeling::Atom*>();
    for (std::vector<MolecularModeling::Atom*>::iterator it2 = SideAtomArm.begin(); it2 != SideAtomArm.end(); it2++)
    {
      MolecularModeling::Atom* side_atom = *it2;
      if (side_atom != NULL)
      {
	MolecularModeling::AtomVector side_atom_neighbors = side_atom->GetNode()->GetNodeNeighbors();
	for (MolecularModeling::AtomVector::iterator atom_it = side_atom_neighbors.begin(); atom_it != side_atom_neighbors.end(); atom_it++){
	  MolecularModeling::Atom* neighbor = *atom_it;
	    if (neighbor->GetIsCycle()){
              all_plus_one_side_atoms.push_back(side_atom);
	      break;
	    }
	}
      }
    }

    for (std::vector<MolecularModeling::Atom*>::iterator it3= all_plus_one_side_atoms.begin(); it3 != all_plus_one_side_atoms.end(); it3++ )
    {
      MolecularModeling::Atom* plus_one_side_atom = *it3;
      std::vector<MolecularModeling::Atom*> visited_atoms = std::vector<MolecularModeling::Atom*>();
      if ( this->CheckIfPlusOneSideAtomBelongsToCurrentMonosaccharide(SideAtomArm, cycle_atoms_, plus_one_side_atom) )
      {
        this->SetCompleteSideGroupAtoms(SideAtomArm, plus_one_side_atom, cycle_atoms_ ,visited_atoms);
      }
    }
  }
    // }//for
}//InitiateDetectionOfCompleteSideGroupAtoms

bool Glycan::Monosaccharide::CheckIfPlusOneSideAtomBelongsToCurrentMonosaccharide (std::vector<MolecularModeling::Atom*>& SideAtomArm, std::vector<MolecularModeling::Atom*>& cycle_atoms, MolecularModeling::Atom* working_atom )
{
  bool belongs_to_current_monosaccharide = true;
  std::vector<MolecularModeling::Atom*> working_node_neighbors = working_atom -> GetNode() -> GetNodeNeighbors();
  std::vector<MolecularModeling::Atom*> attached_anomeric_carbons = std::vector<MolecularModeling::Atom*>();
  for (std::vector<MolecularModeling::Atom*>::iterator it = working_node_neighbors.begin(); it != working_node_neighbors.end(); it++)
  {
    if (*it != NULL)
    {
      MolecularModeling::Atom* working_node_neighbor = *it;
      if (working_node_neighbor->GetIsAnomericCarbon())
      {
        attached_anomeric_carbons.push_back(working_node_neighbor);
      }
    }
  }//for

  if (attached_anomeric_carbons.size() == 2)
  {
  //TODO:figure out which branch is shorter
  }

  if (attached_anomeric_carbons.size() == 1)
  {
    if (std::find(cycle_atoms.begin(),cycle_atoms.end(),attached_anomeric_carbons.at(0)) != cycle_atoms.end())
    { //if the anomeric carbon is from the current monosaccharide
      bool sidechain_is_terminal = true;
      std::vector<MolecularModeling::Atom*> cycle_and_visisted_atoms = cycle_atoms;
      this->CheckIfSideChainIsTerminal(working_atom,cycle_and_visisted_atoms,sidechain_is_terminal);
      if (!sidechain_is_terminal )
      {
        SideAtomArm.erase(std::find(SideAtomArm.begin(), SideAtomArm.end(), working_atom));	//Unless it's terminal, this oxygen belongs to the other monosaccharide.
        belongs_to_current_monosaccharide = false;
      }
    }
  }
  //If number of attached anomeric carbon is zero, then this atom belongs to current side chain, and it is already. No need to do anything!
  return belongs_to_current_monosaccharide;
}


void Glycan::Monosaccharide::SetCompleteSideGroupAtoms(std::vector<MolecularModeling::Atom*>& SideAtomArm, MolecularModeling::Atom* working_atom, std::vector<MolecularModeling::Atom*>& cycle_atoms, std::vector<MolecularModeling::Atom*>& visited_atoms)
{
  std::vector<MolecularModeling::Atom*> working_node_neighbors = working_atom -> GetNode() -> GetNodeNeighbors();
  if (working_atom->GetResidue() == cycle_atoms[0]->GetResidue() && !working_atom->GetIsCycle() && std::find(visited_atoms.begin(), visited_atoms.end(), working_atom) == visited_atoms.end()){
      if (std::find(SideAtomArm.begin(), SideAtomArm.end(), working_atom) == SideAtomArm.end()){
          SideAtomArm.push_back(working_atom);
      }
  }
  /*for (std::vector<MolecularModeling::Atom*>::iterator it = working_node_neighbors.begin(); it != working_node_neighbors.end(); it++)
  {  //Detect atoms that should be added to the current side chain arm.
    MolecularModeling::Atom* working_node_neighbor = *it;
    if (working_node_neighbor != NULL)
    {
      if ((!working_node_neighbor->GetIsCycle()) &&
          (std::find(visited_atoms.begin(), visited_atoms.end(), working_node_neighbor) == visited_atoms.end()) &&
          (!working_node_neighbor->GetResidue()->CheckIfProtein()) )
      {//If not part of another side chain && not part of cycle && not visited && is not protein (NLN,LNK etc)
        std::vector<MolecularModeling::Atom*> working_node_neighbor_neighbors = std::vector<MolecularModeling::Atom*>();
        std::string working_neighbor_name = working_node_neighbor-> GetName();
        if (working_neighbor_name.substr(0,1) == "O" ||
            working_neighbor_name.substr(0,1) == "N" ||
            working_neighbor_name.substr(0,1) == "S")
        { //if this neighbor is an O,N or S
          std::vector<MolecularModeling::Atom*> attached_anomeric_carbons = std::vector<MolecularModeling::Atom*>();
          for (std::vector<MolecularModeling::Atom*>::iterator it2 = working_node_neighbor_neighbors.begin(); it2 != working_node_neighbor_neighbors.end(); it2++)
          {
            if (*it2 != NULL )
            {
              if ( (*it2)-> GetIsAnomericCarbon() )
              {
                attached_anomeric_carbons.push_back(*it2); //Check how many anomeric carbons this working neighbor is attacheed to, should be either 0,1,2.
              }
            }
          }
          if (attached_anomeric_carbons.size() == 1)
          {  // if attached to only one anomeric carbon
            if (attached_anomeric_carbons.front() != working_atom)
            {  // if the working atom IS NOT the anomeric carbon the neighbor is attached to,the current side chain gets this oxygen.
              SideAtomArm.push_back(working_node_neighbor);
              working_node_neighbor-> SetIsSideChain(true);
            }
            // if the workin atom IS INDEED the anomeric carbon the neighbor is attached to, the current side chain doesn't get this oxygen, unless this is a terminal hydroxyl.
            else
            {
              std::vector<MolecularModeling::Atom*> cycle_and_visited_atoms = std::vector<MolecularModeling::Atom*>();
              cycle_and_visited_atoms.push_back(attached_anomeric_carbons.front() );
              bool side_chain_is_terminal = true;
              CheckIfSideChainIsTerminal (working_node_neighbor, cycle_and_visited_atoms, side_chain_is_terminal);
              if (side_chain_is_terminal)
              {
                SideAtomArm.push_back(working_node_neighbor);
                working_node_neighbor-> SetIsSideChain(true);
              }
            }
          }
          if (attached_anomeric_carbons.size() == 0)
          {
            SideAtomArm.push_back(working_node_neighbor);
            working_node_neighbor-> SetIsSideChain(true);
          }
        } // if is oxygen,nitrogen or sulfur
        else
        {  // if this neighbor is not oxygen,nitrogen or sulfur, the current side chain should get this atom
          SideAtomArm.push_back(working_node_neighbor);
          working_node_neighbor-> SetIsSideChain(true);
        }
      }//if not cycle
    }
  }*///for
  visited_atoms.push_back(working_atom);
  for (std::vector<MolecularModeling::Atom*>::iterator it3 = working_node_neighbors.begin(); it3 != working_node_neighbors.end(); it3++)
  {
    if (*it3 != NULL)
    {
      MolecularModeling::Atom* working_node_neighbor = *it3;
      //if a node neighbor is not part of a ring,not visited, and not part on protein, change working atom to this neighbor, and start a new recursion call.
      if (working_node_neighbor != working_atom &&
          !working_node_neighbor->GetIsCycle() &&
          std::find(visited_atoms.begin(), visited_atoms.end(), working_node_neighbor) == visited_atoms.end() &&
          !working_node_neighbor->GetResidue()->CheckIfProtein())
      {
        working_atom = working_node_neighbor;
        SetCompleteSideGroupAtoms(SideAtomArm, working_atom, cycle_atoms ,visited_atoms);
      }//if
    }//if
  }//for
  return;
}//SetCompleteSideGroupAtoms

void Glycan::Monosaccharide::CheckIfSideChainIsTerminal(MolecularModeling::Atom* starting_atom, std::vector<MolecularModeling::Atom*> & cycle_and_visited_atoms, bool & is_terminal)
{
  cycle_and_visited_atoms.push_back(starting_atom);
  std::vector<MolecularModeling::Atom*> starting_atom_node_neighbors = starting_atom->GetNode()->GetNodeNeighbors();
  for (std::vector<MolecularModeling::Atom*>::iterator it = starting_atom_node_neighbors.begin(); it != starting_atom_node_neighbors.end(); it++)
  {
    if (*it != NULL)
    {
      MolecularModeling::Atom* current_atom_node_neighbor = *it;
      if (std::find(cycle_and_visited_atoms.begin(), cycle_and_visited_atoms.end(), current_atom_node_neighbor) == cycle_and_visited_atoms.end() )
      {
        if (current_atom_node_neighbor->GetIsCycle() || current_atom_node_neighbor->GetResidue()->CheckIfProtein() )
        {
          is_terminal = false;
        }
        if(std::find(cycle_and_visited_atoms.begin(), cycle_and_visited_atoms.end(), *it) == cycle_and_visited_atoms.end() &&
          !(*it)->GetResidue()->CheckIfProtein()  && is_terminal )
        {
          starting_atom = *it;
          CheckIfSideChainIsTerminal(starting_atom,cycle_and_visited_atoms,is_terminal);
        }//if
      }
    }
  }//for
  return;
} //CheckIfSideChainIsTerminal

void Glycan::Monosaccharide::ExtractDerivatives(MolecularModeling::Assembly* this_assembly)
{
  int local_debug = -1;
  if (local_debug > 0)
  {
    std::stringstream debugStr;
    debugStr << "On monosaccharide: " << mono_id_;
    gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
  }
  for(std::vector<MolecularModeling::Atom*>::iterator it = cycle_atoms_.begin(); it != cycle_atoms_.end() - 1; it++) ///iterate on cycle atoms except the oxygen of the ring
  {
    int index = distance(cycle_atoms_.begin(), it);
    MolecularModeling::Atom* target = (*it);
    std::string key = "";
    std::string value = "";
    std::vector<MolecularModeling::Atom*> t_neighbors = target->GetNode()->GetNodeNeighbors();
    if (local_debug > 0)
    {
      std::stringstream debugStr;
      debugStr << "On atom: " << target->GetName();
      gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
    }
    for(std::vector<MolecularModeling::Atom*>::iterator it1 = t_neighbors.begin(); it1 != t_neighbors.end(); it1++)
    {
      std::vector<MolecularModeling::Atom*> pattern_atoms = std::vector<MolecularModeling::Atom*>();
      MolecularModeling::Atom* t_neighbor = (*it1);
      if(local_debug > 0)
      {
        std::stringstream debugStr;
        debugStr << "On neighbor: " << t_neighbor->GetName();
        gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
      }
      if(t_neighbor->GetName().at(0) == 'N' && cycle_atoms_str_.find(t_neighbor->GetId()) == std::string::npos)///check formulas with nitrogen
      {
        if((value = this_assembly->CheckxC_N(target, cycle_atoms_str_/*, pattern_atoms*/)).compare("") != 0)///xCH-N
          break;
        if((value = this_assembly->CheckxC_NxO_CO_C(target, cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH3
          break;
        if((value = this_assembly->CheckxC_NxO_CO_CO(target, cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH2OH
          break;
        if((value = this_assembly->CheckxC_NxO_SO3(target, cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-SO3
          break;
        if((value = this_assembly->CheckxC_NxO_PO3(target, cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-PO3
          break;
        if((value = this_assembly->CheckxC_NxO_C(target, cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-CH3
          break;
      }
      if(t_neighbor->GetName().at(0) == 'O' && cycle_atoms_str_.find(t_neighbor->GetId()) == std::string::npos)///check formulas with oxygen
      {
        if((value = this_assembly->CheckxC_NxO_CO_C(target, cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH3
          break;
        if((value = this_assembly->CheckxC_NxO_CO_CO(target, cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH2OH
          break;
        if((value = this_assembly->CheckxC_NxO_SO3(target, cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-SO3
          break;
        if((value = this_assembly->CheckxC_NxO_PO3(target, cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-PO3
          break;
        if((value = this_assembly->CheckxC_NxO_C(target, cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-CH3
          break;
        if((value = this_assembly->CheckxCOO(target, cycle_atoms_str_/*, pattern_atoms*/)).compare("") != 0)///xC-(O,O) and xC-(O,OH)
          break;
      }
    }
    if(value.compare("") != 0)///if any pattern matched add it to the index-derivative map
    {
      if(index == 0)
        key = "a";
      else
        key = gmml::ConvertT(index + 1);
      derivatives_map_.push_back({key, value});
    }
    else if(index != 0)
    {
      key = gmml::ConvertT(index + 1);

      //TODO add function to add formula of unknown derivative
      value = GetFormula(target);
      if(value.compare("") != 0)
      {
        if(((index == 3) && (this->sugar_name_.ring_type_.compare("F") == 0)) ||
           ((index == 4) && (this->sugar_name_.ring_type_.compare("P") == 0)))
        {
          if ((value != "C1O1") && (this->sugar_name_.ring_type_.compare("P") == 0))
          {
            // gmml::log(__LINE__, __FILE__, gmml::INF, key);
            // gmml::log(__LINE__, __FILE__, gmml::INF, value);
            unknown_derivatives_.push_back({key, value});
            derivatives_map_.push_back({key, ""});
          }
          else if ((value != "C2O2") && (this->sugar_name_.ring_type_.compare("F") == 0))
          {
            unknown_derivatives_.push_back({key, value});
            derivatives_map_.push_back({key, ""});
          }
        }
        else if(index != 4 && index !=3 && index != 0)
        {
          gmml::log(__LINE__, __FILE__, gmml::INF, key);
          gmml::log(__LINE__, __FILE__, gmml::INF, std::to_string(index));
          gmml::log(__LINE__, __FILE__, gmml::INF, value);
          unknown_derivatives_.push_back({key, value});
          derivatives_map_.push_back({key, ""});
        }
      }
    }
  }
  for(std::vector<std::vector<MolecularModeling::Atom*>>::iterator it = side_atoms_.begin(); it != side_atoms_.end(); it++) ///iterate on side atoms
  {
    int index = distance(side_atoms_.begin(), it);
    int side_branch_last_carbon_index = 0;
    std::vector<MolecularModeling::Atom*> sides = (*it);
    MolecularModeling::Atom* target = NULL;
    std::string key = "";
    std::string value = "";
    // if(it == side_atoms_.begin())///side atoms of anomeric carbon
    // {
    //     if(sides.at(0) != NULL)
    //         target = sides.at(0);///first index of each side is for carbon atoms in the std::vector<std::vector<MolecularModeling::Atom*>> structure
    // }
    // //WHERE ARE ALL OF THE OTHER SIDE ATOMS?!?!?!?!?!?!?!?!?!??!
    if(it == side_atoms_.end() - 1)//side atoms of last carbon of the ring
    {
        if(sides.at(0) != NULL)
        {
            for(side_branch_last_carbon_index = sides.size() - 1; sides.at(side_branch_last_carbon_index) == NULL; side_branch_last_carbon_index-- )
            {
    //           target = sides.at(side_branch_last_carbon_index);
            }

        }
    }
    // else
    // {
    if(sides.at(0) != NULL)
      target = sides.at(0);
    // }
    if ((local_debug > 0) && (target != NULL))
    {
      std::stringstream debugStr;
      debugStr << "On side atom: " << target->GetName();
      gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
    }
    if(target != NULL)
    {
      std::vector<MolecularModeling::Atom*> t_neighbors = target->GetNode()->GetNodeNeighbors();
      for(std::vector<MolecularModeling::Atom*>::iterator it1 = t_neighbors.begin(); it1 != t_neighbors.end(); it1++)
      {
        std::vector<MolecularModeling::Atom*> pattern_atoms = std::vector<MolecularModeling::Atom*>();
        MolecularModeling::Atom* t_neighbor = (*it1);
        if (local_debug > 0)
        {
          std::stringstream debugStr;
          debugStr << "On side atom neighbor: " << t_neighbor->GetName();
          gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
        }
        if(t_neighbor->GetName().at(0) == 'N' && cycle_atoms_str_.find(t_neighbor->GetId()) == std::string::npos)///check formulas with nitrogen
        {
          if((value = this_assembly->CheckxC_N(target, cycle_atoms_str_/*, pattern_atoms*/)).compare("") != 0)///xCH-N
            break;
          if((value = this_assembly->CheckxC_NxO_CO_C(target, cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH3
            break;
          if((value = this_assembly->CheckxC_NxO_CO_CO(target, cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH2OH
            break;
          if((value = this_assembly->CheckxC_NxO_SO3(target, cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-SO3
            break;
          if((value = this_assembly->CheckxC_NxO_PO3(target, cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-PO3
            break;
          if((value = this_assembly->CheckxC_NxO_C(target, cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-CH3
            break;
          // if((value = this_assembly->CheckxC_NxO_CH3C_COO(target, cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-CH3C-C-(O,O/OH)
          //   break;
        }
        if(t_neighbor->GetName().at(0) == 'O' && cycle_atoms_str_.find(t_neighbor->GetId()) == std::string::npos)///check formulas with oxygen
        {
          if((value = this_assembly->CheckxC_NxO_CO_C(target, cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH3
            break;
          if((value = this_assembly->CheckxC_NxO_CO_CO(target, cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH2OH
            break;
          if((value = this_assembly->CheckxC_NxO_SO3(target, cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-SO3
            break;
          if((value = this_assembly->CheckxC_NxO_PO3(target, cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-PO3
            break;
          if((value = this_assembly->CheckxC_NxO_C(target, cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-CH3
            break;
          if((value = this_assembly->CheckxCOO(target, cycle_atoms_str_/*, pattern_atoms*/)).compare("") != 0)///xC-(O,O) and xC-(O,OH)
            break;
        }
      }
      if(value.compare("") != 0)///if any pattern matched add it to the index-derivative map
      {
        if(index == 0)
          key = "-1";
        else
        {
          switch (side_branch_last_carbon_index)
          {
            case 0:
              key = "+1";
              break;
            case 1:
              key = "+2";
              break;
            case 2:
              key = "+3";
              break;
          }
        }
        derivatives_map_.push_back({key, value});
      }
      else
      {
        if(index == 0)
          key = "-1";
        else
        {
          switch (side_branch_last_carbon_index)
          {
            case 0:
              key = "+1";
              break;
            case 1:
              key = "+2";
              break;
            case 2:
              key = "+3";
              break;
          }
        }
        value = GetFormula(target);
        if(value.compare("") != 0)
        {
          if(key != "-1")
          {
            if(((index == 3) && (this->sugar_name_.ring_type_.compare("F") == 0)) ||
               ((index == 4) && (this->sugar_name_.ring_type_.compare("P") == 0)))
            {
              if (value != "C1O1")
              {
                unknown_derivatives_.push_back({key, value});
                derivatives_map_.push_back({key, ""});
              }
            }
          }
          else
          {
            unknown_derivatives_.push_back({key, value});
            derivatives_map_.push_back({key, ""});
          }
        }
      }
    }
  }
}

std::string Glycan::Monosaccharide::GetFormula(MolecularModeling::Atom* target)//
{
  std::string thisDerivative = "";
  std::vector<std::pair<std::string, int> > elementVector;
  std::vector<MolecularModeling::Atom*> t_neighbors = target->GetNode()->GetNodeNeighbors();
  target->GetNode()->SetIsVisited(true);
  for(std::vector<MolecularModeling::Atom*>::iterator it1 = t_neighbors.begin(); it1 != t_neighbors.end(); it1++)
  {
    MolecularModeling::Atom* t_neighbor = (*it1);
    if((cycle_atoms_str_.find(t_neighbor->GetId()) == std::string::npos) && (t_neighbor->GetElementSymbol() != "H"))//this neighbor isn't in the ring or the hydrogen
    {
      CountElements(t_neighbor, elementVector);
    }
  }
  for(std::vector<std::pair<std::string, int> >::iterator it = elementVector.begin(); it != elementVector.end(); it++)
  {
    std::string thisElement = (*it).first;
    int thisElementCount = (*it).second;
    if((thisElementCount > 0) && (thisElement != "H"))//Not counting Hydrogens
    {
      thisDerivative = thisDerivative + thisElement + std::to_string(thisElementCount);
    }
  }
  if(thisDerivative != "O1")
  {
    // gmml::log(__LINE__, __FILE__, gmml::INF, thisDerivative);
    return thisDerivative;
  }
  else
  {
    return "";
  }
}

void Glycan::Monosaccharide::CountElements(MolecularModeling::Atom* thisAtom, std::vector<std::pair<std::string, int> >& elementVector)
{
  if(elementVector.size() == 0)
  {
    elementVector = {{"C", 0}, {"H", 0}, {"Ac", 0},{"Ag", 0},{"Al", 0},{"Am", 0},{"Ar", 0},{"As", 0},{"At", 0},{"Au", 0},
                     {"B", 0}, {"Ba", 0},{"Be", 0},{"Bh", 0},{"Bi", 0},{"Bk", 0},{"Br", 0},{"Ca", 0},{"Cd", 0},{"Ce", 0},
                     {"Cf", 0},{"Cl", 0},{"Cm", 0},{"Co", 0},{"Cr", 0},{"Cs", 0},{"Cu", 0},{"Dd", 0},{"Dy", 0},{"Er", 0},
                     {"Es", 0},{"Eu", 0},{"F", 0}, {"Fe", 0},{"Fm", 0},{"Fr", 0},{"Ga", 0},{"Gd", 0},{"Ge", 0},{"He", 0},
                     {"Hf", 0},{"Hg", 0},{"Ho", 0},{"Hs", 0},{"I", 0}, {"In", 0},{"Ir", 0},{"K", 0}, {"Kr", 0},{"La", 0},
                     {"Li", 0},{"Lr", 0},{"Lu", 0},{"Md", 0},{"Mg", 0},{"Mn", 0},{"Mo", 0},{"Mt", 0},{"N", 0}, {"Na", 0},
                     {"Nb", 0},{"Nd", 0},{"Ne", 0},{"Ni", 0},{"No", 0},{"Np", 0},{"O", 0}, {"Os", 0},{"P", 0}, {"Pa", 0},
                     {"Pb", 0},{"Pd", 0},{"Pm", 0},{"Po", 0},{"Pr", 0},{"Pt", 0},{"Pu", 0},{"Ra", 0},{"Rb", 0},{"Re", 0},
                     {"Rf", 0},{"Rh", 0},{"Rn", 0},{"Ru", 0},{"S", 0}, {"Sb", 0},{"Sc", 0},{"Se", 0},{"Sg", 0},{"Si", 0},
                     {"Sm", 0},{"Sn", 0},{"Sr", 0},{"Ta", 0},{"Tb", 0},{"Tc", 0},{"Te", 0},{"Th", 0},{"Ti", 0},{"Tl", 0},
                     {"Tm", 0},{"U", 0}, {"V", 0}, {"W", 0}, {"Xe", 0},{"Y", 0}, {"Yb", 0},{"Zn", 0},{"Zr", 0}};
  }
  std::string this_atom_element = thisAtom->GetElementSymbol();
  for(std::vector<std::pair<std::string, int> >::iterator it = elementVector.begin(); it != elementVector.end(); it++)
  {
    if(this_atom_element == (*it).first)
    {
      (*it).second ++;
      // gmml::log(__LINE__, __FILE__, gmml::INF, thisAtom->GetId());
      break;
    }
  }
  thisAtom->GetNode()->SetIsVisited(true);
  std::vector<MolecularModeling::Atom*> thisAtomNeighbors = thisAtom->GetNode()->GetNodeNeighbors();
  if(thisAtomNeighbors.size() > 1)
  {
    for(std::vector<MolecularModeling::Atom*>::iterator it = thisAtomNeighbors.begin(); it != thisAtomNeighbors.end(); it++)
    {
      MolecularModeling::Atom* thisNeighbor = (*it);
      if (!thisNeighbor->GetNode()->GetIsVisited() && !thisNeighbor->GetIsCycle())
      {
        if(!thisNeighbor->GetResidue()->CheckIfProtein())
          CountElements(thisNeighbor, elementVector);
        else
          gmml::log(__LINE__, __FILE__, gmml::ERR, "Ran into a protein counting elements");
      }
    }
  }
}

std::vector<MolecularModeling::Atom*> Glycan::Monosaccharide::ExtractAdditionalSideAtoms()
{
  std::vector<MolecularModeling::Atom*> plus_sides = std::vector<MolecularModeling::Atom*>();
  if(side_atoms_.at(side_atoms_.size() - 1).at(0) != NULL)///if there exist a +1 carbon atom. in side_atoms_ structure (std::vector<std::vector<MolecularModeling::Atom*>>) the first index of the last element is dedicated to +1 atom
  {
    plus_sides.push_back(side_atoms_.at(side_atoms_.size() - 1).at(0));
    std::vector<MolecularModeling::Atom*> plus_one_atom_neighbors = side_atoms_.at(side_atoms_.size() - 1).at(0)->GetNode()->GetNodeNeighbors();
    for(std::vector<MolecularModeling::Atom*>::iterator it1 = plus_one_atom_neighbors.begin(); it1 != plus_one_atom_neighbors.end(); it1++)
    {
      if((*it1)->GetName().at(0) == 'C' && cycle_atoms_str_.find((*it1)->GetId()) == std::string::npos)///+2 carbon atom found
      {
        MolecularModeling::Atom* plus_two = (*it1);
        plus_sides.push_back(plus_two);
        side_atoms_.at(side_atoms_.size() - 1).at(1) = plus_two;///in side_atoms_ structure (std::vector<std::vector<MolecularModeling::Atom*>>) the second index of the last element is dedicated to +2 atom

        std::vector<MolecularModeling::Atom*> plus_two_atom_neighbors = plus_two->GetNode()->GetNodeNeighbors();
        for(std::vector<MolecularModeling::Atom*>::iterator it2 = plus_two_atom_neighbors.begin(); it2 != plus_two_atom_neighbors.end(); it2++)
        {
          MolecularModeling::Atom* plus_three = (*it2);
          if(plus_three->GetName().at(0) == 'C' && plus_three->GetId().compare(side_atoms_.at(side_atoms_.size() - 1).at(0)->GetId()) != 0)///+3 carbon atom found
          {
            plus_sides.push_back(plus_three);
            side_atoms_.at(side_atoms_.size() - 1).at(2) = plus_three;///in side_atoms_ structure (std::vector<std::vector<MolecularModeling::Atom*>>) the third index of the last element is dedicated to +3 atom
            break;
          }
        }
        break;
      }
    }
  }
  return plus_sides;
}

void Glycan::Monosaccharide::GenerateCompleteName(std::vector<MolecularModeling::Atom*> &plus_sides, Glycan::Monosaccharide* this_mono, MolecularModeling::Assembly* this_assembly)
{
  // std::cout << plus_sides.size() << "plus atoms\n";
  // if (!plus_sides.empty())
  // {
    if( plus_sides.size() <= 1 )
    {
      ///COMPLETE NAME GENERATION BASED ON DERIVATIVE MAP
      this_mono->UpdateComplexSugarChemicalCode();
      Glycan::SugarName updated_name = gmml::ComplexSugarNameLookup(chemical_code_->toString());
      if((updated_name.monosaccharide_name_ != "") && (updated_name.monosaccharide_name_ != sugar_name_.monosaccharide_name_))
      {
        sugar_name_ = updated_name;
      }
      this_mono->UpdatePdbCode();
      this_mono->GenerateCompleteSugarName(this_assembly);
    }
    else
    {
      if( plus_sides.size() == 3 )
      {
        std::vector< std::string >::iterator index_it;
        if( ( index_it = find( this_mono->chemical_code_->right_up_.begin(), this_mono->chemical_code_->right_up_.end(), "+1" ) ) != this_mono->chemical_code_->right_up_.end() )
        {
          ///CHECKING R or S
          std::stringstream plus_one;
          std::string orientation = this_assembly->CalculateRSOrientations( this_mono->cycle_atoms_.at( this_mono->cycle_atoms_.size() - 2 ), plus_sides.at( 0 ), plus_sides.at( 1 ) );
          plus_one << "+1" << orientation;
          ( *index_it ) = plus_one.str();
          std::stringstream plus_two;
          orientation = this_assembly->CalculateRSOrientations( plus_sides.at( 0 ), plus_sides.at( 1 ), plus_sides.at( 2 ) );
          plus_two << "+2" << orientation;
          this_mono->chemical_code_->right_up_.push_back( plus_two.str() );
          this_mono->chemical_code_->right_up_.push_back( "+3" );
        }
        else if( ( index_it = find( this_mono->chemical_code_->right_down_.begin(), this_mono->chemical_code_->right_down_.end(), "+1" ) ) != this_mono->chemical_code_->right_down_.end() )
        {
          ///CHECKING R or S
          std::stringstream plus_one;
          std::string orientation = this_assembly->CalculateRSOrientations( this_mono->cycle_atoms_.at( this_mono->cycle_atoms_.size() - 2 ), plus_sides.at( 0 ), plus_sides.at( 1 ) );
          plus_one << "+1" << orientation;
          ( *index_it ) = plus_one.str();
          std::stringstream plus_two;
          orientation = this_assembly->CalculateRSOrientations( plus_sides.at( 0 ), plus_sides.at( 1 ), plus_sides.at( 2 ) );
          plus_two << "+2" << orientation;
          this_mono->chemical_code_->right_down_.push_back( plus_two.str() );
          this_mono->chemical_code_->right_down_.push_back( "+3" );
        }
      }
      else if( plus_sides.size() == 2 )
      {
        std::vector< std::string >::iterator index_it;
        if( ( index_it = find( this_mono->chemical_code_->right_up_.begin(), this_mono->chemical_code_->right_up_.end(), "+1" ) ) != this_mono->chemical_code_->right_up_.end() )
        {
          ///CHECKING R or S
          std::stringstream plus_one;
          std::string orientation = this_assembly->CalculateRSOrientations( this_mono->cycle_atoms_.at( this_mono->cycle_atoms_.size() - 2 ), plus_sides.at( 0 ), plus_sides.at( 1 ) );
          plus_one << "+1" << orientation;
          (*index_it) = plus_one.str();
          this_mono->chemical_code_->right_up_.push_back("+2");
        }
        else if( ( index_it = find( this_mono->chemical_code_->right_down_.begin(), this_mono->chemical_code_->right_down_.end(), "+1" ) ) != this_mono->chemical_code_->right_down_.end() )
        {
          ///CHECKING R or S
          std::stringstream plus_one;
          std::string orientation = this_assembly->CalculateRSOrientations(this_mono->cycle_atoms_.at( this_mono->cycle_atoms_.size() - 2 ), plus_sides.at( 0 ), plus_sides.at( 1 ) );
          plus_one << "+1" << orientation;
          ( *index_it ) = plus_one.str();
          this_mono->chemical_code_->right_down_.push_back( "+2" );
        }
      }
      ///UPDATING CHEMICAL CODE
      this_mono->UpdateComplexSugarChemicalCode();
      Glycan::SugarName updated_name = gmml::ComplexSugarNameLookup(chemical_code_->toString());
      if((updated_name.monosaccharide_name_ != "") && (updated_name.monosaccharide_name_ != sugar_name_.monosaccharide_name_))
      {
        sugar_name_ = updated_name;
      }
      this_mono->UpdatePdbCode();
      // std::cout << "Complex structure side group atoms: " << std::endl;
      // // gmml::log(__LINE__, __FILE__,  gmml::INF, "Complex structure side group atoms: ");
      // for( std::vector< std::vector<MolecularModeling::Atom*> >::iterator it1 = this_mono->side_atoms_.begin(); it1 != this_mono->side_atoms_.end(); it1++ )
      // {
      //   std::stringstream complex_sugar_side;
      //   std::vector<MolecularModeling::Atom*> sides = ( *it1 );
      //   if( it1 == this_mono->side_atoms_.begin() )
      //   {///side atoms of anomeric carbon
      //     if( sides.at( 0 ) != NULL && sides.at( 1 ) != NULL )
      //     {
      //       complex_sugar_side << "[1] -> " << sides.at( 0 )->GetId() << ", " << sides.at( 1 )->GetId();
      //       std::cout << complex_sugar_side.str() << std::endl;
      //       // gmml::log(__LINE__, __FILE__,  gmml::INF, complex_sugar_side.str());
      //     }
      //     else if( sides.at( 1 ) != NULL )
      //     {
      //       complex_sugar_side << "[1] -> " << sides.at( 1 )->GetId();
      //       std::cout << complex_sugar_side.str() << std::endl;
      //       // gmml::log(__LINE__, __FILE__,  gmml::INF, complex_sugar_side.str());
      //     }
      //     else if( sides.at(0) != NULL )
      //     {
      //       complex_sugar_side << "[1] -> " << sides.at( 0 )->GetId();
      //       std::cout << complex_sugar_side.str() << std::endl;
      //       // gmml::log(__LINE__, __FILE__,  gmml::INF, complex_sugar_side.str());
      //     }
      //   }
      //   else if( it1 == this_mono->side_atoms_.end() - 1 )
      //   {//side atoms of last carbon of the ring
      //     complex_sugar_side << "[" << this_mono->cycle_atoms_.size() - 1 << "]";
      //     for( unsigned int i = 0; i < plus_sides.size() ; i++ )
      //     {
      //       complex_sugar_side << " -> " << sides.at( i )->GetId();
      //     }
      //     std::cout << complex_sugar_side.str() << std::endl;
      //     // gmml::log(__LINE__, __FILE__,  gmml::INF, complex_sugar_side.str());
      //   }
      //   else if( sides.at( 1 ) != NULL )
      //   {
      //     int cycle_atom_index = distance( this_mono->side_atoms_.begin(), it1 );
      //     complex_sugar_side << "[" << cycle_atom_index + 1 << "] -> " << sides.at( 1 )->GetId();
      //     std::cout << complex_sugar_side.str() << std::endl;
      //     // gmml::log(__LINE__, __FILE__,  gmml::INF, complex_sugar_side.str());
      //   }
      // }
      // std::cout << std::endl << "Complex sugar chemical code:" << std::endl;
      // // gmml::log(__LINE__, __FILE__,  gmml::INF, "Complex sugar chemical code:");
      // // gmml::log(__LINE__, __FILE__,  gmml::INF, chemical_code_->toString());
      // chemical_code_->Print( std::cout );
      // std::cout << chemical_code_->toString();
      ///FINDING COMPLEX CHEMICAL CODE IN COMPLEX SUGAR NAME LOOKUP TABLE
      this_mono->sugar_name_ = gmml::ComplexSugarNameLookup( this_mono->chemical_code_->toString() );
      // gmml::log(__LINE__, __FILE__, gmml::INF,  std::to_string(plus_sides.size()));
      if( plus_sides.size() == 2 )
      {
        ///COMPLETE NAME GENERATION BASED ON DERIVATIVE MAP
        this_mono->GenerateCompleteSugarName(this_assembly);
        // gmml::log(__LINE__, __FILE__, gmml::INF, "Generating complete sugar name");
      }
    }
  // }
}

void Glycan::Monosaccharide::GenerateCompleteSugarName(MolecularModeling::Assembly* this_assembly)
{
    std::stringstream in_bracket;
    std::stringstream head;
    std::stringstream tail;
    bool minus_one = false;
    // if(std::find_if( derivatives_map_.begin(), derivatives_map_.end(),
    //   [](const std::pair<std::string, std::string>& element){ return element.first == "-1";} ) == derivatives_map_.end())
    // if(derivatives_map_.find("-1") != derivatives_map_.end())
    for(std::vector<std::string>::iterator it = chemical_code_->right_up_.begin(); it != chemical_code_->right_up_.end(); it++)
    {
      std::string key = (*it);

      if(key == "-1")
      {
        // std::cout << "minus 1\n";
        minus_one = true;
      }
    }
    for(std::vector<std::string>::iterator it = chemical_code_->right_down_.begin(); it != chemical_code_->right_down_.end(); it++)
    {
      std::string key = (*it);

      if(key == "-1")
      {
        // std::cout << "minus 1\n";
        minus_one = true;
      }
    }
    // std::cout << derivatives_map_.size() << "derivatives\n";
    for(std::vector<std::pair<std::string, std::string> >::iterator it1 = derivatives_map_.begin(); it1 != derivatives_map_.end(); it1++)
    {
        std::string key = (*it1).first;
        std::string value = (*it1).second;
        // std::cout << key << ": " << value << "\n";
        std::string long_name_pattern = "";
        std::string cond_name_pattern = "";
        std::string long_name_pattern_at_minus_one = "";
        std::string long_name_pattern_at_plus_one = "";
        std::string pattern = "";
        // gmml::log(__LINE__, __FILE__,  gmml::INF, "Naming by pettern below");
        // gmml::log(__LINE__, __FILE__,  gmml::INF, value);
        std::string unknownDerivativePattern = "", unknownDerivativeKey = "";
        for(std::vector<std::pair<std::string, std::string> >::iterator it = this->unknown_derivatives_.begin(); it != this->unknown_derivatives_.end(); it++)
        {
          std::string thisKey = (*it).first;
          std::string thisPattern = (*it).second;
          if(thisKey == key)
          {
            unknownDerivativeKey = thisKey;
            unknownDerivativePattern = thisPattern;
            break;
          }
        }

        if(value.compare("xCH-N") == 0)
        {
            long_name_pattern = "-osamine";
            cond_name_pattern = "N";
            pattern = "CH-N";
            this_assembly->AddModificationRuleOneInfo(key, pattern, this, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        else if(value.compare("xC-N-C=OCH3") == 0)
        {
            long_name_pattern = "N-acetyl-";
            cond_name_pattern = "NAc";
            pattern = "C-N-C=OCH3";
            this_assembly->AddModificationRuleOneInfo(key, pattern, this, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        else if(value.compare("xC-N-C=OCH2OH") == 0)
        {
            long_name_pattern = "N-glycolyl-";
            cond_name_pattern = "NGc";
            pattern = "C-N-C=OCH2OH";
            this_assembly->AddModificationRuleOneInfo(key, pattern, this, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        else if(value.compare("xC-N-SO3") == 0)
        {
            long_name_pattern = "N-sulfo-";
            cond_name_pattern = "NS";
            pattern = "C-N-SO3";
            this_assembly->AddModificationRuleOneInfo(key, pattern, this, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        else if(value.compare("xC-N-PO3") == 0)
        {
            long_name_pattern = "N-phospho-";
            cond_name_pattern = "NP";
            pattern = "C-N-PO3";
            this_assembly->AddModificationRuleOneInfo(key, pattern, this, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        else if(value.compare("xC-N-CH3") == 0)
        {
            long_name_pattern = "N-methyl-";
            cond_name_pattern = "NMe";
            pattern = "C-N-CH3";
            this_assembly->AddModificationRuleOneInfo(key, pattern, this, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        else if(value.compare("xC-O-C=OCH3") == 0)
        {
            long_name_pattern = "-acetyl-";
            cond_name_pattern = "Ac";
            pattern = "C-O-C=OCH3";
            this_assembly->AddDerivativeRuleInfo(key, pattern, this, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        else if(value.compare("xC-O-C=OCH2OH") == 0)
        {
            long_name_pattern = "-glycolyl-";
            cond_name_pattern = "Gc";
            pattern = "C-O-C=OCH2OH";
            this_assembly->AddDerivativeRuleInfo(key, pattern, this, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        else if(value.compare("xC-O-SO3") == 0)
        {
            long_name_pattern = "-sulfo-";
            cond_name_pattern = "S";
            pattern = "C-O-SO3";
            this_assembly->AddDerivativeRuleInfo(key, pattern, this, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        else if(value.compare("xC-O-PO3") == 0)
        {
            long_name_pattern = "-phospho-";
            cond_name_pattern = "P";
            pattern = "C-O-PO3";
            this_assembly->AddDerivativeRuleInfo(key, pattern, this, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        else if(value.compare("xC-O-CH3") == 0)
        {
            long_name_pattern = "-methyl-";
            cond_name_pattern = "Me";
            pattern = "C-O-CH3";
            this_assembly->AddDerivativeRuleInfo(key, pattern, this, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        else if(value.compare("xC-(O,OH)") == 0)
        {
            long_name_pattern_at_minus_one = "-ulosonic acid";
            long_name_pattern_at_plus_one = "-uronic acid";
            cond_name_pattern = "AH";
            pattern = "C-(O,OH)";
            this_assembly->AddModificationRuleTwoInfo(key, pattern, this, long_name_pattern_at_minus_one, long_name_pattern_at_plus_one, cond_name_pattern, tail, minus_one, in_bracket);
        }
        else if(value.compare("xC-(O,O)") == 0)
        {
            long_name_pattern_at_minus_one = "-ulosonate";
            long_name_pattern_at_plus_one = "-uronate";
            cond_name_pattern = "A";
            pattern = "C-(O,O)";
            this_assembly->AddModificationRuleTwoInfo(key, pattern, this, long_name_pattern_at_minus_one, long_name_pattern_at_plus_one, cond_name_pattern, tail, minus_one, in_bracket);
        }
        else if(value.compare("xCHH") == 0)
        {//Deoxy
            long_name_pattern = "-deoxy-";
            cond_name_pattern = "H";
            pattern = "CHH";
            this_assembly->AddDerivativeRuleInfo(key, pattern, this, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        else if(unknownDerivativePattern != "")
        {
          this->on_R_++;
          std::string Rnum;
          Rnum = std::to_string(on_R_);
          long_name_pattern = "-<" + unknownDerivativePattern + ">-";
          cond_name_pattern = "<R" + Rnum + ">";
          pattern = unknownDerivativePattern;
          // gmml::log(__LINE__, __FILE__, gmml::INF, key);
          if(key != "a")
            this_assembly->AddUnknownDerivativeRuleInfo(key, pattern, this, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
          else
          {
            this_assembly->AddDerivativeRuleInfo(key, pattern, this, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
          }
        }
    }
    if(in_bracket.str().size() != 0)
    {
        std::stringstream short_name;
        std::string sn;
        if(sugar_name_.monosaccharide_short_name_.compare("") != 0)
            sn = sugar_name_.monosaccharide_short_name_;
        else
            sn = sugar_name_.monosaccharide_stereochemistry_short_name_;

        ///moving a, b or x to after the bracket: short-name + [...] + a/b/x and removing ", " from the end of bracket stream
        int condensed_name_size = sn.size();
        std::string condensed_name = sn;
        std::string new_name_part1 = "";
        char new_name_part2 = ' ';
        if(condensed_name_size > 0)
        {
          new_name_part1 = condensed_name.substr(0, (condensed_name_size - 1));///short_name
          new_name_part2 = condensed_name.at(condensed_name_size - 1);///a/b/x
        }
        short_name << new_name_part1 << "[" << in_bracket.str().substr(0, in_bracket.str().size() - 1) << "]" << new_name_part2;

        sugar_name_.monosaccharide_short_name_ = short_name.str();
    }
    else if(sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0 && sugar_name_.monosaccharide_short_name_.compare("") == 0)
    {
        sugar_name_.monosaccharide_short_name_ = sugar_name_.monosaccharide_stereochemistry_short_name_;
    }
    std::stringstream long_name;
    if(sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
    {
      if((sugar_name_.monosaccharide_stereochemistry_name_.substr(sugar_name_.monosaccharide_stereochemistry_name_.size() - 3) ==
        "ose") && (tail.str() == "-osamine"))
      {
        long_name << head.str();
        long_name << sugar_name_.monosaccharide_stereochemistry_name_.substr(0,sugar_name_.monosaccharide_stereochemistry_name_.size() - 3);
        long_name << tail.str().substr(1);
        sugar_name_.monosaccharide_stereochemistry_name_ = long_name.str();
      }
      else
      {
        long_name << head.str() << sugar_name_.monosaccharide_stereochemistry_name_ << tail.str();
        sugar_name_.monosaccharide_stereochemistry_name_ = long_name.str();
      }
    }
}

void Glycan::Monosaccharide::UpdatePdbCode()
{
  std::string code = chemical_code_->toString();
  for( int i = 0; i < gmml::COMPLEXSUGARNAMELOOKUPSIZE; i++ )
  {
    if( code.compare( gmml::COMPLEXSUGARNAMELOOKUP[ i ].chemical_code_string_ ) == 0 )
    {
      sugar_name_.pdb_code_ = gmml::COMPLEXSUGARNAMELOOKUP[ i ].pdb_code_;
      return;
    }
  }
  for( int i = 0; i < gmml::SUGARNAMELOOKUPSIZE; i++ )
  {
    if( code.compare( gmml::SUGARNAMELOOKUP[ i ].chemical_code_string_ ) == 0 )
    {
      sugar_name_.pdb_code_ = gmml::SUGARNAMELOOKUP[ i ].pdb_code_;
      return;
    }
  }
}

void Glycan::Monosaccharide::UpdateComplexSugarChemicalCode()
{
  for( std::vector<std::pair< std::string, std::string> >::iterator it1 = derivatives_map_.begin(); it1 != derivatives_map_.end(); it1++ )
  {
    std::string key = ( *it1 ).first;
    std::string value = ( *it1 ).second;
    if( value.compare( "" ) != 0 )
    {
      std::string code = "";
      if( value.compare( "xCH-N" ) == 0 )
      {
        code = key + "N";
      }
      else if( value.compare( "xC-N-C=OCH3" ) == 0 )
      {
        code = key + "NAc";
      }
      else if( value.compare( "xC-N-C=OCH2OH" ) == 0 )
      {
        code = key + "NGc";
      }
      else if( value.compare( "xC-N-SO3" ) == 0 )
      {
        code = key + "NS";
      }
      else if( value.compare( "xC-N-PO3" ) == 0 )
      {
        code = key + "NP";
      }
      else if( value.compare( "xC-N-CH3" ) == 0 )
      {
        code = key + "NMe";
      }
      else if( value.compare( "xC-O-C=OCH3" ) == 0 )
      {
        code = key + "Ac";
      }
      else if( value.compare( "xC-O-C=OCH3OH" ) == 0 )
      {
        code = key + "Gc";
      }
      else if( value.compare( "xC-O-SO3" ) == 0 )
      {
        code = key + "S";
      }
      else if( value.compare( "xC-O-PO3" ) == 0 )
      {
        code = key + "P";
      }
      else if( value.compare( "xC-O-CH3" ) == 0 )
      {
        code = key + "Me";
      }
      else if( value.compare( "xC-(O,OH)" ) == 0 )
      {
        code = key + "AH";
      }
      else if( value.compare( "xC-(O,O)" ) == 0 )
      {
        code = key + "A";
      }
      // else if( value.compare( "xCHH" ) == 0 )
      // {
      //   code = key + "xCHH";
      // }
      std::vector< std::string >::iterator index_it;
      if( ( key.compare( "a" ) == 0 ) || ( key.compare( "-1" ) == 0 ) || ( key.compare( "+1" ) == 0 ) || ( key.compare( "+2" )== 0 ) || ( key.compare( "+3" ) == 0 ) )
      {
        if( ( index_it = find( chemical_code_->right_up_.begin(), chemical_code_->right_up_.end(), key ) ) != chemical_code_->right_up_.end() )
        {
          ( *index_it ) = code;
        }
        else if( ( index_it = find( chemical_code_->right_down_.begin(), chemical_code_->right_down_.end(), key ) ) != chemical_code_->right_down_.end() )
        {
          ( *index_it ) = code;
        }
      }
      else if( ( key.compare( "2" ) == 0 ) || ( key.compare( "3" ) == 0 ) || ( key.compare( "4" ) == 0 ) )
      {
        if( ( index_it = find( chemical_code_->left_up_.begin(), chemical_code_->left_up_.end(), key ) ) != chemical_code_->left_up_.end() )
        {
          ( *index_it ) = code;
        }
        else if( ( index_it = find( chemical_code_->left_down_.begin(), chemical_code_->left_down_.end(), key ) ) != chemical_code_->left_down_.end() )
        {
          ( *index_it ) = code;
        }
      }
    }
  }
}

void Glycan::Monosaccharide::CheckMonoNaming(std::string original_residue, std::string original_residue_id)
{
  if( sugar_name_.monosaccharide_stereochemistry_name_.compare( "" ) == 0 && sugar_name_.monosaccharide_name_.compare( "" ) == 0 )
    {
      // //Rob wanted closest match removed; this is a quick and dirty fix, TODO we still would like to report close matches.
      // //
      // //FINDING CLOSEST MATCH FOR THE CHEMICAL CODE IN THE LOOKUP TABLE
      // std::vector< Glycan::SugarName > closest_matches = std::vector< Glycan::SugarName >();
      // gmml::ClosestMatchSugarStereoChemistryNameLookup( chemical_code_->toString(), closest_matches );
      // if( sugar_name_.monosaccharide_name_.compare( "" ) == 0)
      // {
      //   sugar_name_.monosaccharide_name_ = sugar_name_.monosaccharide_stereochemistry_name_;
      // }
      // if( sugar_name_.monosaccharide_short_name_.compare( "" ) == 0 )
      // {
      //   sugar_name_.monosaccharide_short_name_ = sugar_name_.monosaccharide_stereochemistry_short_name_;
      // }
      //
      // ///ADDING NOTES/ISSUES OF MONOSACCHARIDE STRUCTURE
      // Glycan::Note* matching_note = new Glycan::Note();
      // matching_note->type_ = Glycan::COMMENT;
      // matching_note->category_ = Glycan::MONOSACCHARIDE;
      //
      // std::stringstream ss;
      // ss << "No match found for "; //<<  THIS NEEDS TO BE RESIDUE NAME sugar_name_.monosaccharide_stereochemistry_short_name_
      // ss << original_residue_id;
      // ss << ". close matches: ";
      // for( std::vector< Glycan::SugarName >::iterator ite = closest_matches.begin(); ite != closest_matches.end(); ite++ )
      // {
      //   Glycan::SugarName sn = ( *ite );
      //   if( ite == closest_matches.end() - 1 ) {
      //     ss << sn.monosaccharide_stereochemistry_short_name_;
      //   }
      //   else {
      //     ss << sn.monosaccharide_stereochemistry_short_name_ << ", ";
      //   }
      // }
      // matching_note->description_ = ss.str();
      // mono_notes_.push_back(matching_note);
      // GetAuthorNaming(amino_lib_files, CCD_Path);
      if( sugar_name_.monosaccharide_stereochemistry_name_.compare( "" ) == 0 )
      {
        //TODO call a formula generating function
        sugar_name_.monosaccharide_stereochemistry_name_ = cycle_atoms_[0]->GetResidue()->GetName();
        sugar_name_.monosaccharide_stereochemistry_short_name_ = cycle_atoms_[0]->GetResidue()->GetName();
        sugar_name_.monosaccharide_name_ = cycle_atoms_[0]->GetResidue()->GetName();
        sugar_name_.monosaccharide_short_name_ = cycle_atoms_[0]->GetResidue()->GetName();
      }

      // else
      // {
      //   std::cout << "No exact match found for the chemical code, the following information comes from one of the closest matches:" << std::endl;
      // }
    }
}

// void Glycan::Monosaccharide::addDeoxyToName(Glycan::SugarName base_name, Glycan::ChemicalCode chemical_code, std::vector<int> deoxy_locations)
// {
//   //Figure out how to modify all parts of the sugar name (base_name) correctly for the deoxy locations specified in deoxy_locations
//
//   //Following the logic in GenerateCompleteSugarName
//
//   //chemical_code_string_ won't change for now at least
//
//   //monosaccharide_stereochemistry_name_
//
//   if(base_name.monosaccharide_stereochemistry_name_ != "")
//   {
//
//   }
//
//   //monosaccharide_stereochemistry_short_name_
//
//   //isomer_ does not change
//
//   //name_
//
//   //ring_type_ doesn not change
//
//   //configuration_ does not change
//
//   //monosaccharide_name_
//
//   //monosaccharide_short_name_
//
//   //pdb_code_
// }

Glycan::ChemicalCode* Glycan::Monosaccharide::BuildChemicalCode(std::vector<std::string> orientations)
{
  Glycan::ChemicalCode* code = new Glycan::ChemicalCode();
  if(orientations.size() == 5 )
    code->base_ = "P";
  else if(orientations.size() == 4 )
    code->base_ = "F";
  else
    code->base_ = "?";

  ///Side atom(s) of anomeric
  ///Has only non-ring oxygen neighbor
  if(orientations.at(0).compare("U") == 0)
    code->right_up_.push_back("a");
  else if(orientations.at(0).compare("D") == 0)
    code->right_down_.push_back("a");

  ///Has non-ring oxygen and carbon neighbors
  else if(orientations.at(0).compare("U-1U") == 0)
  {
    code->right_up_.push_back("a");
    code->right_up_.push_back("-1");
  }
  else if(orientations.at(0).compare("D-1U") == 0)
  {
    code->right_down_.push_back("a");
    code->right_up_.push_back("-1");
  }
  else if(orientations.at(0).compare("U-1D") == 0)
  {
    code->right_up_.push_back("a");
    code->right_down_.push_back("-1");
  }
  else if(orientations.at(0).compare("D-1D") == 0)
  {
    code->right_down_.push_back("a");
    code->right_down_.push_back("-1");
  }

  ///Has only non-ring carbon neighbor
  else if(orientations.at(0).compare("-1U") == 0)
    code->right_up_.push_back("-1");
  else if(orientations.at(0).compare("-1D") == 0)
    code->right_down_.push_back("-1");

  ///Side atom of other carbons of the ring
  for(std::vector<std::string>::iterator it = orientations.begin() + 1; it != orientations.end() - 1; it++)
  {
    std::string orientation = (*it);
    int index = distance(orientations.begin(), it);
    if(orientation.compare("U") == 0)
      code->left_up_.push_back(gmml::ConvertT(index + 1));
    else if(orientation.compare("D") == 0)
      code->left_down_.push_back(gmml::ConvertT(index + 1));
    else if(orientation.compare("N") == 0)
    {
      std::stringstream ss;
      ss << gmml::ConvertT(index + 1) << "d";
      code->left_middle_.push_back(ss.str() );
    }
  }

  ///Side atom(s) of last carbon
  if(orientations.at(orientations.size() - 1).compare("U") == 0)
    code->right_up_.push_back("+1");
  else if(orientations.at(orientations.size() - 1).compare("D") == 0)
    code->right_down_.push_back("+1");
  ///Type Deoxy
  else if(orientations.at(orientations.size() - 1).compare("Ud") == 0)
    code->right_up_.push_back("+1d");
  else if(orientations.at(orientations.size() - 1).compare("Dd") == 0)
    code->right_down_.push_back("+1d");

  return code;
}
