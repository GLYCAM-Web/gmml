//Created 12/6/18
//Dave Montgomery

#include <vector>
#include <string>
#include <sstream>

#include "../../includes/MolecularModeling/assembly.hpp"
#include "../../includes/utils.hpp"
#include "../../includes/MolecularModeling/residue.hpp"
#include "../../includes/Glycan/oligosaccharide.hpp"

using Glycan::Oligosaccharide;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
Oligosaccharide::Oligosaccharide(std::vector<Glycan::Monosaccharide*> monos, gmml::ResidueNameMap dataset_residue_names, MolecularModeling::Assembly* assembly)
{
  int number_of_covalent_links = 0;
  int number_of_probable_non_covalent_complexes = 0;
  assembly_ = assembly;
  createOligosaccharideGraphs(monos, dataset_residue_names, number_of_covalent_links, number_of_probable_non_covalent_complexes);
  // createOligosaccharides(monos);
}
Oligosaccharide::Oligosaccharide(MolecularModeling::Assembly* assembly)
{
  assembly_ = assembly;
}
Oligosaccharide::Oligosaccharide()
{

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
/*! \fn
* A function to print out the oligosacchride name and the linkages between its monosacchrides
* Print out the information in a defined structure
* @param out An output stream, the print result will be written in the given output stream
*/
void Glycan::Oligosaccharide::Print(std::ostream& out)
{
  oligosaccharide_name_ = "";
  oligosaccharide_linkages_ = "";
  oligosaccharide_residue_linkages_ = "";
  int linkage_index = 0;
  std::stringstream ss;
  if(terminal_.compare("") != 0)
  {
    if(std::find_if( root_->derivatives_map_.begin(), root_->derivatives_map_.end(),
      [](const std::pair<std::string, std::string>& element){ return element.first == "-1";}) == root_->derivatives_map_.end() )
    // if(root_->derivatives_map_.find("-1") == root_->derivatives_map_.end())
    {
      ss << "1-" << terminal_;
    }
    else
    {
      ss << "2-" << terminal_;
    }
  }
  oligosaccharide_name_ = ss.str();
  oligosaccharide_terminal_ = oligosaccharide_name_;
  bool is_cycle = false;
  this->GenerateNameLinkage(oligosaccharide_name_, oligosaccharide_linkages_, linkage_index, root_->mono_id_, is_cycle);

  if(is_cycle)
  {
    int end_index = oligosaccharide_name_.find_last_of('-');
    std::vector<std::string> tokens = gmml::Split(oligosaccharide_name_, "-");
    std::string sub_name = oligosaccharide_name_.substr(0, end_index);
    std::stringstream new_name;
    // gmml::log(__LINE__, __FILE__,  gmml::INF, " This.size 1" );
    new_name << "[" << tokens.at(tokens.size() - 1).at(0) << sub_name << "-]";
    oligosaccharide_name_ = new_name.str();
  }
  // gmml::log(__LINE__, __FILE__,  gmml::INF, oligosaccharide_name_ );
  // Example oligo sequence LRhapa1-3LRhapa1-3DGlcpNAcb1-OME
  // Example oligo linkages
  //{2}RAM(401_A)C1-RAE(402_A)C3, Glycosidic linkage: RAE(402_A)O3
  //{1}RAE(402_A)C1-MAG(403_A)C3, Glycosidic linkage: MAG(403_A)O3
  //For each line of linkages the info of the first residue (e.g. RAE(402_A)) will be added to the oligosaccharide_residue_linkages_.
  std::vector<std::string> oligo_linkages_tokens = gmml::Split(oligosaccharide_linkages_, "\n");
  std::vector<std::string> oligo_name_tokens = gmml::Split(oligosaccharide_name_, "-"); ///According to position of brackets in the name sequence,
                                                                            ///the brackets will be added to the oligosaccharide_residue_linkages_
  std::stringstream residue_links_stream;
  std::string full_glycosidic_linkage = "";
  std::vector<std::string> link_tokens = std::vector<std::string>();
  std::string link_left_side = "";
  std::string link_right_side = "";
  size_t end_index = 0;
  std::string first_residue_of_linkage  = "";
  std::string second_residue_of_linkage  = "";
  //TODO fix the double condition for loop
  // gmml::log(__LINE__, __FILE__,  gmml::INF, " This.size double for loop ..." );
  // gmml::log(__LINE__, __FILE__,  gmml::INF, oligosaccharide_linkages_ );
  for(unsigned int i = 0; i < oligo_linkages_tokens.size() && i < oligo_name_tokens.size(); i++) ///Processing linkages line by line
  {
    full_glycosidic_linkage = oligo_linkages_tokens.at(i);
    link_tokens = gmml::Split(full_glycosidic_linkage, "}-,");
    link_left_side = link_tokens.at(1);     ///Getting the first residue of linkage in each line. e.g. {2}RAM(401_A)C1-RAE(402_A)C3, Glycosidic linkage: RAE(402_A)O3
                                            ///-> first_mono == RAM(401_A)C1
    end_index = link_left_side.find_last_of(")");
    first_residue_of_linkage = link_left_side.substr(0, end_index + 1); ///filtering out atom name. e.g. RAM(401_A)
    std::string oligo_temp = oligo_name_tokens.at(i);
    residue_links_stream << updateResidueLink(oligo_temp, first_residue_of_linkage) << "-";
    if ((i >= oligo_linkages_tokens.size() - 1) && (link_tokens.size() > 2))
    {
      link_right_side = link_tokens.at(2);     ///Getting the second residue of linkage in the line. e.g. {1}NAG(1521_A)C1-NAG(1520_A)C4, Glycosidic linkage: NAG(1520_A)O4
                                                  ///-> second_mono == NAG(1520_A)C4
      end_index = link_right_side.find_last_of(")");
      second_residue_of_linkage = link_right_side.substr(0, end_index + 1); ///filtering out atom name. e.g. NAG(1520_A)
      std::string oligo_temp = oligo_name_tokens.at(i+1);
      residue_links_stream << updateResidueLink(oligo_temp, second_residue_of_linkage);
      // gmml::log(__LINE__, __FILE__,  gmml::INF, residue_links_stream.str() );
    }
  }

  std::string root_residue_name = "";
  std::string root_residue_number = "";
  std::vector<MolecularModeling::Atom*> mono_ring_atoms = root_->cycle_atoms_;
  if(mono_ring_atoms.at(0) != NULL)
  {
    std::string atom_id = mono_ring_atoms.at(0)->GetId();
    std::vector<std::string> atom_id_tokens = gmml::Split(atom_id, "_");
    root_residue_name = atom_id_tokens.at(2);
    root_residue_number = atom_id_tokens.at(4);
    // gmml::log(__LINE__, __FILE__,  gmml::INF, " This.size another time ..." );
    if(oligo_linkages_tokens.size() == 0 && oligosaccharide_name_.compare("") != 0)
    {
      if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
      {
        residue_links_stream << root_residue_name << "(" << root_residue_number  << "_" << /*root_->bfmp_ring_conformation_ <<*/ ")" ;
      }
      else
      {
        residue_links_stream << root_residue_name << "(" << root_residue_number <<  "_" << atom_id_tokens.at(3) << "_" << /*root_->bfmp_ring_conformation_ <<*/ ")" ;
      }
    }
  }

  if(terminal_.compare("") != 0)
  {
    std::vector<MolecularModeling::Atom*> root_anomeric_atom = root_->side_atoms_.at(0);
    // @TODO Where should we check for the terminal atom?
    if(root_anomeric_atom.at( 0 ) != NULL || root_anomeric_atom.at( 1 ) != NULL)
    {
      std::string atom_id = "";
      if( root_anomeric_atom.at( 0 ) != NULL )
      {
        atom_id = root_anomeric_atom.at( 0 )->GetId();
      }
      else if( root_anomeric_atom.at( 1 ) != NULL )
      {
        atom_id = root_anomeric_atom.at( 1 )->GetId();
      }
      std::vector<std::string> atom_id_tokens = gmml::Split(atom_id, "_");
      std::stringstream vec_size;
      vec_size << atom_id_tokens.size();
      std::string terminal_residue_name = atom_id_tokens.at(2);
      std::string terminal_residue_number = atom_id_tokens.at(4);
      if (root_residue_name.compare(terminal_residue_name) != 0 && root_residue_number.compare(terminal_residue_number) != 0)
      {
        //Additional underscore(_) symbol is added at the end as terminal doesn't have BFMP
        //It makes it easier to remove/add BFMP conformation from the whole oligo-sequence
        //e.g. FUL(404_D_1C4)-NAG(401_D_2d5)-ASN(80_D_)
        //As we have underscore symbol at the end of terminal, we can remove everything after last _ in each residue of oligo
        //oligo-seq without BFMP would look like: FUL(404_D)-NAG(401_D)-ASN(80_D)
        if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
        {
          residue_links_stream << "-" << terminal_residue_name << "(" << terminal_residue_number << "_" << ")" ;
        }
        else
        {
          residue_links_stream << "-" << terminal_residue_name << "(" << terminal_residue_number <<  "_" << atom_id_tokens.at(3) << "_" << ")" ;
        }
      }
    }
  }
  oligosaccharide_residue_linkages_ = residue_links_stream.str();
  gmml::FindReplaceString(oligosaccharide_residue_linkages_, "-]", "]");


  //I have commented out the lines that print to out, as we decided nothing should print that isn't an error.  This function does other things so this was the quickest workaround. DM 1-15-20
  // out << oligosaccharide_name_;
  // gmml::log(__LINE__, __FILE__,  gmml::INF, oligosaccharide_name_);
  // out << std::endl;
  // gmml::log(__LINE__, __FILE__,  gmml::INF, oligosaccharide_linkages_);
}


/*! \fn
* A function to update residueLinkStream based on current oligosaccharide
* @param oligo_temp The current oligosaccharide of oligo-sequence
* @param residue_of_linkage The PDB-ID of current residue
*/
std::string Glycan::Oligosaccharide::updateResidueLink(std::string oligo_temp, std::string residue_of_linkage)
{
  std::stringstream residue_link;
  for (unsigned int number=0; number < oligo_temp.size(); number++)
  {
    if (isdigit (oligo_temp[number]))
    {
      std::string s(1, oligo_temp[number]);
      gmml::FindReplaceString(oligo_temp, s, "");
    }
  }
  if(oligo_temp.find("[") != std::string::npos && oligo_temp.find("]") != std::string::npos)
  {
    if (oligo_temp.find("[") < oligo_temp.find("]"))
    {
      residue_link << residue_of_linkage;
    }
    else
    {
      if (oligo_temp.find("[") == 1 && oligo_temp.find("]") == 0)
      {
        residue_link << "]-[" << residue_of_linkage;
      }
    }
  }
  else if(oligo_temp.find("]") != std::string::npos)
  {
    residue_link << "]-" << residue_of_linkage;
  }
  else if(oligo_temp.find("[") != std::string::npos)
  {
    residue_link << "[" << residue_of_linkage;
  }
  else
  {
    residue_link << residue_of_linkage;
  }
  return residue_link.str();
}

/*! \fn
* A recursive function in order to generate the name sequence and the linkages between sacchrides of an oligosacchride
* This functions starts from the root sacchride and generates the name sequence and linkages in reverse order
* @param oligosaccharide_name The name sequence of the oligosacchride that has been generated by the function so far
* @param oligosaccharide_linkages The linkages between sacchrides of the oligosacchride which has been generated by the function so far
* @param i The counter of linkages between the sacchrides of the oligosacchride
* @param main_root_id The identifier of the main root in the oligosacchride
* @param is_cycle The boolean variable which represnets whether the sacchrides involved in the oligosacchride forms a cycle
*/
void Glycan::Oligosaccharide::GenerateNameLinkage(std::string& oligosaccharide_name, std::string& oligosaccharide_linkages, int& i, int main_root_id, bool& is_cycle)
{
  std::stringstream name;
  std::stringstream res_linkage;
  std::stringstream res_linkage1;
  if(child_oligos_.size() == 0)
  {
    if(root_->sugar_name_.monosaccharide_short_name_.compare("") != 0)
    {
      name << root_->sugar_name_.monosaccharide_short_name_ << oligosaccharide_name;
    }
    else
    {
      name << root_->sugar_name_.monosaccharide_stereochemistry_short_name_ << oligosaccharide_name;
    }
    oligosaccharide_name = name.str();
  }
  else if(child_oligos_.size() == 1)
  {
    if(root_->sugar_name_.monosaccharide_short_name_.compare("") != 0)
    {
      name << root_->sugar_name_.monosaccharide_short_name_ << oligosaccharide_name;
    }
    else
    {
      name << root_->sugar_name_.monosaccharide_stereochemistry_short_name_ << oligosaccharide_name;
    }
    oligosaccharide_name = name.str();

    Monosaccharide* mono1 = root_;
    Oligosaccharide* child_oligo = child_oligos_.at(0);
    Monosaccharide* mono2 = child_oligo->root_;
    // std::string mono1_bfmp = mono1->bfmp_ring_conformation_;
    // std::string mono2_bfmp = mono2->bfmp_ring_conformation_;
    std::string mono1_bfmp = "";
    std::string mono2_bfmp = "";

    std::vector<std::string> tokens = gmml::Split(child_oligos_linkages_.at(0), "-");

    std::stringstream link1;
    std::string from = tokens.at(0);

    std::string mono1_ring_atoms_str = mono1->cycle_atoms_str_;
    std::vector<MolecularModeling::Atom*> mono1_ring_atoms = mono1->cycle_atoms_;
    std::vector<std::vector<MolecularModeling::Atom*> > mono1_side_atoms = mono1->side_atoms_;
    int mono1_start_index = 1;
    std::vector<MolecularModeling::Atom*> mono1_anomeric_side_atoms = mono1_side_atoms.at(0);
    if(mono1_anomeric_side_atoms.at(0) != NULL)
    {
      mono1_start_index = 2;
    }
    if(mono1_ring_atoms_str.find(from) != std::string::npos)
    {
      for(std::vector<MolecularModeling::Atom*>::iterator it = mono1_ring_atoms.begin(); it != mono1_ring_atoms.end(); it++)
      {
        MolecularModeling::Atom* atom = *it;
        if(atom != NULL && atom->GetId().compare(from) == 0)
        {
          int from_index = std::distance(mono1_ring_atoms.begin(), it) + mono1_start_index;
          link1 << "-" << gmml::ConvertT(from_index);
          std::vector<std::string> atom_id_tokens = gmml::Split(from, "_");
          std::vector<std::string> glycosidic_linkage_id_tokens = gmml::Split(tokens.at(1), "_");
          if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
          {res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) <<  "_" << mono1_bfmp << ")"  << atom_id_tokens.at(0);}
          else
          {
            res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) <<  "_" << atom_id_tokens.at(3) <<  "_" << mono1_bfmp << ")"  << atom_id_tokens.at(0);
          }

          if(glycosidic_linkage_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
          {
            res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << ")" << glycosidic_linkage_id_tokens.at(0) << std::endl;
          }
          else
          {
            res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << "_" << glycosidic_linkage_id_tokens.at(3) << ")"  << glycosidic_linkage_id_tokens.at(0) << std::endl;
          }
          break;
        }
      }
    }
    else
    {
      std::vector<MolecularModeling::Atom*> last_carbon_side_atoms = mono1_side_atoms.at(mono1_side_atoms.size() - 1);
      for(std::vector<MolecularModeling::Atom*>::iterator it = last_carbon_side_atoms.begin(); it != last_carbon_side_atoms.end(); it++)
      {
        MolecularModeling::Atom* atom = *it;
        if(atom != NULL && atom->GetId().compare(from) == 0)
        {
          int from_index = std::distance(last_carbon_side_atoms.begin(), it) + mono1_side_atoms.size() + mono1_start_index;
          link1 << "-" << gmml::ConvertT(from_index);
          std::vector<std::string> atom_id_tokens = gmml::Split(from, "_");
          std::vector<std::string> glycosidic_linkage_id_tokens = gmml::Split(tokens.at(1), "_");
          if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
          {
            res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) <<  "_" << mono1_bfmp << ")" << atom_id_tokens.at(0);
          }
          else
          {
            res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << atom_id_tokens.at(3) <<  "_" << mono1_bfmp << ")" << atom_id_tokens.at(0);
          }

          if(glycosidic_linkage_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
          {
            res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << ")" << glycosidic_linkage_id_tokens.at(0) << std::endl;
          }
          else
          {
            res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << "_" << glycosidic_linkage_id_tokens.at(3) << ")"  << glycosidic_linkage_id_tokens.at(0) << std::endl;
          }
          break;
        }
      }
    }

    std::stringstream link;
    std::string to = tokens.at(2);

    std::string mono2_ring_atoms_str = mono2->cycle_atoms_str_;
    std::vector<MolecularModeling::Atom*> mono2_ring_atoms = mono2->cycle_atoms_;
    std::vector<std::vector<MolecularModeling::Atom*> > mono2_side_atoms = mono2->side_atoms_;
    std::vector<MolecularModeling::Atom*> mono2_anomeric_side_atoms = mono2_side_atoms.at(0);

    int mono2_start_index = 1;
    if(mono2_anomeric_side_atoms.at(0) != NULL)
    {
      mono2_start_index = 2;
    }
    if(mono2_ring_atoms_str.find(to) != std::string::npos)
    {
      for(std::vector<MolecularModeling::Atom*>::iterator it = mono2_ring_atoms.begin(); it != mono2_ring_atoms.end(); it++)
      {
        MolecularModeling::Atom* atom = *it;
        if(atom != NULL && atom->GetId().compare(to) == 0)
        {
          int to_index = std::distance(mono2_ring_atoms.begin(), it) + mono2_start_index;
          link << gmml::ConvertT(to_index);
          std::vector<std::string> atom_id_tokens = gmml::Split(to, "_");
          if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
          {
            res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) <<  "_" << mono2_bfmp << ")" << atom_id_tokens.at(0);
          }
          else
          {
            res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << atom_id_tokens.at(3) <<  "_" << mono2_bfmp << ")" << atom_id_tokens.at(0);
          }
          break;
        }
      }
    }
    else
    {
      std::vector<MolecularModeling::Atom*> last_carbon_side_atoms = mono2_side_atoms.at(mono2_side_atoms.size() - 1);
      for(std::vector<MolecularModeling::Atom*>::iterator it = last_carbon_side_atoms.begin(); it != last_carbon_side_atoms.end(); it++)
      {
        MolecularModeling::Atom* atom = *it;
        if(atom != NULL && atom->GetId().compare(to) == 0)
        {
          int to_index = std::distance(last_carbon_side_atoms.begin(), it) + mono2_side_atoms.size() + mono2_start_index;
          link << gmml::ConvertT(to_index);
          std::vector<std::string> atom_id_tokens = gmml::Split(to, "_");
          if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
          {
            res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) <<  "_" << mono2_bfmp << ")" << atom_id_tokens.at(0);
          }
          else
          {
            res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << atom_id_tokens.at(3) <<  "_" << mono2_bfmp << ")" << atom_id_tokens.at(0);
          }
          break;
        }
      }
    }

    link << link1.str() << oligosaccharide_name;
    oligosaccharide_name = link.str();
    res_linkage << res_linkage1.str() << oligosaccharide_linkages;
    oligosaccharide_linkages = res_linkage.str();

    i++;
    if(main_root_id == child_oligo->root_->mono_id_)
    {
      is_cycle = true;
    }
    child_oligo->GenerateNameLinkage(oligosaccharide_name, oligosaccharide_linkages, i, main_root_id, is_cycle);
  }
  else
  {
    std::map<int, int> root_child_index_map = std::map<int, int>();

    int mono1_start_index = 1;
    if(root_->side_atoms_.at(0).at(0) != NULL)///mono has a carbon at position -1 so the ring atom indices should start from 2
    {
      mono1_start_index = 2;
    }
    for(std::vector<std::string>::iterator it = child_oligos_linkages_.begin(); it != child_oligos_linkages_.end(); it++)
    {
      std::string child_linkage = (*it);
      int linkage_index = std::distance(child_oligos_linkages_.begin(), it);
      std::vector<std::string> child_linkage_tokens = gmml::Split(child_linkage, "-");
      std::string from = child_linkage_tokens.at(0);///the atom id from current mono that is involved in the linkage

      if(root_->cycle_atoms_str_.find(from) != std::string::npos)///if the atom that is involved in the linkage is one of the current mono cycle atoms
      {
        for(std::vector<MolecularModeling::Atom*>::iterator it = root_->cycle_atoms_.begin(); it != root_->cycle_atoms_.end(); it++)
        {
          MolecularModeling::Atom* atom = *it;
          if(atom != NULL && atom->GetId().compare(from) == 0)
          {
              int from_index = std::distance(root_->cycle_atoms_.begin(), it) + mono1_start_index;///index position of the atom in the mono
              root_child_index_map[linkage_index] = from_index;
              break;
          }
        }
      }
      else///the atom that is involved in the linkage is one of the current mono side atoms
      {
        std::vector<MolecularModeling::Atom*> last_carbon_side_atoms = root_->side_atoms_.at(root_->side_atoms_.size() - 1);
        for(std::vector<MolecularModeling::Atom*>::iterator it = last_carbon_side_atoms.begin(); it != last_carbon_side_atoms.end(); it++)
        {
          MolecularModeling::Atom* atom = *it;
          if(atom != NULL && atom->GetId().compare(from) == 0)
          {
            int from_index = std::distance(last_carbon_side_atoms.begin(), it) + root_->side_atoms_.size() + mono1_start_index;///index position of the side atom in the mono
            root_child_index_map[linkage_index] = from_index;
            break;
          }
        }
      }
    }

    bool isFirstChild = true;
    while(root_child_index_map.size() > 0)
    {
      res_linkage.str("");
      res_linkage1.str("");
      int min = 100;
      int index_min = 0;///index of the child oligo with the min children in child_oligos
      for(std::map<int, int>::iterator it0 = root_child_index_map.begin(); it0 != root_child_index_map.end(); it0++)
      {
        int map_key = (*it0).first;
        int map_value = (*it0).second;
        if(map_value < min)
        {
          index_min = map_key;
          min = map_value;
        }
      }
      if(isFirstChild)
      {
        isFirstChild = false;
        if(root_->sugar_name_.monosaccharide_short_name_.compare("") != 0)
        {
          name << root_->sugar_name_.monosaccharide_short_name_ << oligosaccharide_name;
        }
        else
        {
          name << root_->sugar_name_.monosaccharide_stereochemistry_short_name_ << oligosaccharide_name;
        }
        oligosaccharide_name = name.str();
      }

      Monosaccharide* mono1 = root_;
      Oligosaccharide* child_oligo = child_oligos_.at(index_min);
      Monosaccharide* mono2 = child_oligo->root_;

      // std::string mono1_bfmp = mono1->bfmp_ring_conformation_;
      // std::string mono2_bfmp = mono2->bfmp_ring_conformation_;
      std::string mono1_bfmp = "";
      std::string mono2_bfmp = "";

      std::vector<std::string> tokens = gmml::Split(child_oligos_linkages_.at(index_min), "-");

      std::stringstream link1;
      std::string from = tokens.at(0);

      std::string mono1_ring_atoms_str = mono1->cycle_atoms_str_;
      std::vector<MolecularModeling::Atom*> mono1_ring_atoms = mono1->cycle_atoms_;
      std::vector<std::vector<MolecularModeling::Atom*> > mono1_side_atoms = mono1->side_atoms_;
      int mono1_start_index = 1;
      std::vector<MolecularModeling::Atom*> mono1_anomeric_side_atoms = mono1_side_atoms.at(0);
      if(mono1_anomeric_side_atoms.at(0) != NULL)
      {
        mono1_start_index = 2;
      }
      if(mono1_ring_atoms_str.find(from) != std::string::npos)
      {
        for(std::vector<MolecularModeling::Atom*>::iterator it = mono1_ring_atoms.begin(); it != mono1_ring_atoms.end(); it++)
        {
          MolecularModeling::Atom* atom = *it;
          if(atom != NULL && atom->GetId().compare(from) == 0)
          {
            int from_index = std::distance(mono1_ring_atoms.begin(), it) + mono1_start_index;
            link1 << "-" << gmml::ConvertT(from_index);
            std::vector<std::string> atom_id_tokens = gmml::Split(from, "_");
            std::vector<std::string> glycosidic_linkage_id_tokens = gmml::Split(tokens.at(1), "_");
            if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
            {
              res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << mono1_bfmp << ")"  << atom_id_tokens.at(0);
            }
            else
            {
              res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << atom_id_tokens.at(3) << "_" << mono1_bfmp << ")" << atom_id_tokens.at(0);
            }

            if(glycosidic_linkage_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
            {
              res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << ")" << glycosidic_linkage_id_tokens.at(0) << std::endl;
            }
            else
            {
              res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << "_" << glycosidic_linkage_id_tokens.at(3) << ")"  << glycosidic_linkage_id_tokens.at(0) << std::endl;
            }
            break;
          }
        }
      }
      else
      {
        std::vector<MolecularModeling::Atom*> last_carbon_side_atoms = mono1_side_atoms.at(mono1_side_atoms.size() - 1);
        for(std::vector<MolecularModeling::Atom*>::iterator it = last_carbon_side_atoms.begin(); it != last_carbon_side_atoms.end(); it++)
        {
          MolecularModeling::Atom* atom = *it;
          if(atom != NULL && atom->GetId().compare(from) == 0)
          {
            int from_index = std::distance(last_carbon_side_atoms.begin(), it) + mono1_side_atoms.size() + mono1_start_index;
            link1 << "-" << gmml::ConvertT(from_index);
            std::vector<std::string> atom_id_tokens = gmml::Split(from, "_");
            std::vector<std::string> glycosidic_linkage_id_tokens = gmml::Split(tokens.at(1), "_");
            if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
            {
              res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << mono1_bfmp << ")" << atom_id_tokens.at(0);
            }
            else
            {
              res_linkage1 << "-" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << atom_id_tokens.at(3) << "_" << mono1_bfmp << ")" << atom_id_tokens.at(0);
            }

            if(glycosidic_linkage_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
            {
              res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << ")" << glycosidic_linkage_id_tokens.at(0) << std::endl;
            }
            else
            {
              res_linkage1 <<  ", Glycosidic linkage: " << glycosidic_linkage_id_tokens.at(2) << "(" << glycosidic_linkage_id_tokens.at(4) << "_" << glycosidic_linkage_id_tokens.at(3) << ")"  << glycosidic_linkage_id_tokens.at(0) << std::endl;
            }
            break;
          }
        }
      }

      std::stringstream link;
      std::string to = tokens.at(2);

      std::string mono2_ring_atoms_str = mono2->cycle_atoms_str_;
      std::vector<MolecularModeling::Atom*> mono2_ring_atoms = mono2->cycle_atoms_;
      std::vector<std::vector<MolecularModeling::Atom*> > mono2_side_atoms = mono2->side_atoms_;
      std::vector<MolecularModeling::Atom*> mono2_anomeric_side_atoms = mono2_side_atoms.at(0);

      int mono2_start_index = 1;
      if(mono2_anomeric_side_atoms.at(0) != NULL)
      {
        mono2_start_index = 2;
      }
      if(mono2_ring_atoms_str.find(to) != std::string::npos)
      {
        for(std::vector<MolecularModeling::Atom*>::iterator it = mono2_ring_atoms.begin(); it != mono2_ring_atoms.end(); it++)
        {
          MolecularModeling::Atom* atom = *it;
          if(atom != NULL && atom->GetId().compare(to) == 0)
          {
            int to_index = std::distance(mono2_ring_atoms.begin(), it) + mono2_start_index;
            link << gmml::ConvertT(to_index);
            std::vector<std::string> atom_id_tokens = gmml::Split(to, "_");
            if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
            {
              res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << mono2_bfmp << ")" << atom_id_tokens.at(0);
            }
            else
            {
              res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << atom_id_tokens.at(3) << "_" << mono2_bfmp << ")" << atom_id_tokens.at(0);
            }
            break;
          }
        }
      }
      else
      {
        std::vector<MolecularModeling::Atom*> last_carbon_side_atoms = mono2_side_atoms.at(mono2_side_atoms.size() - 1);
        for(std::vector<MolecularModeling::Atom*>::iterator it = last_carbon_side_atoms.begin(); it != last_carbon_side_atoms.end(); it++)
        {
          MolecularModeling::Atom* atom = *it;
          if(atom != NULL && atom->GetId().compare(to) == 0)
          {
            int to_index = std::distance(last_carbon_side_atoms.begin(), it) + mono2_side_atoms.size() + mono2_start_index;
            link << gmml::ConvertT(to_index);
            std::vector<std::string> atom_id_tokens = gmml::Split(to, "_");
            if(atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
            {
              res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << mono2_bfmp << ")" << atom_id_tokens.at(0);
            }
            else
            {
              res_linkage << "{" << i+1 << "}" << atom_id_tokens.at(2) << "(" << atom_id_tokens.at(4) << "_" << atom_id_tokens.at(3) << "_" << mono2_bfmp << ")" << atom_id_tokens.at(0);
            }
            break;
          }
        }
      }

      if(root_child_index_map.size() > 1)
      {
        link << link1.str() << "]" << oligosaccharide_name;
      }
      else
      {
        link << link1.str() << oligosaccharide_name;
      }
      oligosaccharide_name = link.str();
      res_linkage << res_linkage1.str() << oligosaccharide_linkages;
      oligosaccharide_linkages = res_linkage.str();

      i++;
      child_oligo->GenerateNameLinkage(oligosaccharide_name, oligosaccharide_linkages, i, main_root_id, is_cycle);
      std::stringstream name;
      if(root_child_index_map.size() > 1)
      {
        name << "[" << oligosaccharide_name;
      }
      else
      {
        name << oligosaccharide_name;
      }
      oligosaccharide_name = name.str();
      root_child_index_map.erase(index_min);
    }
  }
}

// void Glycan::Oligosaccharide::createOligosaccharideGraphs(std::vector<Glycan::Monosaccharide*> detected_monos)
// {
//   gmml::log(__LINE__, __FILE__, gmml::INF, "Creating Oligo Graph");
//   for(std::vector<Glycan::Monosaccharide*>::iterator thisMono = detected_monos.begin(); thisMono != detected_monos.end(); thisMono++)
//   {
//     Glycan::Monosaccharide* thisMonosaccharide = (*thisMono);
//     if(thisMonosaccharide != NULL)
//     {
//       std::vector<MolecularModeling::Atom*> sideAtoms = (*thisMonosaccharide->side_atoms_.begin());//Checking Anomeric Side atoms
//       for(std::vector<MolecularModeling::Atom*>::iterator thisSideAtom = sideAtoms.begin(); thisSideAtom!= sideAtoms.end(); thisSideAtom++)
//       {
//         if(*thisSideAtom != NULL)
//         {
//           bool foundInOtherMono = false;
//           MolecularModeling::Atom* thisSideAtomPointer = (*thisSideAtom);
//           std::vector<MolecularModeling::Atom*> sideAtomNeighbors = thisSideAtomPointer->GetNode()->GetNodeNeighbors();
//           for(std::vector<MolecularModeling::Atom*>::iterator thisSideAtomNeighbor = sideAtomNeighbors.begin(); thisSideAtomNeighbor!= sideAtomNeighbors.end(); thisSideAtomNeighbor++)
//           {
//             MolecularModeling::Atom* thisNeighbor = *thisSideAtomNeighbor;
//             for(std::vector<Glycan::Monosaccharide*>::iterator otherMono = detected_monos.begin(); otherMono != detected_monos.end(); otherMono++)
//             {
//               if(thisMono != otherMono)///Cheking monos other than the current mono
//               {
//                 Glycan::Monosaccharide* otherMonosaccharide = (*otherMono);
//                 for(std::vector<std::vector<MolecularModeling::Atom*> >::iterator otherMonoSides = otherMonosaccharide->side_atoms_.begin(); otherMonoSides != otherMonosaccharide->side_atoms_.end(); otherMonoSides++)
//                 {
//                   for(std::vector<MolecularModeling::Atom*>::iterator otherSideAtom = (*otherMonoSides).begin(); otherSideAtom!= (*otherMonoSides).end(); otherSideAtom++)
//                   {
//                     if(thisNeighbor == (*otherSideAtom))//this neighbor is a side atom in another mono
//                     {
//                       std::string anomericCarbonID = thisMonosaccharide->cycle_atoms_.at(0)->GetId();
//                       anomericCarbonID = anomericCarbonID[1];
//                       std::string otherCarbonID;
//                       std::vector<MolecularModeling::Atom*> otherMonoNeighbors = (*otherSideAtom)->GetNode()->GetNodeNeighbors();
//                       for (std::vector<MolecularModeling::Atom*>::iterator otherMonoNeighbor = otherMonoNeighbors.begin(); otherMonoNeighbor != otherMonoNeighbors.end(); otherMonoNeighbor++)
//                       {
//                         MolecularModeling::Atom* otherNeighbor = (*otherMonoNeighbor);
//                         if (otherNeighbor->GetIsCycle())
//                         {
//                           otherCarbonID = otherNeighbor->GetId();
//                           otherCarbonID = otherCarbonID[1];
//                         }
//                       }
//                       Glycan::GlycosidicLinkage* thisLinkage = new Glycan::GlycosidicLinkage(thisMonosaccharide, otherMonosaccharide, anomericCarbonID, otherCarbonID);
//                       thisMonosaccharide->mono_neighbors_.insert(std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*>(thisLinkage, otherMonosaccharide));
//                       std::stringstream logss;
//                       logss << thisMonosaccharide->sugar_name_.monosaccharide_stereochemistry_short_name_ << " is connected to " << otherMonosaccharide->sugar_name_.monosaccharide_stereochemistry_short_name_ << " via " << anomericCarbonID << "-" << otherCarbonID;
//                       gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
//                     }
//                   }
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }
// }

void Glycan::Oligosaccharide::createOligosaccharideGraphs(std::vector<Glycan::Monosaccharide*> detected_monos,
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
              //              std::cout << ss.str() << std::endl;
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
                //                std::cout << ss.str() << std::endl;
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
        Glycan::Oligosaccharide* oligo = new Glycan::Oligosaccharide();
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
            Glycan::Oligosaccharide* oligo = new Glycan::Oligosaccharide();
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
    // return oligosaccharides;
}

std::vector<Glycan::Oligosaccharide*> Glycan::Oligosaccharide::createOligosaccharides(std::vector<Glycan::Monosaccharide*> detected_monos)
{
  std::vector<Glycan::Oligosaccharide*> detected_oligos;
  for(std::vector<Glycan::Monosaccharide*>::iterator it = detected_monos.begin(); it != detected_monos.end(); it++)
  {
    Glycan::Monosaccharide* this_mono = *it;
    if(this_mono->mono_neighbors_.empty())
    {
      this_mono->is_root_ = true;
      this_mono->is_visited_ = true;
    }
    else
    {
      this_mono->is_root_ = true;
      for(std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> >::iterator monoNeighbor = this_mono->mono_neighbors_.begin(); monoNeighbor != this_mono->mono_neighbors_.end(); monoNeighbor++)
      {
        Glycan::GlycosidicLinkage* thisLinkage = (*monoNeighbor).first;
        // std::stringstream ss;
        // ss << this_mono->cycle_atoms_[0]->GetResidue()->GetId() << " is being compared to " << thisLinkage->non_reducing_mono_->cycle_atoms_[0]->GetResidue()->GetId() << " and the linkage type is " << thisLinkage->inverse_linkage_type_;
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
    }
    if((this_mono->is_root_) && (!this_mono->is_visited_))
    {
      // gmml::log(__LINE__, __FILE__,  gmml::INF, "This mono is the root");
      Glycan::Oligosaccharide* this_Oligo = new Glycan::Oligosaccharide;
      // gmml::log(__LINE__, __FILE__, gmml::INF, this_mono->sugar_name_.monosaccharide_short_name_);
      traverseGraph(this_mono, this_Oligo);
      reindexRGroups(this_Oligo);
      detected_oligos.push_back(this_Oligo);
      // std::string iupac = "Oligo IUPAC Name: " + this_Oligo->IUPAC_name_;
      // gmml::log(__LINE__, __FILE__, gmml::INF, iupac);
      // std::string oligoname =  "Oligo Name: " +this_Oligo->oligosaccharide_name_;
      // gmml::log(__LINE__, __FILE__, gmml::INF, oligoname);
    }
  }
  return detected_oligos;
}

void Glycan::Oligosaccharide::indexMonosaccharides()
{
  std::stringstream testLog;
  int local_debug = -1;

  //Generate Name Index (IUPAC naming order)
  for(std::vector<Glycan::Monosaccharide*>::iterator mono = this->mono_nodes_.begin(); mono != this->mono_nodes_.end(); mono++)
  {
    this->number_branches_ = 0;
    Glycan::Monosaccharide* thisMono = *mono;
    thisMono->is_counted_ = false;
    thisMono->is_visited_ = false;
    thisMono->is_indexed_ = false;
  }
  this->on_index_ = 0;
  std::vector<Glycan::Monosaccharide*> branchedMonos, testBranchedMonos;
  indexMono(this->root_, 0, branchedMonos);
  if(branchedMonos.size() > 0)
  {
    if(local_debug > 0)
    {
      testLog << "There are " << branchedMonos.size() << " branches";
      gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
      testLog.str(std::string());//clear stringstream
    }
    indexBranches(branchedMonos);
  }
  if(local_debug > 0)
  {
    testLog << "Done with Indexing";
    gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
    testLog.str(std::string());//clear stringstream
    testLog << "This Oligo has " << this->number_branches_ << " branches";
    gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
    testLog.str(std::string());//clear stringstream
  }

  //Generate IUPAC Index (Longest chain, then closest branch to aglycon following IUPAC rules for equal branches & branch ordering)
  for(std::vector<Glycan::Monosaccharide*>::iterator mono = this->mono_nodes_.begin(); mono != this->mono_nodes_.end(); mono++)
  {
    this->number_branches_ = 0;
    Glycan::Monosaccharide* thisMono = *mono;
    thisMono->is_counted_ = false;
    thisMono->is_visited_ = false;
    thisMono->is_indexed_ = false;
  }
  this->on_IUPAC_index_ = 0;
  branchedMonos.clear();
  branchedMonos = IUPACindexMono(this->root_, 0, branchedMonos);
  if(branchedMonos.size() > 0)
  {
    IUPACindexBranches(branchedMonos);
  }
  if(local_debug > 0)
  {
    testLog << "Done with IUPAC Indexing";
    gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
    testLog.str(std::string());//clear stringstream
  }
}

void Glycan::Oligosaccharide::indexMono(Glycan::Monosaccharide* thisMono, int thisIndex, std::vector<Glycan::Monosaccharide*> &branchedMonos)
{//Index in the order they appear in the OligoIUPACname so Query can check indexes for subgraph matching correctly
  int local_debug = -1;
  std::stringstream testLog;
  if(!thisMono->is_indexed_)
  {
    if(local_debug > 0)
    {
      std::stringstream testLog;
      testLog << "Mono: " << thisMono->cycle_atoms_[0]->GetResidue()->GetId();
      testLog << ", Index: " << thisIndex;
      gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
      testLog.str(std::string());//clear stringstream
    }

    thisMono->oligo_parent_->on_index_++;
    thisMono->is_indexed_ = true;
    thisMono->is_visited_ = true;
    thisMono->oligosaccharide_index_ = thisIndex;
    if(thisMono->mono_neighbors_.size() == 0)
    {
      if(local_debug > 0)
      {
        testLog << "This mono has no neighbors";
        gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
        testLog.str(std::string());//clear stringstream
      }
      //No neighbors so nothing to do
    }
    else if(thisMono->mono_neighbors_.size() == 1)
    {
      if(local_debug > 0)
      {
        testLog << "This mono has one neighbor";
        gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
        testLog.str(std::string());//clear stringstream
      }
      if(thisMono->is_root_)
      {
        if(local_debug > 0)
        {
          testLog << "This mono is the root";
          gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        if(!thisMono->mono_neighbors_[0].second->is_indexed_)
        {
          indexMono(thisMono->mono_neighbors_[0].second, thisIndex + 1, branchedMonos);
        }
      }
      else
      {
        if(local_debug > 0)
        {
          testLog << "This mono is not the root";
          gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        //No neighbors besides previous which is already indexed so nothing to do
      }
    }
    else if(thisMono->mono_neighbors_.size() == 2)
    {
      if(local_debug > 0)
      {
        testLog << "This mono has two neighbors";
        gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
        testLog.str(std::string());//clear stringstream
      }

      if(thisMono->is_root_) //BRANCHED
      {
        if(local_debug > 0)
        {
          testLog << "This mono is the root and is branched";
          gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        std::vector<int> branchMaxLengths(thisMono->mono_neighbors_.size(), 1);
        getBranchMaxLengths(thisMono, branchMaxLengths);
        cleanCountedBranches(thisMono);
        int numEqualNeighbors = 0;
        std::vector<int> equalNeighborsLengths;
        for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
        {
          if(local_debug > 0)
          {
            std::stringstream ss;
            ss << i << ", " << branchMaxLengths[i];
            gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
          }

          if(std::find(equalNeighborsLengths.begin(), equalNeighborsLengths.end(), branchMaxLengths[i]) == equalNeighborsLengths.end())
          {//if this length is not already considered
            equalNeighborsLengths.push_back(branchMaxLengths[i]);
            int branchesOfThisLength = std::count(branchMaxLengths.begin(), branchMaxLengths.end(), branchMaxLengths[i]);
            if(branchesOfThisLength > 1)
            {
              numEqualNeighbors += branchesOfThisLength;
            }
          }
        }
        if(numEqualNeighbors == 0)
        {
          //Getting shortest branch since it is closest to root in IUPAC Name
          unsigned int shortestBranchLocation = std::distance(branchMaxLengths.begin(), std::min_element(branchMaxLengths.begin(), branchMaxLengths.end()));
          unsigned int otherBranchLocation;
          //only 2 neighbors so no need for loops

          if(shortestBranchLocation == 0)
          {
            otherBranchLocation = 1;
          }
          else
          {
            otherBranchLocation = 0;
          }
          //Branches get indexed later
          this->number_branches_++;
          thisMono->mono_neighbors_[otherBranchLocation].second->oligo_branch_index_ = this->number_branches_;
          branchedMonos.push_back(thisMono->mono_neighbors_[otherBranchLocation].second);

          //Longest branch is main chain
          indexMono(thisMono->mono_neighbors_[shortestBranchLocation].second, thisIndex + 1, branchedMonos);
        }
        else
        {//The 2 branches are the same length so we need to order by linkage
          std::vector<std::string> linkageStringVector;
          for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
          {
            linkageStringVector.push_back(thisMono->mono_neighbors_[i].first->inverse_linkage_type_);
          }
          //go through vector
          //index lower linkage locant
          //go to next linkage
          //Getting highest linkage first as it is closer to root in IUPAC name
          unsigned int highestIndexLocation = std::distance(linkageStringVector.begin(), std::max_element(linkageStringVector.begin(), linkageStringVector.end()));
          unsigned int otherIndexLocation;
          if(highestIndexLocation == 0)
          {
            otherIndexLocation = 1;
          }
          else
          {
            otherIndexLocation = 0;
          }
          //Branches get indexed later
          this->number_branches_++;
          thisMono->mono_neighbors_[otherIndexLocation].second->oligo_branch_index_ = this->number_branches_;
          branchedMonos.push_back(thisMono->mono_neighbors_[otherIndexLocation].second);

          //If equal branches, lower linkage is main chain so index higher linkage first
          indexMono(thisMono->mono_neighbors_[highestIndexLocation].second, thisIndex + 1, branchedMonos);
        }
      }
      else
      {//not branched
        if(local_debug > 0)
        {
          testLog << "This mono is not branched";
          gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }
        //One of the 2 neighbors has already been indexed; Index the other
        if(thisMono->mono_neighbors_[0].second->is_indexed_)
          indexMono(thisMono->mono_neighbors_[1].second, thisIndex + 1, branchedMonos);
        else
          indexMono(thisMono->mono_neighbors_[0].second, thisIndex + 1, branchedMonos);
      }
    }
    else//more than 2 neighbors so branched
    {
      if(local_debug > 0)
      {
        testLog << "This mono has " << thisMono->mono_neighbors_.size() << " neighbors";
        gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
        testLog.str(std::string());//clear stringstream
        for(int i = 0; i < thisMono->mono_neighbors_.size(); i++)
        {
          testLog << thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId() << " indexed bool is " << thisMono->mono_neighbors_[i].second->is_indexed_;
          gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }
      }

      std::vector<int> branchMaxLengths(thisMono->mono_neighbors_.size(), 1);
      getBranchMaxLengths(thisMono, branchMaxLengths);
      cleanCountedBranches(thisMono);
      int numEqualNeighbors = 0;
      std::vector<int> equalNeighborsLengths;
      for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
      {
        if(local_debug > 0)
        {
          testLog << i << ", " << branchMaxLengths[i];
          gmml::log(__LINE__, __FILE__,  gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        if(std::find(equalNeighborsLengths.begin(), equalNeighborsLengths.end(), branchMaxLengths[i]) == equalNeighborsLengths.end()) //if this length is not already considered
        {
          equalNeighborsLengths.push_back(branchMaxLengths[i]);
          int branchesOfThisLength = std::count(branchMaxLengths.begin(), branchMaxLengths.end(), branchMaxLengths[i]);
          if(thisMono->is_root_)
          {
            if(branchesOfThisLength > 1)
            {
              numEqualNeighbors += branchesOfThisLength;
            }
          }
          else
          {
            if(branchMaxLengths[i] == 1 && branchesOfThisLength > 2)
            {//branch going through already indexed monos gets assigned 1
              numEqualNeighbors += branchesOfThisLength;
            }
            else if(branchMaxLengths[i] != 1 && branchesOfThisLength > 1)
            {
              numEqualNeighbors += branchesOfThisLength;
            }
          }
        }
      }
      if(numEqualNeighbors == 0)
      {
        if(local_debug > 0)
        {
          testLog << "This mono has all different length branches";
          gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        int shortestBranchLocation = std::distance(branchMaxLengths.begin(), std::min_element(branchMaxLengths.begin(), branchMaxLengths.end()));
        if(thisMono->mono_neighbors_[shortestBranchLocation].second->is_indexed_)//if is previous mono
        {
          branchMaxLengths[shortestBranchLocation] = 999;
          shortestBranchLocation = std::distance(branchMaxLengths.begin(), std::min_element(branchMaxLengths.begin(), branchMaxLengths.end()));
        }
        branchMaxLengths[shortestBranchLocation] = 999; //Already have the location and need to reset the max to the next highest for proper ordering
        for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
        {
          int nextShortestBranchLocation = std::distance(branchMaxLengths.begin(), std::min_element(branchMaxLengths.begin(), branchMaxLengths.end()));
          if(local_debug > 0)
          {
            testLog << "This branch length: " << branchMaxLengths[nextShortestBranchLocation];
            gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
            testLog.str(std::string());//clear stringstream
          }

          if((branchMaxLengths[nextShortestBranchLocation] != 999) && (!thisMono->mono_neighbors_[nextShortestBranchLocation].second->is_indexed_))
          {
            if(local_debug > 0)
            {
              testLog << thisMono->mono_neighbors_[nextShortestBranchLocation].second->cycle_atoms_[0]->GetResidue()->GetId();
              gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
              testLog.str(std::string());//clear stringstream
            }
            this->number_branches_++;
            thisMono->mono_neighbors_[nextShortestBranchLocation].second->oligo_branch_index_ = this->number_branches_;
            branchedMonos.push_back(thisMono->mono_neighbors_[nextShortestBranchLocation].second);
            branchMaxLengths[nextShortestBranchLocation] = 999;
          }
          else
          {
            if(local_debug > 0)
            {
              testLog << "This did not get added to branched monos";
              gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
              testLog.str(std::string());//clear stringstream
            }
            branchMaxLengths[nextShortestBranchLocation] = 999;
          }
        }
        indexMono(thisMono->mono_neighbors_[shortestBranchLocation].second, thisIndex + 1, branchedMonos);
      }
      else
      {
        if(local_debug > 0)
        {
          testLog << "This mono has " << numEqualNeighbors << " equal length branches";
          gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        //handle equal length branch ordering by linkage type
        std::vector<std::string> linkageStringVector;
        for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
        {
          if(thisMono->mono_neighbors_[i].second->is_indexed_)
          {
            if(local_debug > 0)
            {
              testLog << thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId() << "has been indexed";
              gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
              testLog.str(std::string());//clear stringstream
            }
            branchMaxLengths[i] = 999;
          }
        }
        unsigned int shortestBranchLocation = std::distance(branchMaxLengths.begin(), std::min_element(branchMaxLengths.begin(), branchMaxLengths.end()));
        if((std::count(equalNeighborsLengths.begin(), equalNeighborsLengths.end(), branchMaxLengths[shortestBranchLocation]) == 0) &&  (branchMaxLengths[shortestBranchLocation] != 999))//no other branches of this length
        {
          branchMaxLengths[shortestBranchLocation] = 999;//Have location and need to get next shortest
          for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
          {
            int nextShortestBranchLocation = std::distance(branchMaxLengths.begin(), std::min_element(branchMaxLengths.begin(), branchMaxLengths.end()));
            if((branchMaxLengths[nextShortestBranchLocation] != 999) && (!thisMono->mono_neighbors_[nextShortestBranchLocation].second->is_indexed_))
            {
              this->number_branches_++;
              thisMono->mono_neighbors_[nextShortestBranchLocation].second->oligo_branch_index_ = this->number_branches_;
              branchedMonos.push_back(thisMono->mono_neighbors_[nextShortestBranchLocation].second);
              branchMaxLengths[nextShortestBranchLocation] = 999;
            }
            else if (thisMono->mono_neighbors_[nextShortestBranchLocation].second->is_indexed_)
            {
              if(local_debug > 0)
              {
                testLog << "This mono has been indexed";
                gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
                testLog.str(std::string());//clear stringstream
              }
            }
            else
            {
              if(local_debug > 0)
                {
                  testLog << "This mono branch length: " << branchMaxLengths[nextShortestBranchLocation];
                  gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
                  testLog.str(std::string());//clear stringstream
                }
            }
          }
          if(local_debug > 0)
          {
            testLog << "This branch does not have any equal length branches";
            gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
            testLog.str(std::string());//clear stringstream
          }
          indexMono(thisMono->mono_neighbors_[shortestBranchLocation].second, thisIndex + 1, branchedMonos);
        }
        else if(branchMaxLengths[shortestBranchLocation] != 999)
        {//equal branches, but not already considered
          int equalBranchLength = branchMaxLengths[shortestBranchLocation];
          branchMaxLengths[shortestBranchLocation] = 999;
          if(local_debug > 0)
          {
            testLog << "This branch does have equal length branches of length " << equalBranchLength;
            gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
            testLog.str(std::string());//clear stringstream
            testLog << "The shortest branch is at: " << shortestBranchLocation;
            gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
            testLog.str(std::string());//clear stringstream
          }

          for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
          {
            if(local_debug > 0)
            {
              testLog << thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId() << "This branch is length: " << branchMaxLengths[i];
              gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
              testLog.str(std::string());//clear stringstream
              testLog <<  "On neighbor: " << i;
              gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
              testLog.str(std::string());//clear stringstream
            }

            if((branchMaxLengths[i] == equalBranchLength) && (!thisMono->mono_neighbors_[i].second->is_indexed_))
            {//This branch is not yet considered and is equal length
              if(local_debug > 0)
              {
                testLog << thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId() << " This branch is the same length as the shortest";
                gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
                testLog.str(std::string());//clear stringstream
              }
              linkageStringVector.push_back(thisMono->mono_neighbors_[i].first->inverse_linkage_type_);
            }
            else if(thisMono->mono_neighbors_[i].second->is_indexed_)
            {//Has already been considered
              if(local_debug > 0)
              {
                testLog << thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId() << "This branch is indexed";
                gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
                testLog.str(std::string());//clear stringstream
              }
              linkageStringVector.push_back("-999"); //no linkage starts with -999; will be lower than all others and we traverse highest -> lowest linkage
            }
            else if (i == shortestBranchLocation)
            {
              if(local_debug > 0)
              {
                testLog << thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId() << " This is the shortest branch";
                gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
                testLog.str(std::string());//clear stringstream
              }
              linkageStringVector.push_back(thisMono->mono_neighbors_[i].first->inverse_linkage_type_);
              // equalBranches.push_back(thisMono->mono_neighbors_[i]);
            }
            else
            {
              if(local_debug > 0)
              {
                testLog << thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId() << " This branch is length: " << branchMaxLengths[i];
                gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
                testLog.str(std::string());//clear stringstream
              }
              linkageStringVector.push_back("-999");
            }
          }
        }
        //go through vector
        //index higher linkage locant
        if(local_debug > 0)
        {
          testLog << "Linkage vector is " << linkageStringVector.size() << " long";
          gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        int highestIndexLocation = std::distance(linkageStringVector.begin(), std::max_element(linkageStringVector.begin(), linkageStringVector.end()));
        // if(thisMono->mono_neighbors_[lowestIndexLocation].second->is_indexed_)//if is previous mono
        // {
        //   linkageStringVector[lowestIndexLocation] = "999";
        //   lowestIndexLocation = std::distance(linkageStringVector.begin(), std::min_element(linkageStringVector.begin(), linkageStringVector.end()));
        // }
        linkageStringVector[highestIndexLocation] = "-999";
        for(unsigned int i = 0; i < linkageStringVector.size(); i++)
        {
          int nextHighestIndexLocation = std::distance(linkageStringVector.begin(), std::max_element(linkageStringVector.begin(), linkageStringVector.end()));
          if((!thisMono->mono_neighbors_[nextHighestIndexLocation].second->is_indexed_) && (linkageStringVector[nextHighestIndexLocation] != "-999"))//if isnt previous mono and not already considered
          {
            this->number_branches_++;
            thisMono->mono_neighbors_[nextHighestIndexLocation].second->oligo_branch_index_ = this->number_branches_;
            branchedMonos.push_back(thisMono->mono_neighbors_[nextHighestIndexLocation].second);
            linkageStringVector[nextHighestIndexLocation] = "-999";
          }
        }
        indexMono(thisMono->mono_neighbors_[highestIndexLocation].second, thisIndex + 1, branchedMonos);
      }
    }
  }
}

void Glycan::Oligosaccharide::indexBranches(std::vector<Glycan::Monosaccharide*> branchedMonos)
{
  int local_debug = -1;
  if(local_debug > 0)
  {
    std::stringstream debugSS;
    debugSS << "On Oligo index: " << this->on_index_;
    gmml::log(__LINE__, __FILE__, gmml::INF, debugSS.str());
  }
  for(std::vector<Glycan::Monosaccharide*>::iterator mono = branchedMonos.begin(); mono != branchedMonos.end(); mono++ )
  {
    Glycan::Monosaccharide* thisMono = *mono;
    std::vector<Glycan::Monosaccharide*> subBranchedMonos;
    indexMono(thisMono, this->on_index_, subBranchedMonos);
    if(subBranchedMonos.size() > 0)
    {
      indexBranches(subBranchedMonos);
    }
  }
}
std::vector<Glycan::Monosaccharide*> Glycan::Oligosaccharide::IUPACindexMono(Glycan::Monosaccharide* thisMono, int thisIndex, std::vector<Glycan::Monosaccharide*> &branchedMonos)
{
  int local_debug = -1;
  std::stringstream testLog;
  if(!thisMono->is_indexed_)
  {
    thisMono->IUPAC_index_ = thisIndex;
    thisMono->oligo_parent_->on_IUPAC_index_++;
    thisMono->is_indexed_ = true;
    thisMono->is_visited_ = true;

    if(local_debug > 0)
    {
      testLog << "Mono: " << thisMono->cycle_atoms_[0]->GetResidue()->GetId();
      testLog << ", Index: " << thisIndex;
      gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
      testLog.str(std::string());//clear stringstream
    }

    if(thisMono->mono_neighbors_.size() == 0)
    {
      if(local_debug > 0)
      {
        testLog << "This mono has no neighbors";
        gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
        testLog.str(std::string());//clear stringstream
      }
      // No need to do anything if no neighbors
    }
    else if(thisMono->mono_neighbors_.size() == 1)
    {
      if(local_debug > 0)
      {
        testLog << "This mono has one neighbor";
        gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
        testLog.str(std::string());//clear stringstream
      }

      if(thisMono->is_root_)
      {
        if(local_debug > 0)
        {
          testLog << "This mono is the root";
          gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        if(!thisMono->mono_neighbors_[0].second->is_indexed_)
        {
          IUPACindexMono(thisMono->mono_neighbors_[0].second, thisIndex + 1, branchedMonos);
        }
      }
      else
      {
        if(local_debug > 0)
        {
          testLog << "This mono is not the root";
          gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        //No need to do anything; the previous mono called this function and there is no next mono
      }
    }
    else if(thisMono->mono_neighbors_.size() == 2)
    {
      if(local_debug > 0)
      {
        testLog << "This mono has two neighbors";
        gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
        testLog.str(std::string());//clear stringstream
      }
      if(thisMono->is_root_) //BRANCHED root
      {
        if(local_debug > 0)
        {
          testLog << "This mono is the root and is branched";
          gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        std::vector<int> branchMaxLengths(thisMono->mono_neighbors_.size(), 1);
        getBranchMaxLengths(thisMono, branchMaxLengths);
        cleanCountedBranches(thisMono);
        int numEqualNeighbors = 0;
        std::vector<int> equalNeighborsLengths;
        for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
        {
          if(local_debug > 0)
          {
            testLog << i << ", " << branchMaxLengths[i];
            gmml::log(__LINE__, __FILE__,  gmml::INF, testLog.str());
            testLog.str(std::string());//clear stringstream
          }

          if(std::find(equalNeighborsLengths.begin(), equalNeighborsLengths.end(), branchMaxLengths[i]) == equalNeighborsLengths.end())
          {//if this length is not already considered
            equalNeighborsLengths.push_back(branchMaxLengths[i]);
            int branchesOfThisLength = std::count(branchMaxLengths.begin(), branchMaxLengths.end(), branchMaxLengths[i]);
            if(branchesOfThisLength > 1)
            {
              numEqualNeighbors += branchesOfThisLength;
            }
          }
        }
        if(numEqualNeighbors == 0)
        {
          unsigned int longestBranchLocation = std::distance(branchMaxLengths.begin(), std::max_element(branchMaxLengths.begin(), branchMaxLengths.end()));
          unsigned int otherBranchLocation;
          //only 2 neighbors so no need for loops

          if(longestBranchLocation == 0)
          {
            otherBranchLocation = 1;
          }
          else
          {
            otherBranchLocation = 0;
          }
          //Branches get indexed later
          this->number_branches_++;
          thisMono->mono_neighbors_[otherBranchLocation].second->oligo_branch_index_ = this->number_branches_;
          branchedMonos.push_back(thisMono->mono_neighbors_[otherBranchLocation].second);

          //Longest branch is main chain
          IUPACindexMono(thisMono->mono_neighbors_[longestBranchLocation].second, thisIndex + 1, branchedMonos);
        }
        else
        {//The 2 branches are the same length so we need to order by linkage
          std::vector<std::string> linkageStringVector;
          for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
          {
            linkageStringVector.push_back(thisMono->mono_neighbors_[i].first->inverse_linkage_type_);
          }
          //go through vector
          //index lower linkage locant
          //go to next linkage
          unsigned int lowestIndexLocation = std::distance(linkageStringVector.begin(), std::min_element(linkageStringVector.begin(), linkageStringVector.end()));
          unsigned int otherIndexLocation;
          if(lowestIndexLocation == 0)
          {
            otherIndexLocation = 1;
          }
          else
          {
            otherIndexLocation = 0;
          }
          //Branches get indexed later
          this->number_branches_++;
          thisMono->mono_neighbors_[otherIndexLocation].second->oligo_branch_index_ = this->number_branches_;
          branchedMonos.push_back(thisMono->mono_neighbors_[otherIndexLocation].second);

          //If equal branches, lower linkage is main chain
          IUPACindexMono(thisMono->mono_neighbors_[lowestIndexLocation].second, thisIndex + 1, branchedMonos);
        }
      }
      else
      {
        if(local_debug > 0)
        {
          testLog << "This mono is not branched";
          gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        //One of the neighbors is already indexed, move to the other
        if(thisMono->mono_neighbors_[0].second->is_indexed_)
          IUPACindexMono(thisMono->mono_neighbors_[1].second, thisIndex + 1, branchedMonos);
        else
          IUPACindexMono(thisMono->mono_neighbors_[0].second, thisIndex + 1, branchedMonos);
      }
    }
    else//more than 2 neighbors so branched
    {
      if(local_debug > 0)
      {
        testLog << "This mono has " << thisMono->mono_neighbors_.size() << " neighbors";
        gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
        testLog.str(std::string());//clear stringstream
      }

      std::vector<int> branchMaxLengths(thisMono->mono_neighbors_.size(), 1);
      getBranchMaxLengths(thisMono, branchMaxLengths);
      cleanCountedBranches(thisMono);
      int numEqualNeighbors = 0;
      std::vector<int> equalNeighborsLengths;
      for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
      {
        if(local_debug > 0)
        {
          testLog << i << ", " << branchMaxLengths[i];
          gmml::log(__LINE__, __FILE__,  gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        if(std::find(equalNeighborsLengths.begin(), equalNeighborsLengths.end(), branchMaxLengths[i]) == equalNeighborsLengths.end())
        {//if this length is not already considered
          equalNeighborsLengths.push_back(branchMaxLengths[i]);
          int branchesOfThisLength = std::count(branchMaxLengths.begin(), branchMaxLengths.end(), branchMaxLengths[i]);
          if(thisMono->is_root_)
          {
            if(branchesOfThisLength > 1)
            {
              numEqualNeighbors += branchesOfThisLength;
            }
          }
          else
          {
            if(branchMaxLengths[i] == 1 && branchesOfThisLength > 2)
            {
              numEqualNeighbors += branchesOfThisLength;
            }
            else if(branchMaxLengths[i] != 1 && branchesOfThisLength > 1)
            {
              numEqualNeighbors += branchesOfThisLength;
            }
          }
        }
      }
      if(numEqualNeighbors == 0)
      {
        if(local_debug > 0)
        {
          testLog << "This mono has all different length branches";
          gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        int longestBranchLocation = std::distance(branchMaxLengths.begin(), std::max_element(branchMaxLengths.begin(), branchMaxLengths.end()));
        if(thisMono->mono_neighbors_[longestBranchLocation].second->is_indexed_)
        {//if is previous mono, clear and find next longest branch
          branchMaxLengths[longestBranchLocation] = -1;
          longestBranchLocation = std::distance(branchMaxLengths.begin(), std::max_element(branchMaxLengths.begin(), branchMaxLengths.end()));
        }
        branchMaxLengths[longestBranchLocation] = -1; //Already have the location and need to reset the max to the next highest for proper ordering
        for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
        {
          int nextLongestBranchLocation = std::distance(branchMaxLengths.begin(), std::max_element(branchMaxLengths.begin(), branchMaxLengths.end()));
          if(local_debug > 0)
          {
            testLog << "This branch length: " << branchMaxLengths[nextLongestBranchLocation];
            gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
            testLog.str(std::string());//clear stringstream
          }

          if((branchMaxLengths[nextLongestBranchLocation] != -1) && (!thisMono->mono_neighbors_[nextLongestBranchLocation].second->is_indexed_))
          {
            if(local_debug > 0)
            {
              testLog << thisMono->mono_neighbors_[nextLongestBranchLocation].second->cycle_atoms_[0]->GetResidue()->GetId();
              gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
              testLog.str(std::string());//clear stringstream
            }
            this->number_branches_++;
            thisMono->mono_neighbors_[nextLongestBranchLocation].second->oligo_branch_index_ = this->number_branches_;
            branchedMonos.push_back(thisMono->mono_neighbors_[nextLongestBranchLocation].second);
            branchMaxLengths[nextLongestBranchLocation] = -1;
          }
          else
          {
            if(local_debug > 0)
            {
              testLog << "This did not get added to branched monos";
              gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
              testLog.str(std::string());//clear stringstream
            }
            branchMaxLengths[nextLongestBranchLocation] = -1;
          }
        }
        IUPACindexMono(thisMono->mono_neighbors_[longestBranchLocation].second, thisIndex + 1, branchedMonos);
      }
      else
      {
        if(local_debug > 0)
        {
          testLog << "This mono has " << numEqualNeighbors << " equal length branches";
          gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        std::vector<std::string> linkageStringVector;
        for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
        {
          if(thisMono->mono_neighbors_[i].second->is_indexed_)
          {
            if(local_debug > 0)
            {
              testLog << thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId() << "has been indexed";
              gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
              testLog.str(std::string());//clear stringstream
            }
            branchMaxLengths[i] = -1;
          }
        }
        unsigned int longestBranchLocation = std::distance(branchMaxLengths.begin(), std::max_element(branchMaxLengths.begin(), branchMaxLengths.end()));

        if((std::count(equalNeighborsLengths.begin(), equalNeighborsLengths.end(), branchMaxLengths[longestBranchLocation]) == 0) &&  (branchMaxLengths[longestBranchLocation] != -1))//no other branches of this length and this hasn't been considered
        {
          branchMaxLengths[longestBranchLocation] = -1; //Have location and need to get next longest
          for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
          {
            int nextLongestBranchLocation = std::distance(branchMaxLengths.begin(), std::max_element(branchMaxLengths.begin(), branchMaxLengths.end()));
            if((branchMaxLengths[nextLongestBranchLocation] != -1) && (!thisMono->mono_neighbors_[nextLongestBranchLocation].second->is_indexed_))
            {
              if(local_debug > 0)
              {
                testLog << "This should work...";
                gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
                testLog.str(std::string());//clear stringstream
              }
              this->number_branches_++;
              thisMono->mono_neighbors_[nextLongestBranchLocation].second->oligo_branch_index_ = this->number_branches_;
              branchedMonos.push_back(thisMono->mono_neighbors_[nextLongestBranchLocation].second);
              branchMaxLengths[nextLongestBranchLocation] = -1;
            }
            else if (thisMono->mono_neighbors_[nextLongestBranchLocation].second->is_indexed_)
            {
              if(local_debug > 0)
              {
                testLog << "This mono has been indexed";
                gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
                testLog.str(std::string());//clear stringstream
              }
            }
            else
            {
              if(local_debug > 0)
              {
                testLog << "This mono branch length: " << branchMaxLengths[nextLongestBranchLocation];
                gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
                testLog.str(std::string());//clear stringstream
              }
            }
          }
          if(local_debug > 0)
          {
            testLog << "This branch does not have any equal length branches";
            gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
            testLog.str(std::string());//clear stringstream
          }
          IUPACindexMono(thisMono->mono_neighbors_[longestBranchLocation].second, thisIndex + 1, branchedMonos);
        }
        else if(branchMaxLengths[longestBranchLocation] != -1)
        {//equal branches, but not already considered
          int equalBranchLength = branchMaxLengths[longestBranchLocation];
          branchMaxLengths[longestBranchLocation] = -1;
          if(local_debug > 0)
          {
            testLog << "This branch does have equal length branches of length " << equalBranchLength;
            gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
            testLog.str(std::string());//clear stringstream
            testLog << "The longest branch is at: " << longestBranchLocation;
            gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
            testLog.str(std::string());//clear stringstream
          }

          for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
          {
            if(local_debug > 0)
            {
              testLog << thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId() << "This branch is length: " << branchMaxLengths[i];
              gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
              testLog.str(std::string());//clear stringstream
              testLog <<  "On neighbor: " << i;
              gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
              testLog.str(std::string());//clear stringstream
            }

            if((branchMaxLengths[i] == equalBranchLength) && (!thisMono->mono_neighbors_[i].second->is_indexed_))
            {//This branch is not yet considered and is equal length
              if(local_debug > 0)
              {
                testLog << thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId() << " This branch is the same length as the longest";
                gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
                testLog.str(std::string());//clear stringstream
              }

              linkageStringVector.push_back(thisMono->mono_neighbors_[i].first->inverse_linkage_type_);
            }
            else if(thisMono->mono_neighbors_[i].second->is_indexed_)
            {//Already considered
              if(local_debug > 0)
              {
                testLog << thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId() << "This branch is indexed";
                gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
                testLog.str(std::string());//clear stringstream
              }

              linkageStringVector.push_back("999"); //no linkage starts with 999; will be higher than all others and we traverse lowest to highest linkage
            }
            else if (i == longestBranchLocation)
            {
              if(local_debug > 0)
              {
                testLog << thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId() << " This is the longest branch";
                gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
                testLog.str(std::string());//clear stringstream
              }

              linkageStringVector.push_back(thisMono->mono_neighbors_[i].first->inverse_linkage_type_);
            }
            else
            {
              if(local_debug > 0)
              {
                testLog << thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId() << " This branch is length: " << branchMaxLengths[i];
                gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
                testLog.str(std::string());//clear stringstream
              }

              linkageStringVector.push_back("999");
            }
          }
        }
        //go through vector
        //index lower linkage locant
        if(local_debug > 0)
        {
          testLog << "Linkage vector is " << linkageStringVector.size() << " long";
          gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        int lowestIndexLocation = std::distance(linkageStringVector.begin(), std::min_element(linkageStringVector.begin(), linkageStringVector.end()));
        linkageStringVector[lowestIndexLocation] = "999";
        for(unsigned int i = 0; i < linkageStringVector.size(); i++)
        {
          int nextLowestIndexLocation = std::distance(linkageStringVector.begin(), std::min_element(linkageStringVector.begin(), linkageStringVector.end()));
          if((!thisMono->mono_neighbors_[nextLowestIndexLocation].second->is_indexed_) && (linkageStringVector[nextLowestIndexLocation] != "999"))//if isnt previous mono and not considered
          {
            this->number_branches_++;
            thisMono->mono_neighbors_[nextLowestIndexLocation].second->oligo_branch_index_ = this->number_branches_;
            branchedMonos.push_back(thisMono->mono_neighbors_[nextLowestIndexLocation].second);
            linkageStringVector[nextLowestIndexLocation] = "999";
          }
        }
        IUPACindexMono(thisMono->mono_neighbors_[lowestIndexLocation].second, thisIndex + 1, branchedMonos);
      }
    }
  }
  return branchedMonos;
}

void Glycan::Oligosaccharide::IUPACindexBranches(std::vector<Glycan::Monosaccharide*> branchedMonos)
{
  for(std::vector<Glycan::Monosaccharide*>::iterator mono = branchedMonos.begin(); mono != branchedMonos.end(); mono++ )
  {
    Glycan::Monosaccharide* thisMono = *mono;
    std::vector<Glycan::Monosaccharide*> subBranchedMonos;
    subBranchedMonos = IUPACindexMono(thisMono, this->on_IUPAC_index_, subBranchedMonos);
    if(subBranchedMonos.size() > 0)
    {
      IUPACindexBranches(subBranchedMonos);
    }
  }
}
void Glycan::Oligosaccharide::traverseGraph(Glycan::Monosaccharide* thisMono, Glycan::Oligosaccharide* thisOligo)
{
  int local_debug = -1;
  std::stringstream testLog;
  thisMono->oligo_parent_ = thisOligo;
  std::string thisMonoName = thisMono->sugar_name_.monosaccharide_short_name_;
  std::string thisMonoAuthorName = thisMono->author_sugar_name_.monosaccharide_short_name_;
  bool countedBackwards = false;
  if(thisMonoName.size() == 0)
  {
    thisMonoName = thisMono->sugar_name_.monosaccharide_stereochemistry_short_name_;
  }
  if(thisMonoName.size() == 0) // still
  {
    thisMonoName = thisMono->sugar_name_.monosaccharide_stereochemistry_name_;
  }
  if(thisMonoName.size() == 0) // still
  {

  }
  if(thisMonoAuthorName.size() == 0)
  {
    thisMonoAuthorName = thisMono->author_sugar_name_.monosaccharide_stereochemistry_short_name_;
  }
  if(thisMonoAuthorName.size() == 0) // still
  {
    thisMonoAuthorName = thisMono->author_sugar_name_.monosaccharide_stereochemistry_name_;
  }
  if(local_debug > 0)
  {
    testLog << thisMono->cycle_atoms_[0]->GetResidue()->GetId() << " has " << thisMono->mono_neighbors_.size() << " neighbors.";
    gmml::log(__LINE__, __FILE__,  gmml::INF, testLog.str());
    testLog.str(std::string());//clear stringstream
  }
   //Name terminal and monosaccharide if no neighbors
  if(thisMono->is_root_)
  {
    thisOligo->root_ = thisMono;
    thisOligo->IUPAC_name_ = "";
    thisOligo->author_IUPAC_name_ = "";
    countedBackwards = true;
    std::vector<MolecularModeling::Atom*> terminal_atoms = std::vector<MolecularModeling::Atom*> ();
    std::string terminal;
    MolecularModeling::Atom* anomeric_o = NULL;
    if(thisMono->side_atoms_.size() != 0)
    {
      if(thisMono->side_atoms_.at(0).at(1) != NULL)
      {
        anomeric_o = thisMono->side_atoms_.at(0).at(1);
      }
    }
    if(anomeric_o != NULL)
    {
      terminal = " " + CheckTerminals(anomeric_o, terminal_atoms);//nameTerminal(thisMono);

      if(terminal == " Unknown")
      {
        if(anomeric_o->GetResidue()->CheckIfProtein())
        {
          terminal = anomeric_o->GetResidue()->GetName();
        }
      }
      else
      {
        terminal = gmml::Trim(terminal);
      }
    }
    if(terminal == " Unknown")
    {
      if(thisMono->anomeric_carbon_pointer_ != NULL)
      {
        terminal = thisMono->GetFormula(thisMono->anomeric_carbon_pointer_);
        //add terminal as derivative to this mono at anomeric carbon
        for(std::vector<std::string>::iterator it = thisMono->chemical_code_->right_down_.begin(); it != thisMono->chemical_code_->right_down_.end(); it++)
        {
          if((*it).find("a")!= std::string::npos)
          {
            thisMono->unknown_derivatives_.push_back(std::make_pair("a", terminal));
            thisMono->derivatives_map_.push_back(std::make_pair("a", ""));
          }
        }
        for(std::vector<std::string>::iterator it = thisMono->chemical_code_->right_up_.begin(); it != thisMono->chemical_code_->right_up_.end(); it++)
        {
          if((*it).find("a")!= std::string::npos)
          {
            thisMono->unknown_derivatives_.push_back(std::make_pair("a", terminal));
            thisMono->derivatives_map_.push_back(std::make_pair("a", ""));
          }
        }
        // for(std::vector<std::pair<std::string, std::string> >::iterator it = thisMono->derivatives_map_.begin(); it != thisMono->derivatives_map_.end(); it++)
        // {
        //   std::string key = (*it).first();
        //   std::string value = (*it).second();
        //   if(key == "a")
        //   {
        //     value = terminal;
        //   }
        // }
        //Rename sugar
        thisMono->GenerateCompleteSugarName(assembly_);
        if(local_debug > 0)
        {
          gmml::log(__LINE__, __FILE__,  gmml::INF, thisMono->sugar_name_.monosaccharide_name_);
        }
      }
      if(local_debug > 0)
      {
        gmml::log(__LINE__, __FILE__,  gmml::INF, "Named terminal by Formula");
      }
    }
    std::string anomeric_carbon_id;
    if(thisMono->anomeric_carbon_pointer_ != NULL)
    {
      anomeric_carbon_id = thisMono->anomeric_carbon_pointer_->GetId();
      anomeric_carbon_id = anomeric_carbon_id[1];
    }
    else
    {
      anomeric_carbon_id  = "?";
    }
    terminal = anomeric_carbon_id + "-" + terminal;
    thisOligo->mono_nodes_.push_back(thisMono);
    // thisMono->oligosaccharide_index_ = 0;
    if(thisMono->mono_neighbors_.size() == 0) //Just a monosaccharide
    {
      thisOligo->IUPAC_name_ = thisMonoName + terminal;
      thisOligo->author_IUPAC_name_ = thisMonoAuthorName + terminal;
    }
    else if((thisMono->mono_neighbors_.size() == 1) && (thisMono->mono_neighbors_[0].first->anomeric_anomeric_linkage_ == false))
    {
      if(local_debug > 0)
      {
        gmml::log(__LINE__, __FILE__,  gmml::INF, "Not anomeric");
      }
      thisOligo->IUPAC_name_ = terminal;
      thisOligo->author_IUPAC_name_ = terminal;
    }
    else if (thisMono->mono_neighbors_.size() > 1)
    {
      bool isAnomeric = false;
      for (std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> >::iterator it = thisMono->mono_neighbors_.begin(); it != thisMono->mono_neighbors_.end(); it++)
      {
        Glycan::GlycosidicLinkage* thisLinkage = (*it).first;
        if(thisLinkage->anomeric_anomeric_linkage_ == true)
        {
          isAnomeric = true;
        }
      }
      if(isAnomeric == false)
      {
        if(local_debug > 0)
        {
          gmml::log(__LINE__, __FILE__,  gmml::INF, "Not anomeric");
        }
        thisOligo->IUPAC_name_ = terminal;
        thisOligo->author_IUPAC_name_ = terminal;
      }
    }
  }


  if(thisMono->mono_neighbors_.size() == 0)
  {

  }
  else if(thisMono->mono_neighbors_.size() == 1) //no branching
  {
    if(local_debug > 0)
    {
      testLog << "This mono has one neighbor";
      gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
      testLog.str(std::string());//clear stringstream
    }
    if(thisMono->is_root_)
    {
      if(local_debug > 0)
      {
        testLog << "This mono is the root";
        gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
        testLog.str(std::string());//clear stringstream
      }
      if(!thisMono->mono_neighbors_[0].second->is_visited_)//temp fix to handle anomeric anomeric infinite loop
      {
        Glycan::Monosaccharide* thisMonoNeighbor = thisMono->mono_neighbors_[0].second;
        thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[0].first->inverse_linkage_type_ +
                                 thisMonoName + thisOligo->IUPAC_name_;
        thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[0].first->inverse_linkage_type_ +
                                        thisMonoAuthorName + thisOligo->author_IUPAC_name_;
        thisMono->is_visited_ = true;
        // thisMonoNeighbor->is_visited_ = true;
        // thisMonoNeighbor->oligosaccharide_index_ = 1;
        thisOligo->mono_nodes_.push_back(thisMonoNeighbor);
        if(local_debug > 0)
        {
          gmml::log(__LINE__, __FILE__,  gmml::INF, thisOligo->IUPAC_name_);
        }
        traverseGraph(thisMonoNeighbor, thisOligo);
      }
    }
    else//no more neighbors (one is attached already, and this isn't root)
    {
      thisMono->is_visited_ = true;
      if(local_debug > 0)
      {
        testLog << "This mono is not the root";
        gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
        testLog.str(std::string());//clear stringstream
      }
      thisOligo->IUPAC_name_ = thisMonoName + thisOligo->IUPAC_name_;
      thisOligo->author_IUPAC_name_ = thisMonoAuthorName + thisOligo->author_IUPAC_name_;
      if(local_debug > 0)
      {
        gmml::log(__LINE__, __FILE__,  gmml::INF, thisOligo->IUPAC_name_);
      }
    }
  }
  else if(thisMono->mono_neighbors_.size() == 2)
  {
    if(local_debug > 0)
    {
      testLog << "This mono has two neighbors";
      gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
      testLog.str(std::string());//clear stringstream
    }

    if(thisMono->is_root_) //BRANCHED
    {
      if(local_debug > 0)
      {
        testLog << "This mono is the root and is branched";
        gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
        testLog.str(std::string());//clear stringstream
      }

      thisOligo->IUPAC_name_ = thisMonoName + thisOligo->IUPAC_name_;
      thisOligo->author_IUPAC_name_ = thisMonoAuthorName + thisOligo->author_IUPAC_name_;
      std::vector<int> branchMaxLengths(thisMono->mono_neighbors_.size(), 1);
      getBranchMaxLengths(thisMono, branchMaxLengths);
      cleanCountedBranches(thisMono);
      int numEqualNeighbors = 0;
      std::vector<int> equalNeighborsLengths;
      for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
      {
        if(local_debug > 0)
        {
          testLog << i << ", " << branchMaxLengths[i];
          gmml::log(__LINE__, __FILE__,  gmml::INF, testLog.str());
          testLog.str(std::string());//clear stringstream
        }

        if(std::find(equalNeighborsLengths.begin(), equalNeighborsLengths.end(), branchMaxLengths[i]) == equalNeighborsLengths.end()) //if this length is not already considered
        {
          equalNeighborsLengths.push_back(branchMaxLengths[i]);
          int branchesOfThisLength = std::count(branchMaxLengths.begin(), branchMaxLengths.end(), branchMaxLengths[i]);
          if(branchesOfThisLength > 1)
          {
            numEqualNeighbors += branchesOfThisLength;
          }
        }
      }
      if(numEqualNeighbors == 0)
      {
        int shortestBranchLocation = std::distance(branchMaxLengths.begin(), std::min_element(branchMaxLengths.begin(), branchMaxLengths.end()));
        for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
        {
          shortestBranchLocation = std::distance(branchMaxLengths.begin(), std::min_element(branchMaxLengths.begin(), branchMaxLengths.end()));
          if(i < thisMono->mono_neighbors_.size() - 1) //no countedBackwards flag as there is nothing backwards from root
          {
            thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ + "]" +
                                    thisOligo->IUPAC_name_;
            thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ + "]" +
                                            thisOligo->author_IUPAC_name_;
          }
          else//last branch doesnt get brackets
          {
            thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ +
                                    thisOligo->IUPAC_name_;
            thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ +
                                    thisOligo->author_IUPAC_name_;
          }
          thisMono->is_visited_ = true;
          thisOligo->mono_nodes_.push_back(thisMono->mono_neighbors_[shortestBranchLocation].second);
          traverseGraph(thisMono->mono_neighbors_[shortestBranchLocation].second, thisOligo);
          if(i < thisMono->mono_neighbors_.size() - 1)
          {
            if(local_debug > 0)
            {
              gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
            }

            thisOligo->IUPAC_name_ = "[" + thisOligo->IUPAC_name_;
            thisOligo->author_IUPAC_name_ = "[" + thisOligo->author_IUPAC_name_;
          }
          branchMaxLengths[shortestBranchLocation] = 9999;
        }
      }
      else
      {
        std::vector<std::string> linkageStringVector;
        int shortestBranchLocation = std::distance(branchMaxLengths.begin(), std::min_element(branchMaxLengths.begin(), branchMaxLengths.end()));
        for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
        {
          shortestBranchLocation = std::distance(branchMaxLengths.begin(), std::min_element(branchMaxLengths.begin(), branchMaxLengths.end()));
          if(std::count(equalNeighborsLengths.begin(), equalNeighborsLengths.end(), branchMaxLengths[shortestBranchLocation]) == 0) //no other branches of this length
          {
            linkageStringVector.push_back("0");
            if((i < thisMono->mono_neighbors_.size() - 1) && (!thisMono->mono_neighbors_[shortestBranchLocation].second->is_visited_))
            {
              thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ + "]" +
                                       thisOligo->IUPAC_name_;
              thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ + "]" +
                                       thisOligo->author_IUPAC_name_;
            }
            else if (!thisMono->mono_neighbors_[shortestBranchLocation].second->is_visited_)
            {
              thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ +
                                      thisOligo->IUPAC_name_;
              thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ +
                                      thisOligo->author_IUPAC_name_;
            }
            thisMono->is_visited_ = true;
            thisOligo->mono_nodes_.push_back(thisMono->mono_neighbors_[shortestBranchLocation].second);
            traverseGraph(thisMono->mono_neighbors_[shortestBranchLocation].second, thisOligo);
            if(i < thisMono->mono_neighbors_.size() - 1)
            {
              if(local_debug > 0)
              {
                gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
              }
              thisOligo->IUPAC_name_ = "[" + thisOligo->IUPAC_name_;
              thisOligo->author_IUPAC_name_ = "[" + thisOligo->author_IUPAC_name_;
            }
            branchMaxLengths[shortestBranchLocation] = 9999; //arbitralily large number so it moves to the next shortest branch
          }
          else
          {
            int equalBranchLength = branchMaxLengths[shortestBranchLocation];
            if((branchMaxLengths[i] == equalBranchLength) && (!thisMono->mono_neighbors_[i].second->is_visited_))
            {
              if(local_debug > 0)
              {
                gmml::log(__LINE__, __FILE__,  gmml::INF, thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId());
              }

              linkageStringVector.push_back(thisMono->mono_neighbors_[i].first->inverse_linkage_type_);
              // equalBranches.push_back(thisMono->mono_neighbors_[i]);
            }
            else if((branchMaxLengths[i] == equalBranchLength)&& (thisMono->mono_neighbors_[i].second->is_visited_))
            {
              if(local_debug > 0)
              {
                gmml::log(__LINE__, __FILE__,  gmml::INF, thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId());
              }

              linkageStringVector.push_back("0");
            }
            else
            {
              linkageStringVector.push_back("0");
            }
            branchMaxLengths[shortestBranchLocation] = 9999;
          }
        }
        for(unsigned int i = 0; i < linkageStringVector.size(); i++)
        {
          int highestIndexLocation = std::distance(linkageStringVector.begin(), std::max_element(linkageStringVector.begin(), linkageStringVector.end()));
          if(!thisMono->mono_neighbors_[highestIndexLocation].second->is_visited_)
          {
            if(i < linkageStringVector.size() - 1)
            {
              thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ + "]" +
                                        thisOligo->IUPAC_name_;
              thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ + "]" +
                                              thisOligo->author_IUPAC_name_;
            }
            else
            {
              if(local_debug > 0)
              {
                gmml::log(__LINE__, __FILE__,  gmml::INF, thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_);
              }

              thisOligo->IUPAC_name_ = /*thisMono->mono_neighbors_[highestIndexLocation].second->sugar_name_.monosaccharide_short_name_ +*/
                                       thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ +
                                       thisOligo->IUPAC_name_;
              thisOligo->author_IUPAC_name_ = /*thisMono->mono_neighbors_[highestIndexLocation].second->sugar_name_.monosaccharide_short_name_ +*/
                                       thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ +
                                       thisOligo->author_IUPAC_name_;
            }
            thisMono->is_visited_ = true;
            thisOligo->mono_nodes_.push_back(thisMono->mono_neighbors_[highestIndexLocation].second);
            traverseGraph(thisMono->mono_neighbors_[highestIndexLocation].second, thisOligo);
            if(i < linkageStringVector.size() - 1)
            {
              thisOligo->IUPAC_name_ = "[" + thisOligo->IUPAC_name_;
              thisOligo->author_IUPAC_name_ = "[" + thisOligo->author_IUPAC_name_;
            }
            linkageStringVector[highestIndexLocation] = "0"; //small number string so it moves to the next higest linkage
          }
        }
      }
    }
    else //not root with 2 neighbors
    {
      if(local_debug > 0)
      {
        testLog << "This mono is not branched";
        gmml::log(__LINE__, __FILE__, gmml::INF, testLog.str());
        testLog.str(std::string());//clear stringstream
      }

      int neighborNum = 0;
      Glycan::Monosaccharide* thisMonoNeighbor = thisMono->mono_neighbors_[0].second;
      if(thisMonoNeighbor->is_visited_)
      {
        thisMonoNeighbor = thisMono->mono_neighbors_[1].second;
        neighborNum = 1;
      }
      if(!thisMonoNeighbor->is_visited_)
      {
        thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[neighborNum].first->inverse_linkage_type_ +
                                thisMonoName + thisOligo->IUPAC_name_;
        thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[neighborNum].first->inverse_linkage_type_ +
                                thisMonoAuthorName + thisOligo->author_IUPAC_name_;
        thisMono->is_visited_ = true;
        thisOligo->mono_nodes_.push_back(thisMonoNeighbor);
        if(local_debug > 0)
        {
          gmml::log(__LINE__, __FILE__,  gmml::INF, thisOligo->IUPAC_name_);
        }

        traverseGraph(thisMonoNeighbor, thisOligo);
      }
    }
  }
  else //>2 neighbors so branching
  {
    thisOligo->IUPAC_name_ = thisMonoName + thisOligo->IUPAC_name_;
    thisOligo->author_IUPAC_name_ = thisMonoAuthorName + thisOligo->author_IUPAC_name_;
    std::vector<int> branchMaxLengths(thisMono->mono_neighbors_.size(), 1);
    getBranchMaxLengths(thisMono, branchMaxLengths);
    cleanCountedBranches(thisMono);
    int numEqualNeighbors = 0;
    std::vector<int> equalNeighborsLengths;
    for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
    {
      if(local_debug > 0)
      {
        testLog << "Branch " << i << " has a max length of  " << branchMaxLengths[i];
        gmml::log(__LINE__, __FILE__,  gmml::INF, testLog.str());
        testLog.str(std::string());//clear stringstream
      }

      if(std::find(equalNeighborsLengths.begin(), equalNeighborsLengths.end(), branchMaxLengths[i]) == equalNeighborsLengths.end()) //if this length is not already considered
      {
        // equalNeighborsLengths.push_back(branchMaxLengths[i]);
        int branchesOfThisLength = std::count(branchMaxLengths.begin(), branchMaxLengths.end(), branchMaxLengths[i]);
        if(branchMaxLengths[i] == 1 && branchesOfThisLength > 2)
        {
          numEqualNeighbors += branchesOfThisLength;
          equalNeighborsLengths.push_back(branchMaxLengths[i]);
        }
        else if(branchMaxLengths[i] != 1 && branchesOfThisLength > 1)
        {
          numEqualNeighbors += branchesOfThisLength;
          equalNeighborsLengths.push_back(branchMaxLengths[i]);
        }
        // else if(thisMono->is_root_ && branchMaxLengths[i] ==1 && branchesOfThisLength > 1)
        // {
        //   numEqualNeighbors += branchesOfThisLength;
        // }
      }
    }
    if(numEqualNeighbors == 0)
    {
      if(local_debug > 0)
      {
        testLog << "There are no equal branches.";
        gmml::log(__LINE__, __FILE__,  gmml::INF, testLog.str());
        testLog.str(std::string());//clear stringstream
      }
      for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
      {
        int shortestBranchLocation = std::distance(branchMaxLengths.begin(), std::min_element(branchMaxLengths.begin(), branchMaxLengths.end()));
        if(!thisMono->mono_neighbors_[shortestBranchLocation].second->is_visited_)
        {
          if((i < thisMono->mono_neighbors_.size() - 1) && (countedBackwards))
          {
            thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ + "]" +
                                    thisOligo->IUPAC_name_;
            thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ + "]" +
                                            thisOligo->author_IUPAC_name_;
          }
          else if((i < thisMono->mono_neighbors_.size() - 2) && (!countedBackwards))
          {
            thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ + "]" +
                                    thisOligo->IUPAC_name_;
            thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ + "]" +
                                     thisOligo->author_IUPAC_name_;
          }
          else//last branch doesnt get brackets
          {
            thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ +
                                    thisOligo->IUPAC_name_;
            thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ +
                                    thisOligo->author_IUPAC_name_;
          }
          thisMono->is_visited_ = true;
          if(local_debug > 0)
          {
            gmml::log(__LINE__, __FILE__,  gmml::INF, thisOligo->IUPAC_name_);
          }

          thisOligo->mono_nodes_.push_back(thisMono->mono_neighbors_[shortestBranchLocation].second);
          traverseGraph(thisMono->mono_neighbors_[shortestBranchLocation].second, thisOligo);
          if((i < thisMono->mono_neighbors_.size() - 1) && (countedBackwards))
          {
            thisOligo->IUPAC_name_ = "[" + thisOligo->IUPAC_name_;
            thisOligo->author_IUPAC_name_ = "[" + thisOligo->author_IUPAC_name_;
            if(local_debug > 0)
            {
              gmml::log(__LINE__, __FILE__,  gmml::INF, thisOligo->IUPAC_name_);
            }
          }
          else if((i < thisMono->mono_neighbors_.size() - 2) && (!countedBackwards))
          {
            thisOligo->IUPAC_name_ = "[" + thisOligo->IUPAC_name_;
            thisOligo->author_IUPAC_name_ = "[" + thisOligo->author_IUPAC_name_;
            if(local_debug > 0)
            {
              gmml::log(__LINE__, __FILE__,  gmml::INF, thisOligo->IUPAC_name_);
            }
          }
        }
        else
        {
          countedBackwards = true;
        }
        branchMaxLengths[shortestBranchLocation] = 9999; //arbitralily large number so it moves to the next shortest branch
      }
    }
    else
    {
      if(local_debug > 0)
      {
        testLog << "There are equal branches.";
        gmml::log(__LINE__, __FILE__,  gmml::INF, testLog.str());
        testLog.str(std::string());//clear stringstream
      }
      //TODO:  fix the logic of this function; the loop has to be within the if statement
      //handle branch ordering
      std::vector<std::string> linkageStringVector;
      for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
      {
        linkageStringVector.push_back(thisMono->mono_neighbors_[i].first->inverse_linkage_type_);
      }
      std::vector<bool> isNeighborNamed(thisMono->mono_neighbors_.size(), false);
      int shortestBranchLocation = std::distance(branchMaxLengths.begin(), std::min_element(branchMaxLengths.begin(), branchMaxLengths.end()));
      countedBackwards = false;
      for(unsigned int i = 0; i < thisMono->mono_neighbors_.size(); i++)
      {
        shortestBranchLocation = std::distance(branchMaxLengths.begin(), std::min_element(branchMaxLengths.begin(), branchMaxLengths.end()));
        if(isNeighborNamed[shortestBranchLocation])
        {
          branchMaxLengths[shortestBranchLocation] = 9999;
          shortestBranchLocation = std::distance(branchMaxLengths.begin(), std::min_element(branchMaxLengths.begin(), branchMaxLengths.end()));
        }
        if(std::count(equalNeighborsLengths.begin(), equalNeighborsLengths.end(), branchMaxLengths[shortestBranchLocation]) == 0) //no other branches of this length
        {
          if(local_debug > 0)
          {
            testLog << "No other branches of length: " << branchMaxLengths[shortestBranchLocation];
            gmml::log(__LINE__, __FILE__,  gmml::INF, testLog.str());
            testLog.str(std::string());//clear stringstream
          }
          linkageStringVector.push_back("0");
          if(!thisMono->mono_neighbors_[shortestBranchLocation].second->is_visited_)
          {
            if((i < thisMono->mono_neighbors_.size() - 1) && (countedBackwards))
            {
              thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ + "]" +
                                       thisOligo->IUPAC_name_;
              thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ + "]" +
                                       thisOligo->author_IUPAC_name_;
            }
            else if((i < thisMono->mono_neighbors_.size() - 2) && (!countedBackwards))
            {
              thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ + "]" +
                                      thisOligo->IUPAC_name_;
              thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ + "]" +
                                      thisOligo->author_IUPAC_name_;
            }
            else
            {
              thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ +
                                      thisOligo->IUPAC_name_;
              thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[shortestBranchLocation].first->inverse_linkage_type_ +
                                      thisOligo->author_IUPAC_name_;
            }
            thisMono->is_visited_ = true;
            thisOligo->mono_nodes_.push_back(thisMono->mono_neighbors_[shortestBranchLocation].second);
            traverseGraph(thisMono->mono_neighbors_[shortestBranchLocation].second, thisOligo);
            if((i < thisMono->mono_neighbors_.size() - 1) && (countedBackwards))
            {
              thisOligo->IUPAC_name_ = "[" + thisOligo->IUPAC_name_;
              thisOligo->author_IUPAC_name_ = "[" + thisOligo->author_IUPAC_name_;
              if(local_debug > 0)
              {
                gmml::log(__LINE__, __FILE__,  gmml::INF, thisOligo->IUPAC_name_);
              }
            }
            else if((i < thisMono->mono_neighbors_.size() - 2) && (!countedBackwards))
            {
              thisOligo->IUPAC_name_ = "[" + thisOligo->IUPAC_name_;
              thisOligo->author_IUPAC_name_ = "[" + thisOligo->author_IUPAC_name_;
              if(local_debug > 0)
              {
                gmml::log(__LINE__, __FILE__,  gmml::INF, thisOligo->IUPAC_name_);
              }
            }
            isNeighborNamed[shortestBranchLocation] = true;
          }
          else
          {
            countedBackwards = true;
          }
          branchMaxLengths[shortestBranchLocation] = 9999; //arbitralily large number so it moves to the next shortest branch
        }
        else
        {
          if(local_debug > 0)
          {
            testLog << "The equal branches are of length: " << branchMaxLengths[shortestBranchLocation];
            gmml::log(__LINE__, __FILE__,  gmml::INF, testLog.str());
            testLog.str(std::string());//clear stringstream
          }
          std::vector<std::string> equalLinkageStringVector;
          for(unsigned int j = 0; j < thisMono->mono_neighbors_.size(); j++)
          {
            if(branchMaxLengths[j] == branchMaxLengths[shortestBranchLocation])
            {
              equalLinkageStringVector.push_back(linkageStringVector[j]);
              if(local_debug > 0)
              {
                testLog << "This branch has index: " << linkageStringVector[j];
                gmml::log(__LINE__, __FILE__,  gmml::INF, testLog.str());
                testLog.str(std::string());//clear stringstream
              }
            }
            else
            {
              equalLinkageStringVector.push_back("0");
            }
          }
          for(unsigned int j = 0; j < equalLinkageStringVector.size(); j++)
          {

            int highestIndexLocation = std::distance(equalLinkageStringVector.begin(), std::max_element(equalLinkageStringVector.begin(), equalLinkageStringVector.end()));
            if(local_debug > 0)
            {
              testLog << "The highest index is: " << equalLinkageStringVector[highestIndexLocation];
              gmml::log(__LINE__, __FILE__,  gmml::INF, testLog.str());
              testLog.str(std::string());//clear stringstream
            }
            if(!thisMono->mono_neighbors_[highestIndexLocation].second->is_visited_)
            {
              if(thisMono->is_root_ && i < thisMono->mono_neighbors_.size() - 1)
              {
                thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ + "]" +
                                          thisOligo->IUPAC_name_;
                thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ + "]" +
                                                thisOligo->author_IUPAC_name_;
              }
              else if((i < thisMono->mono_neighbors_.size() - 1) && (countedBackwards))
              {
                if(local_debug > 0)
                {
                  gmml::log(__LINE__, __FILE__,  gmml::INF, "Counted Backwards");
                }
                thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ + "]" +
                                          thisOligo->IUPAC_name_;
                thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ + "]" +
                                                thisOligo->author_IUPAC_name_;
              }
              else if((i < thisMono->mono_neighbors_.size() - 2) && (!countedBackwards))
              {
                if(local_debug > 0)
                {
                  gmml::log(__LINE__, __FILE__,  gmml::INF, "Didn't count Backwards");
                }
                thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ + "]" +
                                         thisOligo->IUPAC_name_;
                thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ + "]" +
                                         thisOligo->author_IUPAC_name_;
              }
              else
              {
                thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ +
                                         thisOligo->IUPAC_name_;
                thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ +
                                         thisOligo->author_IUPAC_name_;
              }
              if(local_debug > 0)
              {
                gmml::log(__LINE__, __FILE__,  gmml::INF, thisOligo->IUPAC_name_);
              }

              thisMono->is_visited_ = true;
              thisOligo->mono_nodes_.push_back(thisMono->mono_neighbors_[highestIndexLocation].second);
              traverseGraph(thisMono->mono_neighbors_[highestIndexLocation].second, thisOligo);
              if(thisMono->is_root_ && i < thisMono->mono_neighbors_.size() - 1)
              {
                thisOligo->IUPAC_name_ = "[" + thisOligo->IUPAC_name_;
                thisOligo->author_IUPAC_name_ = "[" + thisOligo->author_IUPAC_name_;
              }
              else if((i < thisMono->mono_neighbors_.size() - 1) && (countedBackwards))
              {
                thisOligo->IUPAC_name_ = "[" + thisOligo->IUPAC_name_;
                thisOligo->author_IUPAC_name_ = "[" + thisOligo->author_IUPAC_name_;
                if(local_debug > 0)
                {
                  gmml::log(__LINE__, __FILE__,  gmml::INF, thisOligo->IUPAC_name_);
                }
              }
              else if((i < thisMono->mono_neighbors_.size() - 2) && (!countedBackwards))
              {
                thisOligo->IUPAC_name_ = "[" + thisOligo->IUPAC_name_;
                thisOligo->author_IUPAC_name_ = "[" + thisOligo->author_IUPAC_name_;
                if(local_debug > 0)
                {
                  gmml::log(__LINE__, __FILE__,  gmml::INF, thisOligo->IUPAC_name_);
                }
              }
            }
            else
            {
              countedBackwards = true;
            }
            equalLinkageStringVector[highestIndexLocation] = "0"; //small number string so it moves to the next lowest linkage
            i++;
            isNeighborNamed[highestIndexLocation] = true;
          }
        }
      //     int equalBranchLength = branchMaxLengths[shortestBranchLocation];
      //     if((branchMaxLengths[i] == equalBranchLength) && (!thisMono->mono_neighbors_[i].second->is_visited_))
      //     {
      //       // gmml::log(__LINE__, __FILE__,  gmml::INF, thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId());
      //       linkageStringVector.push_back(thisMono->mono_neighbors_[i].first->inverse_linkage_type_);
      //     }
      //     else if((branchMaxLengths[i] == equalBranchLength)&& (thisMono->mono_neighbors_[i].second->is_visited_))
      //     {
      //       // gmml::log(__LINE__, __FILE__,  gmml::INF, thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId());
      //       linkageStringVector.push_back("0");
      //     }
      //     else
      //     {
      //       // gmml::log(__LINE__, __FILE__,  gmml::INF, thisMono->mono_neighbors_[i].second->cycle_atoms_[0]->GetResidue()->GetId());
      //       linkageStringVector.push_back("0");
      //     }
      //     branchMaxLengths[shortestBranchLocation] = 9999;
      //   }
      // }
      // //go through vector
      // // gmml::log(__LINE__, __FILE__,  gmml::INF, std::to_string(linkageStringVector.size()));
      // countedBackwards = false;
      // for(unsigned int i = 0; i < linkageStringVector.size(); i++)
      // {
      //   // gmml::log(__LINE__, __FILE__,  gmml::INF, linkageStringVector[i]);
      //   int highestIndexLocation = std::distance(linkageStringVector.begin(), std::max_element(linkageStringVector.begin(), linkageStringVector.end()));
      //   if(!thisMono->mono_neighbors_[highestIndexLocation].second->is_visited_)
      //   {
      //     if((i < thisMono->mono_neighbors_.size() - 1) && (countedBackwards))
      //     {
      //       // gmml::log(__LINE__, __FILE__,  gmml::INF, "Counted Backwards");
      //       thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ + "]" +
      //                                 thisOligo->IUPAC_name_;
      //       thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ + "]" +
      //                                       thisOligo->author_IUPAC_name_;
      //     }
      //     else if((i < thisMono->mono_neighbors_.size() - 2) && (!countedBackwards))
      //     {
      //       // gmml::log(__LINE__, __FILE__,  gmml::INF, "Didn't count Backwards");
      //       thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ + "]" +
      //                                thisOligo->IUPAC_name_;
      //       thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ + "]" +
      //                                thisOligo->author_IUPAC_name_;
      //     }
      //     else
      //     {
      //       thisOligo->IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ +
      //                                thisOligo->IUPAC_name_;
      //       thisOligo->author_IUPAC_name_ = thisMono->mono_neighbors_[highestIndexLocation].first->inverse_linkage_type_ +
      //                                thisOligo->author_IUPAC_name_;
      //     }
      //     // gmml::log(__LINE__, __FILE__,  gmml::INF, thisOligo->IUPAC_name_);
      //     thisMono->is_visited_ = true;
      //     thisOligo->mono_nodes_.push_back(thisMono->mono_neighbors_[highestIndexLocation].second);
      //     traverseGraph(thisMono->mono_neighbors_[highestIndexLocation].second, thisOligo);
      //     if((i < thisMono->mono_neighbors_.size() - 1) && (countedBackwards))
      //     {
      //       thisOligo->IUPAC_name_ = "[" + thisOligo->IUPAC_name_;
      //       thisOligo->author_IUPAC_name_ = "[" + thisOligo->author_IUPAC_name_;
      //       // gmml::log(__LINE__, __FILE__,  gmml::INF, thisOligo->IUPAC_name_);
      //     }
      //     else if((i < thisMono->mono_neighbors_.size() - 2) && (!countedBackwards))
      //     {
      //       thisOligo->IUPAC_name_ = "[" + thisOligo->IUPAC_name_;
      //       thisOligo->author_IUPAC_name_ = "[" + thisOligo->author_IUPAC_name_;
      //       // gmml::log(__LINE__, __FILE__,  gmml::INF, thisOligo->IUPAC_name_);
      //     }
      //   }
      //   else
      //   {
      //     countedBackwards = true;
      //   }
      //   linkageStringVector[highestIndexLocation] = "0"; //small number string so it moves to the next lowest linkage
      }
    }
  }
}

void Glycan::Oligosaccharide::cleanCountedBranches(Glycan::Monosaccharide* this_mono)
{
  this_mono->is_counted_ = false;
  // gmml::log(__LINE__, __FILE__,  gmml::INF, "Cleaning branches count flag");
  for(std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> >::iterator monoNeighbor = this_mono->mono_neighbors_.begin(); monoNeighbor != this_mono->mono_neighbors_.end(); monoNeighbor++)
  {
    std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> monoNeighborPair = *monoNeighbor;
    Glycan::Monosaccharide* this_mono_neighbor = monoNeighborPair.second;
  //   Glycan::Monosaccharide* this_mono_neighbor_holder;

    // this_mono_neighbor->is_counted_ = false;
    if(/*(!this_mono_neighbor->is_visited_) &&*/ (this_mono_neighbor->is_counted_))
    {
      cleanCountedBranches(this_mono_neighbor);
      // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
    }
  }

}

void Glycan::Oligosaccharide::getBranchMaxLengths(Glycan::Monosaccharide* this_mono, std::vector<int> &branchLengths)
{
  int local_debug = -1;
  int branchIndex = -1;
  if(local_debug > 0)
  {
    std::stringstream testLog;
    testLog << this_mono->cycle_atoms_[0]->GetResidue()->GetId() << "Has " << std::to_string(this_mono->mono_neighbors_.size()) << " neighbors.";
    gmml::log(__LINE__, __FILE__,  gmml::INF, testLog.str());
  }
  for(std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> >::iterator monoNeighbor = this_mono->mono_neighbors_.begin(); monoNeighbor != this_mono->mono_neighbors_.end(); monoNeighbor++)
  {
    std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> monoNeighborPair = *monoNeighbor;
    Glycan::Monosaccharide* this_mono_neighbor = monoNeighborPair.second;
    branchIndex++;
    this_mono->is_counted_ = true;
    if(!this_mono_neighbor->is_visited_)
    {
      if(!this_mono_neighbor->is_counted_)
      {
        this_mono_neighbor->is_counted_ = true;
        if(local_debug > 0)
        {
          gmml::log(__LINE__, __FILE__,  gmml::INF, this_mono_neighbor->cycle_atoms_[0]->GetResidue()->GetId());
        }
        if(this_mono_neighbor->mono_neighbors_.size() > 1)
        {
          while(this_mono_neighbor->mono_neighbors_.size() > 1)//At least one more mono on the branch (if only 1, it's the previous mono so stop counting)
          {
            unsigned int countedNeighbors = 0;
            if(local_debug > 0)
            {
              gmml::log(__LINE__, __FILE__,  gmml::INF, std::to_string(this_mono_neighbor->mono_neighbors_.size()));
            }
            if(this_mono_neighbor->mono_neighbors_.size() == 2)//Not branched
            {
              if(local_debug > 0)
              {
                gmml::log(__LINE__, __FILE__,  gmml::INF, "This neighbot is not branched");
              }
              for(std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> >::iterator monoNeighborNeighbor = this_mono_neighbor->mono_neighbors_.begin(); monoNeighborNeighbor != this_mono_neighbor->mono_neighbors_.end(); monoNeighborNeighbor++)
              {
                std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> monoNeighborNeighborPair = *monoNeighborNeighbor;
                Glycan::Monosaccharide* this_neighbor_neighbor = monoNeighborNeighborPair.second;
                if(local_debug > 0)
                {
                  std::stringstream ss;
                  ss << this_mono_neighbor->cycle_atoms_[0]->GetResidue()->GetId() << " " << this_neighbor_neighbor->cycle_atoms_[0]->GetResidue()->GetId();
                  gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
                }
                if(!this_neighbor_neighbor->is_counted_) // if it isn't the previous mono
                {
                  if(local_debug > 0)
                  {
                    std::stringstream ss;
                    ss << this_mono_neighbor->sugar_name_.monosaccharide_short_name_ << " " << this_neighbor_neighbor->sugar_name_.monosaccharide_short_name_;
                    gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
                  }
                  this_neighbor_neighbor->is_counted_ = true;
                  branchLengths[branchIndex] = branchLengths[branchIndex] + 1;
                  this_mono_neighbor = this_neighbor_neighbor;//Move to next mono
                  if(local_debug > 0)
                  {
                    std::stringstream ss;
                    ss << "Branch " << branchIndex << " has " << branchLengths[branchIndex] + 1 << " neighbors";
                    gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
                  }
                  break;
                }
                else
                {
                  // goto exit;
                  countedNeighbors++;
                }
              }
              if(countedNeighbors == this_mono_neighbor->mono_neighbors_.size())
              {
                if(local_debug > 0)
                {
                  gmml::log(__LINE__, __FILE__,  gmml::INF, "About to break");
                }
                break;
              }
            }
            else //Branched neighbor
            {
              std::vector<int> subBranchLengths(this_mono_neighbor->mono_neighbors_.size(), 1); //an array the size of the number of branches, with initial lengths as 1
              if(local_debug > 0)
              {
                gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
              }
              getBranchMaxLengths(this_mono_neighbor, subBranchLengths);
              branchLengths[branchIndex] = branchLengths[branchIndex] + *std::max_element(subBranchLengths.begin(), subBranchLengths.end());
              break;
            }
          }
        }
      }
      else if(!this_mono_neighbor->is_counted_)
      {
        this_mono_neighbor->is_counted_ = true;
      }
    }
  }
}

void Glycan::Oligosaccharide::reindexRGroups(Glycan::Oligosaccharide* this_Oligo)
{
  int numR = 0;
  for(std::vector<Glycan::Monosaccharide*>::reverse_iterator rit = this_Oligo->mono_nodes_.rbegin(); rit != this_Oligo->mono_nodes_.rend(); rit++)
  {
    Glycan::Monosaccharide* thisMono = *rit;
    if(thisMono->on_R_ > 0)//There is an R group
    {
      std::string logstatement = thisMono->sugar_name_.monosaccharide_short_name_ + std::to_string(thisMono->on_R_);
      // gmml::log(__LINE__, __FILE__,  gmml::INF, logstatement);
      for(int i = 0; i < thisMono->on_R_; i++)
      {
        numR++;
        std::string thisR = "R"+std::to_string(i+1);
        // std::size_t found3 = thisMono->sugar_name_.monosaccharide_name_.find(thisR);
        std::size_t found4 = thisMono->sugar_name_.monosaccharide_short_name_.find(thisR);
        if(found4 == std::string::npos)
        {
          numR--;
          continue;
        }
        std::size_t found5 = thisMono->SNFG_name_.find(thisR);
        std::size_t found6 = thisMono->author_SNFG_name_.find(thisR);
        if(numR < 10)
        {
          if(found4 != std::string::npos)
          {
            found4 += 1;
            thisMono->sugar_name_.monosaccharide_short_name_[found4] = std::to_string(numR)[0];
          }
          if(found5 != std::string::npos)
          {
            found5 += 1;
            thisMono->SNFG_name_[found5] = std::to_string(numR)[0];
          }
          if(found6 != std::string::npos)
          {
            found6 += 1;
            thisMono->author_SNFG_name_[found6] = std::to_string(numR)[0];
          }
        }
        else
        {
          if(found4 != std::string::npos)
          {
            found4 += 1;
            thisMono->sugar_name_.monosaccharide_short_name_.insert(found4 + 1, std::to_string(numR));
            thisMono->sugar_name_.monosaccharide_short_name_.erase(found4, 1);
          }
          if(found5 != std::string::npos)
          {
            found5 += 1;
            thisMono->SNFG_name_.insert(found5 + 1, std::to_string(numR));
            thisMono->SNFG_name_.erase(found5, 1);
          }
          if(found6 != std::string::npos)
          {
            found6 += 1;
            thisMono->author_SNFG_name_.insert(found6 + 1, std::to_string(numR));
            thisMono->author_SNFG_name_.erase(found6, 1);
          }
        }
      }
    }
  }
  std::size_t offset = 0;
  // gmml::log(__LINE__, __FILE__,  gmml::INF, std::to_string(numR));
  for(int i = 0; i < numR; i++)
  {
    std::size_t found = this_Oligo->IUPAC_name_.find("<", offset);
    if(found != std::string::npos)
    {
      found += 2;
      if(i < 10)
        this_Oligo->IUPAC_name_[found] = std::to_string(i+1)[0];
      else
      {
        this_Oligo->IUPAC_name_.insert(found, std::to_string(i));
        this_Oligo->IUPAC_name_.erase(found, 1);
      }
      offset = found;
    }
    else
      break;
  }
}

std::string Glycan::Oligosaccharide::CheckOMETerminal(MolecularModeling::Atom* target, std::vector<MolecularModeling::Atom*> & terminal_atoms)
{
    terminal_atoms = std::vector<MolecularModeling::Atom*> ();
    std::vector<MolecularModeling::Atom*>  atoms_1 = std::vector<MolecularModeling::Atom*> ();
    std::vector<MolecularModeling::Atom*>  atoms_2 = std::vector<MolecularModeling::Atom*> ();
    std::stringstream pattern;
    pattern << "O";
    atoms_1.push_back(target);
    atoms_2.push_back(target);
    MolecularModeling::Atom* C1 = NULL;
    MolecularModeling::Atom* C2 = NULL;
    std::vector<MolecularModeling::Atom*>  o_neighbors = target->GetNode()->GetNodeNeighbors();
    for(std::vector<MolecularModeling::Atom*> ::iterator it = o_neighbors.begin(); it != o_neighbors.end(); it++)
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
        std::vector<MolecularModeling::Atom*>  c1_neighbors = C1->GetNode()->GetNodeNeighbors();
        for(std::vector<MolecularModeling::Atom*> ::iterator it = c1_neighbors.begin(); it != c1_neighbors.end(); it++)
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
        std::vector<MolecularModeling::Atom*>  c2_neighbors = C2->GetNode()->GetNodeNeighbors();
        for(std::vector<MolecularModeling::Atom*> ::iterator it = c2_neighbors.begin(); it != c2_neighbors.end(); it++)
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

std::string Glycan::Oligosaccharide::CheckROHTerminal(MolecularModeling::Atom* target, std::vector<MolecularModeling::Atom*> & terminal_atoms)
{
    terminal_atoms = std::vector<MolecularModeling::Atom*> ();
    std::vector<MolecularModeling::Atom*>  o_neighbors = target->GetNode()->GetNodeNeighbors();
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

std::string Glycan::Oligosaccharide::CheckTBTTerminal(MolecularModeling::Atom *target, std::vector<MolecularModeling::Atom*> & terminal_atoms)
{
    terminal_atoms = std::vector<MolecularModeling::Atom*> ();
    std::vector<MolecularModeling::Atom*>  atoms_1 = std::vector<MolecularModeling::Atom*> ();
    std::vector<MolecularModeling::Atom*>  atoms_2 = std::vector<MolecularModeling::Atom*> ();
    std::stringstream pattern;
    pattern << "O";
    atoms_1.push_back(target);
    atoms_2.push_back(target);
    MolecularModeling::Atom* C1 = NULL;
    MolecularModeling::Atom* C2 = NULL;
    std::vector<MolecularModeling::Atom*>  o_neighbors = target->GetNode()->GetNodeNeighbors();
    for(std::vector<MolecularModeling::Atom*> ::iterator it = o_neighbors.begin(); it != o_neighbors.end(); it++)
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
        std::vector<MolecularModeling::Atom*>  c1_neighbors = C1->GetNode()->GetNodeNeighbors();
        for(std::vector<MolecularModeling::Atom*> ::iterator it = c1_neighbors.begin(); it != c1_neighbors.end(); it++)
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
            std::vector<MolecularModeling::Atom*>  c1c1_neighbors = C1C1->GetNode()->GetNodeNeighbors();
            for(std::vector<MolecularModeling::Atom*> ::iterator it = c1c1_neighbors.begin(); it != c1c1_neighbors.end(); it++)
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
            std::vector<MolecularModeling::Atom*>  c1c2_neighbors = C1C2->GetNode()->GetNodeNeighbors();
            for(std::vector<MolecularModeling::Atom*> ::iterator it = c1c2_neighbors.begin(); it != c1c2_neighbors.end(); it++)
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
            std::vector<MolecularModeling::Atom*>  c1c3_neighbors = C1C3->GetNode()->GetNodeNeighbors();
            for(std::vector<MolecularModeling::Atom*> ::iterator it = c1c3_neighbors.begin(); it != c1c3_neighbors.end(); it++)
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
        std::vector<MolecularModeling::Atom*>  c2_neighbors = C2->GetNode()->GetNodeNeighbors();
        for(std::vector<MolecularModeling::Atom*> ::iterator it = c2_neighbors.begin(); it != c2_neighbors.end(); it++)
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
            std::vector<MolecularModeling::Atom*>  c2c1_neighbors = C2C1->GetNode()->GetNodeNeighbors();
            for(std::vector<MolecularModeling::Atom*> ::iterator it = c2c1_neighbors.begin(); it != c2c1_neighbors.end(); it++)
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
            std::vector<MolecularModeling::Atom*>  c2c2_neighbors = C2C2->GetNode()->GetNodeNeighbors();
            for(std::vector<MolecularModeling::Atom*> ::iterator it = c2c2_neighbors.begin(); it != c2c2_neighbors.end(); it++)
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
            std::vector<MolecularModeling::Atom*>  c2c3_neighbors = C2C3->GetNode()->GetNodeNeighbors();
            for(std::vector<MolecularModeling::Atom*> ::iterator it = c2c3_neighbors.begin(); it != c2c3_neighbors.end(); it++)
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

std::string Glycan::Oligosaccharide::CheckTerminals(MolecularModeling::Atom* target, std::vector<MolecularModeling::Atom*>& terminal_atoms)
{
  if(target !=NULL)
  {
    std::vector<MolecularModeling::Atom*> o_neighbors = target->GetNode()->GetNodeNeighbors();
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
        MolecularModeling::Atom* target_o_neighbor = NULL;
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
            std::vector<MolecularModeling::Residue*> residues = this->assembly_->GetAllResiduesOfAssembly();
            MolecularModeling::Residue* target_residue = NULL;
            for(std::vector<MolecularModeling::Residue*>::iterator it = residues.begin(); it != residues.end(); it++)
            {
                MolecularModeling::Residue* residue = *it;
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
                // std::cout << "This return." << std::endl;
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

void Glycan::Oligosaccharide::CheckLinkageNote(Glycan::Monosaccharide* mono1, Glycan::Monosaccharide* mono2, std::string linkage, std::vector<std::string>& checked_linkages)
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
                    oligo_notes_.push_back(linkage_note);
        			}
			}
		}
	}
    }
}

void Glycan::Oligosaccharide::BuildOligosaccharideTreeStructure(Glycan::Monosaccharide *key, std::vector<Glycan::Monosaccharide*> values, Glycan::Oligosaccharide *oligo,
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
                    std::cout << "HELP US!" << key_mono_linkages.at( i ) << std::endl;
                }*/
                std::string link = key_mono_linkages.at(it_index);
                std::stringstream reverse_link;
                reverse_link << gmml::Split(link, "-").at(2) << "-" << gmml::Split(link, "-").at(1) << "-" << gmml::Split(link, "-").at(0);
                if(find(visited_linkages.begin(), visited_linkages.end(), link) == visited_linkages.end() &&
                        find(visited_linkages.begin(), visited_linkages.end(), reverse_link.str()) == visited_linkages.end())
                {
                    //std::cout << "key id " << key->mono_id_  << ", value id " << value_mono->mono_id_ << std::endl;
                    Glycan::Oligosaccharide* child_oligo = new Glycan::Oligosaccharide();
                    // CalculateOligosaccharideBFactor(child_oligo, values);
                    std::vector<Glycan::Monosaccharide*> value_mono_values = monos_table[value_mono];
                    visited_linkages.push_back(link);
                    //std::cout << "call " << value_mono->mono_id_ << std::endl;
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

void Glycan::Oligosaccharide::CalculateOligosaccharideBFactor(Glycan::Oligosaccharide* oligo, std::vector<Glycan::Monosaccharide*> monos)
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
