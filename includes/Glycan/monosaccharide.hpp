#ifndef MONOSACCHARIDE_HPP
#define MONOSACCHARIDE_HPP

#include <string>
#include <vector>
#include <map>

#include "../MolecularModeling/atom.hpp"
#include "../MolecularModeling/atomnode.hpp"
#include "../MolecularModeling/assembly.hpp"
#include "chemicalcode.hpp"
#include "sugarname.hpp"
#include "note.hpp"

namespace MolecularModeling{class Assembly;}

namespace Glycan
{
  using MolecularModeling::Assembly;
  class GlycosidicLinkage;
  class Oligosaccharide;
  class Monosaccharide
  {
    public:
      int mono_id_;                                                        /*!< The unique identifier of a monosacchride >*/
      std::vector<std::vector<MolecularModeling::Atom*> > side_atoms_;    /*!< The list of side atoms of the ring of the monosacchride >*/
      std::vector<MolecularModeling::Atom*> cycle_atoms_;                 /*!< The list of ring atoms of the ring of the monosacchride >*/
      ChemicalCode* chemical_code_;                                       /*!< The chemical code structure of the monosacchride (Glycode: http://glycam.org/docs/gmml/2016/03/31/glycode-internal-monosaccharide-representation/)>*/
      ChemicalCode* author_chemical_code_;                                        /*!< The chemical code structure of the monosacchride based on the author's residue naming (Glycode: http://glycam.org/docs/gmml/2016/03/31/glycode-internal-monosaccharide-representation/)>*/
      SugarName sugar_name_;                                              /*!< The sugar name object assigned to the monosacchride >*/
      SugarName author_sugar_name_;                                       /*!< The sugar name object assigned to the residue name given by the author >*/
      std::vector<std::pair<std::string, std::string> > derivatives_map_; /*!< A mapping between the monosacchride's atom position/index and the derivative/modification attached to it >*/
      std::vector<std::pair<std::string, std::string> > unknown_derivatives_; /*!< Added so the unknown derivates don't mess up naming based on derivatives_map_ >*/
      std::string cycle_atoms_str_;                                       /*!< The string version of atom identifiers of the ring of the monosacchride >*/
      std::string anomeric_status_;                                       /*!< The detection status of the anomeric carbon >*/
      MolecularModeling::Atom* anomeric_carbon_pointer_ = NULL;
      std::string bfmp_ring_conformation_;                                /*!< The ring conformation of the monosaccharide, currently detected by an external program (BFMP) >*/
      float b_factor_;                                                    /*!< The Average B Factor of the monosaccharide >*/
      std::string SNFG_name_;
      std::string author_SNFG_name_;
      std::vector<Note*> mono_notes_;
      std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> > mono_neighbors_;
      bool is_visited_;
      bool is_root_;
      bool is_counted_;
      bool is_indexed_;
      int oligosaccharide_index_;
      int IUPAC_index_;
      int oligo_branch_index_;
      Glycan::Oligosaccharide* oligo_parent_ = NULL;
      std::string residue_name_;
      int on_R_ = 0;
      MolecularModeling::Assembly* assembly_;
      /*! \fn
        * Default constructor
        */
      Monosaccharide();
      Monosaccharide(std::string* cycle_atoms_str, std::vector<MolecularModeling::Atom*>& cycle_atoms, MolecularModeling::Assembly* this_assembly, std::string CCD_Path);
      Monosaccharide(const Monosaccharide &mono);

      //////////////////////////////////////////////////////////
      //                       FUNCTIONS                      //
      //////////////////////////////////////////////////////////

      void createSNFGname();
      void createAuthorSNFGname();
      MolecularModeling::Atom* FindAnomericCarbon( Glycan::Note* anomeric_note, std::vector< std::string > & anomeric_carbons_status, std::vector<MolecularModeling::Atom*> cycle, std::string cycle_atoms_str ) ;
      std::vector<std::string> GetSideGroupOrientations(MolecularModeling::Assembly* this_assembly);
      void InitiateDetectionOfCompleteSideGroupAtoms ();
      bool CheckIfPlusOneSideAtomBelongsToCurrentMonosaccharide (std::vector<MolecularModeling::Atom*>& SideAtomArm, std::vector<MolecularModeling::Atom*>& cycle_atoms, MolecularModeling::Atom* working_atom );
      void SetCompleteSideGroupAtoms(std::vector<MolecularModeling::Atom*>& SideAtomArm, MolecularModeling::Atom* working_atom, std::vector<MolecularModeling::Atom*>& cycle_atoms, std::vector<MolecularModeling::Atom*>& visited_atoms);
      void CheckIfSideChainIsTerminal(MolecularModeling::Atom* starting_atom, std::vector<MolecularModeling::Atom*> & cycle_and_visited_atoms, bool & is_terminal);
      void ExtractDerivatives(MolecularModeling::Assembly* this_assembly);
      std::string GetFormula(MolecularModeling::Atom* target);
      void CountElements(MolecularModeling::Atom* thisAtom, std::vector<std::pair<std::string, int> >& elementVector);
      std::vector<MolecularModeling::Atom*> ExtractAdditionalSideAtoms();
      void GenerateCompleteName(std::vector<MolecularModeling::Atom*> &plus_sides, Glycan::Monosaccharide* this_mono, MolecularModeling::Assembly* this_assembly);
      void GenerateCompleteSugarName(MolecularModeling::Assembly* this_assembly);
      void UpdatePdbCode();
      void UpdateComplexSugarChemicalCode() ;
      void CheckMonoNaming(std::string original_residue, std::string original_residue_id);
      // void addDeoxyToName(Glycan::SugarName base_name, Glycan::ChemicalCode chemical_code, std::vector<int> deoxy_locations);
      Glycan::ChemicalCode* BuildChemicalCode(std::vector<std::string> orientations);


  } ;
}

#endif // MONOSACCHARIDE_HPP
