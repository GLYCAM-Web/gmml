#ifndef OLIGOSACCHARIDE_HPP
#define OLIGOSACCHARIDE_HPP

#include <vector>
#include <string>
#include <sstream>

#include "../MolecularModeling/assembly.hpp"
#include "../utils.hpp"
#include "../MolecularModeling/residue.hpp"
#include "monosaccharide.hpp"
#include "glycosidiclinkage.hpp"
#include "note.hpp"

namespace Glycan
{
  class Monosaccharide;
  class GlycosidicLinkage;
  class Oligosaccharide
  {
  public:
    Monosaccharide* root_;                              /*!< The root sacchride of the oligosacchride >*/
    std::vector<Oligosaccharide*> child_oligos_;        /*!< The oligosacchrides that are attached to the current oligosaccharide >*/
    std::vector<std::string> child_oligos_linkages_;    /*!< The linkages between the current oligosaccharide and oligosacchrides that are attached to it. i.e. C4_13_4GB_?_2_?_?_1-O4_23_4GB_?_2_?_?_1-C1_24_0MA_?_3_?_?_1 >*/
    std::string oligosaccharide_name_;                  /*!< The complete name sequence of the oligosacchride >*/
    std::string IUPAC_name_;
    std::string author_IUPAC_name_;
    std::vector<std::string> sorted_monos;              /*!< The monosaccharides sorted in the order they are printed in the sequence >*/
    std::string oligosaccharide_linkages_;              /*!< The complete sequence of the linkages of the oligosacchride >*/
    std::string oligosaccharide_residue_linkages_;      /*!< The complete sequence of the residue linkages of the oligosacchride >*/
    std::string terminal_;                              /*!< The terminal residue name of the oligosacchride >*/
    std::string oligosaccharide_terminal_;              /*!< The terminal residue name of the oligosacchride with linkage >*/
    float oligosaccharide_b_factor_;                     /*!< The b_factor for the oligosaccharide >*/
    std::vector<Glycan::Monosaccharide*> mono_nodes_;
    int number_branches_;
    MolecularModeling::Assembly* assembly_;
    std::vector<Note*> oligo_notes_;
    Oligosaccharide();
    Oligosaccharide(MolecularModeling::Assembly* assembly);
    Oligosaccharide(std::vector<Glycan::Monosaccharide*> monos, gmml::ResidueNameMap dataset_residue_names, MolecularModeling::Assembly* assembly);


    //FUNCTIONS
    /*! \fn
              * A function to print out the oligosacchride name and the linkages between its monosacchrides
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
    void Print(std::ostream& out = std::cout);
    /*! \fn
              * A function to update residueLinkStream based on current oligosaccharide
              * @param oligo_temp The current oligosaccharide of oligo-sequence
              * @param residue_of_linkage The PDB-ID of current residue
              */
    std::string updateResidueLink(std::string oligo_temp, std::string residue_of_linkage);
    /*! \fn
              * A recursive function in order to generate the name sequence and the linkages between sacchrides of an oligosacchride
              * This functions starts from the root sacchride and generates the name sequence and linkages in reverse order
              * @param oligosaccharide_name The name sequence of the oligosacchride that has been generated by the function so far
              * @param oligosaccharide_linkages The linkages between sacchrides of the oligosacchride which has been generated by the function so far
              * @param i The counter of linkages between the sacchrides of the oligosacchride
              * @param main_root_id The identifier of the main root in the oligosacchride
              * @param is_cycle The boolean variable which represnets whether the sacchrides involved in the oligosacchride forms a cycle
              */
    void GenerateNameLinkage(std::string& oligosaccharide_name, std::string& oligosaccharide_linkages, int& i, int main_root_id, bool& is_cycle);

    void createOligosaccharideGraphs(std::vector<Glycan::Monosaccharide*> detected_monos, gmml::ResidueNameMap dataset_residue_names,
                                      int& number_of_covalent_links, int& number_of_probable_non_covalent_complexes);
    std::vector<Glycan::Oligosaccharide*> createOligosaccharides(std::vector<Glycan::Monosaccharide*> detected_monos);
    void traverseGraph(Glycan::Monosaccharide* thisMono, Glycan::Oligosaccharide* thisOligo);
    void getBranchMaxLengths(Glycan::Monosaccharide* this_mono, std::vector<int> &branchLengths);
    void cleanCountedBranches(Glycan::Monosaccharide* this_mono);
    std::string CheckTerminals(MolecularModeling::Atom* target, std::vector<MolecularModeling::Atom*>& terminal_atoms);
    void CheckLinkageNote(Glycan::Monosaccharide* mono1, Glycan::Monosaccharide* mono2, std::string linkage, std::vector<std::string>& checked_linkages);
    void BuildOligosaccharideTreeStructure(Glycan::Monosaccharide *key, std::vector<Glycan::Monosaccharide*> values, Glycan::Oligosaccharide *oligo,
                                          std::vector<int>& visited_monos, std::map<Glycan::Monosaccharide*, std::vector<Glycan::Monosaccharide*> > monos_table,
                                          std::map<Glycan::Monosaccharide*, std::vector<std::string> > monos_table_linkages, std::vector<std::string>& visited_linkages);
    void CalculateOligosaccharideBFactor(Glycan::Oligosaccharide* oligo, std::vector<Glycan::Monosaccharide*> monos);
    std::string CheckOMETerminal(MolecularModeling::Atom* target, std::vector<MolecularModeling::Atom*> & terminal_atoms);
    std::string CheckROHTerminal(MolecularModeling::Atom* target, std::vector<MolecularModeling::Atom*> & terminal_atoms);
    std::string CheckTBTTerminal(MolecularModeling::Atom* target, std::vector<MolecularModeling::Atom*> & terminal_atoms);

  };
}

#endif // OLIGOSACCHARIDE_HPP
