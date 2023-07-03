#ifndef GLYCOSIDICLINKAGE_HPP
#define GLYCOSIDICLINKAGE_HPP

#include <vector>
#include <string>
#include <sstream>
#include <list>

#include "../MolecularModeling/assembly.hpp"
#include "../utils.hpp"
#include "../MolecularModeling/residue.hpp"

namespace Glycan
{
    class Monosaccharide;
    class Oligosaccharide;

    class GlycosidicLinkage
    {
        /*TODO Follow C++ Standards
            Things like Getters/Setters, private variables, const correctness, etc.
        */
        /*TODO Fix initialization of Glycosidic Linkage
            I have a strong feeling this makes 2 glycosidic
            bond objects per linkage.
            It is also possible that this can be called by the wrong
            monosaccharide or the (non_)reducing_mono_ logic is
            wonky.
        */
      public:
        Oligosaccharide* parent_oligo_       = NULL;
        Monosaccharide* reducing_mono_       = NULL;
        Monosaccharide* non_reducing_mono_   = NULL;
        Monosaccharide* non_reducing_mono_2_ = NULL; // In case of anomeric-anomeric linkage
        bool anomeric_anomeric_linkage_      = false;
        std::string residue_linkage_name_;     // ie "GAL_F_2_?_?_1 (1-3) A2G_F_1_?_?_1"
        std::string linkage_name_;             // ie "DManpa1-4DGlcpNAcb"
        std::string linkage_type_;             // ie "1-4"
        std::string inverse_linkage_type_;     // ie "4-1" (needed for writing as you travel in reverse)
        std::string anomeric_configuration_;   // can be "a"(α) or "b"(β)
        std::string anomeric_configuration_2_; // For Anomeric-Anomeric

        // hydroxyl_configuration_ is for the hydroxyl on the reducing_mono_
        // TODO figure out how to handle this with ring flipping
        std::string hydroxyl_configuration_; // can be "axial" or "equatorial"

        // TODO index linkages for better error reporting on the frontend
        int index_in_oligo_ = -1;

        /*TODO Add documentation
          What is the expected range of angles? 0-360? Etc.
          also maybe define gg/gt/tg?
          This will make accessors for radians and whatever other functions
          that will use these easier to write
        */
        double phi_angle_       = -9999;
        double phi_prime_angle_ = -9999;
        double psi_angle_       = -9999;
        double omega_angle_     = -9999;

        double phi_CHI_Energy_   = -1;
        double psi_CHI_Energy_   = -1;
        double omega_CHI_Energy_ = -1;
        /*No phi_prime_CHI_Energy_
          We have not investigated the CHI Energy landscape
          of anomeric-anomeric linkages, so there is no
          anomeric-anomeric CHI Energy function.
        */
        double total_CHI_Energy_ = -1;
        int phi_CHI_function_    = -1; /*!<The number of the CHI energy functuion used for the phi angle; Needed for
                                          properly selecting energy plots>*/
        int psi_CHI_function_    = -1; /*!<The number of the CHI energy functuion used for the psi angle; Needed for
                                          properly selecting energy plots>*/
        int omega_CHI_function_  = -1; /*!<The number of the CHI energy functuion used for the omega angle; Needed for
                                          properly selecting energy plots>*/

        /*TODO Update Atom Handling
            This needs to include the start and end nodes (anomeric C and Cx)
            and all bridge nodes in the linkage, assuming it isn't
              mono - something crazy - another mono.
            I'm thinking a vector of atom pointers?
              [0]       Anomeric C
              [1...n-1] Bridge atoms
              [n]       Cx
            So a typical linkage would be:
              [0] Anomeric C
              [1] Glycosidic O (Ox)
              [2] Cx
            And if it were anomeric-anomeric, Cx would be the anomeric C of
            the other monosaccharide

            I also think there should be a vector of atoms needed for torsion angle
            calculations, and potentially torsion angle manipulation in the future
            DM 1/22
        */
        MolecularModeling::Atom* reducing_mono_carbon_       = NULL;
        MolecularModeling::Atom* non_reducing_mono_carbon_   = NULL;
        MolecularModeling::Atom* non_reducing_mono_2_carbon_ = NULL;
        MolecularModeling::Atom* glycosidic_oxygen_          = NULL;

        // CONSTRUCTOR
        GlycosidicLinkage(Monosaccharide* sourceMono, Monosaccharide* targetMono, std::string source_carbon_ID,
                          std::string target_carbon_ID);

        // DESTRUCTOR
        /*TODO Add Destructor*/

        // FUNCTIONS
        double CalculatePhiAngle(std::vector<MolecularModeling::Atom*> linkage_Atoms);
        double CalculatePsiAngle(std::vector<MolecularModeling::Atom*> linkage_Atoms);
        double CalculateOmegaAngle(std::vector<MolecularModeling::Atom*> linkage_Atoms);
        double CalculatePhiChiEnergy();
        double CalculatePsiChiEnergy();
        int pickPsiChiEnergyFunction();
        double CalculateOmegaChiEnergy();

        /** \brief Gets the configuration of the glycosidic Oxygen

            This function gets the orientation of the reducing sugar's oxygen
            (or other substituent) at the carbon involved in the linkage.
            IE in DManpa1-3DManpb1-OH, it would be the oxygen attached to C3 of
            DManpb.
            For x-6 linkages, it returns the orientation of the O4 hydroxyl, as
            that is what affects the linkage energy
        **/
        std::string determineLinkageConfiguration();
    };
} // namespace Glycan

#endif // GLYCOSIDICLINKAGE_HPP
