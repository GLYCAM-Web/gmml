#include "../../../includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataContainer;
using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector;

//////////////////////////////////////////////////////////
//                      QUERY FUNCTIONS                 //
//////////////////////////////////////////////////////////
// Pass in the two atoms on either side the residue-residue linkage


// I think this can go away:
//bool OverWriteEntry(gmml::MolecularMetadata::GLYCAM::DihedralAngleData entry1, gmml::MolecularMetadata::GLYCAM::DihedralAngleData entry2)
//{
//   return (entry1.number_of_bonds_from_anomeric_carbon_ == entry2.number_of_bonds_from_anomeric_carbon_
//            && entry1.index_ == entry2.index_ );
//}

DihedralAngleDataVector DihedralAngleDataContainer::GetEntriesForLinkage( MolecularModeling::Atom* linking_atom1, MolecularModeling::Atom* linking_atom2)
{
    DihedralAngleDataVector matChing_entries;
    Glycam06NamesToTypesLookupContainer metadata_residueNamesToTypes;
    // Go through each entry in the metadata
    for (const auto& entry : dihedralAngleDataVector_)
    {
        // Create a regex of each entry's linking_atom1_ and 2_. These are regex queries.
        //std::cout << "Compare entry " << entry.linking_atom1_ << "-" << entry.linking_atom2_ << " : " << linking_atom1->GetName() << "-" << linking_atom2->GetName() <<"\n";
        std::regex regex1(entry.linking_atom1_, std::regex_constants::ECMAScript);
        std::regex regex2(entry.linking_atom2_, std::regex_constants::ECMAScript);
        // If metadata entry matches (regex query) to the two linking atom names
        if ( (std::regex_search(linking_atom1->GetName(), regex1)) && (std::regex_search(linking_atom2->GetName(), regex2)) )
        {
            // Some entries have conditions for the residue, that they have certain tags. Make sure any conditions are met:
            std::vector<std::string> residue1_types = metadata_residueNamesToTypes.GetTypesForResidue(linking_atom1->GetResidue()->GetName());
            std::vector<std::string> residue2_types = metadata_residueNamesToTypes.GetTypesForResidue(linking_atom2->GetResidue()->GetName());
            if ( (checkIfResidueConditionsAreSatisfied(residue1_types, entry.residue1_conditions_))
                 && (checkIfResidueConditionsAreSatisfied(residue2_types, entry.residue2_conditions_)) )
            {
                //Always add a later entry, but remove earlier match if number_of_bonds_from_anomeric_carbon_ AND index number are the same.
                // I've overloaded the == and != operators in the DihedralAngleData struct to evaluate those.
                // This next line removes any elements of matChing_entries that match "entry", then the line after adds entry.
                matChing_entries.erase(std::remove(matChing_entries.begin(), matChing_entries.end(), entry), matChing_entries.end());
                matChing_entries.push_back(entry);
            }
        }
    }
    return matChing_entries;
}
//////////////////////////////////////////////////////////
//                    PRIVATE FUNCTIONS                 //
//////////////////////////////////////////////////////////
// Some entries have conditions for the first or second residue to have a particular type (aka tag).
// Most entries have "none" for condition. This checks first if condition is "none", and therefore satisfied.
// Otherwise (else if) it checks if any of the residue_types match the condition for the entry, e.g. gauche_effect=galacto.
bool DihedralAngleDataContainer::checkIfResidueConditionsAreSatisfied(std::vector<std::string> residue_types, std::vector<std::string> entry_conditions)
{
    for (const auto& entry_condition : entry_conditions)
    {   // If no condition, return true. If can't find the condition in the list return false, otherwise, having found the condition(s), return true.
       // std::cout << "Condition: " << entry_condition << std::endl;
        if (entry_condition.compare("none")==0)
        {
        //    std::cout << "Returning true as conditions are none" << std::endl;
            return true;
        }
        else if (!(std::find(residue_types.begin(), residue_types.end(), entry_condition) != residue_types.end()))
        {
        //    std::cout << "Returning false as did not find the condition in residue tags" << std::endl;
            return false; //If any condition isn't satisified. return false.
        }
    }
  //  std::cout << "Found all conditions in dihedralangledata.hpp::checkIfResidueConditionsAreSatisfied" << std::endl;
    return true;
}

//////////////////////////////////////////////////////////
//                    INITIALIZER                       //
//////////////////////////////////////////////////////////

// Struct is copied here for reference.
//struct DihedralAngleData
//{
//    std::string linking_atom1_ ;
//    std::string linking_atom2_ ;
//    std::string dihedral_angle_name_ ;
//    double default_angle_value_ ;
//    double lower_deviation_ ;
//    double upper_deviation_ ;
//    double weight_;
//    std::string rotamer_type_ ; // permutation or conformer
//    std::string rotamer_name_ ;
//    int number_of_bonds_from_anomeric_carbon_;
//    int index_ ; // Used to indicate whether multiple entries are meant to overwrite each other or generate an additional angle
//    StringVector residue1_conditions_ ;
//    StringVector residue2_conditions_ ;
//    std::string atom1_ ;
//    std::string atom2_ ;
//    std::string atom3_ ;
//    std::string atom4_ ;
//} ;

/*
 * A note on index_.
 * number_of_bonds_from_anomeric_carbon_. So Phi is 1, Psi is 2, Omg is 3.
 * The index refers the rotamer number. If there are two Phi rotamers, they will have an index of 1 and 2.
 * Chi angle index numbering varies depending on side chain length, so in ASN the Chi1 is 4 bonds away from the sugar, so would be 4
 * If two entries have matChing regex and have the same index_ number (e.g. 1), the first will be overwritten.
 * If two entries have matChing regex and different index_ numbers (e.g. 1,2,3) they will all be used to create multiple rotamers/conformers
 * If two entries have different regex, but apply to the same dihedral angle (e.g. Phi), give them the same index_ number (e.g. 1). Later will overwrite earlier?
 */

DihedralAngleDataContainer::DihedralAngleDataContainer()
{   // const AmberAtomTypeInfo Glycam06j1AtomTypes[] =
    dihedralAngleDataVector_ =
    { //Regex1, Regex2   , Name   , Angle  , Upper  , Lower  , Weight, Entry Type    , Tag  , D , R , Residue 1 type , Residue 2 type,           , Atom names                                                               // Atom names this applies to
        { "C1"   , "O[1-9]" , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "t"  , 1 , 1 , {"aldose"}     , {"none"}                  , "C2" , "C1" , "O." , "C."  }, // Phi should be C2-C1(ano)-Ox-Cx, or C1-C2(ano)-Ox-Cx
        { "C2"   , "O[1-9]" , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "t"  , 1 , 1 , {"none"}       , {"none"}                  , "C3" , "C2" , "O." , "C."  }, // Phi should be C2-C1(ano)-Ox-Cx, or C1-C2(ano)-Ox-Cx
        { "C2"   , "O[1-9]" , "Phi"  , -60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "-g" , 1 , 2 , {"ulosonate", "alpha"}  , {"none"}         , "C3" , "C2" , "O." , "C."  },

        { "C."   , "O[1-5]" , "Psi"  ,   0.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 2 , 1 , {"none"}       , {"none"}                  , "C." , "O." , "C." , "H."  }, // Psi should be C(ano)-Ox-Cx-Hx, if Cx is ring, otherwise, C(ano)-Ox-Cx-C(x-1)
        { "C."   , "O[6-9]" , "Psi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "t"  , 2 , 1 , {"none"}       , {"none"}                  , "C." , "O." , "C." , "C."  },

        { "C."   , "O6"     , "Omega", -60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "gg" , 3 , 1 , {"none"}       , {"none"}                  , "O6" , "C6" , "C5" , "O5"  }, // omg is O6-C5-C5-O5
        { "C."   , "O6"     , "Omega",  60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "gt" , 3 , 2 , {"none"}       , {"none"}                  , "O6" , "C6" , "C5" , "O5"  },
        { "C."   , "O6"     , "Omega", 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "tg" , 3 , 3 , {"none"}       , {"gauche-effect=galacto"} , "O6" , "C6" , "C5" , "O5"  },
      //  Need to check weights when adding dihedral. All possible can return low weight structures.
//      { "C."   , "O6"     , "Omega", 180.0  ,  20.0  ,  20.0  , 0.001 , "permutation" , "tg" , 3 , 3 , {"none"}       , {"gauche-effect=gluco"}   , "O6" , "C6" , "C5" , "O5"  },

         // 2-8 linkages
      //{ "C2", "O8"   , "omg7" ,    };
      // Common sugar derivatives
      // Phosphate/sulfate
        { "[SP]1", "O[1-9]" , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 1 , 1 , {"none"}       , {"none"}                  , "O." , ".1" , "O." , "C."  },
        { "[SP]1", "O[1-5]" , "Psi"  ,   0.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 2 , 1 , {"none"}       , {"none"}                  , ".1" , "O." , "C." , "H."  },
        { "[SP]1", "O[6-9]" , "Psi"  , -60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 2 , 1 , {"none"}       , {"none"}                  , ".1" , "O." , "C." , "C."  },
        { "[SP]1", "O6"     , "Omega", -60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "gg" , 3 , 1 , {"none"}       , {"none"}                  , "O6" , "C6" , "C5" , "O5"  },
        { "[SP]1", "O6"     , "Omega",  60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "gt" , 3 , 2 , {"none"}       , {"none"}                  , "O6" , "C6" , "C5" , "O5"  },
        { "[SP]1", "O6"     , "Omega", 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "tg" , 3 , 3 , {"none"}       , {"gauche-effect=galacto"} , "O6" , "C6" , "C5" , "O5"  },
      // Ac
        { "C1A"  , "O[1-9]" , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "t"  , 1 , 1 , {"none"}       , {"none"}                  , "C2A", "C1A", "O." , "C."  },
        { "C1A"  , "O[1-5]" , "Psi"  ,   0.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "c"  , 2 , 1 , {"none"}       , {"none"}                  , "C1A", "O." , "C." , "H."  },
        { "C1A"  , "O[6-9]" , "Psi"  , -60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , "-g" , 2 , 1 , {"none"}       , {"none"}                  , "C1A", "O." , "C." , "C."  },
      // Me
        { "CH3"  , "O[1-5]" , "Psi"  ,   0.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 2 , 1 , {"none"}       , {"none"}                  , "CH3", "O." , "C." , "H."  },
        { "CH3"  , "O[6-9]" , "Psi"  , -60.0  ,  20.0  ,  20.0  , 1.0   , "permutation" , ""   , 2 , 1 , {"none"}       , {"none"}                  , "CH3", "O." , "C." , "C."  },


        // Protein linkages
        // ASN // Values are from Petrescu et al 2004.
        { "C."   , "ND2"    , "Chi1" , 191.6  ,  14.4  ,  14.4  , 0.497 , "conformer"   , "A"  , 4 , 1 , {"none"}       , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
        { "C."   , "ND2"    , "Chi2" , 177.6  ,  43.0  ,  43.0  , 0.497 , "conformer"   , "A"  , 3 , 1 , {"none"}       , {"amino-acid"}            , "ND2", "CG" , "CB" , "CA"  },
        { "C."   , "ND2"    , "Psi"  , 177.3  ,  12.3  ,  12.3  , 0.497 , "conformer"   , "A"  , 2 , 1 , {"none"}       , {"amino-acid"}            , "C." , "ND2", "CG" , "CB"  },
        { "C1"   , "ND2"    , "Phi"  , 261.0  ,  21.3  ,  21.3  , 0.497 , "conformer"   , "A"  , 1 , 1 , {"none"}       , {"amino-acid"}            , "C." , "C." , "ND2", "CG"  },

        { "C."   , "ND2"    , "Chi1" ,  63.6  ,   8.9  ,   8.9  , 0.178 , "conformer"   , "B"  , 4 , 2 , {"none"}       , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
        { "C."   , "ND2"    , "Chi2" , 191.1  ,  31.6  ,  31.6  , 0.178 , "conformer"   , "B"  , 3 , 2 , {"none"}       , {"amino-acid"}            , "ND2", "CG" , "CB" , "CA"  },
        { "C."   , "ND2"    , "Psi"  , 178.5  ,  13.9  ,  13.9  , 0.178 , "conformer"   , "B"  , 2 , 2 , {"none"}       , {"amino-acid"}            , "C." , "ND2", "CG" , "CB"  },
        { "C1"   , "ND2"    , "Phi"  , 253.7  ,  21.5  ,  21.5  , 0.178 , "conformer"   , "B"  , 1 , 2 , {"none"}       , {"amino-acid"}            , "C." , "C." , "ND2", "CG"  },

        { "C."   , "ND2"    , "Chi1" , 290.6  ,  12.7  ,  12.7  , 0.235 , "conformer"   , "C"  , 4 , 3 , {"none"}       , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
        { "C."   , "ND2"    , "Chi2" , 152.9  ,  23.9  ,  23.9  , 0.235 , "conformer"   , "C"  , 3 , 3 , {"none"}       , {"amino-acid"}            , "ND2", "CG" , "CB" , "CA"  },
        { "C."   , "ND2"    , "Psi"  , 173.1  ,  12.2  ,  12.2  , 0.235 , "conformer"   , "C"  , 2 , 3 , {"none"}       , {"amino-acid"}            , "C." , "ND2", "CG" , "CB"  },
        { "C1"   , "ND2"    , "Phi"  , 268.0  ,  20.3  ,  20.3  , 0.235 , "conformer"   , "C"  , 1 , 3 , {"none"}       , {"amino-acid"}            , "C." , "C." , "ND2", "CG"  },

        { "C."   , "ND2"    , "Chi1" , 302.3  ,  11.5  ,  11.5  , 0.090 , "conformer"   , "D"  , 4 , 4 , {"none"}       , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
        { "C."   , "ND2"    , "Chi2" , 255.0  ,  28.8  ,  28.8  , 0.090 , "conformer"   , "D"  , 3 , 4 , {"none"}       , {"amino-acid"}            , "ND2", "CG" , "CB" , "CA"  },
        { "C."   , "ND2"    , "Psi"  , 178.1  ,  11.5  ,  11.5  , 0.090 , "conformer"   , "D"  , 2 , 4 , {"none"}       , {"amino-acid"}            , "C." , "ND2", "CG" , "CB"  },
        { "C1"   , "ND2"    , "Phi"  , 267.5  ,  23.9  ,  23.9  , 0.090 , "conformer"   , "D"  , 1 , 4 , {"none"}       , {"amino-acid"}            , "C." , "C." , "ND2", "CG"  },


        // THR // Values are from Lovell et al "PENULTIMATE ROTAMER LIBRARY"
        { "C."   , "OG1"    , "Chi1" ,  59.0  ,  10.0  ,  10.0  , 0.49  , "permutation" , "g"  , 3 , 1 , {"none"}       , {"amino-acid"}            , "OG1", "CB" , "CA" , "N"   },
        { "C."   , "OG1"    , "Chi1" ,-171.0  ,   6.0  ,   6.0  , 0.07  , "permutation" , "t"  , 3 , 2 , {"none"}       , {"amino-acid"}            , "OG1", "CB" , "CA" , "N"   },
        { "C."   , "OG1"    , "Chi1" , -61.0  ,   7.0  ,   7.0  , 0.43  , "permutation" , "-g" , 3 , 3 , {"none"}       , {"amino-acid"}            , "OG1", "CB" , "CA" , "N"   },

        { "C."   , "OG1"    , "Psi"  , -60.0  ,  60.0  ,  60.0  , 1.000 , "permutation" , "-g" , 2 , 1 , {"none"}       , {"amino-acid"}            , "C." , "OG1", "CB" , "CA"  },
        { "C."   , "OG1"    , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "t"  , 1 , 1 , {"none"}       , {"amino-acid"}            , "C." , "C." , "OG1", "CB"  },

         // SER // Values not checked
        { "C."   , "OG"     , "Chi1" ,  64.0  ,  10.0  ,  10.0  , 0.48  , "permutation" , "g"  , 3 , 1 , {"none"}       , {"amino-acid"}            , "OG" , "CB" , "CA" , "N"   },
        { "C."   , "OG"     , "Chi1" , 178.0  ,  11.0  ,  11.0  , 0.22  , "permutation" , "t"  , 3 , 2 , {"none"}       , {"amino-acid"}            , "OG" , "CB" , "CA" , "N"   },
        { "C."   , "OG"     , "Chi1" , -65.0  ,   9.0  ,   9.0  , 0.29  , "permutation" , "-g" , 3 , 3 , {"none"}       , {"amino-acid"}            , "OG" , "CB" , "CA" , "N"   },

        { "C."   , "OG"     , "Psi"  , -60.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "-g" , 2 , 1 , {"none"}       , {"amino-acid"}            , "C." , "OG" , "CB" , "CA"  },
        { "C."   , "OG"     , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "t"  , 1 , 1 , {"none"}       , {"amino-acid"}            , "C." , "C." , "OG1", "CB"  },
         // TYR // Values not checked
        { "C."   , "OH"     , "Chi1" , -60.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "-g" , 7 , 1 , {"none"}       , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
        { "C."   , "OH"     , "Chi1" ,  60.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "g"  , 7 , 2 , {"none"}       , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
        { "C."   , "OH"     , "Chi1" , 180.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "t"  , 7 , 3 , {"none"}       , {"amino-acid"}            , "CG" , "CB" , "CA" , "N"   },
        { "C."   , "OH"     , "Chi2" , -60.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "-g" , 6 , 1 , {"none"}       , {"amino-acid"}            , "CD1", "CG" , "CB" , "CA"  },
        { "C."   , "OH"     , "Psi"  , -60.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "-g" , 2 , 1 , {"none"}       , {"amino-acid"}            , "C." , "OH" , "CZ" , "CE1" },
        { "C."   , "OH"     , "Phi"  , 180.0  ,  20.0  ,  20.0  , 1.000 , "permutation" , "t"  , 1 , 1 , {"none"}       , {"amino-acid"}            , "C." , "C." , "OH ", "CZ"  },

    };
}

//    Statistical analysis of the protein environment of N-glycosylation sites: implications for occupancy, structure, and folding
//    Andrei-J. Petrescu  Adina-L. Milac  Stefana M. Petrescu  Raymond A. Dwek Mark R. Wormald
//    Glycobiology, Volume 14, Issue 2, 1 February 2004, Pages 103–114,
