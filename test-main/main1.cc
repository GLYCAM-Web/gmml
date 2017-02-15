#include <iostream>
#include <sstream>
#include <string>
#include "includes/gmml.hpp"

#define PdbPreprocessingTest

#ifdef ParameterFileTest
using namespace ParameterFileSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        ParameterFile* temp = new ParameterFile(argv[1], 1);
        temp->Print(std::cout);
        std::stringstream name;
        name << argv[1] << "_MOD";
        temp->Write(name.str());
    }
    return 0;
}
#endif

#ifdef PrepFileTest
using namespace PrepFileSpace;
int main(int argc, char *argv[])
{    
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        PrepFile* temp = new PrepFile(argv[1]);

        // Print out in a customized format
                temp->Print(std::cout);

        // Get all residue names in the given prep file
//        std::vector<std::string> residues = temp->GetAllResidueNames();
//        for(std::vector<std::string>::iterator it = residues.begin(); it != residues.end(); it++)
//        {
//            std::cout << (*it) << std::endl;
//        }
    }
    return 0;
}
#endif

#ifdef LibFileTest
using namespace LibraryFileSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        LibraryFile* temp = new LibraryFile(argv[1]);
        // Print out in a customized format
        //        temp->Print(std::cout);

        // Get all residue names in the given lib file
        std::vector<std::string> residues = temp->GetAllResidueNames();
        for(std::vector<std::string>::iterator it = residues.begin(); it != residues.end(); it++)
        {
            std::cout << (*it) << std::endl;
        }
    }
    return 0;
}
#endif

#ifdef CrdFileTest
using namespace CoordinateFileSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        CoordinateFile* temp = new CoordinateFile(argv[1]);
        temp->Print(std::cout);
    }
    return 0;
}
#endif

#ifdef PdbFileTest

using namespace PdbFileSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        PdbFile* temp = new PdbFile(argv[1]);
        // Print out in a customized format
//        temp->Print(std::cout);

//        Get all residue names of the given file
//                PdbFileSpace::PdbFile::PdbPairVectorAtomNamePositionFlag residues = temp->GetAllResidueNames();
//        for(PdbFileSpace::PdbFile::PdbPairVectorAtomNamePositionFlag::iterator it = residues.begin(); it != residues.end(); it++)
//        {
//            std::cout << (*it).first << std::endl;
//        }

        // Get all residues of a given file
//        PdbFile::PdbResidueVector residues = temp->GetAllResidues();
//        for(PdbFile::PdbResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
//        {
//            PdbFileSpace::PdbResidue* residue = (*it);
//            residue->Print(std::cout);
//        }

        // Remove a reside from pdb file
//        temp->DeleteResidue(residues.at(0));

        // Remove an atom from pdb file
//        std::string line = "ATOM      9  N   LYS A   2     136.231  18.765  -7.391  1.00 63.59           N  \n";
//        PdbAtom* atom = new PdbAtom(line);
//        temp->DeleteAtom(atom);

        // Split Atom Card in a Model Card
//        temp->SplitAtomCardOfModelCard('A', 25);

        // Insert residue
//        PdbAtomCard* new_atom_card = new PdbAtomCard();
//        new_atom_card->SetRecordName("ATOM");
//        PdbAtomCard::PdbAtomMap new_atom_map;
//        PdbAtomCard::PdbAtomMap atom_card = (*(((((*(temp->GetModels()->GetModels()).begin()).second)->GetModelResidueSet())->GetAtoms()).begin()))->GetAtoms();
//        int serial_number = 1;
//        for(PdbAtomCard::PdbAtomMap::iterator it = atom_card.begin(); it != atom_card.end(); it++)
//        {
//            PdbAtom* atom = (*it).second;
//            if(atom->GetAtomResidueSequenceNumber() == 1 && atom->GetAtomChainId() == 'A')
//            {
//                new_atom_map[serial_number] = new PdbAtom(serial_number,atom->GetAtomName(), ' ', atom->GetAtomResidueName(), 'A', 34,
//                                                          ' ', atom->GetAtomOrthogonalCoordinate(), atom->GetAtomOccupancy(), atom->GetAtomTempretureFactor(),
//                                                          atom->GetAtomElementSymbol(), atom->GetAtomCharge());
//                serial_number++;
//            }
//        }
//        new_atom_card->SetAtoms(new_atom_map);
//        temp->InsertResidueAfter(new_atom_card);

        // Writer
        std::string name = temp->GetPath();
        int index = name.find_last_of('.',name.length()-1);
        name = name.substr(0,index);
        std::stringstream ss;
        ss << name << "_updated.pdb";
        temp->Write(ss.str());

        // Get all atoms of a residue in the pdb file
//        PdbFile::PdbAtomVector atoms = temp->GetAllAtomsOfResidue(residues.at(0));
//        for(PdbFile::PdbAtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
//        {
//            PdbFileSpace::PdbAtom* atom = (*it);
//            atom->Print(std::cout);
//        }

        // Get atom by its name
//        PdbAtom* atom = temp->GetAtomOfResidueByName(residues.at(0), "CA");
//        atom->Print(std::cout);

//        PdbAtom* atom1 = temp->GetAtomOfResidueByName(residues.at(0), "CB");
//        std::cout << atom->GetAtomOrthogonalCoordinate().Distance(atom1->GetAtomOrthogonalCoordinate()) << std::endl;
    }
    return 0;
}
#endif

#ifdef PdbPreprocessorTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;
int main(int argc, char *argv[])
{
    if(argc < 3)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        PdbFile* temp = new PdbFile(argv[1]);
        PdbFile* temp1 = new PdbFile(argv[2]);

        // Get recognized and unrecognized residue names
        PdbFileSpace::PdbFile::PdbPairVectorAtomNamePositionFlag pdb_residue_names = temp->GetAllResidueNames();
        std::vector<std::string> dataset_residue_names = temp1->GetAllResidueNames();

        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        std::vector<std::string> recognized_residue_names = pdb_preprocessor->GetRecognizedResidueNames(pdb_residue_names, dataset_residue_names);
        std::vector<std::string> unrecognized_residue_names = pdb_preprocessor->GetUnrecognizedResidueNames(pdb_residue_names, dataset_residue_names);

        std::cout << "********************* Recognized Residue Names ******************" << std::endl;
        for(std::vector<std::string>::iterator it = recognized_residue_names.begin(); it != recognized_residue_names.end(); it++)
        {
            std::cout << *it << std::endl;
        }

        std::cout << "********************* Unrecognized Residue Names ******************" << std::endl;
        for(std::vector<std::string>::iterator it = unrecognized_residue_names.begin(); it != unrecognized_residue_names.end(); it++)
        {
            std::cout << *it << std::endl;
        }

        // Get recognized and unrecognized residues
        std::vector<PdbResidue*> pdb_residues = temp->GetAllResidues();

        std::vector<PdbResidue*> recognized_residues = pdb_preprocessor->GetRecognizedResidues(pdb_residues, recognized_residue_names);
        std::vector<PdbResidue*> unrecognized_residues = pdb_preprocessor->GetUnrecognizedResidues(pdb_residues, unrecognized_residue_names);

        std::cout << "********************* Recognized Residues ******************" << std::endl;
        for(std::vector<PdbResidue*>::iterator it = recognized_residues.begin(); it != recognized_residues.end(); it++)
        {
            (*it)->Print(std::cout);
        }

        std::cout << "********************* Unrecognized Residues ******************" << std::endl;
        for(std::vector<PdbResidue*>::iterator it = unrecognized_residues.begin(); it != unrecognized_residues.end(); it++)
        {
            (*it)->Print(std::cout);
        }
    }
    return 0;
}
#endif

#ifdef PdbUnrecognizedResidueTest
#include "includes/utils.hpp"
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;

int main(int argc, char *argv[])
{
    if(argc < 5)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        std::string pdb_file_path = argv[1];
        int number_of_lib_files = gmml::ConvertString<int>(argv[2]);
        int number_of_prep_files = gmml::ConvertString<int>(argv[3]);

        std::vector<std::string> lib_files_paths;
        for(int i = 0; i < number_of_lib_files; i++)
        {
            lib_files_paths.push_back(argv[4+i]);
        }
        std::vector<std::string> prep_files_paths;
        for(int i = 0; i < number_of_prep_files; i++)
        {
            prep_files_paths.push_back(argv[4+number_of_lib_files+i]);
        }

        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractUnrecognizedResidues(pdb_file_path, lib_files_paths, prep_files_paths);

        PdbPreprocessor::PdbPreprocessorUnrecognizedResidueVector unrecognized_residues = pdb_preprocessor->GetUnrecognizedResidues();
        for(PdbPreprocessor::PdbPreprocessorUnrecognizedResidueVector::iterator it = unrecognized_residues.begin(); it != unrecognized_residues.end(); it++)
        {
            (*it)->Print(std::cout);
        }
    }
    return 0;
}
#endif

#ifdef PdbRecognizedResidueTest
#include "includes/utils.hpp"
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;

int main(int argc, char *argv[])
{
    if(argc < 5)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        std::string pdb_file_path = argv[1];
        int number_of_lib_files = gmml::ConvertString<int>(argv[2]);
        int number_of_prep_files = gmml::ConvertString<int>(argv[3]);

        std::vector<std::string> lib_files_paths;
        for(int i = 0; i < number_of_lib_files; i++)
        {
            lib_files_paths.push_back(argv[4+i]);
        }
        std::vector<std::string> prep_files_paths;
        for(int i = 0; i < number_of_prep_files; i++)
        {
            prep_files_paths.push_back(argv[4+number_of_lib_files+i]);
        }

        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractRecognizedResidues(pdb_file_path, lib_files_paths, prep_files_paths);

        PdbPreprocessor::PdbPreprocessorRecognizedResidueVector recognized_residues = pdb_preprocessor->GetRecognizedResidues();
        for(PdbPreprocessor::PdbPreprocessorRecognizedResidueVector::iterator it = recognized_residues.begin(); it != recognized_residues.end(); it++)
        {
            (*it)->Print(std::cout);
        }
    }
    return 0;
}
#endif

#ifdef PdbCYSResidueTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractCYSResidues(argv[1]);
        PdbPreprocessor::PdbPreprocessorDisulfideBondVector cys_residues = pdb_preprocessor->GetDisulfideBonds();
        for(PdbPreprocessor::PdbPreprocessorDisulfideBondVector::iterator it = cys_residues.begin(); it != cys_residues.end(); it++)
        {
            (*it)->Print(std::cout);
        }
    }
    return 0;
}
#endif

#ifdef PdbHISResidueTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractHISResidues(argv[1]);
        PdbPreprocessor::PdbPreprocessorHistidineMappingVector his_residues = pdb_preprocessor->GetHistidineMappings();
        for(PdbPreprocessor::PdbPreprocessorHistidineMappingVector::iterator it = his_residues.begin(); it != his_residues.end(); it++)
        {
            (*it)->Print(std::cout);
        }

    }
    return 0;
}
#endif

#ifdef PdbUnknownHeavyAtomsTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;
int main(int argc, char *argv[])
{
    if(argc < 5)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        std::string pdb_file_path = argv[1];
        int number_of_lib_files = gmml::ConvertString<int>(argv[2]);
        int number_of_prep_files = gmml::ConvertString<int>(argv[3]);

        std::vector<std::string> lib_files_paths;
        for(int i = 0; i < number_of_lib_files; i++)
        {
            lib_files_paths.push_back(argv[4+i]);
        }
        std::vector<std::string> prep_files_paths;
        for(int i = 0; i < number_of_prep_files; i++)
        {
            prep_files_paths.push_back(argv[4+number_of_lib_files+i]);
        }

        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractUnknownHeavyAtoms(pdb_file_path, lib_files_paths, prep_files_paths);

        PdbPreprocessor::PdbPreprocessorUnrecognizedHeavyAtomVector unknown_heavy_atoms = pdb_preprocessor->GetUnrecognizedHeavyAtoms();
        for(PdbPreprocessor::PdbPreprocessorUnrecognizedHeavyAtomVector::iterator it = unknown_heavy_atoms.begin(); it != unknown_heavy_atoms.end(); it++)
        {
            (*it)->Print(std::cout);
        }
    }
    return 0;
}
#endif

#ifdef PdbRemovedHydrogensTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;
int main(int argc, char *argv[])
{
    if(argc < 5)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        std::string pdb_file_path = argv[1];
        int number_of_lib_files = gmml::ConvertString<int>(argv[2]);
        int number_of_prep_files = gmml::ConvertString<int>(argv[3]);

        std::vector<std::string> lib_files_paths;
        for(int i = 0; i < number_of_lib_files; i++)
        {
            lib_files_paths.push_back(argv[4+i]);
        }
        std::vector<std::string> prep_files_paths;
        for(int i = 0; i < number_of_prep_files; i++)
        {
            prep_files_paths.push_back(argv[4+number_of_lib_files+i]);
        }

        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractRemovedHydrogens(pdb_file_path, lib_files_paths, prep_files_paths);

        PdbPreprocessor::PdbPreprocessorReplacedHydrogenVector removed_hydrogen_atoms = pdb_preprocessor->GetReplacedHydrogens();
        for(PdbPreprocessor::PdbPreprocessorReplacedHydrogenVector::iterator it = removed_hydrogen_atoms.begin(); it != removed_hydrogen_atoms.end(); it++)
        {
            (*it)->Print(std::cout);
        }
    }
    return 0;
}
#endif

#ifdef PdbChainTerminationTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractAminoAcidChains(argv[1]);
        PdbPreprocessor::PdbPreprocessorChainTerminationVector chains = pdb_preprocessor->GetChainTerminations();
        for(PdbPreprocessor::PdbPreprocessorChainTerminationVector::iterator it = chains.begin(); it != chains.end(); it++)
        {
            (*it)->Print(std::cout);
        }

    }
    return 0;
}
#endif

#ifdef PdbGapTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractGapsInAminoAcidChains(argv[1]);
        PdbPreprocessor::PdbPreprocessorMissingResidueVector gaps = pdb_preprocessor->GetMissingResidues();
        for(PdbPreprocessor::PdbPreprocessorMissingResidueVector::iterator it = gaps.begin(); it != gaps.end(); it++)
        {
            (*it)->Print(std::cout);
        }

    }
    return 0;
}
#endif

#ifdef PdbAlternateResidueTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractAlternateResidue(argv[1]);
        PdbPreprocessor::PdbPreprocessorAlternateResidueMap alternate_residues = pdb_preprocessor->GetAlternateResidueMap();
        for(PdbPreprocessor::PdbPreprocessorAlternateResidueMap::iterator it = alternate_residues.begin(); it != alternate_residues.end(); it++)
        {
            PdbPreprocessorAlternateResidue* alternate_residue = (*it).second;
            alternate_residue->Print(std::cout);
        }
    }
    return 0;
}
#endif

#ifdef PdbRemoveUnrecognizedResidueTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;

int main(int argc, char *argv[])
{
    if(argc < 5)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        std::string pdb_file_path = argv[1];
        int number_of_lib_files = gmml::ConvertString<int>(argv[2]);
        int number_of_prep_files = gmml::ConvertString<int>(argv[3]);

        std::vector<std::string> lib_files_paths;
        for(int i = 0; i < number_of_lib_files; i++)
        {
            lib_files_paths.push_back(argv[4+i]);
        }
        std::vector<std::string> prep_files_paths;
        for(int i = 0; i < number_of_prep_files; i++)
        {
            prep_files_paths.push_back(argv[4+number_of_lib_files+i]);
        }

        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractUnrecognizedResidues(pdb_file_path, lib_files_paths, prep_files_paths);

        PdbPreprocessor::PdbPreprocessorUnrecognizedResidueVector unrecognized_residues = pdb_preprocessor->GetUnrecognizedResidues();
        for(PdbPreprocessor::PdbPreprocessorUnrecognizedResidueVector::iterator it = unrecognized_residues.begin(); it != unrecognized_residues.end(); it++)
        {
            (*it)->Print(std::cout);
        }
        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        pdb_preprocessor->RemoveUnrecognizedResidues(pdb_file,pdb_preprocessor->GetUnrecognizedResidues());

        // Write into the file
        std::string name = pdb_file->GetPath();
        int index = name.find_last_of('.',name.length()-1);
        name = name.substr(0,index);
        std::stringstream ss;
        ss << name << "_updated.pdb";
        pdb_file->Write(ss.str());
    }
    return 0;
}
#endif


#ifdef PdbUpdateCYSResidueTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractCYSResidues(argv[1]);
        PdbPreprocessor::PdbPreprocessorDisulfideBondVector cys_residues = pdb_preprocessor->GetDisulfideBonds();
        for(PdbPreprocessor::PdbPreprocessorDisulfideBondVector::iterator it = cys_residues.begin(); it != cys_residues.end(); it++)
        {
            (*it)->Print(std::cout);
        }
        PdbFile* pdb_file = new PdbFile(argv[1]);
//        pdb_preprocessor->GetDisulfideBonds().at(0)->SetIsBonded(false);
        pdb_preprocessor->UpdateCYSResidues(pdb_file, pdb_preprocessor->GetDisulfideBonds());

        // Write into the file
        std::string name = pdb_file->GetPath();
        int index = name.find_last_of('.',name.length()-1);
        name = name.substr(0,index);
        std::stringstream ss;
        ss << name << "_updated.pdb";
        pdb_file->Write(ss.str());
    }
    return 0;
}
#endif

#ifdef PdbUpdateHISResidueTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractHISResidues(argv[1]);
        PdbPreprocessor::PdbPreprocessorHistidineMappingVector his_residues = pdb_preprocessor->GetHistidineMappings();
        for(PdbPreprocessor::PdbPreprocessorHistidineMappingVector::iterator it = his_residues.begin(); it != his_residues.end(); it++)
        {
            (*it)->Print(std::cout);
        }

        PdbFile* pdb_file = new PdbFile(argv[1]);
        pdb_preprocessor->UpdateHISMapping(pdb_file, pdb_preprocessor->GetHistidineMappings());

        // Write into the file
        std::string name = pdb_file->GetPath();
        int index = name.find_last_of('.',name.length()-1);
        name = name.substr(0,index);
        std::stringstream ss;
        ss << name << "_updated.pdb";
        pdb_file->Write(ss.str());
    }
    return 0;
}
#endif

#ifdef PdbRemoveUnknownHeavyAtomsTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;
int main(int argc, char *argv[])
{
    if(argc < 5)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        std::string pdb_file_path = argv[1];
        int number_of_lib_files = gmml::ConvertString<int>(argv[2]);
        int number_of_prep_files = gmml::ConvertString<int>(argv[3]);

        std::vector<std::string> lib_files_paths;
        for(int i = 0; i < number_of_lib_files; i++)
        {
            lib_files_paths.push_back(argv[4+i]);
        }
        std::vector<std::string> prep_files_paths;
        for(int i = 0; i < number_of_prep_files; i++)
        {
            prep_files_paths.push_back(argv[4+number_of_lib_files+i]);
        }

        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractUnknownHeavyAtoms(pdb_file_path, lib_files_paths, prep_files_paths);

        PdbPreprocessor::PdbPreprocessorUnrecognizedHeavyAtomVector unknown_heavy_atoms = pdb_preprocessor->GetUnrecognizedHeavyAtoms();
        for(PdbPreprocessor::PdbPreprocessorUnrecognizedHeavyAtomVector::iterator it = unknown_heavy_atoms.begin(); it != unknown_heavy_atoms.end(); it++)
        {
            (*it)->Print(std::cout);
        }

        PdbFile* pdb_file = new PdbFile(argv[1]);
        pdb_preprocessor->RemoveResiduesOfUnknownHeavyAtoms(pdb_file, pdb_preprocessor->GetUnrecognizedHeavyAtoms());

        // Write into the file
        std::string name = pdb_file->GetPath();
        int index = name.find_last_of('.',name.length()-1);
        name = name.substr(0,index);
        std::stringstream ss;
        ss << name << "_updated.pdb";
        pdb_file->Write(ss.str());
    }
    return 0;
}
#endif

#ifdef PdbRemoveRemovedHydrogensTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;
int main(int argc, char *argv[])
{
    if(argc < 5)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        std::string pdb_file_path = argv[1];
        int number_of_lib_files = gmml::ConvertString<int>(argv[2]);
        int number_of_prep_files = gmml::ConvertString<int>(argv[3]);

        std::vector<std::string> lib_files_paths;
        for(int i = 0; i < number_of_lib_files; i++)
        {
            lib_files_paths.push_back(argv[4+i]);
        }
        std::vector<std::string> prep_files_paths;
        for(int i = 0; i < number_of_prep_files; i++)
        {
            prep_files_paths.push_back(argv[4+number_of_lib_files+i]);
        }

        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractRemovedHydrogens(pdb_file_path, lib_files_paths, prep_files_paths);

        PdbPreprocessor::PdbPreprocessorReplacedHydrogenVector removed_hydrogen_atoms = pdb_preprocessor->GetReplacedHydrogens();
        for(PdbPreprocessor::PdbPreprocessorReplacedHydrogenVector::iterator it = removed_hydrogen_atoms.begin(); it != removed_hydrogen_atoms.end(); it++)
        {
            (*it)->Print(std::cout);
        }

        PdbFile* pdb_file = new PdbFile(argv[1]);
        pdb_preprocessor->RemoveRemovedHydrogens(pdb_file,pdb_preprocessor->GetReplacedHydrogens());

        // Write into the file
        std::string name = pdb_file->GetPath();
        int index = name.find_last_of('.',name.length()-1);
        name = name.substr(0,index);
        std::stringstream ss;
        ss << name << "_updated.pdb";
        pdb_file->Write(ss.str());
    }
    return 0;
}
#endif

#ifdef PdbUpdateChainTerminationTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;

int main(int argc, char *argv[])
{
    if(argc < 3)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        std::string pdb_file_path = argv[1];
        int number_of_lib_files = gmml::ConvertString<int>(argv[2]);

        std::vector<std::string> lib_files_paths;
        for(int i = 0; i < number_of_lib_files; i++)
        {
            lib_files_paths.push_back(argv[3+i]);
        }

        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractAminoAcidChains(pdb_file_path);
        pdb_preprocessor->GetChainTerminations().at(0)->SetSelectedNTermination(gmml::COCH3);
        pdb_preprocessor->GetChainTerminations().at(1)->SetSelectedCTermination(gmml::NH2);
        pdb_preprocessor->GetChainTerminations().at(2)->SetSelectedNTermination(gmml::COCH3);
        pdb_preprocessor->GetChainTerminations().at(2)->SetSelectedCTermination(gmml::NH2);

        PdbPreprocessor::PdbPreprocessorChainTerminationVector chains = pdb_preprocessor->GetChainTerminations();
        for(PdbPreprocessor::PdbPreprocessorChainTerminationVector::iterator it = chains.begin(); it != chains.end(); it++)
        {
            (*it)->Print(std::cout);
        }

        PdbFile* pdb_file = new PdbFile(argv[1]);
        pdb_preprocessor->UpdateAminoAcidChains(pdb_file, lib_files_paths, pdb_preprocessor->GetChainTerminations());

        // Write into the file
        std::string name = pdb_file->GetPath();
        int index = name.find_last_of('.',name.length()-1);
        name = name.substr(0,index);
        std::stringstream ss;
        ss << name << "_updated.pdb";
        pdb_file->Write(ss.str());

    }
    return 0;
}
#endif

#ifdef PdbUpdateGapTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;
int main(int argc, char *argv[])
{
    if(argc < 3)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        std::string pdb_file_path = argv[1];
        int number_of_lib_files = gmml::ConvertString<int>(argv[2]);

        std::vector<std::string> lib_files_paths;
        for(int i = 0; i < number_of_lib_files; i++)
        {
            lib_files_paths.push_back(argv[3+i]);
        }

        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractGapsInAminoAcidChains(pdb_file_path);
        pdb_preprocessor->GetMissingResidues().at(0)->SetSelectedNTermination(gmml::NH3);
        pdb_preprocessor->GetMissingResidues().at(0)->SetSelectedCTermination(gmml::CO2);
        PdbPreprocessor::PdbPreprocessorMissingResidueVector gaps = pdb_preprocessor->GetMissingResidues();
        for(PdbPreprocessor::PdbPreprocessorMissingResidueVector::iterator it = gaps.begin(); it != gaps.end(); it++)
        {
            (*it)->Print(std::cout);
        }

        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        pdb_preprocessor->UpdateGapsInAminoAcidChains(pdb_file, lib_files_paths, pdb_preprocessor->GetMissingResidues());

        // Write into the file
        std::string name = pdb_file->GetPath();
        int index = name.find_last_of('.',name.length()-1);
        name = name.substr(0,index);
        std::stringstream ss;
        ss << name << "_updated.pdb";
        pdb_file->Write(ss.str());
    }
    return 0;
}
#endif

#ifdef PdbRemoveUnselectedAlternateResiduesTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
        pdb_preprocessor->ExtractAlternateResidue(argv[1]);
        PdbPreprocessor::PdbPreprocessorAlternateResidueMap alternate_residues = pdb_preprocessor->GetAlternateResidueMap();
        for(PdbPreprocessor::PdbPreprocessorAlternateResidueMap::iterator it = alternate_residues.begin(); it != alternate_residues.end(); it++)
        {
            PdbPreprocessorAlternateResidue* alternate_residue = (*it).second;
            alternate_residue->Print(std::cout);
        }

        PdbFile* pdb_file = new PdbFile(argv[1]);
        pdb_preprocessor->RemoveUnselectedAlternateResidues(pdb_file, pdb_preprocessor->GetAlternateResidueMap());

        // Write into the file
        std::string name = pdb_file->GetPath();
        int index = name.find_last_of('.',name.length()-1);
        name = name.substr(0,index);
        std::stringstream ss;
        ss << name << "_updated.pdb";
        pdb_file->Write(ss.str());
    }
    return 0;
}
#endif

#ifdef PdbPreprocessingTest
using namespace PdbFileSpace;
using namespace PdbPreprocessorSpace;
int main(int argc, char *argv[])
{
    if(argc < 7)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        std::string pdb_file_path = argv[1];
        int number_of_amino_lib_files = gmml::ConvertString<int>(argv[2]);
        int number_of_glycam_lib_files = gmml::ConvertString<int>(argv[3]);
        int number_of_other_lib_files = gmml::ConvertString<int>(argv[4]);
        int number_of_prep_files = gmml::ConvertString<int>(argv[5]);

        std::vector<std::string> amino_lib_files_paths = std::vector<std::string>();
        for(int i = 0; i < number_of_amino_lib_files; i++)
        {
            amino_lib_files_paths.push_back(argv[6+i]);
        }
        std::vector<std::string> glycam_lib_files_paths = std::vector<std::string>();
        for(int i = 0; i < number_of_glycam_lib_files; i++)
        {
            glycam_lib_files_paths.push_back(argv[6+number_of_amino_lib_files+i]);
        }
        std::vector<std::string> other_lib_files_paths = std::vector<std::string>();
        for(int i = 0; i < number_of_other_lib_files; i++)
        {
            other_lib_files_paths.push_back(argv[6+number_of_amino_lib_files+number_of_glycam_lib_files+i]);
        }
        std::vector<std::string> prep_files_paths = std::vector<std::string>();
        for(int i = 0; i < number_of_prep_files; i++)
        {
            prep_files_paths.push_back(argv[6+number_of_amino_lib_files+number_of_glycam_lib_files+number_of_other_lib_files+i]);
        }

        PdbFile* pdb_file = new PdbFile(pdb_file_path);
        PdbPreprocessor* pdb_preprocessor = new PdbPreprocessor();
//        pdb_preprocessor->ExtractGapsInAminoAcidChains(pdb_file, lib_files_paths);
        pdb_preprocessor->Preprocess(pdb_file, amino_lib_files_paths, glycam_lib_files_paths, other_lib_files_paths, prep_files_paths);

        /////////////////////////////////////////////////////////////////////////
        //                              CODE HERE                              //
        /////////////////////////////////////////////////////////////////////////
        // Apply user changes to pdb_preprocessor attributes to be reflected in the output files

//        pdb_preprocessor->Print(std::cout);
//        pdb_preprocessor->GetChainTerminations().at(0)->SetSelectedNTermination(gmml::COCH3);
//        pdb_preprocessor->GetChainTerminations().at(0)->SetSelectedCTermination(gmml::NH2);
        pdb_preprocessor->ApplyPreprocessing(pdb_file, amino_lib_files_paths, glycam_lib_files_paths, prep_files_paths);
//        pdb_preprocessor->UpdateGapsInAminoAcidChains(pdb_file, lib_files_paths, pdb_preprocessor->GetMissingResidues());

        // Write into the file
        std::string name = pdb_file->GetPath();
        int index = name.find_last_of('.',name.length()-1);
        name = name.substr(0,index);
        std::stringstream ss;
        ss << name << "_updated.pdb";
        pdb_file->WriteWithTheGivenModelNumber(ss.str());
    }
    return 0;
}
#endif

#ifdef TopologyFileTest
#include "includes/FileSet/TopologyFileSpace/topologyfile.hpp"
using namespace TopologyFileSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        TopologyFile* temp = new TopologyFile(argv[1]);
//        temp->Print(std::cout);

        // Write into the file
        std::string name = temp->GetPath();
        int index = name.find_last_of('.',name.length()-1);
        name = name.substr(0,index);
        std::stringstream ss;
        ss << name << "_updated.parm7";
        temp->Write(ss.str());
    }
    return 0;
}
#endif

#ifdef PdbqtFileTest
#include "includes/FileSet/PdbqtFileSpace/pdbqtfile.hpp"
using namespace PdbqtFileSpace;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        PdbqtFile* temp = new PdbqtFile(argv[1]);
        temp->Print(std::cout);

        // Write into the file
//        std::string name = temp->GetPath();
//        int index = name.find_last_of('.',name.length()-1);
//        name = name.substr(0,index);
//        std::stringstream ss;
//        ss << name << "_updated.pdbqt";
//        temp->Write(ss.str());
    }
    return 0;
}
#endif

#ifdef BuildAssemblyTest
using namespace MolecularModeling;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
//        std::vector<std::string> file_path = gmml::Split(argv[2],"#");
//        gmml::InputFileType type = gmml::ConvertString2AssemblyInputFileType(argv[1]);
//        gmml::BuildingStructureOption building_type = gmml::ConvertString2AssemblyBuildingStructureOption(argv[3]);
//        Assembly* temp = new Assembly(file_path, type);


//        std::vector<std::string> options = std::vector<std::string>();
//        options = gmml::Split(argv[4], "#");
//        std::vector<std::string> database_files = std::vector<std::string>();
//        database_files = gmml::Split(argv[5], "#");
//        temp->BuildStructure(building_type, options, database_files);
//        temp->BuildStructure(building_type, options, file_path);
//        temp->Print(std::cout);

//        temp->BuildAssemblyFromLibraryFile(file_path.at(0));
//        temp->BuildStructureByLIBFileInformation();
//        temp->Print(std::cout);
//        LibraryFileSpace::LibraryFile* library_file = temp->BuildLibraryFileStructureFromAssembly();
//        library_file->Write("library.text");

//        CoordinateFileSpace::CoordinateFile* coordinate_file = temp->BuildCoordinateFileStructureFromAssembly();
//        coordinate_file->Write("coordinate.text");

//        Assembly* temp = new Assembly();
//        temp->BuildAssemblyFromTopologyCoordinateFile(argv[1], argv[2], argv[3]);
//        temp->SetSourceFile(argv[1]);
//        temp->BuildStructureByTOPFileInformation();
//        temp->Print(std::cout);
//        TopologyFileSpace::TopologyFile* top_file = temp->BuildTopologyFileStructureFromAssembly(argv[3]);
//        top_file->Write("top_file.parm7");
//        CoordinateFileSpace::CoordinateFile* crd_file = temp->BuildCoordinateFileStructureFromAssembly();
//        crd_file->Write("crd_file.rst7");
//        PdbFileSpace::PdbFile* pdb_file = temp->BuildPdbFileStructureFromAssembly();
//        pdb_file->Write("pdb_file.pdb");
//        PrepFileSpace::PrepFile* prep_file = temp->BuildPrepFileStructureFromAssembly(argv[3]);
//        prep_file->Write("prep_file.prep");

//        Assembly* temp = new Assembly();
//        temp->BuildAssemblyFromPdbFile(argv[1]);
//        temp->BuildStructureByDistance();
//        temp->Print(std::cout);
//        PdbFileSpace::PdbFile* pdb_file = temp->BuildPdbFileStructureFromAssembly();
//        pdb_file->Write("pdb_file.pdb");
//        temp->Select("1.*:MAN,GAL@C1,C2:ALA@O1;1.1.2,1.1.3:NAG@N2");

//        Assembly* temp = new Assembly();
//        temp->BuildAssemblyFromPrepFile(argv[1], argv[2]);
//        temp->SetSourceFile(argv[1]);
//        temp->BuildStructureByPrepFileInformation();
//        temp->Print(std::cout);
//        PrepFileSpace::PrepFile* prep_file = temp->BuildPrepFileStructureFromAssembly(argv[2]);
//        prep_file->Write("prep_file.prep");
//        PdbFileSpace::PdbFile* pdb_file = temp->BuildPdbFileStructureFromAssembly();
//        pdb_file->Write("pdb_file.pdb");
//        TopologyFileSpace::TopologyFile* topology_file = temp->BuildTopologyFileStructureFromAssembly(argv[2]);
//        topology_file->Write("top_file.parm7");
//        CoordinateFileSpace::CoordinateFile* coordinate_file = temp->BuildCoordinateFileStructureFromAssembly();
//        coordinate_file->Write("crd_file.rst7");

//        temp->BuildStructureByDistance();
//        temp->Print(std::cout);
//        PdbFileSpace::PdbFile* pdb_file = temp->BuildPdbFileStructureFromAssembly();
//        pdb_file->Write("pdbqt_file.pdb");

        Assembly* temp = new Assembly();
        CondensedSequenceSpace::CondensedSequence::CondensedSequenceAmberPrepResidueTree prep_residues =
                CondensedSequenceSpace::CondensedSequence::CondensedSequenceAmberPrepResidueTree();
        if(temp->CheckCondensedSequenceSanity(argv[1], prep_residues))
            temp->BuildAssemblyFromCondensedSequence(argv[1], argv[2], argv[3], true);
        temp->SetSourceFile(argv[2]);
//        temp->BuildStructureByPrepFileInformation();
//        temp->BuildStructureByDistance();
//        TopologyFileSpace::TopologyFile* top_file = temp->BuildTopologyFileStructureFromAssembly(argv[3]);
//        top_file->Write("top_file.parm7");
//        CoordinateFileSpace::CoordinateFile* crd_file = temp->BuildCoordinateFileStructureFromAssembly();
//        crd_file->Write("crd_file.rst7");
        PdbFileSpace::PdbFile* pdb_file = temp->BuildPdbFileStructureFromAssembly();
        pdb_file->Write("pdb_file.pdb");
//        TopologyFileSpace::TopologyFile* top_file = temp->BuildTopologyFileStructureFromAssembly(argv[3]);
//        top_file->Write("top_file.parm7");
//        CoordinateFileSpace::CoordinateFile* crd_file = temp->BuildCoordinateFileStructureFromAssembly();
//        crd_file->Write("crd_file.rst7");
//        temp = new Assembly();
//        std::vector<std::string> prep_files = std::vector<std::string>();
//        prep_files.push_back(argv[2]);
//        temp->BuildAssemblyFromPdbFile("pdb_file.pdb", std::vector<std::string>(), std::vector<std::string>(), std::vector<std::string>(),
//                                       prep_files, argv[3]);
//        temp->BuildStructureByDistance();
//        PdbFileSpace::PdbFile* pdb_file_1 = temp->BuildPdbFileStructureFromAssembly();
//        pdb_file_1->Write("pdb_file_1.pdb");
//        PrepFileSpace::PrepFile* prep_file = temp->BuildPrepFileStructureFromAssembly(argv[3]);
//        prep_file->Write("prep_file.prep");
//        temp->Print();

    }
    return 0;
}
#endif

#ifdef BuildAssemblyTest1
using namespace MolecularModeling;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        Assembly* temp = new Assembly();
        CondensedSequenceSpace::CondensedSequence::CondensedSequenceAmberPrepResidueTree prep_residues =
                CondensedSequenceSpace::CondensedSequence::CondensedSequenceAmberPrepResidueTree();
        if(temp->CheckCondensedSequenceSanity(argv[1], prep_residues))
        {
            temp->BuildAssemblyFromCondensedSequence(argv[1], argv[2], argv[3], true);
            temp->SetSourceFile(argv[2]);
            PdbFileSpace::PdbFile* pdb_file = temp->BuildPdbFileStructureFromAssembly();
            pdb_file->Write("pdb_file.pdb");

            CondensedSequenceSpace::CondensedSequence* condensed_sequence = new CondensedSequenceSpace::CondensedSequence(argv[1]);
            CondensedSequenceSpace::CondensedSequence::CondensedSequenceRotamersAndGlycosidicAnglesInfo rotamers_glycosidic_angles_info =
                    condensed_sequence->GetCondensedSequenceRotamersAndGlycosidicAnglesInfo(condensed_sequence->GetCondensedSequenceResidueTree());
            /*for(unsigned int i = 0; i < rotamers_glycosidic_angles_info.size(); i++)
                std::cout << rotamers_glycosidic_angles_info.at(i).first << std::endl;
            std::cout << condensed_sequence->CountAllPossible28LinkagesRotamers(rotamers_glycosidic_angles_info) << std::endl;*/
            CondensedSequenceSpace::CondensedSequence::IndexNameMap names = CondensedSequenceSpace::CondensedSequence::IndexNameMap();
            Assembly::AssemblyVector structures = temp->BuildAllRotamersFromCondensedSequence(argv[1], argv[2],
                    argv[3], rotamers_glycosidic_angles_info, names);

            for(int i = 0; i < structures.size(); i++)
            {
                PdbFileSpace::PdbFile* pdb_file = structures.at(i)->BuildPdbFileStructureFromAssembly();
                pdb_file->Write("pdb_file_" + names[i] + ".pdb");
            }
        }
    }
    return 0;
}
#endif

#ifdef BuildPdbFromSugarIdentifier
using namespace MolecularModeling;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        std::vector<std::string> amino_lib_files = gmml::Split(argv[2], ",");
        Assembly* temp = new Assembly();
        temp->BuildAssemblyFromPdbFile(argv[1]);
        temp->BuildStructureByDistance(3);
        std::vector<Glycan::Oligosaccharide*> oligos = temp->ExtractSugars(amino_lib_files);
        gmml::GlycamResidueNamingMap res_name_map = temp->ExtractResidueGlycamNamingMap(oligos);
        temp->UpdateResidueName2GlycamName(res_name_map, argv[3]);
//	temp->Print();


//        Assembly* assembly = new Assembly();
//        for(std::vector<Glycan::Oligosaccharide*>::iterator it = oligos.begin(); it != oligos.end(); it++)
//        {
//            Glycan::Oligosaccharide* oligo = *it;
//            Assembly* olig = new Assembly();
//            olig->BuildAssemblyFromCondensedSequence(oligo->oligosaccharide_name_, argv[3], argv[4]);
//            PdbFileSpace::PdbFile* pdb_file = olig->BuildPdbFileStructureFromAssembly();
//            pdb_file->Write(oligo->oligosaccharide_name_ + ".pdb");
//            assembly->AddAssembly(olig);

//        }
//        temp->Print();
        PdbFileSpace::PdbFile* pdb_file = temp->BuildPdbFileStructureFromAssembly(1);
        pdb_file->Write("pdb_file.pdb");

    }
}
#endif

#ifdef MultipleAssemblyTest
using namespace MolecularModeling;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        int number_of_files = gmml::ConvertString<int>(argv[1]);
        Assembly* assembly = new Assembly();
        Assembly* assembly1 = new Assembly();
        for(int i = 0; i < number_of_files; i++)
        {
            std::vector<std::string> file_path = gmml::Split(argv[(i+1)*2+1],"#");
            gmml::InputFileType type = gmml::ConvertString2AssemblyInputFileType(argv[(i+1)*2]);
            Assembly* temp = new Assembly(file_path, type);
            assembly1->AddAssembly(temp);
        }
        Assembly* assembly2 = new Assembly();
        for(int i = 0; i < number_of_files; i++)
        {
            std::vector<std::string> file_path = gmml::Split(argv[(i+1)*2+1],"#");
            gmml::InputFileType type = gmml::ConvertString2AssemblyInputFileType(argv[(i+1)*2]);
            Assembly* temp = new Assembly(file_path, type);
            assembly2->AddAssembly(temp);
        }
        assembly->AddAssembly(assembly1);
//        assembly->AddAssembly(assembly2);
        Assembly::AssemblyVector a1_assemblies = assembly1->GetAssemblies();
        for(Assembly::AssemblyVector::iterator it = a1_assemblies.begin(); it != a1_assemblies.end(); it++)
        {
            Assembly* temp = (*it);
            std::cout << temp->GetId() << std::endl;
        }
        std::cout << assembly1->GetId() << std::endl;
        Assembly::AssemblyVector a2_assemblies = assembly2->GetAssemblies();
        for(Assembly::AssemblyVector::iterator it = a2_assemblies.begin(); it != a2_assemblies.end(); it++)
        {
            Assembly* temp = (*it);
            std::cout << temp->GetId() << std::endl;
        }
        std::cout << assembly2->GetId() << std::endl;
        std::cout << assembly->GetId() << std::endl;
        assembly->BuildStructureByDistance();
        assembly->Print(std::cout);
//        assembly->Print(std::cout);
//        Assembly::AtomVector atoms = assembly1->Select("1.*:#520,MAN,GAL@#3740-3750,^C:ALA@O1;1.1.2,1.1.3:NAG@O$");
//        for(Assembly::AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
//        {
//            Atom* atom = *it;
//            std::cout << atom->GetResidue()->GetId() << ":" << atom->GetId() << std::endl;
//        }
//        assembly1->BuildStructureByDistance();
//        assembly1->Print(std::cout);
    }
    return 0;
}
#endif

#ifdef IonizingAssemblyTest
using namespace MolecularModeling;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        Assembly* temp = new Assembly();
//        CoordinateFileSpace::CoordinateFile* coordinate_file = temp->BuildCoordinateFileStructureFromAssembly();
//        coordinate_file->Write("coordinate.text");

//        temp->BuildAssemblyFromTopologyCoordinateFile(argv[1], argv[2], argv[3]);
//        temp->BuildStructureByDistance();
//        temp->SetSourceFile(argv[1]);
//        temp->BuildStructureByTOPFileInformation();
//        temp->Print();
//        PrepFileSpace::PrepFile* prep_file = temp->BuildPrepFileStructureFromAssembly(argv[3]);
//        prep_file->Write("prep_file.prep");
//        PdbFileSpace::PdbFile* pdb_file = temp->BuildPdbFileStructureFromAssembly();
//        pdb_file->Write("pdb_file.pdb");
//        CoordinateFileSpace::CoordinateFile* crd_file = temp->BuildCoordinateFileStructureFromAssembly();
//        crd_file->Write("crd_file.rst7");
//        LibraryFileSpace::LibraryFile* lib_file = temp->BuildLibraryFileStructureFromAssembly();
//        lib_file->Write("lib_file.lib");
//        TopologyFileSpace::TopologyFile* top_file = temp->BuildTopologyFileStructureFromAssembly(argv[3]);
//        top_file->Write("top_file.parm7");
//        temp->Ionizing(argv[4], argv[5], argv[6], gmml::ConvertString<int>(argv[7]));
//        PrepFileSpace::PrepFile* prep_file_ion = temp->BuildPrepFileStructureFromAssembly(argv[3]);
//        prep_file_ion->Write("prep_file_ion.prep");
//        PdbFileSpace::PdbFile* pdb_file_ion = temp->BuildPdbFileStructureFromAssembly();
//        pdb_file_ion->Write("pdb_file_ion.pdb");
//        CoordinateFileSpace::CoordinateFile* crd_file_ion = temp->BuildCoordinateFileStructureFromAssembly();
//        crd_file_ion->Write("crd_file_ion.crd");
//        TopologyFileSpace::TopologyFile* top_file_ion = temp->BuildTopologyFileStructureFromAssembly(argv[3], argv[6]);
//        top_file_ion->Write("top_file_ion.top");
//        LibraryFileSpace::LibraryFile* lib_file_ion = temp->BuildLibraryFileStructureFromAssembly();
//        lib_file_ion->Write("lib_file_ion.lib");


        std::vector<std::string> amino_libs = gmml::Split(argv[2], ",");
        std::vector<std::string> glycam_libs = std::vector<std::string>();
        std::vector<std::string> other_libs = std::vector<std::string>();
        std::vector<std::string> prep = gmml::Split(argv[3], ",");
        temp->BuildAssemblyFromPdbFile(argv[1], amino_libs, glycam_libs, other_libs, prep, argv[4]);
        temp->BuildStructureByDistance();
//        PrepFileSpace::PrepFile* prep_file = temp->BuildPrepFileStructureFromAssembly(argv[4]);
//        prep_file->Write("prep_file.prep");
//        PdbFileSpace::PdbFile* pdb_file = temp->BuildPdbFileStructureFromAssembly();
//        pdb_file->Write("pdb_file.pdb");
        CoordinateFileSpace::CoordinateFile* crd_file = temp->BuildCoordinateFileStructureFromAssembly();
        crd_file->Write("crd_file.rst7");
        TopologyFileSpace::TopologyFile* top_file = temp->BuildTopologyFileStructureFromAssembly(argv[4]);
        top_file->Write("top_file.parm7");
//        LibraryFileSpace::LibraryFile* lib_file = temp->BuildLibraryFileStructureFromAssembly();
//        lib_file->Write("lib_file.lib");
//        temp->Ionizing(argv[5], argv[6], argv[7], gmml::ConvertString<int>(argv[8]));
//        PdbFileSpace::PdbFile* pdb_file_ion = temp->BuildPdbFileStructureFromAssembly();
//        pdb_file_ion->Write("pdb_file_ion.pdb");
//        CoordinateFileSpace::CoordinateFile* crd_file_ion = temp->BuildCoordinateFileStructureFromAssembly();
//        crd_file_ion->Write("crd_file_ion.crd");
//        LibraryFileSpace::LibraryFile* lib_file_ion = temp->BuildLibraryFileStructureFromAssembly();
//        lib_file_ion->Write("lib_file_ion.lib");
//        TopologyFileSpace::TopologyFile* top_file_ion = temp->BuildTopologyFileStructureFromAssembly(argv[4], argv[7]);
//        top_file_ion->Write("top_file_ion.top");
//        PrepFileSpace::PrepFile* prep_file_ion = temp->BuildPrepFileStructureFromAssembly(argv[4]);
//        prep_file_ion->Write("prep_file_ion.prep");

//        temp->BuildAssemblyFromPrepFile(argv[1], argv[2]);
//        temp->SetSourceFile(argv[1]);
//        temp->BuildStructureByPrepFileInformation();
//        temp->BuildStructureByDistance();
//        temp->Print(std::cout);
//        PrepFileSpace::PrepFile* prep_file = temp->BuildPrepFileStructureFromAssembly(argv[2]);
//        prep_file->Write("prep_file.prep");
//        PdbFileSpace::PdbFile* pdb_file = temp->BuildPdbFileStructureFromAssembly();
//        pdb_file->Write("pdb_file.pdb");
//        CoordinateFileSpace::CoordinateFile* crd_file = temp->BuildCoordinateFileStructureFromAssembly();
//        crd_file->Write("crd_file.rst7");
//        TopologyFileSpace::TopologyFile* top_file = temp->BuildTopologyFileStructureFromAssembly(argv[2]);
//        top_file->Write("top_file.parm7");
//        LibraryFileSpace::LibraryFile* lib_file = temp->BuildLibraryFileStructureFromAssembly();
//        lib_file->Write("lib_file.lib");
//        temp->Ionizing(argv[3], argv[4], argv[5], gmml::ConvertString<int>(argv[6]));
//        PrepFileSpace::PrepFile* prep_file_ion = temp->BuildPrepFileStructureFromAssembly(argv[2]);
//        prep_file_ion->Write("prep_file_ion.prep");
//        PdbFileSpace::PdbFile* pdb_file_ion = temp->BuildPdbFileStructureFromAssembly();
//        pdb_file_ion->Write("pdb_file_ion.pdb");
//        CoordinateFileSpace::CoordinateFile* crd_file_ion = temp->BuildCoordinateFileStructureFromAssembly();
//        crd_file_ion->Write("crd_file_ion.crd");
//        TopologyFileSpace::TopologyFile* top_file_ion = temp->BuildTopologyFileStructureFromAssembly(argv[2], argv[5]);
//        top_file_ion->Write("top_file_ion.top");
//        LibraryFileSpace::LibraryFile* lib_file_ion = temp->BuildLibraryFileStructureFromAssembly();
//        lib_file_ion->Write("lib_file_ion.lib");

//        temp->BuildAssemblyFromPdbqtFile(argv[1]);
//        temp->BuildStructureByDistance();
//        temp->Print(std::cout);
//        PrepFileSpace::PrepFile* prep_file = temp->BuildPrepFileStructureFromAssembly(argv[2]);
//        prep_file->Write("prep_file.prep");
//        PdbFileSpace::PdbFile* pdb_file = temp->BuildPdbFileStructureFromAssembly();
//        pdb_file->Write("pdb_file.pdb");
//        CoordinateFileSpace::CoordinateFile* crd_file = temp->BuildCoordinateFileStructureFromAssembly();
//        crd_file->Write("crd_file.crd");
//        TopologyFileSpace::TopologyFile* top_file = temp->BuildTopologyFileStructureFromAssembly(argv[2], argv[5]);
//        top_file->Write("top_file.top");
//        LibraryFileSpace::LibraryFile* lib_file = temp->BuildLibraryFileStructureFromAssembly();
//        lib_file->Write("lib_file.lib");
//        temp->Ionizing(argv[3], argv[4], argv[5], gmml::ConvertString<int>(argv[6]));
//        PrepFileSpace::PrepFile* prep_file_ion = temp->BuildPrepFileStructureFromAssembly(argv[2]);
//        prep_file_ion->Write("prep_file_ion.prep");
//        PdbFileSpace::PdbFile* pdb_file_ion = temp->BuildPdbFileStructureFromAssembly();
//        pdb_file_ion->Write("pdb_file_ion.pdb");
//        CoordinateFileSpace::CoordinateFile* crd_file_ion = temp->BuildCoordinateFileStructureFromAssembly();
//        crd_file_ion->Write("crd_file_ion.crd");
//        TopologyFileSpace::TopologyFile* top_file_ion = temp->BuildTopologyFileStructureFromAssembly(argv[2], argv[5]);
//        top_file_ion->Write("top_file_ion.top");
//        LibraryFileSpace::LibraryFile* lib_file_ion = temp->BuildLibraryFileStructureFromAssembly();
//        lib_file_ion->Write("lib_file_ion.lib");

//        temp->BuildAssemblyFromLibraryFile(argv[1]);
//        temp->SetSourceFile(argv[1]);
//        temp->BuildStructureByLIBFileInformation();
//        temp->Print(std::cout);
//        PrepFileSpace::PrepFile* prep_file = temp->BuildPrepFileStructureFromAssembly(argv[2]);
//        prep_file->Write("prep_file.prep");
//        PdbFileSpace::PdbFile* pdb_file = temp->BuildPdbFileStructureFromAssembly();
//        pdb_file->Write("pdb_file.pdb");
//        CoordinateFileSpace::CoordinateFile* crd_file = temp->BuildCoordinateFileStructureFromAssembly();
//        crd_file->Write("crd_file.rst7");
//        TopologyFileSpace::TopologyFile* top_file = temp->BuildTopologyFileStructureFromAssembly(argv[2]);
//        top_file->Write("top_file.parm7");
//        LibraryFileSpace::LibraryFile* lib_file = temp->BuildLibraryFileStructureFromAssembly();
//        lib_file->Write("lib_file.lib");
//        temp->Ionizing(argv[3], argv[4], argv[5], gmml::ConvertString<int>(argv[6]));
//        PrepFileSpace::PrepFile* prep_file_ion = temp->BuildPrepFileStructureFromAssembly(argv[2]);
//        prep_file_ion->Write("prep_file_ion.prep");
//        PdbFileSpace::PdbFile* pdb_file_ion = temp->BuildPdbFileStructureFromAssembly();
//        pdb_file_ion->Write("pdb_file_ion.pdb");
//        CoordinateFileSpace::CoordinateFile* crd_file_ion = temp->BuildCoordinateFileStructureFromAssembly();
//        crd_file_ion->Write("crd_file_ion.crd");
//        TopologyFileSpace::TopologyFile* top_file_ion = temp->BuildTopologyFileStructureFromAssembly(argv[2], argv[5]);
//        top_file_ion->Write("top_file_ion.top");
//        LibraryFileSpace::LibraryFile* lib_file_ion = temp->BuildLibraryFileStructureFromAssembly();
//        lib_file_ion->Write("lib_file_ion.lib");

    }
    return 0;
}
#endif


#ifdef SolvationTest
using namespace MolecularModeling;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        Assembly* temp = new Assembly();
//        temp->BuildAssemblyFromLibraryFile(argv[1], argv[3]);
//        temp->SetSourceFile(argv[1]);
//        temp->BuildStructureByLIBFileInformation();

        temp->BuildAssemblyFromTopologyCoordinateFile(argv[1], argv[2], argv[3]);
        temp->SetSourceFile(argv[1]);
        temp->BuildStructureByTOPFileInformation();

//        temp->BuildAssemblyFromPdbFile(argv[1]);
//        temp->BuildStructureByDistance(5);

//        temp->BuildAssemblyFromPrepFile(argv[1], argv[3]);
//        temp->SetSourceFile(argv[1]);
//        temp->BuildStructureByDistance();
        temp->Solvation(10  ,3,argv[4]);
//        temp->Print();
//        std::cout << solvated->GetAssemblies().at(0)->GetId() << " " << solvated->GetAssemblies().at(1)->GetId() << std::endl;




        PdbFileSpace::PdbFile* pdb_file = temp->BuildPdbFileStructureFromAssembly(1,0);
        pdb_file->Write("pdb_file.pdb");
        TopologyFileSpace::TopologyFile* top_file = temp->BuildTopologyFileStructureFromAssembly(argv[3]);
        top_file->Write("top_file.parm7");
        CoordinateFileSpace::CoordinateFile* crd_file = temp->BuildCoordinateFileStructureFromAssembly();
        crd_file->Write("crd_file.rst7");

    }
}
#endif

#ifdef SplitAssemblyTest
using namespace MolecularModeling;
int main(int argc, char *argv[])
{
    if(argc < 2)
        std::cout << "Not enough parameters" << std::endl;
    else
    {
        Assembly* temp = new Assembly();
        temp->BuildAssemblyFromPdbFile(argv[1]);
        temp->BuildStructureByDistance();
        temp->Print(std::cout);
        PdbFileSpace::PdbFile* pdb_file = temp->BuildPdbFileStructureFromAssembly();
        pdb_file->Write("pdb_file.pdb");

        Assembly* solvent = new Assembly();
        Assembly* solute = new Assembly();
        temp->SplitSolvent(solvent, solute);
        PdbFileSpace::PdbFile* pdb_solvent = solvent->BuildPdbFileStructureFromAssembly();
        pdb_solvent->Write("solvent.pdb");
        PdbFileSpace::PdbFile* pdb_solute = solute->BuildPdbFileStructureFromAssembly();
        pdb_solute->Write("solute.pdb");
    }
    return 0;
}
#endif
