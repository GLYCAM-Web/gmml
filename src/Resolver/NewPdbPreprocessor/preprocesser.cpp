//#include <vector>
//#include "includes/Resolver/NewPdbPreprocessor/preprocesser.hpp"
//#include "includes/InputSet/PdbFile/pdbResidue.hpp"
//#include "includes/CodeUtils/logging.hpp"
//
//using pdb::Preprocesser;
//
//Preprocesser::Preprocesser(PdbFile &pdbFile)
//{
//    // CYS
//    std::vector<pdb::PdbResidue> cysResidues = pdbFile.GetResidues("CYS");
//    for (auto &cysRes : cysResidues)
//    {
//        if (residue->GetAtom("SG").Distance(cysRes->GetAtom("SG")) < gmml::CysBondLength)
//        {
//            residue->SetName("CYX");
//            cysRes->SetName("CYX");
//            this->AddConnection(residue->GetAtom("SG").GetIndex(), cysRes->GetAtom("SG").GetIndex());
//        }
//    }
//    // HIS
//    std::vector<pdb::PdbResidue> hisResidues = this->GetCoordinateSection().GetResidues("HIS");
//
//    // Chain terminations
//    std::vector<pdb::PdbResidue> residues = this->GetCoordinateSection().GetResidues();
//    for (auto &residue : residues)
//    {
//        if (residue->GetName() == "HIS")
//        {
//            // Check what hydrogen atoms it contains. Decide on name.
//        }
//        else if (residue->GetName() == "CYS")
//        {
//
//            cysResidues.push_back(residue);
//        }
//
//    }
//    for (auto &cysResidue : cysResidues)
//    {
//
//    }
//
//}
