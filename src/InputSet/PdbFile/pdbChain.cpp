//#include "includes/InputSet/PdbFile/pdbChain.hpp"
//#include "includes/CodeUtils/logging.hpp"
//
//using pdb::PdbChain;
////////////////////////////////////////////////////////////
////                       CONSTRUCTOR                    //
////////////////////////////////////////////////////////////
//PdbChain::PdbChain(PdbResidue* pdbResidue)
//{
//    this->AddResidue(pdbResidue);
//}
//PdbChain::PdbChain(std::vector<PdbResidue*> pdbResidues)
//{
//    pdbResidues_ = pdbResidues;
//}
////////////////////////////////////////////////////////////
////                         ACCESSOR                     //
////////////////////////////////////////////////////////////
//pdb::PdbResidue* PdbChain::GetFirstResidue() const
//{
//    return pdbResidues_.front();
//}
//const std::string& PdbChain::GetChainId() const
//{
//    return this->GetFirstResidue()->GetChainId();
//}
////////////////////////////////////////////////////////////
////                    FUNCTIONS                         //
////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////
////                          MUTATOR                     //
////////////////////////////////////////////////////////////
//void PdbChain::AddResidue(PdbResidue* residue)
//{
//    pdbResidues_.push_back(residue);
//    return;
//}
////////////////////////////////////////////////////////////
////                      DISPLAY FUNCTION                //
////////////////////////////////////////////////////////////
//void PdbChain::Print(std::ostream &out) const
//{
//    out << "Printing chain " << this->GetChainId() << "\n";
//    for(auto &residue : pdbResidues_)
//    {
//        residue->Print();
//    }
//}
