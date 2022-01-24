//#include "includes/InputSet/PdbFile/pdbResidue.hpp"
//
//using pdb::Residue;
////////////////////////////////////////////////////////////
////                       CONSTRUCTOR                    //
////////////////////////////////////////////////////////////
//Residue::Residue(AtomRecord* atomRecord)
//{
//    this->AddAtom(atomRecord);
//}
//Residue::Residue(std::vector<AtomRecord*> atomRecords)
//{
//    atomRecords_ = atomRecords;
//}
////////////////////////////////////////////////////////////
////                         ACCESSOR                     //
////////////////////////////////////////////////////////////
//pdb::AtomRecord* Residue::GetFirstAtom() const
//{
//    return atomRecords_.front();
//}
//const std::string& Residue::GetChainId() const
//{
//    return this->GetFirstAtom()->GetChainId();
//}
//const std::string& Residue::GetName() const
//{
//    return this->GetFirstAtom()->GetResidueName();
//}
//const int& Residue::GetSequenceNumber() const
//{
//    return this->GetFirstAtom()->GetResidueSequenceNumber();
//}
//const std::string& Residue::GetInsertionCode() const
//{
//    return this->GetFirstAtom()->GetInsertionCode();
//}
//std::string Residue::GetId() const
//{
//    return this->GetFirstAtom()->GetResidueId();
//}
////////////////////////////////////////////////////////////
////                          MUTATOR                     //
////////////////////////////////////////////////////////////
//void Residue::AddAtom(AtomRecord* atomRecord)
//{
//    std::cout << "Adding atom with residue id " << atomRecord->GetResidueId() << " to residue.\n";
//    atomRecords_.push_back(atomRecord);
//}
////////////////////////////////////////////////////////////
////                      DISPLAY FUNCTION                //
////////////////////////////////////////////////////////////
//void Residue::Print(std::ostream &out) const
//{
//    out << "pdb::Residue : " << this->GetId() << std::endl;
//}
