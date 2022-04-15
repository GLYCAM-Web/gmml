#include "../../../includes/InputSet/PdbFileSpace/pdbcispeptidesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbcispeptidecard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbCISPeptideSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbCISPeptideSection::PdbCISPeptideSection() {}
PdbCISPeptideSection::PdbCISPeptideSection(std::stringstream &stream_block)
{
    std::string line;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        PdbCISPeptideCard* cis_peptide = new PdbCISPeptideCard(line);
        AddCISPeptideCards(cis_peptide);
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbCISPeptideSection::CISPeptideCardVector PdbCISPeptideSection::GetCISPeptideCards()
{
    return cis_peptide_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbCISPeptideSection::SetCISPeptideCards(CISPeptideCardVector cis_peptide)
{
    cis_peptide_.clear();
    for(CISPeptideCardVector::iterator it = cis_peptide.begin(); it != cis_peptide.end(); it++)
    {
        cis_peptide.push_back(*it);
    }
}

void PdbCISPeptideSection::AddCISPeptideCards(PdbCISPeptideCard *cis_peptide)
{
    cis_peptide_.push_back(cis_peptide);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbCISPeptideSection::Print(std::ostream &out)
{
    for(CISPeptideCardVector::iterator it = cis_peptide_.begin(); it != cis_peptide_.end(); it++)
            (*it)->Print(out);
    out << std::endl;
}
