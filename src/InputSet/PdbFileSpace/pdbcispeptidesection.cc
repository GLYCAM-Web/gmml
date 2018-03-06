#include "../../../includes/InputSet/PdbFileSpace/pdbcispeptidesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbcispeptidecard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbCISPeptideSection::PdbCISPeptideSection() {}
PdbCISPeptideSection::PdbCISPeptideSection(stringstream &stream_block)
{
    string line;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
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
void PdbCISPeptideSection::Print(ostream &out)
{
    for(CISPeptideCardVector::iterator it = cis_peptide_.begin(); it != cis_peptide_.end(); it++)
            (*it)->Print(out);
    out << endl;
}
