#include "../../../includes/FileSet/PdbFileSpace/pdbcompoundspecification.h"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbCompoundSpecification::PdbCompoundSpecification() {}

PdbCompoundSpecification::PdbCompoundSpecification(string molecule_id, string molecule_name) : molecule_id_(molecule_id), molecule_name_(molecule_name){}

PdbCompoundSpecification::PdbCompoundSpecification(string molecule_id, string molecule_name, vector<string> &chain_ids, string fragment, vector<string> &molecule_synonyms,
                                                   vector<int> &enzyme_commission_numbers, bool is_engineered, bool has_mutation, string comments) :
    molecule_id_(molecule_id), molecule_name_(molecule_name), chain_ids_(chain_ids), fragment_(fragment), molecule_synonyms_(molecule_synonyms),
    enzyme_commission_numbers_(enzyme_commission_numbers), is_engineered_(is_engineered), has_mutation_(has_mutation), comments_(comments){}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbCompoundSpecification::GetMoleculeId()
{
    return molecule_id_;
}

string PdbCompoundSpecification::GetMoleculeName()
{
    return molecule_name_;
}

vector<string> PdbCompoundSpecification::GetChainIds()
{
    return chain_ids_;
}

string PdbCompoundSpecification::GetFragment()
{
    return fragment;
}

vector<string> PdbCompoundSpecification::GetFragment()
{
    return fragment_;
}

vector<int> PdbCompoundSpecification::GetEnzymeCommissionNumbers()
{
    return enzyme_commission_numbers_;
}

bool PdbCompoundSpecification::GetIsEngineered()
{
    return is_engineered_;
}

bool PdbCompoundSpecification::GetHasMutation()
{
    return has_mutation_;
}

string PdbCompoundSpecification::GetComments()
{
    return comments;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbCompoundSpecification::SetMoleculeId(string molecule_id)
{
    molecule_id_ = molecule_id;
}

void PdbCompoundSpecification::SetMoleculeName(string molecule_name)
{
    molecule_name_ = molecule_name;
}

void PdbCompoundSpecification::SetChainIds(vector<string> chain_ids)
{
    chain_ids_.clear();
    for(vector<string>::const_iterator it = chain_ids.begin(); it != chain_ids.end(); it++)
    {
        chain_ids_.push_back(*it);
    }
}

void PdbCompoundSpecification::AddChainId(string chain_id)
{
    chain_ids_.push_back(chain_id);
}

void PdbCompoundSpecification::SetFragment(string fragment)
{
    fragment_ = fragment;
}

void PdbCompoundSpecification::SetEnzymeCommissionNumbers(vector<int> enzyme_commission_numbers)
{
    enzyme_commission_numbers_.clear();
    for(vector<int>::const_iterator it = enzyme_commission_numbers.begin(); it != enzyme_commission_numbers.end(); it++)
    {
        enzyme_commission_numbers_.push_back(*it);
    }
}

void PdbCompoundSpecification::AddEnzymeCommissionNumber(int enzyme_commission_number)
{
    enzyme_commission_numbers_.push_back(enzyme_commission_number);
}

void PdbCompoundSpecification::SetIsEngineered(bool is_engineered)
{
    is_engineered_ = is_engineered;
}

void PdbCompoundSpecification::SetHasMutation(bool has_mutation)
{
    has_mutation_ = has_mutation;
}

void PdbCompoundSpecification::setComments(string comments)
{
    comments_ = comments;
}
