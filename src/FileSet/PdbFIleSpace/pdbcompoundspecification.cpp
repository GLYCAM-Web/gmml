#include "../../../includes/FileSet/PdbFileSpace/pdbcompoundspecification.h"

using namespace std;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbCompoundSpecification::PdbCompoundSpecification() {}

PdbCompoundSpecification::PdbCompoundSpecification(const string& molecule_id, const string& molecule_name) : molecule_id_(molecule_id), molecule_name_(molecule_name){}

PdbCompoundSpecification::PdbCompoundSpecification(const string &molecule_id, const string &molecule_name, const vector<string> &chain_ids, const string &fragment,
                                                   const vector<string> &molecule_synonyms, vector<int> &enzyme_commission_numbers, bool is_engineered, bool has_mutation, const string& comments) :
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
    return fragment_;
}

vector<string> PdbCompoundSpecification::GetMoleculeSynonyms()
{
    return molecule_synonyms_;
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
    return comments_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbCompoundSpecification::SetMoleculeId(const string molecule_id)
{
    molecule_id_ = molecule_id;
}

void PdbCompoundSpecification::SetMoleculeName(const string molecule_name)
{
    molecule_name_ = molecule_name;
}

void PdbCompoundSpecification::SetChainIds(const vector<string> chain_ids)
{
    chain_ids_.clear();
    for(vector<string>::const_iterator it = chain_ids.begin(); it != chain_ids.end(); it++)
    {
        chain_ids_.push_back(*it);
    }
}

void PdbCompoundSpecification::AddChainId(const string chain_id)
{
    chain_ids_.push_back(chain_id);
}

void PdbCompoundSpecification::SetFragment(const string fragment)
{
    fragment_ = fragment;
}

void PdbCompoundSpecification::SetEnzymeCommissionNumbers(const vector<int> enzyme_commission_numbers)
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

void PdbCompoundSpecification::setComments(const string comments)
{
    comments_ = comments;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
