#include <iostream>

#include "../../../includes/FileSet/PdbFileSpace/pdbcompoundspecification.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbCompoundSpecification::PdbCompoundSpecification() : molecule_id_(""), molecule_name_(""), chain_ids_(), fragment_(""),
    molecule_synonyms_(), enzyme_commission_numbers_(), is_engineered_(false), has_mutation_(false), comments_(""){}

PdbCompoundSpecification::PdbCompoundSpecification(const string& molecule_id, const string& molecule_name) : molecule_id_(molecule_id), molecule_name_(molecule_name){}

PdbCompoundSpecification::PdbCompoundSpecification(const string &molecule_id, const string &molecule_name, const vector<string> &chain_ids, const string &fragment,
                                                   const vector<string> &molecule_synonyms, vector<string> &enzyme_commission_numbers, bool is_engineered, bool has_mutation, const string& comments) :
    molecule_id_(molecule_id), molecule_name_(molecule_name), chain_ids_(chain_ids), fragment_(fragment), molecule_synonyms_(molecule_synonyms),
    enzyme_commission_numbers_(enzyme_commission_numbers), is_engineered_(is_engineered), has_mutation_(has_mutation), comments_(comments){}

PdbCompoundSpecification::PdbCompoundSpecification(stringstream& specification_block) : molecule_id_(""), molecule_name_(""), chain_ids_(), fragment_(""),
    molecule_synonyms_(), enzyme_commission_numbers_(), is_engineered_(false), has_mutation_(false), comments_("")
{    
    string line;
    getline(specification_block, line);
    string temp = line;
    int flag = 0;
    stringstream molecule_name, chain_id, fragment, molecule_synonyms, enzyme_commission_numbers, comments;
    while(!Trim(temp).empty())
    {
        vector<string> tokens = Split(line,":;");

        string token_name = Trim(tokens.at(0));

        if(token_name == "MOL_ID")
        {
            molecule_id_ = tokens.at(1);
            Trim(molecule_id_);
        }
        else if(token_name == "MOLECULE")
        {
            string s = tokens.at(1);
            molecule_name << s;
            flag = 1;
        }
        else if(token_name == "CHAIN")
        {
            string s = tokens.at(1);
            chain_id << s;
            flag = 2;
        }
        else if(token_name == "FRAGMENT")
        {
            string s = tokens.at(1);
            fragment << s;
            flag = 3;
        }
        else if(token_name == "SYNONYM")
        {
            string s = tokens.at(1);
            molecule_synonyms << s;
            flag = 4;
        }
        else if(token_name == "EC")
        {
            string s = tokens.at(1);
            enzyme_commission_numbers << s;
            flag = 5;
        }
        else if(token_name == "ENGINEERED")
        {
            string status = Trim(tokens.at(1));
            if(status == "YES")
                is_engineered_ = true;
            else
                is_engineered_ = false;
        }
        else if(token_name == "MUTATION")
        {
            string status = Trim(tokens.at(1));
            if(status == "YES")
                has_mutation_ = true;
            else
                has_mutation_ = false;
        }
        else if(token_name == "OTHER_DETAILS")
        {
            string s = tokens.at(1);
            comments << s;
            flag = 6;
        }
        else
        {
            switch (flag)
            {
                case 1:
                    molecule_name << line;
                    break;
                case 2:
                    chain_id << line;
                    break;
                case 3:
                    fragment << line;
                    break;
                case 4:
                    molecule_synonyms << line;
                    break;
                case 5:
                    enzyme_commission_numbers << line;
                    break;
                case 6:
                    comments << line;
                    break;
            }
        }
        getline(specification_block, line);
        temp = line;
    }
    if(molecule_name.str().length() > 0)
    {
        molecule_name_ = Split(molecule_name.str(),";").at(0);
        Trim(molecule_name_);
    }

    if(chain_id.str().length() > 0)
    {
        string chain_ids = chain_id.str();
        Trim(chain_ids);
        chain_ids_ = Split(chain_ids,",;");
//        for(vector<string>::iterator it = chain_ids_.begin(); it != chain_ids_.end(); it++)
//        {
//            Trim(*it);
//        }
    }

    if(fragment.str().length() > 0)
    {
        fragment_ = fragment.str();
        Trim(fragment_);
        fragment_ = Split(fragment_, ";").at(0);
        Trim(fragment_);
    }

    if(molecule_synonyms.str().length() > 0)
    {
        string synonyms = molecule_synonyms.str();
        Trim(synonyms);
        molecule_synonyms_ = Split(synonyms, ",;");
//        for(vector<string>::iterator it = molecule_synonyms_.begin(); it != molecule_synonyms_.end(); it++)
//        {
//            Trim(*it);
//        }
    }

    if(enzyme_commission_numbers.str().length() > 0)
    {
        string commission_numbers = enzyme_commission_numbers.str();
        Trim(commission_numbers);
        enzyme_commission_numbers_ = Split(commission_numbers, ",;");
//        for(vector<string>::iterator it = enzyme_commission_numbers_.begin(); it != enzyme_commission_numbers_.end(); it++)
//        {
//            Trim(*it);
//        }
    }

    if(comments.str().length() > 0)
    {
        comments_ = Split(comments.str(),";").at(0);
        Trim(comments_);
    }
}

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

vector<string> PdbCompoundSpecification::GetEnzymeCommissionNumbers()
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

void PdbCompoundSpecification::SetMoleculeSynonyms(std::vector<std::string> molecule_synonyms)
{
    molecule_synonyms_.clear();
    for(vector<string>::const_iterator it = molecule_synonyms.begin(); it != molecule_synonyms.end(); it++)
    {
        molecule_synonyms_.push_back(*it);
    }
}

void PdbCompoundSpecification::AddMoleculeSynonym(const std::string molecule_synonym)
{
    molecule_synonyms_.push_back(molecule_synonym);
}

void PdbCompoundSpecification::SetEnzymeCommissionNumbers(const vector<string> enzyme_commission_numbers)
{
    enzyme_commission_numbers_.clear();
    for(vector<string>::const_iterator it = enzyme_commission_numbers.begin(); it != enzyme_commission_numbers.end(); it++)
    {
        enzyme_commission_numbers_.push_back(*it);
    }
}

void PdbCompoundSpecification::AddEnzymeCommissionNumber(string enzyme_commission_number)
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
void PdbCompoundSpecification::Print(ostream &out)
{
    out << "Molecule ID: " << molecule_id_ << ", Molecule Name: " << molecule_name_ << endl;
    out << "Chain IDs: ";
    for(vector<string>::iterator it = chain_ids_.begin(); it != chain_ids_.end(); it++)
    {
        out << (*it) << ", ";
    }
    out << endl << "Fragment: " << fragment_ << endl << "Molecule Synonyms: ";
    for(vector<string>::iterator it = molecule_synonyms_.begin(); it != molecule_synonyms_.end(); it++)
    {
        out << (*it) << ", ";
    }
    out << endl << "Enzyme Commission Numbers: ";
    for(vector<string>::iterator it = enzyme_commission_numbers_.begin(); it != enzyme_commission_numbers_.end(); it++)
    {
        out << (*it) << ", ";
    }
    out << endl << "Engineered: ";
    if(is_engineered_)
        out << "YES" << endl;
    else
        out << "NO" << endl;
    out << "Mutation: ";
    if(has_mutation_)
        out << "YES" << endl;
    else
        out << "NO" << endl;
    out << "Comments: " << comments_ << endl;
}
