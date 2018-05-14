#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/pdbcompoundspecification.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbCompoundSpecification;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbCompoundSpecification::PdbCompoundSpecification() : molecule_id_(""), molecule_name_(""), chain_ids_(), fragment_(""),
    molecule_synonyms_(), enzyme_commission_numbers_(), is_engineered_(false), has_mutation_(false), comments_(""){}

PdbCompoundSpecification::PdbCompoundSpecification(const std::string& molecule_id, const std::string& molecule_name) : molecule_id_(molecule_id), molecule_name_(molecule_name){}

PdbCompoundSpecification::PdbCompoundSpecification(const std::string &molecule_id, const std::string &molecule_name, const std::vector<std::string> &chain_ids, const std::string &fragment,
                                                   const std::vector<std::string> &molecule_synonyms, std::vector<std::string> &enzyme_commission_numbers, bool is_engineered, bool has_mutation, const std::string& comments) :
    molecule_id_(molecule_id), molecule_name_(molecule_name), chain_ids_(chain_ids), fragment_(fragment), molecule_synonyms_(molecule_synonyms),
    enzyme_commission_numbers_(enzyme_commission_numbers), is_engineered_(is_engineered), has_mutation_(has_mutation), comments_(comments){}

PdbCompoundSpecification::PdbCompoundSpecification(std::stringstream& specification_block) : molecule_id_(""), molecule_name_(""), chain_ids_(), fragment_(""),
    molecule_synonyms_(), enzyme_commission_numbers_(), is_engineered_(false), has_mutation_(false), comments_("")
{
    std::string line;
    getline(specification_block, line);
    std::string temp = line;
    int flag = 0;
    std::stringstream molecule_name, chain_id, fragment, molecule_synonyms, enzyme_commission_numbers, comments;
    while(!gmml::Trim(temp).empty())
    {
        std::vector<std::string> tokens = gmml::Split(line,":;");

        std::string token_name = gmml::Trim(tokens.at(0));

        if(token_name == "MOL_ID")
        {
            molecule_id_ = tokens.at(1);
            gmml::Trim(molecule_id_);
        }
        else if(token_name == "MOLECULE")
        {
            std::string s = tokens.at(1);
            molecule_name << s;
            flag = 1;
        }
        else if(token_name == "CHAIN")
        {
            std::string s = tokens.at(1);
            chain_id << s;
            flag = 2;
        }
        else if(token_name == "FRAGMENT")
        {
            std::string s = tokens.at(1);
            fragment << s;
            flag = 3;
        }
        else if(token_name == "SYNONYM")
        {
            std::string s = tokens.at(1);
            molecule_synonyms << s;
            flag = 4;
        }
        else if(token_name == "EC")
        {
            std::string s = tokens.at(1);
            enzyme_commission_numbers << s;
            flag = 5;
        }
        else if(token_name == "ENGINEERED")
        {
            std::string status = gmml::Trim(tokens.at(1));
            if(status == "YES")
                is_engineered_ = true;
            else
                is_engineered_ = false;
        }
        else if(token_name == "MUTATION")
        {
            std::string status = gmml::Trim(tokens.at(1));
            if(status == "YES")
                has_mutation_ = true;
            else
                has_mutation_ = false;
        }
        else if(token_name == "OTHER_DETAILS")
        {
            std::string s = tokens.at(1);
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
        molecule_name_ = gmml::Split(molecule_name.str(),";").at(0);
        gmml::Trim(molecule_name_);
    }

    if(chain_id.str().length() > 0)
    {
        std::string chain_ids = chain_id.str();
        gmml::Trim(chain_ids);
        chain_ids_ = gmml::Split(chain_ids,",;");
//        for(std::vector<std::string>::iterator it = chain_ids_.begin(); it != chain_ids_.end(); it++)
//        {
//            gmml::Trim(*it);
//        }
    }

    if(fragment.str().length() > 0)
    {
        fragment_ = fragment.str();
        gmml::Trim(fragment_);
        fragment_ = gmml::Split(fragment_, ";").at(0);
        gmml::Trim(fragment_);
    }

    if(molecule_synonyms.str().length() > 0)
    {
        std::string synonyms = molecule_synonyms.str();
        gmml::Trim(synonyms);
        molecule_synonyms_ = gmml::Split(synonyms, ",;");
//        for(std::vector<std::string>::iterator it = molecule_synonyms_.begin(); it != molecule_synonyms_.end(); it++)
//        {
//            gmml::Trim(*it);
//        }
    }

    if(enzyme_commission_numbers.str().length() > 0)
    {
        std::string commission_numbers = enzyme_commission_numbers.str();
        gmml::Trim(commission_numbers);
        enzyme_commission_numbers_ = gmml::Split(commission_numbers, ",;");
//        for(std::vector<std::string>::iterator it = enzyme_commission_numbers_.begin(); it != enzyme_commission_numbers_.end(); it++)
//        {
//            gmml::Trim(*it);
//        }
    }

    if(comments.str().length() > 0)
    {
        comments_ = gmml::Split(comments.str(),";").at(0);
        gmml::Trim(comments_);
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbCompoundSpecification::GetMoleculeId()
{
    return molecule_id_;
}

std::string PdbCompoundSpecification::GetMoleculeName()
{
    return molecule_name_;
}

std::vector<std::string> PdbCompoundSpecification::GetChainIds()
{
    return chain_ids_;
}

std::string PdbCompoundSpecification::GetFragment()
{
    return fragment_;
}

std::vector<std::string> PdbCompoundSpecification::GetMoleculeSynonyms()
{
    return molecule_synonyms_;
}

std::vector<std::string> PdbCompoundSpecification::GetEnzymeCommissionNumbers()
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

std::string PdbCompoundSpecification::GetComments()
{
    return comments_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbCompoundSpecification::SetMoleculeId(const std::string molecule_id)
{
    molecule_id_ = molecule_id;
}

void PdbCompoundSpecification::SetMoleculeName(const std::string molecule_name)
{
    molecule_name_ = molecule_name;
}

void PdbCompoundSpecification::SetChainIds(const std::vector<std::string> chain_ids)
{
    chain_ids_.clear();
    for(std::vector<std::string>::const_iterator it = chain_ids.begin(); it != chain_ids.end(); it++)
    {
        chain_ids_.push_back(*it);
    }
}

void PdbCompoundSpecification::AddChainId(const std::string chain_id)
{
    chain_ids_.push_back(chain_id);
}

void PdbCompoundSpecification::SetFragment(const std::string fragment)
{
    fragment_ = fragment;
}

void PdbCompoundSpecification::SetMoleculeSynonyms(std::vector<std::string> molecule_synonyms)
{
    molecule_synonyms_.clear();
    for(std::vector<std::string>::const_iterator it = molecule_synonyms.begin(); it != molecule_synonyms.end(); it++)
    {
        molecule_synonyms_.push_back(*it);
    }
}

void PdbCompoundSpecification::AddMoleculeSynonym(const std::string molecule_synonym)
{
    molecule_synonyms_.push_back(molecule_synonym);
}

void PdbCompoundSpecification::SetEnzymeCommissionNumbers(const std::vector<std::string> enzyme_commission_numbers)
{
    enzyme_commission_numbers_.clear();
    for(std::vector<std::string>::const_iterator it = enzyme_commission_numbers.begin(); it != enzyme_commission_numbers.end(); it++)
    {
        enzyme_commission_numbers_.push_back(*it);
    }
}

void PdbCompoundSpecification::AddEnzymeCommissionNumber(std::string enzyme_commission_number)
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

void PdbCompoundSpecification::setComments(const std::string comments)
{
    comments_ = comments;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbCompoundSpecification::Print(std::ostream &out)
{
    out << "Molecule ID: " << molecule_id_ << ", Molecule Name: " << molecule_name_ << std::endl;
    out << "Chain IDs: ";
    for(std::vector<std::string>::iterator it = chain_ids_.begin(); it != chain_ids_.end(); it++)
    {
        out << (*it) << ", ";
    }
    out << std::endl << "Fragment: " << fragment_ << std::endl << "Molecule Synonyms: ";
    for(std::vector<std::string>::iterator it = molecule_synonyms_.begin(); it != molecule_synonyms_.end(); it++)
    {
        out << (*it) << ", ";
    }
    out << std::endl << "Enzyme Commission Numbers: ";
    for(std::vector<std::string>::iterator it = enzyme_commission_numbers_.begin(); it != enzyme_commission_numbers_.end(); it++)
    {
        out << (*it) << ", ";
    }
    out << std::endl << "Engineered: ";
    if(is_engineered_)
        out << "YES" << std::endl;
    else
        out << "NO" << std::endl;
    out << "Mutation: ";
    if(has_mutation_)
        out << "YES" << std::endl;
    else
        out << "NO" << std::endl;
    out << "Comments: " << comments_ << std::endl;
}
