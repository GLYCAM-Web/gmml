#ifndef PDBCOMPOUNDSPECIFICATION_HPP
#define PDBCOMPOUNDSPECIFICATION_HPP

#include <string>
#include <vector>

namespace PdbFileSpace
{
    class PdbCompoundSpecification
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbCompoundSpecification();
            PdbCompoundSpecification(const std::string& molecule_id, const std::string& molecule_name);
            PdbCompoundSpecification(const std::string& molecule_id, const std::string& molecule_name, const std::vector<std::string>& chain_ids,
                                     const std::string& fragment, const std::vector<std::string>& molecule_synonyms, std::vector<int>& enzyme_commission_numbers,
                                     bool is_engineered, bool has_mutation, const std::string& comments);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetMoleculeId();
            std::string GetMoleculeName();
            std::vector<std::string> GetChainIds();
            std::string GetFragment();
            std::vector<std::string> GetMoleculeSynonyms();
            std::vector<int> GetEnzymeCommissionNumbers();
            bool GetIsEngineered();
            bool GetHasMutation();
            std::string GetComments();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetMoleculeId(const std::string molecule_id);
            void SetMoleculeName(const std::string molecule_name);
            void SetChainIds(const std::vector<std::string> chain_ids);
            void AddChainId(const std::string chain_id);
            void SetFragment(const std::string fragment);
            void SetMoleculeSynonyms(const std::vector<std::string> molecule_synonyms);
            void AddMoleculeSynonym(const std::string molecule_synonym);
            void SetEnzymeCommissionNumbers(const std::vector<int> enzyme_commission_numbers);
            void AddEnzymeCommissionNumber(int enzyme_commission_number);
            void SetIsEngineered(bool is_engineered);
            void SetHasMutation(bool has_mutation);
            void setComments(const std::string comments);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string molecule_id_;
            std::string molecule_name_;
            std::vector<std::string> chain_ids_;
            std::string fragment_;
            std::vector<std::string> molecule_synonyms_;
            std::vector<int> enzyme_commission_numbers_;
            bool is_engineered_;
            bool has_mutation_;
            std::string comments_;
    };
}
#endif // PDBCOMPOUNDSPECIFICATION_HPP
