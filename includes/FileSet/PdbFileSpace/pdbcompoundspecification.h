#ifndef PDBCOMPOUNDSPECIFICATION_H
#define PDBCOMPOUNDSPECIFICATION_H

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
            PdbCompoundSpecification(std::string molecule_id, std::string molecule_name);
            PdbCompoundSpecification(std::string molecule_id, std::string molecule_name, std::vector<std::string>& chain_ids,
                                     std::string fragment, std::vector<std::string>& molecule_synonyms, std::vector<int>& enzyme_commission_numbers,
                                     bool is_engineered, bool has_mutation, std::string comments);

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
            void SetMoleculeId(std::string molecule_id);
            void SetMoleculeName(std::string molecule_name);
            void SetChainIds(std::vector<std::string> chain_ids);
            void AddChainId(std::string chain_id);
            void SetFragment(std::string fragment);
            void SetMoleculeSynonyms(std::vector<std::string> molecule_synonyms);
            void AddMoleculeSynonym(std::string molecule_synonym);
            void SetEnzymeCommissionNumbers(std::vector<int> enzyme_commission_numbers);
            void AddEnzymeCommissionNumber(int enzyme_commission_number);
            void SetIsEngineered(bool is_engineered);
            void SetHasMutation(bool has_mutation);
            void setComments(std::string comments);

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
#endif // PDBCOMPOUNDSPECIFICATION_H
