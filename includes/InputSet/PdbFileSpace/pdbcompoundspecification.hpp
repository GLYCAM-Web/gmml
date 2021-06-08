// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBCOMPOUNDSPECIFICATION_HPP
#define PDBCOMPOUNDSPECIFICATION_HPP

#include <string>
#include <vector>
#include <iostream>

namespace PdbFileSpace
{
    class PdbCompoundSpecification
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbCompoundSpecification();
            /*! \fn
              * Constructor with required parameters
              * @param molecule_id
              * @param molecule_name
              */
            PdbCompoundSpecification(const std::string& molecule_id, const std::string& molecule_name);
            /*! \fn
              * Constructor with required parameters
              * @param molecule_id Molecule identifier
              * @param molecule_name Molecule name
              * @param chain_ids Chain identifiers of the molecule
              * @param fragment Fragment
              * @param molecule_synonyms Molecule synonyms
              * @param enzyme_commission_numbers Enzyme omission numbers of the molecule
              * @param is_engineered Indicates that the molecule is engineered or not
              * @param has_mutation Indicates that the molecule has been mutated or not
              * @param comments Comments
              */
            PdbCompoundSpecification(const std::string& molecule_id, const std::string& molecule_name, const std::vector<std::string>& chain_ids,
                                     const std::string& fragment, const std::vector<std::string>& molecule_synonyms, std::vector<std::string>& enzyme_commission_numbers,
                                     bool is_engineered, bool has_mutation, const std::string& comments);
            /*! \fn
              * A constructor that get a stream block of compound specifications and parse the whole block to fill the related fields
              * @param stream_block A whole block of compound specifications of a specific molecule id in a pdb file
              */
            PdbCompoundSpecification(std::stringstream& specification_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the molecule id in compound specification
              * @return molecule_id_ attribute of the current object of this class
              */
            std::string GetMoleculeId();
            /*! \fn
              * An accessor function in order to access to the molecule name in compound specification
              * @return molecule_name_ attribute of the current object of this class
              */
            std::string GetMoleculeName();
            /*! \fn
              * An accessor function in order to access to the chain ids in compound specification
              * @return chain_ids_ of the current object of this class
              */
            std::vector<std::string> GetChainIds();
            /*! \fn
              * An accessor function in order to access to the fragment in compound specification
              * @return fragment_ attribute of the current object of this class
              */
            std::string GetFragment();
            /*! \fn
              * An accessor function in order to access to the molecule synonyms in compound specification
              * @return molecule_synonyms_ of the current object of this class
              */
            std::vector<std::string> GetMoleculeSynonyms();
            /*! \fn
              * An accessor function in order to access to the enzyme commission numbers in compound specification
              * @return enzyme_commission_numbers_ of the current object of this class
              */
            std::vector<std::string> GetEnzymeCommissionNumbers();
            /*! \fn
              * An accessor function in order to access to the is engineered attribute in compound specification
              * @return is_engineered_ attribute of the current object of this class
              */
            bool GetIsEngineered();
            /*! \fn
              * An accessor function in order to access to the has mutation attribute in compound specification
              * @return has_mutation_ attribute of the current object of this class
              */
            bool GetHasMutation();
            /*! \fn
              * An accessor function in order to access to the comments in compound specification
              * @return comments_ attribute of the current object of this class
              */
            std::string GetComments();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the molecule id of the current object
              * Set the molecule_id_ attribute of the compound specification
              * @param molecule_id The molecule id of the current object
              */
            void SetMoleculeId(const std::string molecule_id);
            /*! \fn
              * A mutator function in order to set the molecule name of the current object
              * Set the molecule_name_ attribute of the compound specification
              * @param molecule_name The molecule name of the current object
              */
            void SetMoleculeName(const std::string molecule_name);
            /*! \fn
              * A mutator function in order to set the list of chain ids of the current object
              * Set the chain_ids_ of the compound specification
              * @param chain_ids The chain ids of the current object
              */
            void SetChainIds(const std::vector<std::string> chain_ids);
            /*! \fn
              * A function in order to add chain id to the current object
              * @param chain_id The chain id of the current object
              */
            void AddChainId(const std::string chain_id);
            /*! \fn
              * A mutator function in order to set the fragment of the current object
              * Set the fragment_ attribute of the compound specification
              * @param fragment The fragment of the current object
              */
            void SetFragment(const std::string fragment);
            /*! \fn
              * A mutator function in order to set the list of molecule synonyms of the current object
              * Set the molecule_synonyms_ of the compound specification
              * @param molecule_synonyms The molecule synonyms of the current object
              */
            void SetMoleculeSynonyms(std::vector<std::string> molecule_synonyms);
            /*! \fn
              * A function in order to add the molecule synonym to the current object
              * Set the molecule_synonym_ attribute of the compound specification
              * @param molecule_synonym The molecule synonym of the current object
              */
            void AddMoleculeSynonym(const std::string molecule_synonym);
            /*! \fn
              * A mutator function in order to set the enzyme commission numbers of the current object
              * Set the enzyme_commission_numbers_ attribute of the compound specification
              * @param enzyme_commission_numbers The enzyme commission numbers of the current object
              */
            void SetEnzymeCommissionNumbers(std::vector<std::string> enzyme_commission_numbers);
            /*! \fn
              * A function in order to add the enzyme commission number to the current object
              * Set the enzyme_commission_number_ attribute of the compound specification
              * @param enzyme_commission_number The enzyme commission number of the current object
              */
            void AddEnzymeCommissionNumber(std::string enzyme_commission_number);
            /*! \fn
              * A mutator function in order to set the is engineered attribute of the current object
              * Set the is_engineered_ attribute of the compound specification
              * @param is_engineered The is engineered attribute of the current object
              */
            void SetIsEngineered(bool is_engineered);
            /*! \fn
              * A mutator function in order to set thehas mutation attribute of the current object
              * Set the has_mutation_ attribute of the compound specification
              * @param has_mutation The has mutation attribute of the current object
              */
            void SetHasMutation(bool has_mutation);
            /*! \fn
              * A mutator function in order to set the comments of the current object
              * Set the comments_ attribute of the compound specification
              * @param comments The comments of the current object
              */
            void setComments(const std::string comments);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the compund specification contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string molecule_id_;                               /*!< Molecule identifier >*/
            std::string molecule_name_;                             /*!< Molecule name >*/
            std::vector<std::string> chain_ids_;                    /*!< Chain identifiers involving in the molecule >*/
            std::string fragment_;                                  /*!< Fragment specification >*/
            std::vector<std::string> molecule_synonyms_;            /*!< Synonyms of the molecule >*/
            std::vector<std::string> enzyme_commission_numbers_;    /*!< List of enzyme commission numbers >*/
            bool is_engineered_;                                    /*!< Indicates that the molecule is engineered or not >*/
            bool has_mutation_;                                     /*!< Indicates that the molecule has been mutated or not >*/
            std::string comments_;                                  /*!< Comments >*/
    };
}
#endif // PDBCOMPOUNDSPECIFICATION_HPP
