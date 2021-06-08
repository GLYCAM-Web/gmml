#ifndef PDBQTMODEL_HPP
#define PDBQTMODEL_HPP

#include <string>
#include <sstream>
#include <iostream>
#include <vector>

namespace PdbqtFileSpace
{
    class PdbqtRemarkCard;
    class PdbqtTorsionalDoFCard;
    class PdbqtModelResidueSet;
    class PdbqtCompoundCard;
    class PdbqtModel
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * A vector of pdbqt remark cards
              */
            typedef std::vector<PdbqtRemarkCard*> RemarkCardVector;
            /*! \typedef
              * A vector of pdbqt torsional degree of freedom cards
              */
            typedef std::vector<PdbqtTorsionalDoFCard*> TorsionalDoFCardVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbqtModel();
            /*! \fn
              * Constructor with required parameters
              * @param model_block
              */
            PdbqtModel(std::ifstream& model_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the model serial number in a pdbqt model
              * @return model_serial_number_ attribute of the current object of this class
              */
            int GetModelSerialNumber();
            /*! \fn
              * An accessor function in order to access to the model compound card in a pdbqt model
              * @return model_compound_card_ attribute of the current object of this class
              */
            PdbqtCompoundCard* GetModelCompoundCard();
            /*! \fn
              * An accessor function in order to access to the model residue set in a pdbqt model
              * @return model_residue_set_ attribute of the current object of this class
              */
            PdbqtModelResidueSet* GetModelResidueSet();
            /*! \fn
              * An accessor function in order to access to the remark cards
              * @return remarks_ attribute of the current object of this class
              */
            RemarkCardVector GetRemarks();
            /*! \fn
              * An accessor function in order to access to the torsional degree of freedom cards
              * @return torsional_dof_cards_ attribute of the current object of this class
              */
            TorsionalDoFCardVector GetTorsionalDoFCards();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the model serial number of the current object
              * Set the model_serial_number_ attribute of the current pdbqt model
              * @param model_serial_number The model serial number attribute of the current object
              */
            void SetModelSerialNumber(int model_serial_number);
            /*! \fn
              * A mutator function in order to set the model compound card of the current object
              * Set the compound_card_ attribute of the current pdbqt model
              * @param model_compound_card The model compound card attribute of the current object
              */
            void SetModelCompundCard(PdbqtCompoundCard* model_compound_card);
            /*! \fn
              * A mutator function in order to set the model residue set of the current object
              * Set the model_residue_set_ attribute of the current pdbqt model
              * @param model_residue_set The model residue set attribute of the current object
              */
            void SetModelResidueSet(PdbqtModelResidueSet* model_residue_set);
            /*! \fn
              * A mutator function in order to set remark cards of the current object
              * Set the remarks_ attribute of the current pdbqt model
              * @param remarks The remarks attribute of the current object
              */
            void SetRemarks(RemarkCardVector remarks);
            /*! \fn
              * A function in order to add the remarks to the current object
              * Set the remakrs_ attribute of the current pdbqt model
              * @param remark The remark of the current object
              */
            void AddRemark(PdbqtRemarkCard* remark);
            /*! \fn
              * A mutator function in order to set torsional degree of freedom cards of the current object
              * Set the torsional_dof_cards_ attribute of the current pdbqt model
              * @param torsional_dof_cards The tosional degree of freedom attribute of the current object
              */
            void SetTorsionalDofCards(TorsionalDoFCardVector torsional_dof_cards);
            /*! \fn
              * A function in order to add the torsional degree of freedom card to the current object
              * Set the torsional_dof_cards_ attribute of the current pdbqt model
              * @param torsional_dof_card_ The remark of the current object
              */
            void AddTorsionalDoFCard(PdbqtTorsionalDoFCard* torsional_dof_cards);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the model contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            int model_serial_number_;                    /*!< Serial number of a model >*/
            PdbqtCompoundCard* model_compound_card_;/*!< Compund card of a model >*/
            PdbqtModelResidueSet* model_residue_set_;    /*!< Residue sets involving in a model >*/
            RemarkCardVector remarks_;                   /*!< Remark cards of a model >*/
            TorsionalDoFCardVector torsional_dof_cards_; /*!< Torsional degree of freedom cards of a model >*/
    };
}

#endif // PDBQTMODEL_HPP
