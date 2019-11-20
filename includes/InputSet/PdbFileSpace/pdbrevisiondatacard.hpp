// Created by: Dave Montgomery

#ifndef PDBREVISIONDATACARD_HPP
#define PDBREVISIONDATACARD_HPP

#include <string>
#include <iostream>


namespace PdbFileSpace
{
    class PdbRevisionDataCard
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbRevisionDataCard();
            /*! \fn
              * Constructor with required parameters
              * @param line
              */
            PdbRevisionDataCard(std::string& line);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a revision data
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the modification number in a revision data
              * @return mod_num_ attribute of the current object of this class
              */
            int GetModificationNumber();
            /*! \fn
              * An accessor function in order to access to the modification date in a revision data
              * @return mod_date_ attribute of the current object of this class
              */
            std::string GetModificationDate();
            /*! \fn
              * An accessor function in order to access to the modification id in a revision data
              * @return mod_id_ attribute of the current object of this class
              */
            std::string GetModificationID();
            /*! \fn
              * An accessor function in order to access to the modification type attribute in a revision data
              * @return mod_type_ attribute of the current object of this class
              */
            int GetModificationType();
            /*! \fn
              * An accessor function in order to access to the modifcation detail attribute in a revision data
              * @return mod_detail_ attribute of the current object of this class
              */
            std::string GetModificationDetails();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current revision data
              * @param record_name The record name attribute of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the Modiciation Number attribute of the current object
              * Set the mod_num_ attribute of the current revision data
              * @param mod_num The Token attribute of the current object
              */
            void SetModificationNumber(int mod_num_);
            /*! \fn
              * A mutator function in order to set the modification date attribute of the current object
              * Set the mod_date_ attribute of the current revision data
              * @param mod_date The Value attribute of the current object
              */
            void SetModificationDate(std::string mod_date);
            /*! \fn
              * A mutator function in order to set the modification ID attribute of the current object
              * Set the mod_id_ attribute of the current revision data
              * @param mod_id The Value attribute of the current object
              */
            void SetModificationID(std::string mod_id);
            /*! \fn
              * A mutator function in order to set the modification type attribute of the current object
              * Set the mod_type_ attribute of the current revision data
              * @param mod_type The Value attribute of the current object
              */
            void SetModificationType(int mod_type);
            /*! \fn
              * A mutator function in order to set the modification detail attribute of the current object
              * Set the mod_detail_ attribute of the current revision data
              * @param mod_detail The Value attribute of the current object
              */
            void SetModificationDetails(std::string mod_detail);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the revision data contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;                       /*!< Record name of revision data card which is the first column of each line of the card >*/
            int mod_num_;                                   /*!< Number of revision (1 for initial entry) >*/
            std::string mod_date_;                          /*!< Date of revision >*/
            std::string mod_id_;                            /*!< PDB ID of revision >*/
            int mod_type_;                                  /*!< Type of revision (0 os initial entry and 1 as revised) >*/
            std::string mod_detail_;                        /*!< Details of the revision >*/
    };
}

#endif // PDBREVISIONDATACARD_HPP
