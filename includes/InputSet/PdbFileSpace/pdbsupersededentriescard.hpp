// Created by: Dave Montgomery

#ifndef PDBSUPERSEDEDENTRIESCARD_HPP
#define PDBSUPERSEDEDENTRIESCARD_HPP

#include <string>
#include <iostream>


namespace PdbFileSpace
{
    class PdbSupersededEntriesCard
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbSupersededEntriesCard();
            /*! \fn
              * Constructor with required parameters
              * @param line
              */
            PdbSupersededEntriesCard(std::string& line);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a superseded entries
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the superseded date in a superseded entries
              * @return superseded_date_ attribute of the current object of this class
              */
            std::string GetSupersededDate();
            /*! \fn
              * An accessor function in order to access to the superseded id in a superseded entries
              * @return superseded_id_ attribute of the current object of this class
              */
            std::string GetSupersededID();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current superseded entries
              * @param record_name The record name attribute of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the superseded date attribute of the current object
              * Set the superseded_date_ attribute of the current superseded entries
              * @param superseded_date_ The Value attribute of the current object
              */
            void SetSupersededDate(std::string superseded_date);
            /*! \fn
              * A mutator function in order to set the superseded ID attribute of the current object
              * Set the superseded_id_ attribute of the current superseded entries
              * @param superseded_id_ The Value attribute of the current object
              */
            void SetSupersededID(std::string superseded_id);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the superseded entries contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;                       /*!< Record name of superseded entries card which is the first column of each line of the card >*/
            std::string superseded_date_;                    /*!< Date of superseding >*/
            std::string superseded_id_;                     /*!< PDB ID of entries made obsolete by this one >*/
    };
}

#endif // PDBSUPERSEDEDENTRIESCARD_HPP
