#ifndef TOPOLOGYANGLETYPE_HPP
#define TOPOLOGYANGLETYPE_HPP

#include <string>
#include <iostream>
#include <vector>

namespace TopologyFileSpace
{
    class TopologyAngleType
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyAngleType();
            /*! \fn
              * Constructor with required parameters
              * @param index Index of appearing of the angle type in a topology file
              * @param force_constant Force constant value of the angle appears in a topology file
              * @param equilibrium_value Equilibrium value of the angle appears in a topology file
              */
            TopologyAngleType(int index, double force_constant, double equilibrium_value);
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the angle type index
              * @return index_ attribute of the current object of this class
              */
            int GetIndex();
            /*! \fn
              * An accessor function in order to access to the angle type force constant
              * @return force_constant_ attribute of the current object of this class
              */
            double GetForceConstant();
            /*! \fn
              * An accessor function in order to access to the angle type equilibrium value
              * @return equilibrium_value_ attribute of the current object of this class
              */
            double GetEquilibriumValue();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the index of the current object
              * Set the index_ attribute of the current topology angle type
              * @param index The index attribute of the current object
              */
            void SetIndex(int index);
            /*! \fn
              * A mutator function in order to set the force constant of the current object
              * Set the force_constant_ attribute of the current topology angle type
              * @param force_constant The index attribute of the current object
              */
            void SetForceConstant(double force_constant);
            /*! \fn
              * A mutator function in order to set the force constant of the current object
              * Set the force_constant_ attribute of the current topology angle type
              * @param force_constant The force constant attribute of the current object
              */
            void SetEquilibriumValue(double equilibrium_value);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the angle type contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            int index_;                 /*!< Index of appearing of an angle type in a topology file >*/
            double force_constant_;     /*!< Force constant value of the angle appears in a topology file >*/
            double equilibrium_value_;  /*!< Equilibrium value of the angle appears in a topology file >*/

    };
}

#endif // TOPOLOGYANGLETYPE_HPP
