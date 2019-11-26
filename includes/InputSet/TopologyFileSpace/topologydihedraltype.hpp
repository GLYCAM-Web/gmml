#ifndef TOPOLOGYDIHEDRALTYPE_HPP
#define TOPOLOGYDIHEDRALTYPE_HPP

#include <string>
#include <iostream>
#include <vector>

namespace TopologyFileSpace
{
    class TopologyDihedralType
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyDihedralType();
            /*! \fn
              * Constructor with required parameters
              * @param index Index of the dihedral in the topology file
              * @param force_constant Force constant value of the dihedral
              * @param periodicity Periodicity value of the dihedral
              * @param phase Phase value of the dihedral
              * @param scee SCEE coefficient of the dihedral
              * @param scnb SCNB coefficient of the dihedral
              */
            TopologyDihedralType(int index, double force_constant, double periodicity, double phase, double scee, double scnb);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the dihedral type index
              * @return index_ attribute of the current object of this class
              */
            int GetIndex();
            /*! \fn
              * An accessor function in order to access to the dihedral type periodicity
              * @return periodicity_ attribute of the current object of this class
              */
            double GetPeriodicity();
            /*! \fn
              * An accessor function in order to access to the dihedral type phase
              * @return phase_ attribute of the current object of this class
              */
            double GetPhase();
            /*! \fn
              * An accessor function in order to access to the dihedral scee
              * @return scee_ attribute of the current object of this class
              */
            double GetScee();
            /*! \fn
              * An accessor function in order to access to the dihedral type scnb
              * @return scnb_ attribute of the current object of this class
              */
            double GetScnb();
            /*! \fn
              * An accessor function in order to access to the dihedral type force constant
              * @return force_constant_ attribute of the current object of this class
              */
            double GetForceConstant();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the index of the current object
              * Set the index_ attribute of the current topology dihedral type
              * @param index The index attribute of the current object
              */
            void SetIndex(int index);
            /*! \fn
              * A mutator function in order to set the periodicity of the current object
              * Set the periodicity_ attribute of the current topology dihedral type
              * @param periodicity The periodicity attribute of the current object
              */
            void SetPeriodicity(double periodicity);
            /*! \fn
              * A mutator function in order to set the phase of the current object
              * Set the phase_ attribute of the current topology dihedral type
              * @param phase The phase attribute of the current object
              */
            void SetPhase(double phase);
            /*! \fn
              * A mutator function in order to set the scee of the current object
              * Set the scee_ attribute of the current topology dihedral type
              * @param scee The scee attribute of the current object
              */
            void SetScee(double scee);
            /*! \fn
              * A mutator function in order to set the scnb of the current object
              * Set the scnb_ attribute of the current topology dihedral type
              * @param scnb The scnb attribute of the current object
              */
            void SetScnb(double scnb);
            /*! \fn
              * A mutator function in order to set the force constant of the current object
              * Set the force_constant_ attribute of the current topology dihedral type
              * @param force_constant The force constant attribute of the current object
              */
            void SetForceConstant(double force_constant);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the dihedral type contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            int index_;                     /*!< Index of dihedral type in a topology file >*/
            double periodicity_;            /*!< Periodicity of a dihedral >*/
            double phase_;                  /*!< Phase of a dihedral >*/
            double scee_;                   /*!< SCEE coefficient of a dihedral >*/
            double scnb_;                   /*!< SCNB coefficient of a dihedral >*/
            double force_constant_;         /*!< Force constant of dihedral >*/

    };
}

#endif // TOPOLOGYANGLETYPE_HPP
