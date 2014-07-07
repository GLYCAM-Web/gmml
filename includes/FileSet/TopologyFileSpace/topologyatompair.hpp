#ifndef TOPOLOGYATOMPAIR_HPP
#define TOPOLOGYATOMPAIR_HPP

#include <string>
#include <iostream>
#include <map>

namespace TopologyFileSpace
{
    class TopologyAtomPair
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            typedef std::map<std::string, double> TopologyCoefficientMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyAtomPair();

            TopologyAtomPair(std::string pair_type, double coefficient_a, double coefficient_b);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////           
            /*! \fn
              * An accessor function in order to access to the pair type
              * @return pair_type_ attribute of the current object of this class
              */
            std::string GetPairType();
            /*! \fn
              * An accessor function in order to access to the coefficient a
              * @return coefficient_a_ attribute of the current object of this class
              */
            double GetCoefficientA();
            /*! \fn
              * An accessor function in order to access to the coefficient b
              * @return coefficient_b_ attribute of the current object of this class
              */
            double GetCoefficientB();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the pair type of the current object
              * Set the pair_type_ attribute of the current topology atom type
              * @param pair_type The pair type attribute of the current object
              */
            void SetPairType(std::string pair_type);
            /*! \fn
              * A mutator function in order to set the coefficient a of the current object
              * Set the coefficient_a_ attribute of the current topology atom type
              * @param coefficient_a The coefficient a attribute of the current object
              */
            void SetCoefficientA(double coefficient_a);
            /*! \fn
              * A mutator function in order to set the coefficient b of the current object
              * Set the coefficient_b_ attribute of the current topology atom type
              * @param coefficient_b The coefficient b attribute of the current object
              */
            void SetCoefficientB(double coefficient_b);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the atom type contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string pair_type_;
            double coefficient_a_;
            double coefficient_b_;

    };
}

#endif // TOPOLOGYATOMPAIR_HPP
