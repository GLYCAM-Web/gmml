// Created by: Delaram Rahbarinia
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBORIGINXN_HPP
#define PDBORIGINXN_HPP

#include <string>
#include <sstream>
#include <iostream>

#include "../../Geometry/coordinate.hpp"

namespace PdbFileSpace
{
    class PdbOriginXn
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbOriginXn();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbOriginXn(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the origin X N attribute in a origin X N
              * @return origin_x_n_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the N attribute in a origin X N
              * @return n_ attribute of the current object of this class
              */
            int GetN();
            /*! \fn
              * An accessor function in order to access to the origin N attribute in a origin X N
              * @return origin_n_ attribute of the current object of this class
              */
            Geometry::Coordinate GetOrigin();
            /*! \fn
              * An accessor function in order to access to the T attribute in a origin X N
              * @return t_ attribute of the current object of this class
              */
            double GetT();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current origin X N
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the N attribute of the current object
              * Set the n_ attribute of the current origin X N
              * @param n The N of the current object
              */
            void SetN(int n);
            /*! \fn
              * A mutator function in order to set the origin of the current object
              * Set the origin attribute of the current origin X N
              * @param origin The origin of the current object
              */
            void SetOrigin(Geometry::Coordinate origin);
            /*! \fn
              * A mutator function in order to set the T attribute of the current object
              * Set the t_ attribute of the current origin X N
              * @param t The T attribute of the current object
              */
            void SetT(double t);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the origin xn contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            int n_;
            Geometry::Coordinate origin_;
            double t_;
    };
}

#endif // PDBORIGINXN_HPP
