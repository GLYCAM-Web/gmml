#ifndef PDBMATRIXN_HPP
#define PDBMATRIXN_HPP

#include <string>
#include <iostream>

#include "../../../includes/Geometry/coordinate.hpp"

namespace PdbFileSpace
{
    class PdbMatrixN
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            PdbMatrixN();
            PdbMatrixN(std::string& line);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            int GetN();
            int GetSerialNumber();
            Geometry::Coordinate GetTransformationVector();
            double GetV();
            int GetIGiven();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetRecordName(const std::string record_name);
            void SetN(int n);
            void SetSerialNumber(int serial_number);
            void SetTransformationVector(Geometry::Coordinate transfomration_vector);
            void SetV(double v);
            void SetIGiven(int i_given);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            int n_;
            int serial_number_;
            Geometry::Coordinate transfomration_vector_;
            double v_;
            int i_given_;
    };
}

#endif // PDBMATRIXN_HPP
