#ifndef PDBSCALEN_HPP
#define PDBSCALEN_HPP

#include <string>
#include <sstream>
#include "../../../includes/Geometry/coordinate.hpp"

namespace PdbFileSpace
{
    class PdbScaleN
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            PdbScaleN();
            PdbScaleN(std::istringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            int GetN();
            Geometry::Coordinate GetScaleVector();
            double GetU();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetRecordName(const std::string record_name);
            void SetN(int n);
            void SetScaleVector(Geometry::Coordinate scale_vector);
            void SetU(double u);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////



        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            int n_;
            Geometry::Coordinate scale_vector_;
            double u_;
    };
}

#endif // PDBSCALEN_HPP
