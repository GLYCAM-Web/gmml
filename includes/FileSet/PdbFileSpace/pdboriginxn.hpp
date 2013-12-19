#ifndef PDBORIGINXN_HPP
#define PDBORIGINXN_HPP

#include <string>
#include "../../../includes/Geometry/coordinate.hpp"

namespace PdbFileSpace
{
    class PdbOriginXn
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            PdbOriginXn();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            int GetN();
            Geometry::Coordinate GetOrigin();
            double GetT();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetRecordName(const std::string record_name);
            void SetN(int n);
            void SetOrigin(Geometry::Coordinate origin);
            void SetT(double t);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////



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
