#ifndef PDBCRYSTALLOGRAPHICCARD_HPP
#define PDBCRYSTALLOGRAPHICCARD_HPP

#include <string>
#include <sstream>

namespace PdbFileSpace
{
    class PdbCrystallographicCard
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            PdbCrystallographicCard();
            PdbCrystallographicCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            double GetA();
            double GetB();
            double GetC();
            double GetAlpha();
            double GetBeta();
            double GetGamma();
            std::string GetSpaceGroup();
            int GetZValue();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetRecordName(const std::string record_name);
            void SetA(double a);
            void SetB(double b);
            void SetC(double c);
            void SetAlpha(double alpha);
            void SetBeta(double beta);
            void SetGamma(double gamma);
            void SetSpaceGroup(const std::string space_group);
            void SetZValue(int z_value);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////



        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            double a_;
            double b_;
            double c_;
            double alpha_;
            double beta_;
            double gamma_;
            std::string space_group_;
            int z_value_;

    };
}


#endif // PDBCRYSTALLOGRAPHICCARD_HPP
