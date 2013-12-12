#ifndef PDBCOMPOUNDCARD_H
#define PDBCOMPOUNDCARD_H

#include <string>
#include <map>

namespace PdbFileSpace
{
    class PdbCompoundSpecification;

    class PdbCompoundCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            typedef std::map<std::string, PdbCompoundSpecification*> PdbCompoundSpecificationMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbCompoundCard();
            PdbCompoundCard(std::string record_name);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetRecordName();
            PdbCompoundSpecificationMap GetCompoundSpecifications();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetRecordName(std::string record_name);


            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            PdbCompoundSpecificationMap compound_specifications_;

    };
}

#endif // PDBCOMPOUNDCARD_H
