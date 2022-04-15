#include <iostream>

#include "../../../includes/InputSet/PdbFileSpace/pdbcompoundsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbcompoundspecification.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbCompoundSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbCompoundSection::PdbCompoundSection() : record_name_("COMPND"){}
PdbCompoundSection::PdbCompoundSection(const std::string& record_name) : record_name_(record_name) {}

PdbCompoundSection::PdbCompoundSection(std::stringstream& stream_block)
{
    std::string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_record_name_set){
            record_name_ = line.substr(0,6);
            gmml::Trim(record_name_);
            is_record_name_set=true;
        }
        std::stringstream ss, specification_block;
        ss << line.substr(10,70);
        specification_block << ss.str() << std::endl;

        getline(stream_block, line);
        temp = line;
        while (line.find("MOL_ID") == std::string::npos && !gmml::Trim(temp).empty()){
            std::stringstream sss;
            sss << line.substr(10,70);

            specification_block << sss.str() << std::endl;
            getline(stream_block, line);
            temp = line;
        }
        PdbCompoundSpecification* compound_specification = new PdbCompoundSpecification(specification_block);
        compound_specifications_[compound_specification->GetMoleculeId()] = compound_specification;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbCompoundSection::GetRecordName()
{
    return record_name_;
}

PdbCompoundSection::PdbCompoundSpecificationMap PdbCompoundSection::GetCompoundSpecifications()
{
    return compound_specifications_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbCompoundSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbCompoundSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl;
    out << "============= Compound Specification =============" << std::endl;
    for(PdbCompoundSection::PdbCompoundSpecificationMap::iterator it = compound_specifications_.begin(); it != compound_specifications_.end(); it++)
    {
        out << "Molecule ID: " << (it)->first << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
}
