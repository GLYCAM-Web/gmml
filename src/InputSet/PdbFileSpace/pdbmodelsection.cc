#include "../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelsection.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbModelSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbModelSection::PdbModelSection() : record_name_("MODEL") {}

PdbModelSection::PdbModelSection(std::stringstream &stream_block)
{
    std::string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(line.find("MODEL") != std::string::npos)
        {
            if(!is_record_name_set){
                record_name_ = line.substr(0,6);
                gmml::Trim(record_name_);
                is_record_name_set=true;
            }
            std::stringstream model_block;
            while(line.find("MODEL") != std::string::npos || line.find("ATOM") != std::string::npos || line.find("ANISOU") != std::string::npos
                    || line.find("TER") != std::string::npos || line.find("HETATM") != std::string::npos || line.find("ENDMDL") != std::string::npos)
            {
                model_block << line << std::endl;
                getline(stream_block,line);
                temp = line;
            }
            PdbModelCard* pdb_model = new PdbModelCard(model_block);
            models_[pdb_model->GetModelSerialNumber()] = pdb_model;
        }
        else
        {
            if(!is_record_name_set){
                record_name_ = "MODEL ";
                gmml::Trim(record_name_);
                is_record_name_set = true;
            }
            std::stringstream model_block;
            while(line.find("MODEL") != std::string::npos || line.find("ATOM") != std::string::npos || line.find("ANISOU") != std::string::npos
                    || line.find("TER") != std::string::npos || line.find("HETATM") != std::string::npos || line.find("ENDMDL") != std::string::npos)
            {
                model_block << line << std::endl;
                getline(stream_block,line);
                temp = line;
            }
//            cout << model_block.str() << std::endl;
            PdbModelCard* pdb_model = new PdbModelCard(model_block);
            models_[pdb_model->GetModelSerialNumber()] = pdb_model;
        }
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

std::string PdbModelSection::GetRecordName(){
    return record_name_;
}

PdbModelSection::PdbModelCardMap PdbModelSection::GetModels(){
    return models_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbModelSection::SetRecordName(const std::string record_name){
    record_name_ = record_name;
}

void PdbModelSection::SetModels(PdbModelCardMap models){
    models_ = models;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbModelSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "================= Models =================" << std::endl;
    for(PdbModelSection::PdbModelCardMap::iterator it = models_.begin(); it != models_.end(); it++)
    {
        out << "Model Serial Number: ";
        if((it)->first != gmml::iNotSet)
            out << (it)->first << std::endl;
        else
            out << " " << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
}
