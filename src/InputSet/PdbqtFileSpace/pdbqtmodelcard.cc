#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbqtFileSpace::PdbqtModelCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtModelCard::PdbqtModelCard() : record_name_("MODEL"){}

PdbqtModelCard::PdbqtModelCard(std::stringstream &stream_block)
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
            while(line.find("MODEL") != std::string::npos || line.find("COMPND") != std::string::npos || line.find("REMARK") != std::string::npos
                    || line.find("ROOT") != std::string::npos || line.find("ATOM") != std::string::npos || line.find("ENDROOT") != std::string::npos
                    || line.find("BRANCH") != std::string::npos || line.find("ENDBRANCH") != std::string::npos || line.find("HETATM") != std::string::npos
                    || line.find("TORSDOF") != std::string::npos || line.find("ENDMDL") != std::string::npos)
            {
                model_block << line << std::endl;
                if(line.find("ENDMDL") != std::string::npos)
                {
                    PdbqtModel* pdbqt_model = new PdbqtModel(model_block);
                    models_[pdbqt_model->GetModelSerialNumber()] = pdbqt_model;
                    model_block.str("");
                }
                getline(stream_block,line);
                temp = line;
            }
        }
        else
        {
            if(!is_record_name_set){
                record_name_ = "MODEL ";
                gmml::Trim(record_name_);
                is_record_name_set = true;
            }
            std::stringstream model_block;
            while(line.find("MODEL") != std::string::npos || line.find("COMPND") != std::string::npos || line.find("REMARK") != std::string::npos
                    || line.find("ROOT") != std::string::npos || line.find("ATOM") != std::string::npos || line.find("ENDROOT") != std::string::npos
                    || line.find("BRANCH") != std::string::npos || line.find("ENDBRANCH") != std::string::npos || line.find("HETATM") != std::string::npos
                    || line.find("TORSDOF") != std::string::npos || line.find("ENDMDL") != std::string::npos)
            {
                model_block << line << std::endl;
                if(line.find("ENDMDL") != std::string::npos)
                {
                    PdbqtModel* pdbqt_model = new PdbqtModel(model_block);
                    models_[pdbqt_model->GetModelSerialNumber()] = pdbqt_model;
                    model_block.str("");
                }
                getline(stream_block,line);
                temp = line;
            }
        }
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

std::string PdbqtModelCard::GetRecordName(){
    return record_name_;
}

PdbqtModelCard::PdbqtModelMap PdbqtModelCard::GetModels(){
    return models_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbqtModelCard::SetRecordName(const std::string record_name){
    record_name_ = record_name;
}

void PdbqtModelCard::SetModels(PdbqtModelMap models){
    models_ = models;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtModelCard::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "================= Models =================" << std::endl;
    for(PdbqtModelCard::PdbqtModelMap::iterator it = models_.begin(); it != models_.end(); it++)
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
