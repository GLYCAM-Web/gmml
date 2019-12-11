#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbqtFileSpace::PdbqtModelCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtModelCard::PdbqtModelCard() : record_name_("MODEL"){}

PdbqtModelCard::PdbqtModelCard(std::ifstream &stream_block)
{
    std::string line;

    record_name_ = "MODEL";

    while(getline(stream_block, line)){
	if (line.find("MODEL") != std::string::npos || line.find("ROOT") != std::string::npos || line.find("BRANCH") != std::string::npos ||
	    line.find("ATOM") != std::string::npos || line.find("HETATM") != std::string::npos){
	    int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
            stream_block.seekg(offset, stream_block.cur);//Go back one line
            PdbqtModel* pdbqt_model = new PdbqtModel(stream_block);
            models_[pdbqt_model->GetModelSerialNumber()] = pdbqt_model;
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
