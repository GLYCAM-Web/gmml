#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbqtFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtModelCard::PdbqtModelCard() : record_name_("MODEL"){}

PdbqtModelCard::PdbqtModelCard(stringstream &stream_block)
{
    string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
    {
        if(line.find("MODEL") != string::npos)
        {
            if(!is_record_name_set){
                record_name_ = line.substr(0,6);
                Trim(record_name_);
                is_record_name_set=true;
            }
            stringstream model_block;
            while(line.find("MODEL") != string::npos || line.find("COMPND") != string::npos || line.find("REMARK") != string::npos
                    || line.find("ROOT") != string::npos || line.find("ATOM") != string::npos || line.find("ENDROOT") != string::npos
                    || line.find("BRANCH") != string::npos || line.find("ENDBRANCH") != string::npos || line.find("HETATM") != string::npos
                    || line.find("TORSDOF") != string::npos || line.find("ENDMDL") != string::npos)
            {
                model_block << line << endl;
                if(line.find("ENDMDL") != string::npos)
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
                Trim(record_name_);
                is_record_name_set = true;
            }
            stringstream model_block;
            while(line.find("MODEL") != string::npos || line.find("COMPND") != string::npos || line.find("REMARK") != string::npos
                    || line.find("ROOT") != string::npos || line.find("ATOM") != string::npos || line.find("ENDROOT") != string::npos
                    || line.find("BRANCH") != string::npos || line.find("ENDBRANCH") != string::npos || line.find("HETATM") != string::npos
                    || line.find("TORSDOF") != string::npos || line.find("ENDMDL") != string::npos)
            {
                model_block << line << endl;
                if(line.find("ENDMDL") != string::npos)
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

string PdbqtModelCard::GetRecordName(){
    return record_name_;
}

PdbqtModelCard::PdbqtModelMap PdbqtModelCard::GetModels(){
    return models_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbqtModelCard::SetRecordName(const string record_name){
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
void PdbqtModelCard::Print(ostream &out)
{
    out << "Record Name: " << record_name_ << endl <<
           "================= Models =================" << endl;
    for(PdbqtModelCard::PdbqtModelMap::iterator it = models_.begin(); it != models_.end(); it++)
    {
        out << "Model Serial Number: ";
        if((it)->first != iNotSet)
            out << (it)->first << endl;
        else
            out << " " << endl;
        (it)->second->Print();
        out << endl;
    }
}

