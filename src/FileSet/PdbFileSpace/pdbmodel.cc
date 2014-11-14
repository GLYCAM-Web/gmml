#include "../../../includes/FileSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbmodel.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbModel::PdbModel() {}

PdbModel::PdbModel(stringstream &model_block)
{
    string line;
    stringstream residue_set_block;
    getline(model_block, line);
    if(line.find("MODEL") != string::npos)
    {
        if(line.substr(10, 4) == "    ")
            model_serial_number_ = iNotSet;
        else
            model_serial_number_ = ConvertString<int>(line.substr(10,4));
        getline(model_block,line);
        string temp = line;
        while(line.find("ATOM") != string::npos || line.find("ANISOU") != string::npos
              || line.find("TER") != string::npos || line.find("HETATM") != string::npos)
        {
            residue_set_block << line << endl;
            getline(model_block, line);
            temp = line;
        }
        model_residue_set_ = new PdbModelResidueSet(residue_set_block);
    }
    else
    {
        model_serial_number_ = 1;
        string temp = line;
        while(line.find("ATOM") != string::npos || line.find("ANISOU") != string::npos
              || line.find("TER") != string::npos || line.find("HETATM") != string::npos)
        {
            residue_set_block << line << endl;
            getline(model_block, line);
            temp = line;
        }
        model_residue_set_ = new PdbModelResidueSet(residue_set_block);
    }

}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

int PdbModel::GetModelSerialNumber(){
    return model_serial_number_;
}

PdbModelResidueSet* PdbModel::GetModelResidueSet(){
    return model_residue_set_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbModel::SetModelSerialNumber(int model_serial_number){
    model_serial_number_ = model_serial_number;
}

void PdbModel::SetModelResidueSet(PdbModelResidueSet* model_residue_set){
    model_residue_set_ = model_residue_set;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbModel::Print(ostream &out)
{
    out << "Model Serial Number: ";
    if(model_serial_number_ != iNotSet)
        out << model_serial_number_;
    else
        out << " ";
    out << endl
        << "====================== Residue Set =====================" << endl;
    model_residue_set_->Print(out);
    out << endl;
}
