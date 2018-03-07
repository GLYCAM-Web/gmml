#include "../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbModelCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbModelCard::PdbModelCard() {}

PdbModelCard::PdbModelCard(std::stringstream &model_block)
{
    std::string line;
    std::stringstream residue_set_block;
    getline(model_block, line);
    if(line.find("MODEL") != std::string::npos)
    {
        if(line.substr(10, 4) == "    ")
            model_serial_number_ = gmml::iNotSet;
        else
            model_serial_number_ = gmml::ConvertString<int>(line.substr(10,4));
        getline(model_block,line);
        std::string temp = line;
        while(line.find("ATOM") != std::string::npos || line.find("ANISOU") != std::string::npos
              || line.find("TER") != std::string::npos || line.find("HETATM") != std::string::npos)
        {
            residue_set_block << line << std::endl;
            getline(model_block, line);
            temp = line;
        }
        model_residue_set_ = new PdbFileSpace::PdbModelResidueSet(residue_set_block);
    }
    else
    {
        model_serial_number_ = 1;
        std::string temp = line;
        while(line.find("ATOM") != std::string::npos || line.find("ANISOU") != std::string::npos
              || line.find("TER") != std::string::npos || line.find("HETATM") != std::string::npos)
        {
            residue_set_block << line << std::endl;
            getline(model_block, line);
            temp = line;
        }
//        cout << residue_set_block.str() << std::endl;
        model_residue_set_ = new PdbFileSpace::PdbModelResidueSet(residue_set_block);
    }

}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

int PdbModelCard::GetModelSerialNumber(){
    return model_serial_number_;
}

PdbFileSpace::PdbModelResidueSet* PdbModelCard::GetModelResidueSet(){
    return model_residue_set_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void PdbModelCard::SetModelSerialNumber(int model_serial_number){
    model_serial_number_ = model_serial_number;
}

void PdbModelCard::SetModelResidueSet(PdbFileSpace::PdbModelResidueSet* model_residue_set){
    model_residue_set_ = model_residue_set;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbModelCard::Print(std::ostream &out)
{
    out << "Model Serial Number: ";
    if(model_serial_number_ != gmml::iNotSet)
        out << model_serial_number_;
    else
        out << " ";
    out << std::endl
        << "====================== Residue Set =====================" << std::endl;
    model_residue_set_->Print(out);
    out << std::endl;
}
