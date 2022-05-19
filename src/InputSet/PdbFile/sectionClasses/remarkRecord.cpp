#include "includes/InputSet/PdbFile/SectionClasses/remarkRecord.hpp"
#include "includes/utils.hpp"
#include "includes/common.hpp"
#include "includes/CodeUtils/logging.hpp"

using pdb::RemarkRecord;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
RemarkRecord::RemarkRecord()
{
    this->SetResolution(gmml::dNotSet);
    this->SetBFactor(gmml::dNotSet);
}
RemarkRecord::RemarkRecord(std::stringstream &stream_block)
{
    this->SetResolution(gmml::dNotSet);
    this->SetBFactor(gmml::dNotSet);
    std::string line;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(line.find("REMARK") != std::string::npos)
        {
            if (line.find("2 RESOLUTION.") != std::string::npos)
            {
                std::string tmp_resolution = line.substr(23,7);
                gmml::Trim(tmp_resolution);
                try
                {
                    this->SetResolution( std::stof( tmp_resolution ) );
                }
                catch (const std::invalid_argument& error)
                {
                    gmml::log(__LINE__, __FILE__, gmml::ERR, "RESOLUTION is not a valid float value. Value:\t" + tmp_resolution);
                }
            }
            if (line.find("MEAN B VALUE") != std::string::npos)
            {
                int start = line.find(":") + 1;
                std::string tmp_b_factor = line.substr(start,80-start);
                gmml::Trim( tmp_b_factor );
                try
                {
                    this->SetBFactor( std::stof( tmp_b_factor ) );
                }
                catch(const std::invalid_argument& error)
                {
                    gmml::log(__LINE__, __FILE__, gmml::ERR, "MEAN B VALUE is not a valid float value. Value:\t" + tmp_b_factor);
                }
            }
        }
        getline(stream_block,line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void RemarkRecord::SetResolution(const float resolution)
{
    this->resolution_ = resolution;
}

void RemarkRecord::SetBFactor(const float b_factor)
{
    this->b_factor_ = b_factor;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void RemarkRecord::Print(std::ostream &out) const
{
    out << "Resolution: " << this->GetResolution() << ". BFactor: " << this->GetBFactor() << "\n";
}

void RemarkRecord::Write(std::ostream& stream) const
{
    stream << "REMARK   2\n";
    stream << "REMARK   2 RESOLUTION.  " << this->GetResolution() << " ANGSTROMS.\n";
    stream << "REMARK   3\n";
    stream << "REMARK   3   MEAN B VALUE      (OVERALL, A**2) : " << this->GetBFactor() << "\n";
}
