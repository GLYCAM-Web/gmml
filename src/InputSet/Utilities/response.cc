#include <string>
#include <vector>
#include <utility>
#include "../../../includes/InputSet/Utilities/response.hpp"
using InputOutput::Response;

//Constructor
Response::Response()
{
}

//Accessor
std::string Response::GetServiceType()
{
    return this->service_type;
}

std::vector<Glycan::Note*> Response::GetNotices()
{
    return this->notices;
}

Response::OutputDict Response::GetOutputs()
{
    return this->output_dict;
}

std::string Response::GetRequestID()
{
    return this->request_id;
}

Response::Tags Response::GetTags()
{
    return this->tags;
}

//Mutator
void Response::SetServiceType(std::string new_serivce_type)
{
    this->service_type = new_serivce_type;
}

void Response::AddNotice(Glycan::Note* new_notice)
{
    this->notices.push_back(new_notice);
}

void Response::AddOutput(std::pair <std::string, std::string> new_output)
{
    this->output_dict.push_back(new_output);
}

void Response::GetRequestID(std::string new_request_id)
{
    this->request_id = request_id;
}

void Response::AddTag(std::pair<std::string, std::string> new_tag)
{
    this->tags.push_back(new_tag);
}



