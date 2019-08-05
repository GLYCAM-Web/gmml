#ifndef RESPONSE_HPP
#define RESPONSE_HPP

#include <string>
#include <vector>
#include <utility>

#include "../../../includes/Glycan/note.hpp"

namespace InputOutput
{
    class Response{
	public:
	    typedef std::vector<std::pair<std::string, std::string> > OutputDict;
	    typedef std::vector<std::pair<std::string, std::string> > Tags;
	    enum NoticeType {NOTE, WARNING, ERROR, EXIT};

	    //Constructor
	    Response();

	    //Accessor
	    std::string GetServiceType();
	    std::vector<Glycan::Note*> GetNotices();
	    OutputDict GetOutputs();
	    std::string GetRequestID();
	    Tags GetTags();

	    //Mutator
	    void SetServiceType(std::string new_serivce_type);
	    void AddNotice(Glycan::Note* new_notice);
	    void AddOutput(std::pair <std::string, std::string> new_output);
	    void GetRequestID(std::string new_request_id);
	    void AddTag(std::pair<std::string, std::string> new_tag);

	private:
	    std::string service_type;
	    std::vector<Glycan::Note*> notices;
	    OutputDict output_dict;
	    std::string request_id;
	    Tags tags;
    };
}
#endif
