#ifndef GMML_INCLUDES_ABSTRACT_BUILDER_HPP
#define GMML_INCLUDES_ABSTRACT_BUILDER_HPP
// My idea is that once there is more than one thing, or the status class needs more functionality
// this will "have a" status that is a separate class, and this may have other things too.
// Update: the status idea is bad, classes that aren't ok should throw. Stop-gap for now, as I can't yet throw across SWIG and into gems.
#include <string>
namespace Abstract
{
	class absBuilder
	{
	public:
        //example types {ERROR, WARNING, INFO, OK}; ENUMS are annoying to convert for JSON, so strings it is.
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        absBuilder() {statusType_ = "OK";}
        absBuilder(std::string t, std::string s) : statusType_(t), statusMessage_(s) {};
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline bool IsStatusOk() {return this->GetStatusType() == "OK";}
        inline std::string GetStatusType() {return statusType_;}
        inline std::string GetStatusMessage() {return statusMessage_;}
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        inline void SetStatusType(std::string s) {statusType_ = s;}
        inline void SetStatusMessage(std::string s) {statusMessage_ = s;}
        inline void SetStatus(std::string type, std::string message) {this->SetStatusType(type); this->SetStatusMessage(message);}
    private:
        //////////////////////////////////////////////////////////
        //                       ATTRRIBUTES                    //
        //////////////////////////////////////////////////////////
        std::string statusType_;  // enum Type. See enum above.
        std::string statusMessage_;
	};
}
#endif
