#ifndef LABELS_HPP
#define LABELS_HPP

#include <string>
#include <vector>

namespace Abstract
{
	class Labels
	{
	public:
		//////////////////////////////////////////////////////////
		//                       CONSTRUCTOR                    //
		//////////////////////////////////////////////////////////
		Labels() {this->AddLabel("");}
		Labels(std::string label) {this->AddLabel(label);}
		Labels(std::vector<std::string> labels) {this->SetLabels(labels);}
		//////////////////////////////////////////////////////////
		//                       ACCESSOR                       //
		//////////////////////////////////////////////////////////
		inline std::vector<std::string> GetLabels() {return labels_;}
		std::string GetLabel();
		//////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
		inline void SetLabels(std::vector<std::string> labels) {labels_ = labels;}
		inline void AddLabel(std::string label) {labels_.push_back(label);}
		//////////////////////////////////////////////////////////
        //                      FUNCTIONS                       //
        //////////////////////////////////////////////////////////
        std::string FindLabelContaining(const std::string query);
		bool CompareLabels(const std::vector<std::string> otherLabels);
	private:
		//////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
		std::vector<std::string> labels_;
	};
}
#endif // LABELS_HPP