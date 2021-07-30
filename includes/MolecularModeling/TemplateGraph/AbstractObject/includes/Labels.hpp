#ifndef ABSTRACTOBJECT_INCLUDES_LABELS_HPP
#define ABSTRACTOBJECT_INCLUDES_LABELS_HPP

#include <regex>
#include <string>
#include <vector>
#include <algorithm>

//TODO: Pick a better namespace name for this, figure out what we should use instead of `unsigned int`
namespace abstrab
{
class Labels
{
public:
	/************************************************
	 *  CONSTRUCTORS/DESTRUCTORS
	 ***********************************************/
	inline Labels() :
			name_m(""), labels_m( { "" } )
	{
	}

	inline Labels(std::string name_t) :
			name_m(name_t), labels_m( { "" } )
	{
	}

	inline Labels(std::vector<std::string> labels_t) :
			name_m(""), labels_m(labels_t)
	{
	}

	inline Labels(std::string name_t, std::string label_t) :
			name_m(name_t), labels_m(
			{ label_t })
	{
	}

	inline Labels(std::string name_t, std::vector<std::string> labels_t) :
			name_m(name_t), labels_m(labels_t)
	{
	}

	inline ~Labels()
	{
	}

	//copy constructor
	inline Labels(const Labels &rhs) :
			name_m(rhs.name_m), labels_m(rhs.labels_m)
	{

	}
	//move constructor
	inline Labels(Labels &&rhs) :
			name_m(rhs.name_m), labels_m(rhs.labels_m)
	{
	}
	//copy assignment
	inline Labels& operator =(const Labels &rhs)
	{
		this->name_m = rhs.name_m;
		this->labels_m = rhs.labels_m;
		return *this;
	}
	//move assignment
	Labels& operator =(Labels &&rhs)
	{
		this->name_m = rhs.name_m;
		this->labels_m = rhs.labels_m;
		return *this;
	}

	/************************************************
	 *  GETTER/SETTER
	 ***********************************************/
	inline std::string getName() const
	{
		return this->name_m;
	}

	inline std::string getLabel() const
	{
		if (this->labels_m.empty())
			return "";
		return this->labels_m.back();
	}

	inline std::vector<std::string> getLabels() const
	{
		return this->labels_m;
	}

	inline void setName(std::string name_t)
	{
		this->name_m = name_t;
	}

	inline void setLabels(std::string label_t)
	{
		this->labels_m.clear();
		this->addLabel(label_t);
	}

	inline void setLabels(std::vector<std::string> labels_t)
	{
		this->labels_m = labels_t;
	}

	/************************************************
	 *  MUTATORS
	 ***********************************************/
	inline void addLabel(std::string label_t)
	{
		this->labels_m.push_back(label_t);
	}

	inline void addLabels(std::vector<std::string> labels)
	{
		this->labels_m.insert(this->labels_m.end(), labels.begin(),
				labels.end());
	}

	/************************************************
	 *  FUNCTIONS
	 ***********************************************/
	inline bool containsLabel(const std::string query)
	{
		for (std::string currLabel : this->labels_m)
		{
			//TODO: Why was regex used in your version?
			if (currLabel == query)
				return true;
		}
		return false;
	}
	bool compareLabels(const std::vector<std::string> otherLabels)
	{
		return this->compareLabels(otherLabels, 1);
	}
	inline bool compareLabels(const std::vector<std::string> otherLabels,
			unsigned int numMatches)
	{
		unsigned int matchCounter = 0;
		for (std::string currLabel : this->labels_m)
		{
			if (std::find(otherLabels.begin(), otherLabels.end(), currLabel)
					!= otherLabels.end())
			{
				matchCounter++;
				if (matchCounter == numMatches)
					return true;
			}
		}
		return false;
	}

	inline std::string findLabelContaining(const std::string query)
	{
		std::regex regexQuery(query, std::regex_constants::ECMAScript);
		for (auto &label : this->getLabels())
		{
			if (std::regex_search(label, regexQuery))
			{
				return label;
			}
		}
		//std::cout << "Found nothing.\n";
		return "";
	}
private:
	/************************************************
	 *  ATTRIBUTES
	 ***********************************************/
	std::string name_m;
	std::vector<std::string> labels_m;

};
}

#endif // LABELS_HPP

