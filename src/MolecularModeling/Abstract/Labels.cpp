#include <algorithm>
#include <iostream>
#include <regex>
#include "includes/MolecularModeling/Abstract/Labels.hpp"

using Abstract::Labels;

std::string Labels::GetLabel()
{ // When requesting singular Label, likely want the latest one.
	if (labels_.empty())
		return "";
	return labels_.back();
}

std::string Labels::FindLabelContaining(const std::string query)
{ // When requesting singular Label, likely want the latest one.
    std::regex regexQuery(query, std::regex_constants::ECMAScript);
	for (auto &label : this->GetLabels())
	{
		if (std::regex_search(label, regexQuery))
		{
			return label;
		}
	}
	std::cout << "Found nothing.\n";
	return "";
}


bool Labels::CompareLabels(const std::vector<std::string> otherLabels)
{ // If any label here matches any in other label, return true
	for(auto &otherLabel : otherLabels)
	{
		auto labels = this->GetLabels();
		if (std::find(labels.begin(), labels.end(), otherLabel ) != labels.end() )
		{
			//std::cout << "Node labels match for " << this->GetLabel() << " & " << otherLabel << "\n";
			return true;
		}
		//std::cout << "Node labels DONT match for " << this->GetLabel() << " & " << otherLabel << "\n";
	}
	return false;
}