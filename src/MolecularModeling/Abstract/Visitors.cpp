#include <algorithm>
#include <string>
#include <vector>
#include "includes/MolecularModeling/Abstract/Visitors.hpp"

using Abstract::Visitors;

bool Visitors::GetIsVisitedBy(std::string visitor)
{
	if (std::find(visitors_.begin(), visitors_.end(), visitor ) != visitors_.end() )
		return true;
	return false;
}

void Visitors::RemoveVisitor(std::string visitor)
{
	visitors_.erase(std::remove(visitors_.begin(), visitors_.end(), visitor), visitors_.end());
	return; 
}
