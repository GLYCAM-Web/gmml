#ifndef RULESET_HPP
#define RULESET_HPP

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <string>

namespace gmml
{
    namespace MolecularMetadata
    {
	class RuleSet
	{
	    namespace pt = boost::property_tree;
	    RuleSet();
	    RuleSet(std::string metadata_json);

	    //Private
	    pt root;
	    
	};
    }
}
#endif
