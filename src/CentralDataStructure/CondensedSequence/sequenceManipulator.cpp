#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Graph.hpp"
#include "includes/CodeUtils/logging.hpp"
#include <sstream>
#include <sys/stat.h> // for checking if file exists
#include <fstream>  // writing outputDotFile

using cdsCondensedSequence::SequenceManipulator;
using cdsCondensedSequence::ParsedResidue;

bool file_exists(const char *filename)
{
    struct stat buffer;
    return (stat (filename, &buffer) == 0);
}

std::string SequenceManipulator::ReorderSequence()
{	// Just doing the default by ascending link number for now.
	for (auto &residue : this->getResidues())
	{
	    ParsedResidue* brian = static_cast<ParsedResidue*>(residue);
	    brian->sortOutEdgesBySourceTObjectComparator();
	}
	return this->Print();
}

std::vector<ParsedResidue*> SequenceManipulator::GetParsedResiduesOrderedByConnectivity()
{
    std::vector<ParsedResidue*> rawResidues;
    // Go via Graph so order decided by connectivity, depth first traversal:
    glygraph::Graph<cds::Residue> sequenceGraph(this->GetTerminal());
    for(auto &node : sequenceGraph.getNodes())
    {
        rawResidues.push_back(static_cast<ParsedResidue*>(node->getDeriviedClass()));
    }
    return rawResidues;
}

void SequenceManipulator::LabelSequence()
{
	this->SetIndexByConnectivity();
	std::stringstream ss;
	for (auto &residue : this->GetParsedResiduesOrderedByConnectivity())
	{
		ss << residue->GetName() << "&Label=residue-" << residue->getIndex() << ";";
		residue->addLabel(ss.str());
		ss.str( std::string() ); ss.clear();  // Must do both of these to clear the stream
		for (auto &linkage : residue->getOutEdges())
		{
			ss << linkage->getLabel() << "&Label=link-" << linkage->getIndex() << ";";
			linkage->addLabel(ss.str());
			ss.str( std::string() ); ss.clear(); // Must do both of these to clear the stream
		}
	}
	return;
}

void SequenceManipulator::SetIndexByConnectivity()
{
	unsigned long long linkIndex = 0; // Convention to start form 0 for linkages.
	unsigned long long residueIndex = 1; // Convention to start from 1 for residues.
	for (auto &residue : this->GetParsedResiduesOrderedByConnectivity())
	{
		residue->setIndex(residueIndex);
        residue->setNumber(residueIndex); // ToDo temporary, switch to using number here. Keep index as a gmml internal thing, never shown to user.
		++residueIndex;
		for (auto &linkage : residue->getInEdges())
		{
			linkage->setIndex(linkIndex);
			++linkIndex;
		}
	}
	return;
}

std::string SequenceManipulator::Print(const bool withLabels)
{
	if (withLabels)
	{
		this->LabelSequence();
	}
	std::vector<std::string> output;
	int branchStackSize = 0;
	this->RecurvePrint(this->GetTerminal(), branchStackSize, output, withLabels);
	std::reverse(output.begin(), output.end()); // Reverse order, as it starts from terminal.
	std::stringstream ss;
	for (auto &label : output)
	{
		ss << label;
	}
	gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
	return ss.str();
}

void SequenceManipulator::RecurvePrint(ParsedResidue* currentResidue, int& branchStackSize, std::vector<std::string>& output, const bool withLabels)
{
	auto neighbors = currentResidue->GetChildren();
	size_t numberOfNeighbors = neighbors.size();
	// Derivatives. E.g. 2S,3Me in DManp[2S,3Me]a1-6DManpa1-OH
	std::string outputResidueString = currentResidue->GetName(withLabels);
	std::vector<std::string> derivatives;
	for (auto &neighbor : neighbors)
	{
		if (neighbor->GetType() == Abstract::ResidueType::Derivative)
		{
			--numberOfNeighbors;
			derivatives.push_back(neighbor->GetLinkageName(withLabels) + neighbor->GetName(withLabels));
			derivatives.push_back(",");
		}
	}
	if (!derivatives.empty())
	{
		derivatives.pop_back(); // Remove the last ","
		std::reverse(derivatives.begin(), derivatives.end()); // order should be 2S,6S, not 6S,2S.
		outputResidueString += "[";
		for (auto &derivative : derivatives)
		{
			outputResidueString += derivative;
		}
		outputResidueString += "]";
	}
	// Output
	outputResidueString += currentResidue->GetLinkageName(withLabels);
	output.push_back(outputResidueString);
	// End of a branch check
	if (numberOfNeighbors == 0 && branchStackSize > 0)
	{
		output.push_back("["); 
		--branchStackSize;
	}
	size_t loopCount = 0;
	for (auto &neighbor : neighbors)
	{
		if (neighbor->GetType() != Abstract::ResidueType::Derivative)
		{
			++loopCount;
			if (loopCount < numberOfNeighbors)
			{
				output.push_back("]");
				++branchStackSize;
			}
			this->RecurvePrint(neighbor, branchStackSize, output, withLabels);
		}
	}
	return;
}

std::string SequenceManipulator::PrintGraphViz(GraphVizDotConfig &configs)
{
	this->SetIndexByConnectivity();
	std::stringstream ss;
	ss << "graph G {graph [splines=false forcelabels=true  dpi=" << configs.dpi_ << "];\n";
	ss << "node [ shape=\"none\" fontname=Helvetica labelfontsize=12 forcelabels=\"true\";\n";
	ss << "label=\"none\" size=50 fixedsize=\"true\" scale=\"true\"];\n";
	ss << "edge [labelfontsize=12 fontname=Helvetica labeldistance=1.2 labelangle = 320.0];\n";
	ss << "rankdir=LR nodesep=\"0.05\" ranksep=\"0.8\";\n";
	for (auto &residue : this->GetParsedResiduesOrderedByConnectivity())
	{
		if (residue->GetType() != Abstract::ResidueType::Derivative)
		{
			ss << this->GetGraphVizLineForResidue(*residue, configs) << "\n";
		}
	}
	ss << "}\n";
	// Open and overwrite.
	std::ofstream outputDotFile(configs.file_name_, std::ios::trunc);
	outputDotFile << ss.str();
	outputDotFile.close();
	gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
	return ss.str();
}

std::string SequenceManipulator::GetGraphVizLineForResidue(ParsedResidue &residue, GraphVizDotConfig &configs)
{
	std::stringstream logss;
    logss << "Getting GraphVizLine for " << residue.GetName() << "\n"; 
    gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    logss.str( std::string() ); logss.clear();  // Must do both of these to clear the stream;
    std::stringstream ss;
    ss << residue.getIndex() << " [";
    // Aglycone
    if (residue.GetType() == Abstract::ResidueType::Aglycone)
    {
        ss << "shape=box label=\"" << residue.GetMonosaccharideName() << "\"]";
        return ss.str();
    }
    // Sugar
    std::string label = "";
    std::string imageFile = configs.svg_directory_path_ + residue.GetMonosaccharideName() + ".svg";
    logss << "Searching for image: " << imageFile << "\n";
    if(file_exists(imageFile.c_str()))
    {
        logss << "FOUND IT\n";
        (residue.GetRingType() == "f") ? label = "f" : label = "";
        ss << "label=\"" << label << "\" height=\"0.7\" image=\"" << imageFile << "\"];\n";
    }
    else
    {
    	logss << "Not image available, using circle\n";
        ss << "shape=circle height=\"0.7\" label=\"" << residue.GetMonosaccharideName() << "\"];\n";
    }
    // Derivatives
    std::string derivativeStr = "";
    for (auto &childLink : residue.GetChildren())
    {
        if (childLink->GetType() == Abstract::ResidueType::Derivative)
        {
            derivativeStr += childLink->GetLinkageName() + childLink->GetName() + " ";
        }
    }
    if (! derivativeStr.empty())
    {
        ss << "\n" << "b" << residue.getIndex();
        ss << "[ shape=\"plaintext\",fontsize=\"12\",forcelabels=\"true\"; height = \"0.3\"; labelloc = b;  label=\""; 
        ss << derivativeStr << "\"];\n";
        ss << "{ rank=\"same\"; b" << residue.getIndex() << " " << residue.getIndex() << "};\n";
        ss << "{nodesep=\"0.2\";b" << residue.getIndex() << ";" << residue.getIndex() << "};\n";
        ss << "b" << residue.getIndex() << "--" << residue.getIndex() << " [style=invis];\n";
    }
    // Linkage
    for (auto &parent : residue.GetParents())
    { // There is either 1 or 0 parents, this covers both cases.
        ss << residue.getIndex() << "--" << parent->getIndex() << "[label=\"";
        if (configs.show_config_labels_)
        {
        	ss << residue.GetConfiguration();
        }
        if (configs.show_position_labels_)
        {
        	ss << residue.GetLinkage();
        }
        ss << "\"];\n";
        if (configs.show_edge_labels_)
        {
	        for (auto &linkage : residue.getOutEdges())
	        {
	            ss << residue.getIndex() << "--" << parent->getIndex();
	            ss << "[taillabel=< <B>" << linkage->getIndex() << "</B>>, ";
	            ss << "labelfontsize = 14, labeldistance = 2.0, labelangle = -35";
	            ss << "];\n";
	        }
    	}
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    return ss.str();
}
