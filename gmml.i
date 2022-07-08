/* File: gmml.i */
%module gmml
%include <std_string.i>
%include <std_iostream.i>
%include<std_map.i>
%include<std_vector.i>

%{
#define SWIG_FILE_WITH_INIT
//#include "/usr/include/sql.h"
//#include "/usr/include/sqlext.h"

#include "includes/gmml.hpp"
#include "includes/common.hpp"
#include "includes/utils.hpp"
#include "includes/generictypedefs.hpp"
#include "includes/CodeUtils/codetests.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/InputSet/CoordinateFileSpace/coordinatefile.hpp"
#include "includes/InputSet/CoordinateFileSpace/coordinatefileprocessingexception.hpp"
#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/GeometryTopology/plane.hpp"
#include "includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "includes/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.hpp"
#include "includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
#include "includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
#include "includes/ParameterSet/ParameterFileSpace/parameterfileprocessingexception.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"

#include "includes/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.hpp"
#include "includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp"
#include "includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "includes/InputSet/CondensedSequenceSpace/sequencestring.hpp"

#include "includes/InputSet/CondensedSequence/graphVizDotConfig.hpp"
#include "includes/InternalPrograms/DrawGlycan/drawGlycan.hpp"

#include "includes/InternalPrograms/Sequence/sequence.hpp"

//#include "includes/InputSet/CifFileSpace/ciffileatom.hpp"
//#include "includes/InputSet/CifFileSpace/ciffile.hpp"
//#include "includes/InputSet/CifFileSpace/ciffileprocessingexception.hpp"

#include "includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbauthorsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbcaveatsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbcispeptidesection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbcispeptidecard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbcompoundsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbcompoundspecification.hpp"
#include "includes/InputSet/PdbFileSpace/pdbconnectsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbcrystallographiccard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbdatabasereferencesection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbdatabasereference.hpp"
#include "includes/InputSet/PdbFileSpace/pdbdisulfidebondsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbdisulfideresidue.hpp"
#include "includes/InputSet/PdbFileSpace/pdbdisulfideresiduebond.hpp"
#include "includes/InputSet/PdbFileSpace/pdbexperimentaldatasection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "includes/InputSet/PdbFileSpace/pdbformulasection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbformulacard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbheadercard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbhelixsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbhelixcard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbhelixresidue.hpp"
#include "includes/InputSet/PdbFileSpace/pdbheterogencard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbheterogensection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbheterogennamesection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbheterogennamecard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbheterogensynonymsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbheterogensynonymcard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbjournalsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdblinksection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbkeywordssection.hpp"
#include "includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp"
#include "includes/InputSet/PdbFileSpace/pdbmastercard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbmatrixnsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbmatrixncard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbmodelsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "includes/InputSet/PdbFileSpace/pdbmodeltypesection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbnummodelcard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbobsoletesection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbobsoletecard.hpp"
#include "includes/InputSet/PdbFileSpace/pdboriginxnsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdboriginxncard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbremarksection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbresidue.hpp"
#include "includes/InputSet/PdbFileSpace/pdbresiduemodificationsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbresiduemodificationcard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbresiduesequencesection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbresiduesequencecard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbrevisiondatasection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbrevisiondatacard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbscalensection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbscalencard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbsequenceadvancedsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbsequenceadvancedcard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbsheetsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbsheetcard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbsheetstrand.hpp"
#include "includes/InputSet/PdbFileSpace/pdbsheetstrandresidue.hpp"
#include "includes/InputSet/PdbFileSpace/pdbsitesection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbsitecard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbsiteresidue.hpp"
#include "includes/InputSet/PdbFileSpace/pdbsourcecard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbsourcesection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbsplitsection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbsupersededentriescard.hpp"
#include "includes/InputSet/PdbFileSpace/pdbsupersededentriessection.hpp"
#include "includes/InputSet/PdbFileSpace/pdbtitlesection.hpp"

#include "includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "includes/InputSet/PdbqtFileSpace/pdbqtbranchcard.hpp"
#include "includes/InputSet/PdbqtFileSpace/pdbqtcompoundcard.hpp"
#include "includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "includes/InputSet/PdbqtFileSpace/pdbqtfileprocessingexception.hpp"
#include "includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "includes/InputSet/PdbqtFileSpace/pdbqtremarkcard.hpp"
#include "includes/InputSet/PdbqtFileSpace/pdbqtrootcard.hpp"
#include "includes/InputSet/PdbqtFileSpace/pdbqttorsionaldofcard.hpp"

#include "includes/Resolver/PdbPreprocessor/pdbpreprocessor.hpp"
#include "includes/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.hpp"
#include "includes/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.hpp"
#include "includes/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.hpp"
#include "includes/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.hpp"
#include "includes/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.hpp"
#include "includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.hpp"
#include "includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.hpp"
#include "includes/Resolver/PdbPreprocessor/pdbpreprocessoralternateresidue.hpp"
#include "includes/Resolver/PdbPreprocessor/pdbpreprocessorresidueinfo.hpp"

#include "includes/GeometryTopology/rotation.hpp"
#include "includes/GeometryTopology/grid.hpp"
#include "includes/GeometryTopology/cell.hpp"
#include "includes/GeometryTopology/ResidueLinkages/residue_linkage.hpp"
#include "includes/GeometryTopology/ResidueLinkages/rotatable_dihedral.hpp"

#include "includes/MolecularMetadata/GLYCAM/amberatomtypeinfo.hpp"
#include "includes/MolecularMetadata/GLYCAM/bondlengthbytypepair.hpp"
#include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp"
#include "includes/MolecularMetadata/AMBER/amberelements.hpp"
#include "includes/MolecularMetadata/element.hpp"
#include "includes/MolecularMetadata/molecularmetadata.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06residuecodes.hpp"

#include "includes/MolecularModeling/dockingatom.hpp"
#include "includes/MolecularModeling/moleculardynamicatom.hpp"
#include "includes/MolecularModeling/quantommechanicatom.hpp"
#include "includes/MolecularModeling/atom.hpp"
#include "includes/MolecularModeling/residue.hpp"
#include "includes/MolecularModeling/atomnode.hpp"
#include "includes/MolecularModeling/assembly.hpp"
#include "includes/MolecularModeling/molecule.hpp"
#include "includes/MolecularModeling/Selections/selections.hpp"

#include "includes/InputSet/TopologyFileSpace/topologyangle.hpp"
#include "includes/InputSet/TopologyFileSpace/topologyangletype.hpp"
#include "includes/InputSet/TopologyFileSpace/topologyassembly.hpp"
#include "includes/InputSet/TopologyFileSpace/topologyatom.hpp"
#include "includes/InputSet/TopologyFileSpace/topologyatompair.hpp"
#include "includes/InputSet/TopologyFileSpace/topologybond.hpp"
#include "includes/InputSet/TopologyFileSpace/topologybondtype.hpp"
#include "includes/InputSet/TopologyFileSpace/topologydihedral.hpp"
#include "includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp"
#include "includes/InputSet/TopologyFileSpace/topologyfile.hpp"
#include "includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "includes/InputSet/TopologyFileSpace/topologyfileprocessingexception.hpp"

#include "includes/Abstract/builder.hpp"
#include "includes/InternalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"

//#include "includes/External_Libraries/json.hpp"

%}

%inline %{
std::ostream & get_cout() { return std::cout; }
%}

//%include "/usr/include/sql.h"
//%include "/usr/include/sqlext.h"

%include "includes/gmml.hpp"
%include "includes/common.hpp"
%include "includes/utils.hpp"
%include "includes/generictypedefs.hpp"
%include "includes/CodeUtils/codetests.hpp"
#include "includes/CodeUtils/logging.hpp"
%include "includes/InputSet/CoordinateFileSpace/coordinatefile.hpp"
%include "includes/InputSet/CoordinateFileSpace/coordinatefileprocessingexception.hpp"
%include "includes/GeometryTopology/coordinate.hpp"
%include "includes/GeometryTopology/plane.hpp"
%include "includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
%include "includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
%include "includes/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.hpp"
%include "includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
%include "includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
%include "includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
%include "includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
%include "includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
%include "includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
%include "includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
%include "includes/ParameterSet/ParameterFileSpace/parameterfileprocessingexception.hpp"
%include "includes/ParameterSet/PrepFileSpace/prepfile.hpp"
%include "includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
%include "includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
%include "includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"

%include "includes/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.hpp"
%include "includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
%include "includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp"
%include "includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
%include "includes/InputSet/CondensedSequenceSpace/sequencestring.hpp"

%include "includes/InputSet/CondensedSequence/graphVizDotConfig.hpp"
%include "includes/InternalPrograms/DrawGlycan/drawGlycan.hpp"

%include "includes/InternalPrograms/Sequence/sequence.hpp"

%include "includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbauthorsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbcaveatsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbcispeptidesection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbcispeptidecard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbcompoundsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbcompoundspecification.hpp"
%include "includes/InputSet/PdbFileSpace/pdbconnectsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbcrystallographiccard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbdatabasereferencesection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbdatabasereference.hpp"
%include "includes/InputSet/PdbFileSpace/pdbdisulfidebondsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbdisulfideresidue.hpp"
%include "includes/InputSet/PdbFileSpace/pdbdisulfideresiduebond.hpp"
%include "includes/InputSet/PdbFileSpace/pdbexperimentaldatasection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbfile.hpp"
%include "includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
%include "includes/InputSet/PdbFileSpace/pdbformulasection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbformulacard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbheadercard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbhelixsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbhelixcard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbhelixresidue.hpp"
%include "includes/InputSet/PdbFileSpace/pdbheterogensection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbheterogencard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbheterogennamesection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbheterogennamecard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbheterogensynonymsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbheterogensynonymcard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbjournalsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbkeywordssection.hpp"
%include "includes/InputSet/PdbFileSpace/pdblinksection.hpp"
%include "includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
%include "includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp"
%include "includes/InputSet/PdbFileSpace/pdbmastercard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbmatrixnsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbmatrixncard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbmodelsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
%include "includes/InputSet/PdbFileSpace/pdbmodeltypesection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbnummodelcard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbobsoletesection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbobsoletecard.hpp"
%include "includes/InputSet/PdbFileSpace/pdboriginxnsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdboriginxncard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbremarksection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbresiduemodificationsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbresiduemodificationcard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbresiduesequencesection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbresiduesequencecard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbrevisiondatasection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbrevisiondatacard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbscalensection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbscalencard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbsequenceadvancedsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbsequenceadvancedcard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbsheetsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbsheetcard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbsheetstrand.hpp"
%include "includes/InputSet/PdbFileSpace/pdbsheetstrandresidue.hpp"
%include "includes/InputSet/PdbFileSpace/pdbsitesection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbsitecard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbsiteresidue.hpp"
%include "includes/InputSet/PdbFileSpace/pdbsourcecard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbsourcesection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbsplitsection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbsupersededentriescard.hpp"
%include "includes/InputSet/PdbFileSpace/pdbsupersededentriessection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbtitlesection.hpp"
%include "includes/InputSet/PdbFileSpace/pdbresidue.hpp"

%include "includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
%include "includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
%include "includes/InputSet/PdbqtFileSpace/pdbqtbranchcard.hpp"
%include "includes/InputSet/PdbqtFileSpace/pdbqtcompoundcard.hpp"
%include "includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
%include "includes/InputSet/PdbqtFileSpace/pdbqtfileprocessingexception.hpp"
%include "includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
%include "includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
%include "includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
%include "includes/InputSet/PdbqtFileSpace/pdbqtremarkcard.hpp"
%include "includes/InputSet/PdbqtFileSpace/pdbqtrootcard.hpp"
%include "includes/InputSet/PdbqtFileSpace/pdbqttorsionaldofcard.hpp"

//%include "includes/InputSet/CifFileSpace/ciffileatom.hpp"
//%include "includes/InputSet/CifFileSpace/ciffile.hpp"
//%include "includes/InputSet/CifFileSpace/ciffileprocessingexception.hpp"

%include "includes/Resolver/PdbPreprocessor/pdbpreprocessor.hpp"
%include "includes/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.hpp"
%include "includes/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.hpp"
%include "includes/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.hpp"
%include "includes/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.hpp"
%include "includes/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.hpp"
%include "includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.hpp"
%include "includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.hpp"
%include "includes/Resolver/PdbPreprocessor/pdbpreprocessoralternateresidue.hpp"
%include "includes/Resolver/PdbPreprocessor/pdbpreprocessorresidueinfo.hpp"

%include "includes/MolecularMetadata/GLYCAM/amberatomtypeinfo.hpp"
%include "includes/MolecularMetadata/GLYCAM/bondlengthbytypepair.hpp"
%include "includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"
%include "includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp"
%include "includes/MolecularMetadata/AMBER/amberelements.hpp"
%include "includes/MolecularMetadata/element.hpp"
%include "includes/MolecularMetadata/molecularmetadata.hpp"
%include "includes/MolecularMetadata/GLYCAM/glycam06residuecodes.hpp"

%include "includes/MolecularModeling/dockingatom.hpp"
%include "includes/MolecularModeling/moleculardynamicatom.hpp"
%include "includes/MolecularModeling/quantommechanicatom.hpp"
%include "includes/MolecularModeling/atom.hpp"
%include "includes/MolecularModeling/residue.hpp"
%include "includes/MolecularModeling/atomnode.hpp"
%include "includes/MolecularModeling/assembly.hpp"
%include "includes/MolecularModeling/molecule.hpp"
%include "includes/MolecularModeling/Selections/selections.hpp"

%include "includes/GeometryTopology/grid.hpp"
%include "includes/GeometryTopology/cell.hpp"
%include "includes/GeometryTopology/ResidueLinkages/residue_linkage.hpp"
%include "includes/GeometryTopology/ResidueLinkages/rotatable_dihedral.hpp"

%include "includes/InputSet/TopologyFileSpace/topologyangle.hpp"
%include "includes/InputSet/TopologyFileSpace/topologyangletype.hpp"
%include "includes/InputSet/TopologyFileSpace/topologyassembly.hpp"
%include "includes/InputSet/TopologyFileSpace/topologyatom.hpp"
%include "includes/InputSet/TopologyFileSpace/topologyatompair.hpp"
%include "includes/InputSet/TopologyFileSpace/topologybond.hpp"
%include "includes/InputSet/TopologyFileSpace/topologybondtype.hpp"
%include "includes/InputSet/TopologyFileSpace/topologydihedral.hpp"
%include "includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp"
%include "includes/InputSet/TopologyFileSpace/topologyfile.hpp"
%include "includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
%include "includes/InputSet/TopologyFileSpace/topologyfileprocessingexception.hpp"

%include "includes/Abstract/builder.hpp"
%include "includes/InternalPrograms/CarbohydrateBuilder/carbohydrateBuilder.hpp"

//%include "includes/External_Libraries/json.hpp"


%template(string_vector) std::vector<std::string>;
%template(int_vector) std::vector<int>;
%template(double_vector) std::vector<double>;
%template(char_vector) std::vector<char>;
%template(bool_vector) std::vector<bool>;
//std::vector<Dihedral> improper_dihedrals_;
%template(dihedral_vector) std::vector<std::vector<std::string> >;

///GeometryTopology///
//std::vector<GeometryTopology::Coordinate*> coordinates_;
%template(coordinate_vector) std::vector<GeometryTopology::Coordinate*>;

//typedef std::vector<GeometryTopology::Cell*> CellVector;
%template(cell_vector) std::vector<GeometryTopology::Cell*>;


///Prep File///
//std::vector<PrepFileSpace::PrepFileAtom*> atoms_;
%template(prepfileatom_vector) std::vector<PrepFileSpace::PrepFileAtom*>;


///Parameter File///
%template() std::pair<std::string, ParameterFileSpace::ParameterFileAtom*>;
%template(atoms_map_parameter_file) std::map<std::string, ParameterFileSpace::ParameterFileAtom*>;

%template() std::pair<std::vector<std::string>, ParameterFileSpace::ParameterFileBond*>;
%template(bonds_map_parameter_file) std::map<std::vector<std::string>, ParameterFileSpace::ParameterFileBond*>;

%template() std::pair<std::vector<std::string>, ParameterFileSpace::ParameterFileAngle*>;
%template(angles_map_parameter_file) std::map<std::vector<std::string>, ParameterFileSpace::ParameterFileAngle*>;

//std::vector<ParameterFileSpace::ParameterFileDihedralTerm> terms_;
%template(dihedral_terms_vector_parameter_file) std::vector<ParameterFileSpace::ParameterFileDihedralTerm>;
%template() std::pair<std::vector<std::string>, ParameterFileSpace::ParameterFileDihedral*>;
%template(dihedrals_map_parameter_file) std::map<std::vector<std::string>, ParameterFileSpace::ParameterFileDihedral*>;


///Library File///
//typedef std::map<std::string, LibraryFileSpace::LibraryFileResidue*> ResidueMap;
%template() std::pair<std::string, LibraryFileSpace::LibraryFileResidue*>;
%template(residue_map_library_file) std::map<std::string, LibraryFileSpace::LibraryFileResidue*>;

//typedef std::map<int, LibraryFileSpace::LibraryFileAtom*> AtomMap;
%template() std::pair<int, LibraryFileSpace::LibraryFileAtom*>;
%template(atom_map_library_file) std::map<int, LibraryFileSpace::LibraryFileAtom*>;


///Prep File///
//typedef std::map<int, int> Loop;
%template() std::pair<int,int>;
%template(loop_map_prep_file) std::map<int, int>;

//typedef std::map< std::string, PrepFileSpace::PrepFileResidue* > ResidueMap;
%template() std::pair<std::string, PrepFileSpace::PrepFileResidue*>;
%template(residue_map_prep_file) std::map<std::string, PrepFileSpace::PrepFileResidue*>;


///PDB file///
//typedef std::map<int, PrepFileSpace::PdbAtomCard*> PdbAtomMap;
%template() std::pair<int, PdbFileSpace::PdbAtomCard*>;
%template(atom_map_pdb_file) std::map<int, PdbFileSpace::PdbAtomCard*>;

//typedef std::map<std::string, PdbFileSpace::PdbCompoundSpecification*> PdbCompoundSpecificationMap;
%template() std::pair<std::string, PdbFileSpace::PdbCompoundSpecification*>;
%template(compound_specification_map_pdb_file) std::map<std::string, PdbFileSpace::PdbCompoundSpecification*>;

//typedef std::map<int, std::vector<int> > BondedAtomsSerialNumbersMap;
%template() std::pair<int, std::vector<int> >;
%template(bonded_atoms_serial_numbers_map_pdb_file) std::map<int, std::vector<int> >;

//typedef std::map<int, PdbFileSpace::PdbDisulfideResidueBond*> DisulfideResidueBondMap;
%template() std::pair<int, PdbFileSpace::PdbDisulfideResidueBond*>;
%template(disulfide_residue_bond_map_pdb_file) std::map<int, PdbFileSpace::PdbDisulfideResidueBond*>;

//typedef std::vector<PdbFileSpace::PdbDisulfideResidue*> DisulfideResidueVector;
%template(pdbdisulfideresidue_vector) std::vector<PdbFileSpace::PdbDisulfideResidue*>;

//typedef std::map<std::string, PdbFileSpace::PdbFormulaCard*> FormulaCardMap;
%template() std::pair<std::string, PdbFileSpace::PdbFormulaCard*>;
%template(formula_map_pdb_file) std::map<std::string, PdbFileSpace::PdbFormulaCard*>;

//typedef std::vector<PdbFileSpace::PdbHelixResidue*> HelixResidueVector;
%template(pdbhelixresidue_vector) std::vector<PdbFileSpace::PdbHelixResidue*>;

//typedef std::vector<PdbFileSpace::PdbObsoleteCard*> ObsoleteCardVector;
%template(pdbobsolete_card_vector) std::vector<PdbFileSpace::PdbObsoleteCard*>;

//typedef std::map<std::string, PdbFileSpace::PdbHelixCard*> HelixCardMap;
%template() std::pair<std::string, PdbFileSpace::PdbHelixCard*>;
%template(helix_map_pdb_file) std::map<std::string, PdbFileSpace::PdbHelixCard*>;

//typedef PdbFileSpace::PdbAtomCard PdbFileSpace::PdbHeterogenAtomSection;
///???

//typedef std::map<int, PdbFileSpace::PdbHeterogenAtom*> PdbHeterogenAtomCardMap;
//%template() std::pair<int, PdbFileSpace::PdbHeterogenAtomSection::PdbHeterogenAtom*>;
//%template(heterogen_atom_map_pdb_file) std::map<int, PdbFileSpace::PdbHeterogenAtomSection::PdbHeterogenAtom*>;

//typedef std::map<std::string, PdbFileSpace::PdbHeterogenCard*> HeterogenCardMap;
%template() std::pair<std::string, PdbFileSpace::PdbHeterogenCard*>;
%template(heterogen_pdb_map_file) std::map<std::string, PdbFileSpace::PdbHeterogenCard*>;

//typedef std::map<std::string, PdbFileSpace::PdbHeterogenNameCard*> HeterogenNameCardMap;
%template() std::pair<std::string, PdbFileSpace::PdbHeterogenNameCard*>;
%template(heterogen_name_map_pdb_file) std::map<std::string, PdbFileSpace::PdbHeterogenNameCard*>;

//typedef std::map<std::string, PdbFileSpace::PdbHeterogenSynonymCard*> HeterogenSynonymCardMap;
%template() std::pair<std::string, PdbFileSpace::PdbHeterogenSynonymCard*>;
%template(heterogen_synonym_map_pdb_file) std::map<std::string, PdbFileSpace::PdbHeterogenSynonymCard*>;

//typedef std::vector< PdbLinkCardResidue* > LinkResidueVector;
%template(pdblinkcardresidue_vector) std::vector<PdbFileSpace::PdbLinkCardResidue*>;

//typedef std::vector< PdbFileSpace::PdbLinkCard* > LinkCardVector;
%template(pdblink_vector) std::vector<PdbFileSpace::PdbLinkCard*>;

//typedef std::vector<PdbFileSpace::PdbMatrixNCard*> MatrixNVector;
%template(pdbmatrixn_vector) std::vector<PdbFileSpace::PdbMatrixNCard*>;

//typedef std::vector<MatrixNVector> MatrixNVectorVector;
%template(pdbmatrixn_vector_vector) std::vector<PdbFileSpace::PdbMatrixNSection::MatrixNVector>;

//typedef std::map<int, PdbFileSpace::PdbModelCard*> PdbModelCardMap;
%template() std::pair<int, PdbFileSpace::PdbModelCard*>;
%template(model_map_pdb_file) std::map<int, PdbFileSpace::PdbModelCard*>;

//typedef std::vector<PdbFileSpace::PdbAtomSection*> AtomSectionVector;
%template(pdbatomsection_vector) std::vector<PdbFileSpace::PdbAtomSection*>;

//typedef std::vector<PdbHeterogenAtomSection*> HeterogenAtomCardVector;
%template(pdbheterogenatomcard_vector) std::vector<PdbFileSpace::PdbHeterogenAtomSection*>;

//typedef std::vector< PdbOriginXnCard* > OriginXnVector;
%template(pdboriginxn_vector) std::vector<PdbFileSpace::PdbOriginXnCard*>;

//typedef std::map<std::string, PdbResidueModificationCard*> ResidueModificationMap;
%template() std::pair<std::string, PdbFileSpace::PdbResidueModificationCard*>;
%template(residue_modification_map_pdb_file) std::map<std::string, PdbFileSpace::PdbResidueModificationCard*>;

//typedef std::map<char, PdbFileSpace::PdbResidueSequenceCard*> ResidueSequenceMap;
%template() std::pair<char, PdbFileSpace::PdbResidueSequenceCard*>;
%template(residue_sequence_map_pdb_file) std::map<char, PdbFileSpace::PdbResidueSequenceCard*>;

//typedef std::vector< PdbFileSpace::PdbScaleNCard* > ScaleNVector;
%template(pdbscalenvector_vector) std::vector<PdbFileSpace::PdbScaleNCard*>;

//typedef std::vector<PdbFileSpace::PdbSheetStrand*> SheetStrandVector;
%template(pdbsheetstrand_vector) std::vector<PdbFileSpace::PdbSheetStrand*>;

//typedef std::map<std::string, PdbFileSpace::PdbSheetCard*> SheetMap;
%template() std::pair<std::string, PdbFileSpace::PdbSheetCard*>;
%template(sheet_pdb_map_file) std::map<std::string, PdbFileSpace::PdbSheetCard*>;

//typedef std::vector<PdbFileSpace::PdbSheetStrandResidue*> SheetStrandResidueVector;
%template(pdbsheetstrandresidue_vector) std::vector<PdbFileSpace::PdbSheetStrandResidue*>;

//typedef std::vector< PdbFileSpace::PdbSiteResidue* > SiteResidueVector;
%template(pdbsiteresidue_vector) std::vector<PdbFileSpace::PdbSiteResidue*>;

//typedef std::map<std::string, PdbFileSpace::PdbSiteCard*> PdbSiteCardMap;
%template() std::pair<std::string, PdbFileSpace::PdbSiteCard*>;
%template(site_map_pdb_file) std::map<std::string, PdbFileSpace::PdbSiteCard*>;

//typedef std::vector<PdbSourceCard*> SourceCardVector;
%template(pdb_source_vector) std::vector<PdbFileSpace::PdbSourceCard*>;

//typedef std::vector<PdbFileSpace::PdbRevisionDataCard*> RevisionDataCardVector;
%template(pdb_revision_vector) std::vector<PdbFileSpace::PdbRevisionDataCard*>;

//typedef std::vector<PdbFileSpace::PdbSupersededEntriesCard*> SupersededEntriesCardVector;
%template(pdb_superseded_vector) std::vector<PdbFileSpace::PdbSupersededEntriesCard*>;

//typedef std::vector<PdbFileSpace::PdbAtomCard*> PdbAtomCardVector;
%template(pdb_atom_vector) std::vector<PdbFileSpace::PdbAtomCard*>;

//typedef std::vector<PdbFileSpace::PdbSequenceAdvancedCard*> SequenceAdvancedCardVector;
%template(pdb_seq_adv_vector) std::vector<PdbFileSpace::PdbSequenceAdvancedCard*>;

//typedef std::vector<PdbFileSpace::PdbCISPeptideCard*> CISPeptideCardVector;
%template(pdb_cis_peptide_vector) std::vector<PdbFileSpace::PdbCISPeptideCard*>;


///PDB Preprocessor///
//typedef std::vector<PdbPreprocessorSpace::PdbPreprocessorDisulfideBond*> PdbPreprocessorDisulfideBondVector;
%template(pdbpreprocessordisulfidebond_vector) std::vector<PdbPreprocessorSpace::PdbPreprocessorDisulfideBond*>;

//typedef std::vector<PdbPreprocessorSpace::PdbPreprocessorChainTermination*> PdbPreprocessorChainTerminationVector;
%template(pdbpreprocessorchaintermination_vector) std::vector<PdbPreprocessorSpace::PdbPreprocessorChainTermination*>;

//typedef std::vector<PdbPreprocessorSpace::PdbPreprocessorHistidineMapping*> PdbPreprocessorHistidineMappingVector;
%template(pdbpreprocessorhistidinemapping_vector) std::vector<PdbPreprocessorSpace::PdbPreprocessorHistidineMapping*>;

//typedef std::vector<PdbPreprocessorSpace::PdbPreprocessorMissingResidue*> PdbPreprocessorMissingResidueVector;
%template(pdbpreprocessormissingresidue_vector) std::vector<PdbPreprocessorSpace::PdbPreprocessorMissingResidue*>;

//typedef std::vector<PdbPreprocessorSpace::PdbPreprocessorUnrecognizedResidue*> PdbPreprocessorUnrecognizedResidueVector;
%template(pdbpreprocessorunrecognizedresidue_vector) std::vector<PdbPreprocessorSpace::PdbPreprocessorUnrecognizedResidue*>;

//typedef std::vector<PdbPreprocessorSpace::PdbPreprocessorUnrecognizedHeavyAtom*> PdbPreprocessorUnrecognizedHeavyAtomVector;
%template(pdbpreprocessorunrecognizedheavyatom_vector) std::vector<PdbPreprocessorSpace::PdbPreprocessorUnrecognizedHeavyAtom*>;

//typedef std::vector<PdbPreprocessorSpace::PdbPreprocessorReplacedHydrogen*> PdbPreprocessorReplacedHydrogenVector;
%template(pdbpreprocessorreplacedhydrogen_vector) std::vector<PdbPreprocessorSpace::PdbPreprocessorReplacedHydrogen*>;

//typedef std::map<std::string, PdbPreprocessorSpace::PdbPreprocessorAlternateResidue*> PdbPreprocessorAlternateResidueMap;
%template() std::pair<std::string, PdbPreprocessorSpace::PdbPreprocessorAlternateResidue*>;
%template(alternate_residue_map_pdbpreprocessor_file) std::map<std::string, PdbPreprocessorSpace::PdbPreprocessorAlternateResidue*>;

//typedef std::map<std::string, PdbPreprocessorSpace::PdbPreprocessorResidueInfo*> PdbPreprocessorResidueInfoMap;
%template() std::pair<std::string, PdbPreprocessorSpace::PdbPreprocessorResidueInfo*>;
%template(residue_info_map_pdbpreprocessor_file) std::map<std::string, PdbPreprocessorSpace::PdbPreprocessorResidueInfo*>;


///MolecularModeling///
//typedef std::vector<MolecularModeling::Assembly*> AssemblyVector;
%template(assembly_vector) std::vector<MolecularModeling::Assembly* >;

//typedef std::vector<MolecularModeling::Residue*> ResidueVector;
%template(residue_vector) std::vector<MolecularModeling::Residue* >;

//typedef std::vector<MolecularModeling::Atom*> AtomVector;
%template(atom_vector) std::vector<MolecularModeling::Atom* >;

//typedef std::vector<AtomVector > AtomVectorVector;
//%template(atom_vector_vector) std::vector<std::vector<MolecularModeling::Atom* > >;

//typedef std::map<std::string, AtomVector> CycleMap;
%template() std::pair<std::string, std::vector<MolecularModeling::Atom* > >;
%template(cycle_map_assembly_file) std::map<std::string, std::vector<MolecularModeling::Atom* > >;

//typedef std::vector<Glycan::Oligosaccharide*> OligosaccharideVector;
%template(oligosaccharide_vector) std::vector<Glycan::Oligosaccharide* >;

//typedef std::vector<MolecularModeling::ResidueNode*> ResidueNodeVector;
%template(residue_node_vector) std::vector<MolecularModeling::ResidueNode* >;

///Topology File///
//typedef std::map<std::string, TopologyFileSpace::TopologyResidue*> TopologyResidueMap;
%template() std::pair<std::string, TopologyFileSpace::TopologyResidue*>;
%template(residue_map_Topology_file) std::map<std::string, TopologyFileSpace::TopologyResidue*>;

//typedef std::map<std::vector<std::string>, TopologyFileSpace::TopologyBond*> TopologyBondMap;
%template() std::pair<std::string, TopologyFileSpace::TopologyBond*>;
%template(bond_map_Topology_file) std::map<std::string, TopologyFileSpace::TopologyBond*>;

//typedef std::map<std::vector<std::string>, TopologyFileSpace::TopologyAngle*> TopologyAngleMap;
%template() std::pair<std::string, TopologyFileSpace::TopologyAngle*>;
%template(angle_map_Topology_file) std::map<std::string, TopologyFileSpace::TopologyAngle*>;

//typedef std::map<std::vector<std::string>, TopologyFileSpace::TopologyDihedral*> TopologyDihedralMap;
%template() std::pair<std::string, TopologyFileSpace::TopologyDihedral*>;
%template(dihedral_map_Topology_file) std::map<std::string, TopologyFileSpace::TopologyDihedral*>;

//typedef std::map<int, double> TopologyCoefficientMap;
%template() std::pair<int, double>;
%template(coefficient_map_Topology_file) std::map<int, double>;

//typedef std::map<int, TopologyFileSpace::TopologyAtomType*> TopologyAtomPairMap;
%template() std::pair<std::string, TopologyFileSpace::TopologyAtomPair*>;
%template(atom_pair_map_Topology_file) std::map<std::string, TopologyFileSpace::TopologyAtomPair*>;

//typedef std::map<int, TopologyFileSpace::TopologyBondType*> TopologyBondTypeMap;
%template() std::pair<int, TopologyFileSpace::TopologyBondType*>;
%template(bond_type_map_Topology_file) std::map<int, TopologyFileSpace::TopologyBondType*>;

//typedef std::map<int, TopologyFileSpace::TopologyAngleType*> TopologyAngleTypeMap;
%template() std::pair<int, TopologyFileSpace::TopologyAngleType*>;
%template(angle_type_map_Topology_file) std::map<int, TopologyFileSpace::TopologyAngleType*>;

//typedef std::map<int, TopologyFileSpace::TopologyDihedralType*> TopologyDihedralTypeMap;
%template() std::pair<int, TopologyFileSpace::TopologyDihedralType*>;
%template(dihedral_type_map_Topology_file) std::map<int, TopologyFileSpace::TopologyDihedralType*>;


///PDBQT file///
//typedef std::map<std::string, PdbFileSpace::PdbFile::PdbAtomCardVector* > PdbResidueAtomsMap;
%template() std::pair<std::string, PdbFileSpace::PdbFile::PdbAtomCardVector* >;
%template(pdb_residue_atom_map) std::map<std::string, PdbFileSpace::PdbFile::PdbAtomCardVector* >;

//typedef std::vector<PdbqtFileSpace::PdbqtAtom*> PdbqtAtomVector;
%template(pdbqt_atom_vector) std::vector<PdbqtFileSpace::PdbqtAtom*>;

//typedef std::map<std::string, PdbqtFileSpace::PdbqtFile::PdbqtAtomVector* > PdbqtResidueAtomsMap;
%template() std::pair<std::string, PdbqtFileSpace::PdbqtFile::PdbqtAtomVector*>;
%template(pdbqt_residue_atom_map) std::map<std::string, PdbqtFileSpace::PdbqtFile::PdbqtAtomVector*>;

//typedef std::map<int, PdbqtFileSpace::PdbqtModel*> PdbqtModelMap;
%template() std::pair<int, PdbqtFileSpace::PdbqtModel*>;
%template(pdbqt_model_map) std::map<int, PdbqtFileSpace::PdbqtModel*>;

//typedef std::vector<PdbqtFileSpace::PdbqtRemarkCard*> RemarkCardVector;
%template(pdbqt_remark_card_vector) std::vector<PdbqtFileSpace::PdbqtRemarkCard*>;

//typedef std::vector<PdbqtFileSpace::PdbqtTorsionalDoFCard*> TorsionalDoFCardVector;
%template(torsional_dof_card_vector) std::vector<PdbqtFileSpace::PdbqtTorsionalDoFCard*>;

//typedef std::map<int, PdbqtFileSpace::PdbqtAtom*> PdbqtAtomMap;
%template() std::pair<int, PdbqtFileSpace::PdbqtAtom* >;
%template(pdbqt_atom_map) std::map<int, PdbqtFileSpace::PdbqtAtom* >;

//typedef std::vector<PdbqtFileSpace::PdbqtBranchCard*> BranchCardVector;
%template(pdbqt_branch_card_vector) std::vector<PdbqtFileSpace::PdbqtBranchCard*>;

///Cif File///
//typedef std::vector<CifFileSpace::CifFileAtom*> CifFileAtomVector;
//%template(cif_atom_vector) std::vector<CifFileSpace::CifFileAtom*>;

///Condensed Sequence///
//typedef std::vector<CondensedSequenceSpace::CondensedSequenceResidue*> CondensedSequenceResidueVector;
%template(condensedsequence_residue_vector) std::vector<CondensedSequenceSpace::CondensedSequenceResidue*>;

//typedef std::vector<gmml::CondensedSequenceTokenType> CondensedSequenceTokenTypeVector;
%template(condensedsequence_token_type_vector) std::vector<gmml::CondensedSequenceTokenType>;

//typedef std::vector<CondensedSequenceSpace::CondensedSequenceResidue*> CondensedSequenceResidueTree;
//%template(condensedsequence_residue_tree) std::vector<CondensedSequenceSpace::CondensedSequenceResidue*>;

//typedef std::vector<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*> CondensedSequenceGlycam06ResidueTree;
%template(condensedsequence_glycam06_residue_tree) std::vector<CondensedSequenceSpace::CondensedSequenceGlycam06Residue*>;

//typedef std::pair<std::string, CondensedSequenceSpace::RotamersAndGlycosidicAnglesInfo*> RotamerNameInfoPair;
%template(rotamer_name_info_pair) std::pair<std::string, CondensedSequenceSpace::RotamersAndGlycosidicAnglesInfo*>;

//typedef std::vector<RotamerNameInfoPair> CondensedSequenceRatomersAndGlycosidicAnglesInfo;
%template(rotamer_angle_info_vector) std::vector<std::pair<std::string, CondensedSequenceSpace::RotamersAndGlycosidicAnglesInfo*> >;

//std::pair<std::string, double>
%template(string_double_pair) std::pair<std::string, double>;

//std::vector<std::pair<std::string, double> >
%template(glycosidic_angle_name_value_pair_vector) std::vector<std::pair<std::string, double> >;

//std::pair<std::string, std::vector<std::string> >
%template(string_vector_string_pair) std::pair<std::string, std::vector<std::string> >;

//std::vector<std::pair<std::string, std::vector<std::string> > >
%template(string_vector_string_pair_vector) std::vector<std::pair<std::string, std::vector<std::string> > >;

%template(vector_vector_int) std::vector<std::vector<int> >;

//std::vector<std::vector<double> >
%template(vector_vector_double) std::vector<std::vector<double> >;

//typedef std::map<int, std::vector<std::vector<double> > > IndexLinkageConfigurationMap;
%template() std::pair<int, std::vector<std::vector<double> > >;
%template(int_vector_vector_double_map) std::map<int, std::vector<std::vector<double> > >;

//typedef std::map<int, std::string> IndexNameMap;
%template() std::pair<int, std::string>;
%template(int_string_map) std::map<int, std::string>;

//typedef std::map<int, std::string> DerivativeMap;
//%template() std::pair<int, std::string >;
//%template(condensedsequence_derivative_map) std::map<int, std::string >;

///Common///
//typedef std::map<std::string, std::string> ResidueNameMap;
%template() std::pair<std::string, std::string>;
%template(residue_name_map) std::map<std::string, std::string>;

///Utils///
//typedef std::map<int, std::vector<Glycan::SugarName> > SugarNameClosestMatchMap;
//%template() std::pair<int, std::vector<Glycan::SugarName> >;
//%template(sugar_name_closest_match_map) std::map<int, std::vector<Glycan::SugarName> >;

///Carbohydrate Builder///
//typedef std::vector<DihedralOptions> DihedralOptionsVector;
%template(dihedral_options_vector) std::vector<CondensedSequence::DihedralOptions>;

//typedef std::vector<LinkageOptions> LinkageOptionsVector;
%template(linkage_options_vector) std::vector<CondensedSequence::LinkageOptions>;

//typedef std::vector<SingleRotamerInfo> SingleRotamerInfoVector;
%template(single_rotamer_info_vector) std::vector<CondensedSequence::SingleRotamerInfo>;

//typedef std::vector<MolecularModeling::Residue*> ResidueVector;
%template(residue_vector) std::vector<MolecularModeling::Residue* >;

//constexpr operator size_t() { return 0; }
