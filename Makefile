#############################################################################
# Makefile for building: libgmml.so.1.0.0
# Project:  gmml
# Template: lib
#############################################################################

####### Compiler, tools and options

CC            = gcc
CXX           = g++
CFLAGS        = -m64 -pipe -O2 -Wall -W -D_REENTRANT -fPIC
CXXFLAGS      = -m64 -pipe -O2 -Wall -W -D_REENTRANT -fPIC
INCPATH       = -I. -Iincludes/ParameterSet/ParameterFileSpace -Iincludes -Iincludes/Geometry -Iincludes/ParameterSet/PrepFileSpace -Iincludes/ParameterSet/LibraryFileSpace -Iincludes/FileSet/CoordinateFileSpace -Iincludes/FileSet/PdbFileSpace -Iincludes/Resolver/PdbPreprocessor -Iincludes/FileSet/TopologyFileSpace -Iincludes/MolecularModeling -Iincludes/Geometry/InternalCoordinate -I.
LINK          = g++
LFLAGS        = -m64 -Wl,-O1 -shared -Wl,-soname,libgmml.so.1
LIBS          = $(SUBLIBS)  -L/usr/lib/x86_64-linux-gnu -lpthread 
AR            = ar cqs
RANLIB        = 
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
SED           = sed
COPY_FILE     = $(COPY)
COPY_DIR      = $(COPY) -r
STRIP         = strip
INSTALL_FILE  = install -m 644 -p
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = install -m 755 -p
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = src/Geometry/coordinate.cc \
		src/Geometry/plane.cc \
		src/MolecularModeling/assembly.cc \
		src/MolecularModeling/atom.cc \
		src/MolecularModeling/atomnode.cc \
		src/MolecularModeling/dockingatom.cc \
		src/MolecularModeling/element.cc \
		src/MolecularModeling/moleculardynamicatom.cc \
		src/MolecularModeling/quantommechanicatom.cc \
		src/MolecularModeling/residue.cc \
		src/FileSet/CoordinateFileSpace/coordinatefile.cc \
		src/FileSet/CoordinateFileSpace/coordinatefileprocessingexception.cc \
		src/FileSet/PdbFileSpace/pdbatom.cc \
		src/FileSet/PdbFileSpace/pdbatomcard.cc \
		src/FileSet/PdbFileSpace/pdbcompoundcard.cc \
		src/FileSet/PdbFileSpace/pdbcompoundspecification.cc \
		src/FileSet/PdbFileSpace/pdbconnectcard.cc \
		src/FileSet/PdbFileSpace/pdbcrystallographiccard.cc \
		src/FileSet/PdbFileSpace/pdbdisulfidebondcard.cc \
		src/FileSet/PdbFileSpace/pdbdisulfideresidue.cc \
		src/FileSet/PdbFileSpace/pdbdisulfideresiduebond.cc \
		src/FileSet/PdbFileSpace/pdbfile.cc \
		src/FileSet/PdbFileSpace/pdbfileprocessingexception.cc \
		src/FileSet/PdbFileSpace/pdbformula.cc \
		src/FileSet/PdbFileSpace/pdbformulacard.cc \
		src/FileSet/PdbFileSpace/pdbheadercard.cc \
		src/FileSet/PdbFileSpace/pdbhelix.cc \
		src/FileSet/PdbFileSpace/pdbhelixcard.cc \
		src/FileSet/PdbFileSpace/pdbhelixresidue.cc \
		src/FileSet/PdbFileSpace/pdbheterogen.cc \
		src/FileSet/PdbFileSpace/pdbheterogenatomcard.cc \
		src/FileSet/PdbFileSpace/pdbheterogencard.cc \
		src/FileSet/PdbFileSpace/pdbheterogenname.cc \
		src/FileSet/PdbFileSpace/pdbheterogennamecard.cc \
		src/FileSet/PdbFileSpace/pdbheterogensynonym.cc \
		src/FileSet/PdbFileSpace/pdbheterogensynonymcard.cc \
		src/FileSet/PdbFileSpace/pdblink.cc \
		src/FileSet/PdbFileSpace/pdblinkcard.cc \
		src/FileSet/PdbFileSpace/pdblinkresidue.cc \
		src/FileSet/PdbFileSpace/pdbmatrixn.cc \
		src/FileSet/PdbFileSpace/pdbmatrixncard.cc \
		src/FileSet/PdbFileSpace/pdbmodel.cc \
		src/FileSet/PdbFileSpace/pdbmodelcard.cc \
		src/FileSet/PdbFileSpace/pdbmodelresidueset.cc \
		src/FileSet/PdbFileSpace/pdbmodeltypecard.cc \
		src/FileSet/PdbFileSpace/pdbnummodelcard.cc \
		src/FileSet/PdbFileSpace/pdboriginxn.cc \
		src/FileSet/PdbFileSpace/pdboriginxncard.cc \
		src/FileSet/PdbFileSpace/pdbresidue.cc \
		src/FileSet/PdbFileSpace/pdbresiduemodification.cc \
		src/FileSet/PdbFileSpace/pdbresiduemodificationcard.cc \
		src/FileSet/PdbFileSpace/pdbresiduesequence.cc \
		src/FileSet/PdbFileSpace/pdbresiduesequencecard.cc \
		src/FileSet/PdbFileSpace/pdbscalen.cc \
		src/FileSet/PdbFileSpace/pdbscalencard.cc \
		src/FileSet/PdbFileSpace/pdbsheet.cc \
		src/FileSet/PdbFileSpace/pdbsheetcard.cc \
		src/FileSet/PdbFileSpace/pdbsheetstrand.cc \
		src/FileSet/PdbFileSpace/pdbsheetstrandresidue.cc \
		src/FileSet/PdbFileSpace/pdbsite.cc \
		src/FileSet/PdbFileSpace/pdbsitecard.cc \
		src/FileSet/PdbFileSpace/pdbsiteresidue.cc \
		src/FileSet/PdbFileSpace/pdbtitlecard.cc \
		src/FileSet/TopologyFileSpace/topologyangle.cc \
		src/FileSet/TopologyFileSpace/topologyangletype.cc \
		src/FileSet/TopologyFileSpace/topologyassembly.cc \
		src/FileSet/TopologyFileSpace/topologyatom.cc \
		src/FileSet/TopologyFileSpace/topologyatompair.cc \
		src/FileSet/TopologyFileSpace/topologybond.cc \
		src/FileSet/TopologyFileSpace/topologybondtype.cc \
		src/FileSet/TopologyFileSpace/topologydihedral.cc \
		src/FileSet/TopologyFileSpace/topologydihedraltype.cc \
		src/FileSet/TopologyFileSpace/topologyfile.cc \
		src/FileSet/TopologyFileSpace/topologyfileprocessingexception.cc \
		src/FileSet/TopologyFileSpace/topologyresidue.cc \
		src/Geometry/InternalCoordinate/angle.cc \
		src/Geometry/InternalCoordinate/dihedral.cc \
		src/Geometry/InternalCoordinate/distance.cc \
		src/ParameterSet/LibraryFileSpace/libraryfile.cc \
		src/ParameterSet/LibraryFileSpace/libraryfileatom.cc \
		src/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.cc \
		src/ParameterSet/LibraryFileSpace/libraryfileresidue.cc \
		src/ParameterSet/ParameterFileSpace/parameterfile.cc \
		src/ParameterSet/ParameterFileSpace/parameterfileangle.cc \
		src/ParameterSet/ParameterFileSpace/parameterfileatom.cc \
		src/ParameterSet/ParameterFileSpace/parameterfilebond.cc \
		src/ParameterSet/ParameterFileSpace/parameterfiledihedral.cc \
		src/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.cc \
		src/ParameterSet/ParameterFileSpace/parameterfileprocessingexception.cc \
		src/ParameterSet/PrepFileSpace/prepfile.cc \
		src/ParameterSet/PrepFileSpace/prepfileatom.cc \
		src/ParameterSet/PrepFileSpace/prepfileprocessingexception.cc \
		src/ParameterSet/PrepFileSpace/prepfileresidue.cc \
		src/Resolver/PdbPreprocessor/pdbpreprocessor.cc \
		src/Resolver/PdbPreprocessor/pdbpreprocessoralternateresidue.cc \
		src/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.cc \
		src/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.cc \
		src/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.cc \
		src/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.cc \
		src/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.cc \
		src/Resolver/PdbPreprocessor/pdbpreprocessorresidueinfo.cc \
		src/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.cc \
		src/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.cc 
OBJECTS       = coordinate.o \
		plane.o \
		assembly.o \
		atom.o \
		atomnode.o \
		dockingatom.o \
		element.o \
		moleculardynamicatom.o \
		quantommechanicatom.o \
		residue.o \
		coordinatefile.o \
		coordinatefileprocessingexception.o \
		pdbatom.o \
		pdbatomcard.o \
		pdbcompoundcard.o \
		pdbcompoundspecification.o \
		pdbconnectcard.o \
		pdbcrystallographiccard.o \
		pdbdisulfidebondcard.o \
		pdbdisulfideresidue.o \
		pdbdisulfideresiduebond.o \
		pdbfile.o \
		pdbfileprocessingexception.o \
		pdbformula.o \
		pdbformulacard.o \
		pdbheadercard.o \
		pdbhelix.o \
		pdbhelixcard.o \
		pdbhelixresidue.o \
		pdbheterogen.o \
		pdbheterogenatomcard.o \
		pdbheterogencard.o \
		pdbheterogenname.o \
		pdbheterogennamecard.o \
		pdbheterogensynonym.o \
		pdbheterogensynonymcard.o \
		pdblink.o \
		pdblinkcard.o \
		pdblinkresidue.o \
		pdbmatrixn.o \
		pdbmatrixncard.o \
		pdbmodel.o \
		pdbmodelcard.o \
		pdbmodelresidueset.o \
		pdbmodeltypecard.o \
		pdbnummodelcard.o \
		pdboriginxn.o \
		pdboriginxncard.o \
		pdbresidue.o \
		pdbresiduemodification.o \
		pdbresiduemodificationcard.o \
		pdbresiduesequence.o \
		pdbresiduesequencecard.o \
		pdbscalen.o \
		pdbscalencard.o \
		pdbsheet.o \
		pdbsheetcard.o \
		pdbsheetstrand.o \
		pdbsheetstrandresidue.o \
		pdbsite.o \
		pdbsitecard.o \
		pdbsiteresidue.o \
		pdbtitlecard.o \
		topologyangle.o \
		topologyangletype.o \
		topologyassembly.o \
		topologyatom.o \
		topologyatompair.o \
		topologybond.o \
		topologybondtype.o \
		topologydihedral.o \
		topologydihedraltype.o \
		topologyfile.o \
		topologyfileprocessingexception.o \
		topologyresidue.o \
		angle.o \
		dihedral.o \
		distance.o \
		libraryfile.o \
		libraryfileatom.o \
		libraryfileprocessingexception.o \
		libraryfileresidue.o \
		parameterfile.o \
		parameterfileangle.o \
		parameterfileatom.o \
		parameterfilebond.o \
		parameterfiledihedral.o \
		parameterfiledihedralterm.o \
		parameterfileprocessingexception.o \
		prepfile.o \
		prepfileatom.o \
		prepfileprocessingexception.o \
		prepfileresidue.o \
		pdbpreprocessor.o \
		pdbpreprocessoralternateresidue.o \
		pdbpreprocessorchaintermination.o \
		pdbpreprocessordisulfidebond.o \
		pdbpreprocessorhistidinemapping.o \
		pdbpreprocessormissingresidue.o \
		pdbpreprocessorreplacedhydrogen.o \
		pdbpreprocessorresidueinfo.o \
		pdbpreprocessorunrecognizedheavyatom.o \
		pdbpreprocessorunrecognizedresidue.o
DESTDIR       = 
TARGET        = libgmml.so.1.0.0
TARGETA       = libgmml.a
TARGETD       = libgmml.so.1.0.0
TARGET0       = libgmml.so
TARGET1       = libgmml.so.1
TARGET2       = libgmml.so.1.0

first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

all: Makefile  $(TARGET)

$(TARGET):  $(OBJECTS) $(SUBLIBS) $(OBJCOMP)  
	-$(DEL_FILE) $(TARGET) $(TARGET0) $(TARGET1) $(TARGET2)
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS) $(OBJCOMP)
	-ln -s $(TARGET) $(TARGET0)
	-ln -s $(TARGET) $(TARGET1)
	-ln -s $(TARGET) $(TARGET2)



staticlib: $(TARGETA)

$(TARGETA):  $(OBJECTS) $(OBJCOMP) 
	-$(DEL_FILE) $(TARGETA) 
	$(AR) $(TARGETA) $(OBJECTS)

clean:compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) $(TARGET0) $(TARGET1) $(TARGET2) $(TARGETA)

check: first

mocclean: compiler_moc_header_clean compiler_moc_source_clean

mocables: compiler_moc_header_make_all compiler_moc_source_make_all

compiler_moc_header_make_all:
compiler_moc_header_clean:
compiler_rcc_make_all:
compiler_rcc_clean:
compiler_image_collection_make_all:
compiler_image_collection_clean:
compiler_moc_source_make_all:
compiler_moc_source_clean:
compiler_uic_make_all:
compiler_uic_clean:
compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: 

####### Compile

coordinate.o: src/Geometry/coordinate.cc includes/Geometry/coordinate.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o coordinate.o src/Geometry/coordinate.cc

plane.o: src/Geometry/plane.cc includes/Geometry/plane.hpp \
		includes/Geometry/coordinate.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o plane.o src/Geometry/plane.cc

assembly.o: src/MolecularModeling/assembly.cc includes/MolecularModeling/assembly.hpp \
		includes/Geometry/coordinate.hpp \
		includes/Geometry/plane.hpp \
		includes/common.hpp \
		includes/FileSet/PdbFileSpace/pdbfile.hpp \
		includes/FileSet/TopologyFileSpace/topologyfile.hpp \
		includes/FileSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/FileSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/FileSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/FileSet/PdbFileSpace/pdbmodel.hpp \
		includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/FileSet/TopologyFileSpace/topologyassembly.hpp \
		includes/FileSet/TopologyFileSpace/topologyresidue.hpp \
		includes/FileSet/TopologyFileSpace/topologyatom.hpp \
		includes/FileSet/TopologyFileSpace/topologybond.hpp \
		includes/FileSet/TopologyFileSpace/topologybondtype.hpp \
		includes/FileSet/TopologyFileSpace/topologyangle.hpp \
		includes/FileSet/TopologyFileSpace/topologyangletype.hpp \
		includes/FileSet/TopologyFileSpace/topologydihedral.hpp \
		includes/FileSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/FileSet/TopologyFileSpace/topologyatompair.hpp \
		includes/FileSet/PdbFileSpace/pdbtitlecard.hpp \
		includes/FileSet/PdbFileSpace/pdbatomcard.hpp \
		includes/FileSet/PdbFileSpace/pdbheterogenatomcard.hpp \
		includes/FileSet/PdbFileSpace/pdbatom.hpp \
		includes/FileSet/PdbFileSpace/pdbconnectcard.hpp \
		includes/FileSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/utils.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o assembly.o src/MolecularModeling/assembly.cc

atom.o: src/MolecularModeling/atom.cc includes/MolecularModeling/atom.hpp \
		includes/Geometry/coordinate.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/atomnode.hpp \
		includes/MolecularModeling/residue.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o atom.o src/MolecularModeling/atom.cc

atomnode.o: src/MolecularModeling/atomnode.cc includes/MolecularModeling/atomnode.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/Geometry/coordinate.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp \
		includes/MolecularModeling/residue.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o atomnode.o src/MolecularModeling/atomnode.cc

dockingatom.o: src/MolecularModeling/dockingatom.cc includes/MolecularModeling/dockingatom.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o dockingatom.o src/MolecularModeling/dockingatom.cc

element.o: src/MolecularModeling/element.cc includes/MolecularModeling/element.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o element.o src/MolecularModeling/element.cc

moleculardynamicatom.o: src/MolecularModeling/moleculardynamicatom.cc includes/MolecularModeling/moleculardynamicatom.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moleculardynamicatom.o src/MolecularModeling/moleculardynamicatom.cc

quantommechanicatom.o: src/MolecularModeling/quantommechanicatom.cc includes/MolecularModeling/quantommechanicatom.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o quantommechanicatom.o src/MolecularModeling/quantommechanicatom.cc

residue.o: src/MolecularModeling/residue.cc includes/MolecularModeling/residue.hpp \
		includes/MolecularModeling/assembly.hpp \
		includes/Geometry/coordinate.hpp \
		includes/common.hpp \
		includes/FileSet/PdbFileSpace/pdbfile.hpp \
		includes/FileSet/TopologyFileSpace/topologyfile.hpp \
		includes/FileSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/FileSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/FileSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/FileSet/PdbFileSpace/pdbmodel.hpp \
		includes/MolecularModeling/atom.hpp \
		includes/MolecularModeling/moleculardynamicatom.hpp \
		includes/MolecularModeling/quantommechanicatom.hpp \
		includes/MolecularModeling/dockingatom.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o residue.o src/MolecularModeling/residue.cc

coordinatefile.o: src/FileSet/CoordinateFileSpace/coordinatefile.cc includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp \
		includes/FileSet/CoordinateFileSpace/coordinatefile.hpp \
		includes/FileSet/CoordinateFileSpace/coordinatefileprocessingexception.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o coordinatefile.o src/FileSet/CoordinateFileSpace/coordinatefile.cc

coordinatefileprocessingexception.o: src/FileSet/CoordinateFileSpace/coordinatefileprocessingexception.cc includes/common.hpp \
		includes/FileSet/CoordinateFileSpace/coordinatefileprocessingexception.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o coordinatefileprocessingexception.o src/FileSet/CoordinateFileSpace/coordinatefileprocessingexception.cc

pdbatom.o: src/FileSet/PdbFileSpace/pdbatom.cc includes/FileSet/PdbFileSpace/pdbatom.hpp \
		includes/Geometry/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbatom.o src/FileSet/PdbFileSpace/pdbatom.cc

pdbatomcard.o: src/FileSet/PdbFileSpace/pdbatomcard.cc includes/FileSet/PdbFileSpace/pdbatomcard.hpp \
		includes/FileSet/PdbFileSpace/pdbatom.hpp \
		includes/Geometry/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbatomcard.o src/FileSet/PdbFileSpace/pdbatomcard.cc

pdbcompoundcard.o: src/FileSet/PdbFileSpace/pdbcompoundcard.cc includes/FileSet/PdbFileSpace/pdbcompoundcard.hpp \
		includes/FileSet/PdbFileSpace/pdbcompoundspecification.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbcompoundcard.o src/FileSet/PdbFileSpace/pdbcompoundcard.cc

pdbcompoundspecification.o: src/FileSet/PdbFileSpace/pdbcompoundspecification.cc includes/FileSet/PdbFileSpace/pdbcompoundspecification.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbcompoundspecification.o src/FileSet/PdbFileSpace/pdbcompoundspecification.cc

pdbconnectcard.o: src/FileSet/PdbFileSpace/pdbconnectcard.cc includes/FileSet/PdbFileSpace/pdbconnectcard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbconnectcard.o src/FileSet/PdbFileSpace/pdbconnectcard.cc

pdbcrystallographiccard.o: src/FileSet/PdbFileSpace/pdbcrystallographiccard.cc includes/FileSet/PdbFileSpace/pdbcrystallographiccard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbcrystallographiccard.o src/FileSet/PdbFileSpace/pdbcrystallographiccard.cc

pdbdisulfidebondcard.o: src/FileSet/PdbFileSpace/pdbdisulfidebondcard.cc includes/FileSet/PdbFileSpace/pdbdisulfidebondcard.hpp \
		includes/FileSet/PdbFileSpace/pdbdisulfideresiduebond.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbdisulfidebondcard.o src/FileSet/PdbFileSpace/pdbdisulfidebondcard.cc

pdbdisulfideresidue.o: src/FileSet/PdbFileSpace/pdbdisulfideresidue.cc includes/FileSet/PdbFileSpace/pdbdisulfideresidue.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbdisulfideresidue.o src/FileSet/PdbFileSpace/pdbdisulfideresidue.cc

pdbdisulfideresiduebond.o: src/FileSet/PdbFileSpace/pdbdisulfideresiduebond.cc includes/FileSet/PdbFileSpace/pdbdisulfideresiduebond.hpp \
		includes/FileSet/PdbFileSpace/pdbdisulfideresidue.hpp \
		includes/common.hpp \
		includes/utils.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbdisulfideresiduebond.o src/FileSet/PdbFileSpace/pdbdisulfideresiduebond.cc

pdbfile.o: src/FileSet/PdbFileSpace/pdbfile.cc includes/FileSet/PdbFileSpace/pdbfile.hpp \
		includes/FileSet/PdbFileSpace/pdbheadercard.hpp \
		includes/FileSet/PdbFileSpace/pdbtitlecard.hpp \
		includes/FileSet/PdbFileSpace/pdbcompoundcard.hpp \
		includes/FileSet/PdbFileSpace/pdbcompoundspecification.hpp \
		includes/FileSet/PdbFileSpace/pdbnummodelcard.hpp \
		includes/FileSet/PdbFileSpace/pdbmodeltypecard.hpp \
		includes/FileSet/PdbFileSpace/pdbresiduesequencecard.hpp \
		includes/FileSet/PdbFileSpace/pdbresiduesequence.hpp \
		includes/FileSet/PdbFileSpace/pdbresiduemodificationcard.hpp \
		includes/FileSet/PdbFileSpace/pdbresiduemodification.hpp \
		includes/FileSet/PdbFileSpace/pdbheterogencard.hpp \
		includes/FileSet/PdbFileSpace/pdbheterogen.hpp \
		includes/FileSet/PdbFileSpace/pdbheterogennamecard.hpp \
		includes/FileSet/PdbFileSpace/pdbheterogenname.hpp \
		includes/FileSet/PdbFileSpace/pdbheterogensynonymcard.hpp \
		includes/FileSet/PdbFileSpace/pdbheterogensynonym.hpp \
		includes/FileSet/PdbFileSpace/pdbformulacard.hpp \
		includes/FileSet/PdbFileSpace/pdbformula.hpp \
		includes/FileSet/PdbFileSpace/pdbhelixcard.hpp \
		includes/FileSet/PdbFileSpace/pdbhelix.hpp \
		includes/FileSet/PdbFileSpace/pdbhelixresidue.hpp \
		includes/FileSet/PdbFileSpace/pdbsheetcard.hpp \
		includes/FileSet/PdbFileSpace/pdbsheet.hpp \
		includes/FileSet/PdbFileSpace/pdbsheetstrand.hpp \
		includes/FileSet/PdbFileSpace/pdbsheetstrandresidue.hpp \
		includes/FileSet/PdbFileSpace/pdbdisulfidebondcard.hpp \
		includes/FileSet/PdbFileSpace/pdbdisulfideresiduebond.hpp \
		includes/FileSet/PdbFileSpace/pdbdisulfideresidue.hpp \
		includes/FileSet/PdbFileSpace/pdblinkcard.hpp \
		includes/FileSet/PdbFileSpace/pdblink.hpp \
		includes/FileSet/PdbFileSpace/pdblinkresidue.hpp \
		includes/FileSet/PdbFileSpace/pdbsitecard.hpp \
		includes/FileSet/PdbFileSpace/pdbsite.hpp \
		includes/FileSet/PdbFileSpace/pdbsiteresidue.hpp \
		includes/FileSet/PdbFileSpace/pdbcrystallographiccard.hpp \
		includes/FileSet/PdbFileSpace/pdboriginxncard.hpp \
		includes/FileSet/PdbFileSpace/pdboriginxn.hpp \
		includes/Geometry/coordinate.hpp \
		includes/FileSet/PdbFileSpace/pdbscalencard.hpp \
		includes/FileSet/PdbFileSpace/pdbscalen.hpp \
		includes/FileSet/PdbFileSpace/pdbmatrixncard.hpp \
		includes/FileSet/PdbFileSpace/pdbmatrixn.hpp \
		includes/FileSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/FileSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/FileSet/PdbFileSpace/pdbmodel.hpp \
		includes/FileSet/PdbFileSpace/pdbatom.hpp \
		includes/FileSet/PdbFileSpace/pdbatomcard.hpp \
		includes/FileSet/PdbFileSpace/pdbheterogenatomcard.hpp \
		includes/FileSet/PdbFileSpace/pdbconnectcard.hpp \
		includes/FileSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/FileSet/PdbFileSpace/pdbresidue.hpp \
		includes/utils.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbfile.o src/FileSet/PdbFileSpace/pdbfile.cc

pdbfileprocessingexception.o: src/FileSet/PdbFileSpace/pdbfileprocessingexception.cc includes/common.hpp \
		includes/FileSet/PdbFileSpace/pdbfileprocessingexception.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbfileprocessingexception.o src/FileSet/PdbFileSpace/pdbfileprocessingexception.cc

pdbformula.o: src/FileSet/PdbFileSpace/pdbformula.cc includes/FileSet/PdbFileSpace/pdbformula.hpp \
		includes/common.hpp \
		includes/utils.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbformula.o src/FileSet/PdbFileSpace/pdbformula.cc

pdbformulacard.o: src/FileSet/PdbFileSpace/pdbformulacard.cc includes/FileSet/PdbFileSpace/pdbformulacard.hpp \
		includes/FileSet/PdbFileSpace/pdbformula.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbformulacard.o src/FileSet/PdbFileSpace/pdbformulacard.cc

pdbheadercard.o: src/FileSet/PdbFileSpace/pdbheadercard.cc includes/FileSet/PdbFileSpace/pdbheadercard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbheadercard.o src/FileSet/PdbFileSpace/pdbheadercard.cc

pdbhelix.o: src/FileSet/PdbFileSpace/pdbhelix.cc includes/FileSet/PdbFileSpace/pdbhelix.hpp \
		includes/FileSet/PdbFileSpace/pdbhelixresidue.hpp \
		includes/common.hpp \
		includes/utils.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbhelix.o src/FileSet/PdbFileSpace/pdbhelix.cc

pdbhelixcard.o: src/FileSet/PdbFileSpace/pdbhelixcard.cc includes/FileSet/PdbFileSpace/pdbhelixcard.hpp \
		includes/FileSet/PdbFileSpace/pdbhelix.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbhelixcard.o src/FileSet/PdbFileSpace/pdbhelixcard.cc

pdbhelixresidue.o: src/FileSet/PdbFileSpace/pdbhelixresidue.cc includes/FileSet/PdbFileSpace/pdbhelixresidue.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbhelixresidue.o src/FileSet/PdbFileSpace/pdbhelixresidue.cc

pdbheterogen.o: src/FileSet/PdbFileSpace/pdbheterogen.cc includes/FileSet/PdbFileSpace/pdbheterogen.hpp \
		includes/common.hpp \
		includes/utils.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbheterogen.o src/FileSet/PdbFileSpace/pdbheterogen.cc

pdbheterogenatomcard.o: src/FileSet/PdbFileSpace/pdbheterogenatomcard.cc includes/FileSet/PdbFileSpace/pdbheterogenatomcard.hpp \
		includes/FileSet/PdbFileSpace/pdbatom.hpp \
		includes/Geometry/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbheterogenatomcard.o src/FileSet/PdbFileSpace/pdbheterogenatomcard.cc

pdbheterogencard.o: src/FileSet/PdbFileSpace/pdbheterogencard.cc includes/FileSet/PdbFileSpace/pdbheterogencard.hpp \
		includes/FileSet/PdbFileSpace/pdbheterogen.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbheterogencard.o src/FileSet/PdbFileSpace/pdbheterogencard.cc

pdbheterogenname.o: src/FileSet/PdbFileSpace/pdbheterogenname.cc includes/FileSet/PdbFileSpace/pdbheterogenname.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbheterogenname.o src/FileSet/PdbFileSpace/pdbheterogenname.cc

pdbheterogennamecard.o: src/FileSet/PdbFileSpace/pdbheterogennamecard.cc includes/FileSet/PdbFileSpace/pdbheterogennamecard.hpp \
		includes/FileSet/PdbFileSpace/pdbheterogenname.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbheterogennamecard.o src/FileSet/PdbFileSpace/pdbheterogennamecard.cc

pdbheterogensynonym.o: src/FileSet/PdbFileSpace/pdbheterogensynonym.cc includes/FileSet/PdbFileSpace/pdbheterogensynonym.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbheterogensynonym.o src/FileSet/PdbFileSpace/pdbheterogensynonym.cc

pdbheterogensynonymcard.o: src/FileSet/PdbFileSpace/pdbheterogensynonymcard.cc includes/FileSet/PdbFileSpace/pdbheterogensynonymcard.hpp \
		includes/FileSet/PdbFileSpace/pdbheterogensynonym.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbheterogensynonymcard.o src/FileSet/PdbFileSpace/pdbheterogensynonymcard.cc

pdblink.o: src/FileSet/PdbFileSpace/pdblink.cc includes/FileSet/PdbFileSpace/pdblink.hpp \
		includes/FileSet/PdbFileSpace/pdblinkresidue.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdblink.o src/FileSet/PdbFileSpace/pdblink.cc

pdblinkcard.o: src/FileSet/PdbFileSpace/pdblinkcard.cc includes/FileSet/PdbFileSpace/pdblinkcard.hpp \
		includes/FileSet/PdbFileSpace/pdblink.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdblinkcard.o src/FileSet/PdbFileSpace/pdblinkcard.cc

pdblinkresidue.o: src/FileSet/PdbFileSpace/pdblinkresidue.cc includes/FileSet/PdbFileSpace/pdblinkresidue.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdblinkresidue.o src/FileSet/PdbFileSpace/pdblinkresidue.cc

pdbmatrixn.o: src/FileSet/PdbFileSpace/pdbmatrixn.cc includes/FileSet/PdbFileSpace/pdbmatrixn.hpp \
		includes/Geometry/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbmatrixn.o src/FileSet/PdbFileSpace/pdbmatrixn.cc

pdbmatrixncard.o: src/FileSet/PdbFileSpace/pdbmatrixncard.cc includes/FileSet/PdbFileSpace/pdbmatrixn.hpp \
		includes/Geometry/coordinate.hpp \
		includes/FileSet/PdbFileSpace/pdbmatrixncard.hpp \
		includes/utils.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbmatrixncard.o src/FileSet/PdbFileSpace/pdbmatrixncard.cc

pdbmodel.o: src/FileSet/PdbFileSpace/pdbmodel.cc includes/FileSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/FileSet/PdbFileSpace/pdbmodel.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbmodel.o src/FileSet/PdbFileSpace/pdbmodel.cc

pdbmodelcard.o: src/FileSet/PdbFileSpace/pdbmodelcard.cc includes/FileSet/PdbFileSpace/pdbmodel.hpp \
		includes/FileSet/PdbFileSpace/pdbmodelcard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbmodelcard.o src/FileSet/PdbFileSpace/pdbmodelcard.cc

pdbmodelresidueset.o: src/FileSet/PdbFileSpace/pdbmodelresidueset.cc includes/FileSet/PdbFileSpace/pdbmodelresidueset.hpp \
		includes/FileSet/PdbFileSpace/pdbatomcard.hpp \
		includes/FileSet/PdbFileSpace/pdbheterogenatomcard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbmodelresidueset.o src/FileSet/PdbFileSpace/pdbmodelresidueset.cc

pdbmodeltypecard.o: src/FileSet/PdbFileSpace/pdbmodeltypecard.cc includes/FileSet/PdbFileSpace/pdbmodeltypecard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbmodeltypecard.o src/FileSet/PdbFileSpace/pdbmodeltypecard.cc

pdbnummodelcard.o: src/FileSet/PdbFileSpace/pdbnummodelcard.cc includes/FileSet/PdbFileSpace/pdbnummodelcard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbnummodelcard.o src/FileSet/PdbFileSpace/pdbnummodelcard.cc

pdboriginxn.o: src/FileSet/PdbFileSpace/pdboriginxn.cc includes/FileSet/PdbFileSpace/pdboriginxn.hpp \
		includes/Geometry/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdboriginxn.o src/FileSet/PdbFileSpace/pdboriginxn.cc

pdboriginxncard.o: src/FileSet/PdbFileSpace/pdboriginxncard.cc includes/FileSet/PdbFileSpace/pdboriginxn.hpp \
		includes/Geometry/coordinate.hpp \
		includes/FileSet/PdbFileSpace/pdboriginxncard.hpp \
		includes/utils.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdboriginxncard.o src/FileSet/PdbFileSpace/pdboriginxncard.cc

pdbresidue.o: src/FileSet/PdbFileSpace/pdbresidue.cc includes/FileSet/PdbFileSpace/pdbresidue.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbresidue.o src/FileSet/PdbFileSpace/pdbresidue.cc

pdbresiduemodification.o: src/FileSet/PdbFileSpace/pdbresiduemodification.cc includes/FileSet/PdbFileSpace/pdbresiduemodification.hpp \
		includes/common.hpp \
		includes/utils.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbresiduemodification.o src/FileSet/PdbFileSpace/pdbresiduemodification.cc

pdbresiduemodificationcard.o: src/FileSet/PdbFileSpace/pdbresiduemodificationcard.cc includes/FileSet/PdbFileSpace/pdbresiduemodificationcard.hpp \
		includes/FileSet/PdbFileSpace/pdbresiduemodification.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbresiduemodificationcard.o src/FileSet/PdbFileSpace/pdbresiduemodificationcard.cc

pdbresiduesequence.o: src/FileSet/PdbFileSpace/pdbresiduesequence.cc includes/FileSet/PdbFileSpace/pdbresiduesequence.hpp \
		includes/common.hpp \
		includes/utils.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbresiduesequence.o src/FileSet/PdbFileSpace/pdbresiduesequence.cc

pdbresiduesequencecard.o: src/FileSet/PdbFileSpace/pdbresiduesequencecard.cc includes/FileSet/PdbFileSpace/pdbresiduesequencecard.hpp \
		includes/FileSet/PdbFileSpace/pdbresiduesequence.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbresiduesequencecard.o src/FileSet/PdbFileSpace/pdbresiduesequencecard.cc

pdbscalen.o: src/FileSet/PdbFileSpace/pdbscalen.cc includes/FileSet/PdbFileSpace/pdbscalen.hpp \
		includes/Geometry/coordinate.hpp \
		includes/utils.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbscalen.o src/FileSet/PdbFileSpace/pdbscalen.cc

pdbscalencard.o: src/FileSet/PdbFileSpace/pdbscalencard.cc includes/FileSet/PdbFileSpace/pdbscalen.hpp \
		includes/Geometry/coordinate.hpp \
		includes/FileSet/PdbFileSpace/pdbscalencard.hpp \
		includes/utils.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbscalencard.o src/FileSet/PdbFileSpace/pdbscalencard.cc

pdbsheet.o: src/FileSet/PdbFileSpace/pdbsheet.cc includes/FileSet/PdbFileSpace/pdbsheet.hpp \
		includes/FileSet/PdbFileSpace/pdbsheetstrand.hpp \
		includes/FileSet/PdbFileSpace/pdbsheetstrandresidue.hpp \
		includes/common.hpp \
		includes/utils.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbsheet.o src/FileSet/PdbFileSpace/pdbsheet.cc

pdbsheetcard.o: src/FileSet/PdbFileSpace/pdbsheetcard.cc includes/FileSet/PdbFileSpace/pdbsheetcard.hpp \
		includes/FileSet/PdbFileSpace/pdbsheet.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbsheetcard.o src/FileSet/PdbFileSpace/pdbsheetcard.cc

pdbsheetstrand.o: src/FileSet/PdbFileSpace/pdbsheetstrand.cc includes/FileSet/PdbFileSpace/pdbsheetstrand.hpp \
		includes/FileSet/PdbFileSpace/pdbsheetstrandresidue.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbsheetstrand.o src/FileSet/PdbFileSpace/pdbsheetstrand.cc

pdbsheetstrandresidue.o: src/FileSet/PdbFileSpace/pdbsheetstrandresidue.cc includes/FileSet/PdbFileSpace/pdbsheetstrandresidue.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbsheetstrandresidue.o src/FileSet/PdbFileSpace/pdbsheetstrandresidue.cc

pdbsite.o: src/FileSet/PdbFileSpace/pdbsite.cc includes/FileSet/PdbFileSpace/pdbsite.hpp \
		includes/FileSet/PdbFileSpace/pdbsiteresidue.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbsite.o src/FileSet/PdbFileSpace/pdbsite.cc

pdbsitecard.o: src/FileSet/PdbFileSpace/pdbsitecard.cc includes/FileSet/PdbFileSpace/pdbsite.hpp \
		includes/FileSet/PdbFileSpace/pdbsitecard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbsitecard.o src/FileSet/PdbFileSpace/pdbsitecard.cc

pdbsiteresidue.o: src/FileSet/PdbFileSpace/pdbsiteresidue.cc includes/FileSet/PdbFileSpace/pdbsiteresidue.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbsiteresidue.o src/FileSet/PdbFileSpace/pdbsiteresidue.cc

pdbtitlecard.o: src/FileSet/PdbFileSpace/pdbtitlecard.cc includes/FileSet/PdbFileSpace/pdbtitlecard.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbtitlecard.o src/FileSet/PdbFileSpace/pdbtitlecard.cc

topologyangle.o: src/FileSet/TopologyFileSpace/topologyangle.cc includes/FileSet/TopologyFileSpace/topologyangle.hpp \
		includes/FileSet/TopologyFileSpace/topologyangletype.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o topologyangle.o src/FileSet/TopologyFileSpace/topologyangle.cc

topologyangletype.o: src/FileSet/TopologyFileSpace/topologyangletype.cc includes/FileSet/TopologyFileSpace/topologyangletype.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o topologyangletype.o src/FileSet/TopologyFileSpace/topologyangletype.cc

topologyassembly.o: src/FileSet/TopologyFileSpace/topologyassembly.cc includes/FileSet/TopologyFileSpace/topologyassembly.hpp \
		includes/FileSet/TopologyFileSpace/topologyresidue.hpp \
		includes/FileSet/TopologyFileSpace/topologyatom.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o topologyassembly.o src/FileSet/TopologyFileSpace/topologyassembly.cc

topologyatom.o: src/FileSet/TopologyFileSpace/topologyatom.cc includes/FileSet/TopologyFileSpace/topologyatom.hpp \
		includes/FileSet/TopologyFileSpace/topologyatompair.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o topologyatom.o src/FileSet/TopologyFileSpace/topologyatom.cc

topologyatompair.o: src/FileSet/TopologyFileSpace/topologyatompair.cc includes/FileSet/TopologyFileSpace/topologyatompair.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o topologyatompair.o src/FileSet/TopologyFileSpace/topologyatompair.cc

topologybond.o: src/FileSet/TopologyFileSpace/topologybond.cc includes/FileSet/TopologyFileSpace/topologybond.hpp \
		includes/FileSet/TopologyFileSpace/topologybondtype.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o topologybond.o src/FileSet/TopologyFileSpace/topologybond.cc

topologybondtype.o: src/FileSet/TopologyFileSpace/topologybondtype.cc includes/FileSet/TopologyFileSpace/topologybondtype.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o topologybondtype.o src/FileSet/TopologyFileSpace/topologybondtype.cc

topologydihedral.o: src/FileSet/TopologyFileSpace/topologydihedral.cc includes/FileSet/TopologyFileSpace/topologydihedral.hpp \
		includes/FileSet/TopologyFileSpace/topologydihedraltype.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o topologydihedral.o src/FileSet/TopologyFileSpace/topologydihedral.cc

topologydihedraltype.o: src/FileSet/TopologyFileSpace/topologydihedraltype.cc includes/FileSet/TopologyFileSpace/topologydihedraltype.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o topologydihedraltype.o src/FileSet/TopologyFileSpace/topologydihedraltype.cc

topologyfile.o: src/FileSet/TopologyFileSpace/topologyfile.cc includes/FileSet/TopologyFileSpace/topologyfile.hpp \
		includes/FileSet/TopologyFileSpace/topologyatompair.hpp \
		includes/FileSet/TopologyFileSpace/topologybondtype.hpp \
		includes/FileSet/TopologyFileSpace/topologyangletype.hpp \
		includes/FileSet/TopologyFileSpace/topologydihedraltype.hpp \
		includes/FileSet/TopologyFileSpace/topologyassembly.hpp \
		includes/FileSet/TopologyFileSpace/topologyatom.hpp \
		includes/FileSet/TopologyFileSpace/topologybond.hpp \
		includes/FileSet/TopologyFileSpace/topologyangle.hpp \
		includes/FileSet/TopologyFileSpace/topologydihedral.hpp \
		includes/FileSet/TopologyFileSpace/topologyresidue.hpp \
		includes/FileSet/TopologyFileSpace/topologyfileprocessingexception.hpp \
		includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o topologyfile.o src/FileSet/TopologyFileSpace/topologyfile.cc

topologyfileprocessingexception.o: src/FileSet/TopologyFileSpace/topologyfileprocessingexception.cc includes/common.hpp \
		includes/FileSet/TopologyFileSpace/topologyfileprocessingexception.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o topologyfileprocessingexception.o src/FileSet/TopologyFileSpace/topologyfileprocessingexception.cc

topologyresidue.o: src/FileSet/TopologyFileSpace/topologyresidue.cc includes/FileSet/TopologyFileSpace/topologyresidue.hpp \
		includes/FileSet/TopologyFileSpace/topologyatom.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o topologyresidue.o src/FileSet/TopologyFileSpace/topologyresidue.cc

angle.o: src/Geometry/InternalCoordinate/angle.cc includes/Geometry/InternalCoordinate/angle.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o angle.o src/Geometry/InternalCoordinate/angle.cc

dihedral.o: src/Geometry/InternalCoordinate/dihedral.cc includes/Geometry/InternalCoordinate/dihedral.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o dihedral.o src/Geometry/InternalCoordinate/dihedral.cc

distance.o: src/Geometry/InternalCoordinate/distance.cc includes/Geometry/InternalCoordinate/distance.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o distance.o src/Geometry/InternalCoordinate/distance.cc

libraryfile.o: src/ParameterSet/LibraryFileSpace/libraryfile.cc includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o libraryfile.o src/ParameterSet/LibraryFileSpace/libraryfile.cc

libraryfileatom.o: src/ParameterSet/LibraryFileSpace/libraryfileatom.cc includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/Geometry/coordinate.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o libraryfileatom.o src/ParameterSet/LibraryFileSpace/libraryfileatom.cc

libraryfileprocessingexception.o: src/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.cc includes/common.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o libraryfileprocessingexception.o src/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.cc

libraryfileresidue.o: src/ParameterSet/LibraryFileSpace/libraryfileresidue.cc includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/Geometry/coordinate.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o libraryfileresidue.o src/ParameterSet/LibraryFileSpace/libraryfileresidue.cc

parameterfile.o: src/ParameterSet/ParameterFileSpace/parameterfile.cc includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfile.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileprocessingexception.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o parameterfile.o src/ParameterSet/ParameterFileSpace/parameterfile.cc

parameterfileangle.o: src/ParameterSet/ParameterFileSpace/parameterfileangle.cc includes/common.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o parameterfileangle.o src/ParameterSet/ParameterFileSpace/parameterfileangle.cc

parameterfileatom.o: src/ParameterSet/ParameterFileSpace/parameterfileatom.cc includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o parameterfileatom.o src/ParameterSet/ParameterFileSpace/parameterfileatom.cc

parameterfilebond.o: src/ParameterSet/ParameterFileSpace/parameterfilebond.cc includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o parameterfilebond.o src/ParameterSet/ParameterFileSpace/parameterfilebond.cc

parameterfiledihedral.o: src/ParameterSet/ParameterFileSpace/parameterfiledihedral.cc includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp \
		includes/common.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o parameterfiledihedral.o src/ParameterSet/ParameterFileSpace/parameterfiledihedral.cc

parameterfiledihedralterm.o: src/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.cc includes/common.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o parameterfiledihedralterm.o src/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.cc

parameterfileprocessingexception.o: src/ParameterSet/ParameterFileSpace/parameterfileprocessingexception.cc includes/common.hpp \
		includes/ParameterSet/ParameterFileSpace/parameterfileprocessingexception.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o parameterfileprocessingexception.o src/ParameterSet/ParameterFileSpace/parameterfileprocessingexception.cc

prepfile.o: src/ParameterSet/PrepFileSpace/prepfile.cc includes/common.hpp \
		includes/utils.hpp \
		includes/Geometry/coordinate.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o prepfile.o src/ParameterSet/PrepFileSpace/prepfile.cc

prepfileatom.o: src/ParameterSet/PrepFileSpace/prepfileatom.cc includes/common.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o prepfileatom.o src/ParameterSet/PrepFileSpace/prepfileatom.cc

prepfileprocessingexception.o: src/ParameterSet/PrepFileSpace/prepfileprocessingexception.cc includes/common.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o prepfileprocessingexception.o src/ParameterSet/PrepFileSpace/prepfileprocessingexception.cc

prepfileresidue.o: src/ParameterSet/PrepFileSpace/prepfileresidue.cc includes/utils.hpp \
		includes/common.hpp \
		includes/Geometry/coordinate.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o prepfileresidue.o src/ParameterSet/PrepFileSpace/prepfileresidue.cc

pdbpreprocessor.o: src/Resolver/PdbPreprocessor/pdbpreprocessor.cc includes/Resolver/PdbPreprocessor/pdbpreprocessor.hpp \
		includes/FileSet/PdbFileSpace/pdbresidue.hpp \
		includes/FileSet/PdbFileSpace/pdbfile.hpp \
		includes/FileSet/PdbFileSpace/pdbatom.hpp \
		includes/Geometry/coordinate.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfile.hpp \
		includes/common.hpp \
		includes/ParameterSet/PrepFileSpace/prepfile.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessoralternateresidue.hpp \
		includes/Resolver/PdbPreprocessor/pdbpreprocessorresidueinfo.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp \
		includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp \
		includes/ParameterSet/PrepFileSpace/prepfileatom.hpp \
		includes/FileSet/PdbFileSpace/pdbatomcard.hpp \
		includes/FileSet/PdbFileSpace/pdbfileprocessingexception.hpp \
		includes/utils.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbpreprocessor.o src/Resolver/PdbPreprocessor/pdbpreprocessor.cc

pdbpreprocessoralternateresidue.o: src/Resolver/PdbPreprocessor/pdbpreprocessoralternateresidue.cc includes/Resolver/PdbPreprocessor/pdbpreprocessoralternateresidue.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbpreprocessoralternateresidue.o src/Resolver/PdbPreprocessor/pdbpreprocessoralternateresidue.cc

pdbpreprocessorchaintermination.o: src/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.cc includes/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbpreprocessorchaintermination.o src/Resolver/PdbPreprocessor/pdbpreprocessorchaintermination.cc

pdbpreprocessordisulfidebond.o: src/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.cc includes/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbpreprocessordisulfidebond.o src/Resolver/PdbPreprocessor/pdbpreprocessordisulfidebond.cc

pdbpreprocessorhistidinemapping.o: src/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.cc includes/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbpreprocessorhistidinemapping.o src/Resolver/PdbPreprocessor/pdbpreprocessorhistidinemapping.cc

pdbpreprocessormissingresidue.o: src/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.cc includes/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.hpp \
		includes/common.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbpreprocessormissingresidue.o src/Resolver/PdbPreprocessor/pdbpreprocessormissingresidue.cc

pdbpreprocessorreplacedhydrogen.o: src/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.cc includes/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbpreprocessorreplacedhydrogen.o src/Resolver/PdbPreprocessor/pdbpreprocessorreplacedhydrogen.cc

pdbpreprocessorresidueinfo.o: src/Resolver/PdbPreprocessor/pdbpreprocessorresidueinfo.cc includes/Resolver/PdbPreprocessor/pdbpreprocessorresidueinfo.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbpreprocessorresidueinfo.o src/Resolver/PdbPreprocessor/pdbpreprocessorresidueinfo.cc

pdbpreprocessorunrecognizedheavyatom.o: src/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.cc includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbpreprocessorunrecognizedheavyatom.o src/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedheavyatom.cc

pdbpreprocessorunrecognizedresidue.o: src/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.cc includes/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pdbpreprocessorunrecognizedresidue.o src/Resolver/PdbPreprocessor/pdbpreprocessorunrecognizedresidue.cc

####### Install

install:   FORCE

uninstall:   FORCE

FORCE:

