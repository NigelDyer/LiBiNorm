################################################################################
#
# Make all or Make debug makes the debug version
# Make release to make the release version
#
#   Background information on the structure of the makefile can be found at
#   http://www2.warwick.ac.uk/fac/sci/systemsbiology/staff/dyer/software/gccmakefiles/
#
################################################################################

# Standard rules
CCC = g++

# Flags required by all stages of C++ compiler
CCCALLFLAGS= -std=gnu++11 -DBAM_LIBRARY

# Directory information
BAMTOOLSDIR = bamtools/src/
BIOINFORMATICSLIBDIR = bioinformaticsLib/
MCMCLIBDIR = mcmcLib/
HISAT2DIR = hisat2Lib/
LIBINORMSRCDIR = LiBiNormSrc

RELDIR = Release
DEBUGDIR = Debug

ifeq "$(findstring release, $(MAKECMDGOALS))" ""
BUILD=$(DEBUGDIR)
CCCALLFLAGS += -g3 -O0 -D_DEBUG -Wall -Wno-unknown-pragmas 
else
BUILD=$(RELDIR)
CCCALLFLAGS += -O3
endif

INCLUDES= -I../$(BAMTOOLSDIR) -I../$(BIOINFORMATICSLIBDIR) -I$(MCMCLIBDIR) -I$(HISAT2DIR) -I$(LIBINORMSRCDIR)

################################################################################
# Outputs of this Makefile

LIBINORM = LiBiNorm

ifeq ($(OS),Windows_NT)
EXE = .exe
endif

LIBINORMEXE = $(BUILD)/$(LIBINORM)$(EXE) 

TARGS =  $(LIBINORMEXE) 

################################################################################
# Libraries to be linked. 

LIBS        = -lz -pthread
LIBPATH     = 

################################################################################
# Source files and the build specific outputs
#  All of the .cpp files in MCMCLIBDIR are included in the build and the objects are in <build>/mcmcLib
#  The files that are from the external directories are handled slightly differently so that their object files
#  are also placed within the $(BUILD) directory and so are deleted with a make clean

LIBINORMSRC = $(LIBINORMSRCDIR)/$(LIBINORM).cpp

LIBINORMSRCEX = $(filter-out $(LIBINORMSRC), $(shell find $(LIBINORMSRCDIR) -name *.cpp) )

MCMCLIBSRC =  $(shell find $(MCMCLIBDIR) -name *.cpp)

BIOLIBSRC = $(shell find ../$(BIOINFORMATICSLIBDIR) -name *.cpp)

BAMTOOLSSRC = $(filter-out ../$(BAMTOOLSDIR)api/internal/io/TcpSocketEngine_win_p.cpp, \
	$(shell find ../$(BAMTOOLSDIR)api -name *.cpp) ) \
	$(addprefix ../$(BAMTOOLSDIR), toolkit/bamtools_sort.cpp utils/bamtools_options.cpp )

COREOBJS :=  $(addprefix $(BUILD)/, $(LIBINORMSRCEX:%.cpp=%.o) $(MCMCLIBSRC:%.cpp=%.o) \
		$(subst ../,,$(BIOLIBSRC:%.cpp=%.o) $(BAMTOOLSSRC:%.cpp=%.o) ) ) 


################################################################################
# For building necessary outout directories.  The sort method is used to remove duplicated entries

DIRS = $(sort $(dir $(COREOBJS)  ) ) 

$(DIRS) :
	mkdir -p $@


##################################################################
#
#	Instructions for building release and debug object files  These are dependant on the Makefile so Makefile changes
#	force a rebuild.  Specific rules for source files in libraries in other directories

$(BUILD)/%.o : %.cpp  Makefile
	$(CCC) -c $(CCCALLFLAGS) $(INCLUDES) -o $@ $<

$(BUILD)/$(BIOINFORMATICSLIBDIR)%.o : ../$(BIOINFORMATICSLIBDIR)%.cpp Makefile
	$(CCC) -c $(CCCALLFLAGS) $(INCLUDES) -o $@ $<

$(BUILD)/$(BAMTOOLSDIR)%.o : ../$(BAMTOOLSDIR)%.cpp  Makefile
	$(CCC) -c $(CCCALLFLAGS) $(INCLUDES) -Wno-sign-compare -o $@ $<

################################################################################
# The main builds

debug : all
	
release : all    

all:  $(DIRS) $(TARGS)
	@echo "%% $(BUILD) LiBiNorm code built"

#	The final make rule
$(LIBINORMEXE) :$(BUILD)/$(LIBINORMSRC:%.cpp=%.o) $(COREOBJS)
	$(CCC)  $^ -o $@ $(CCCALLFLAGS) $(INCLUDES) $(LIBPATH) $(LIBS)

#	All intermediate and final build products go in one directory and its sub directories 
#	for each of the build types, making a make clean very simple
	
clean : 
	rm -r -f $(RELDIR)
	rm -r -f $(DEBUGDIR)

#	do make depend to update dependancies.  This makes a dependancy list that is dynamically dependant 
#	on the build type.  There are three make depends, one for all of the sources within this directory (SOURCES)
#	and one for the sources that are in other library directories (BIOLIBSRC & BAMTOOLSSRC)
#	The dummy that contains XXZZ is prefixed is part of the process of dealing with the fact that the object files
#	are not in the same directory as the source files.  The dependancies work without removing it  
#   but it looks neater if they are removed.

depend :
	makedepend  -Y $(CCCAALLFLAGS) $(INCLUDES) $(LIBINORMSRC) $(MCMCLIBSRC) $(LIBINORMSRCEX) -p'$$(BUILD)/'
	makedepend  -Y -a $(CCCAALLFLAGS) $(INCLUDES) $(BIOLIBSRC) $(BAMTOOLSSRC) -p'$$(BUILD)/XXZZ/'
	sed -i -- 's/\/XXZZ\/..//g' Makefile	

# DO NOT DELETE THIS LINE -- make depend depends on it.

$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/LiBiNorm.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/LiBiNormCore.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/GeneCountData.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/libParser.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: mcmcLib/mcmc.h mcmcLib/params.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/nelderMeadOptimiser.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/Options.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/ModelData.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/LiBiCount.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/FeatureFileEx.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/Regions.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/src/api/BamReader.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/src/api/api_global.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/src/api/BamAux.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/featureFile.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: ../bamtools/src/api/BamWriter.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/LiBiDedup.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/LiBiConv.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/LiBiTools.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/LiBiVariation.h
$(BUILD)/LiBiNormSrc/LiBiNorm.o: LiBiNormSrc/MakeFastq.h
$(BUILD)/mcmcLib/mcmc.o: ../bioinformaticsLib/rand.h
$(BUILD)/mcmcLib/mcmc.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/mcmcLib/mcmc.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/mcmcLib/mcmc.o: ../bioinformaticsLib/printEx.h
$(BUILD)/mcmcLib/mcmc.o: ../bioinformaticsLib/inQuotes.h mcmcLib/mcmc.h
$(BUILD)/mcmcLib/mcmc.o: mcmcLib/params.h
$(BUILD)/mcmcLib/mcmc.o: ../bioinformaticsLib/nelderMeadOptimiser.h
$(BUILD)/mcmcLib/mcmc.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/mcmcLib/params.o: ../bioinformaticsLib/rand.h
$(BUILD)/mcmcLib/params.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/mcmcLib/params.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/mcmcLib/params.o: ../bioinformaticsLib/printEx.h
$(BUILD)/mcmcLib/params.o: ../bioinformaticsLib/inQuotes.h mcmcLib/params.h
$(BUILD)/mcmcLib/params.o: ../bioinformaticsLib/nelderMeadOptimiser.h
$(BUILD)/mcmcLib/params.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: LiBiNormSrc/FeatureFileEx.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: LiBiNormSrc/Regions.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bamtools/src/api/BamReader.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bamtools/src/api/api_global.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bamtools/src/api/BamAux.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: LiBiNormSrc/Options.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bioinformaticsLib/featureFile.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: LiBiNormSrc/GeneCountData.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bioinformaticsLib/libParser.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: mcmcLib/mcmc.h mcmcLib/params.h
$(BUILD)/LiBiNormSrc/FeatureFileEx.o: ../bioinformaticsLib/nelderMeadOptimiser.h
$(BUILD)/LiBiNormSrc/GeneCountData.o: LiBiNormSrc/Options.h
$(BUILD)/LiBiNormSrc/GeneCountData.o: LiBiNormSrc/GeneCountData.h
$(BUILD)/LiBiNormSrc/GeneCountData.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/GeneCountData.o: ../bioinformaticsLib/libParser.h
$(BUILD)/LiBiNormSrc/GeneCountData.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/GeneCountData.o: mcmcLib/mcmc.h mcmcLib/params.h
$(BUILD)/LiBiNormSrc/GeneCountData.o: ../bioinformaticsLib/nelderMeadOptimiser.h
$(BUILD)/LiBiNormSrc/GeneCountData.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/GeneCountData.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/GeneCountData.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/LiBiNormSrc/GeneCountData.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bamtools/src/api/BamReader.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bamtools/src/api/api_global.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bamtools/src/api/BamAux.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: LiBiNormSrc/LiBiConv.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: LiBiNormSrc/FeatureFileEx.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: LiBiNormSrc/Regions.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: LiBiNormSrc/Options.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bioinformaticsLib/featureFile.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: LiBiNormSrc/GeneCountData.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bioinformaticsLib/libParser.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: mcmcLib/mcmc.h mcmcLib/params.h
$(BUILD)/LiBiNormSrc/LiBiConv.o: ../bioinformaticsLib/nelderMeadOptimiser.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/bamAlignmentEx.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: LiBiNormSrc/Regions.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/src/api/BamReader.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/src/api/api_global.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/src/api/BamAux.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: LiBiNormSrc/Options.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/libParser.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: LiBiNormSrc/LiBiCount.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: LiBiNormSrc/FeatureFileEx.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/featureFile.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: LiBiNormSrc/GeneCountData.h mcmcLib/mcmc.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: mcmcLib/params.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bioinformaticsLib/nelderMeadOptimiser.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: ../bamtools/src/api/BamWriter.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: LiBiNormSrc/LiBiNorm.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: LiBiNormSrc/LiBiNormCore.h
$(BUILD)/LiBiNormSrc/LiBiCount.o: LiBiNormSrc/ModelData.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: LiBiNormSrc/LiBiDedup.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/src/api/BamReader.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/src/api/api_global.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/src/api/BamAux.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bamtools/src/api/BamWriter.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: LiBiNormSrc/Regions.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/LiBiDedup.o: LiBiNormSrc/Options.h
$(BUILD)/LiBiNormSrc/LiBiNormCore.o: LiBiNormSrc/LiBiNormCore.h
$(BUILD)/LiBiNormSrc/LiBiNormCore.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/LiBiNormCore.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/LiBiNormCore.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/LiBiNormCore.o: LiBiNormSrc/GeneCountData.h
$(BUILD)/LiBiNormSrc/LiBiNormCore.o: ../bioinformaticsLib/libParser.h
$(BUILD)/LiBiNormSrc/LiBiNormCore.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/LiBiNormCore.o: mcmcLib/mcmc.h mcmcLib/params.h
$(BUILD)/LiBiNormSrc/LiBiNormCore.o: ../bioinformaticsLib/nelderMeadOptimiser.h
$(BUILD)/LiBiNormSrc/LiBiNormCore.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/LiBiNormSrc/LiBiNormCore.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/LiBiNormCore.o: LiBiNormSrc/Options.h
$(BUILD)/LiBiNormSrc/LiBiNormCore.o: LiBiNormSrc/ModelData.h
$(BUILD)/LiBiNormSrc/LiBiNormCore.o: LiBiNormSrc/LiBiOptimiser.h
$(BUILD)/LiBiNormSrc/LiBiOptimiser.o: ../bioinformaticsLib/rand.h
$(BUILD)/LiBiNormSrc/LiBiOptimiser.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/LiBiNormSrc/LiBiOptimiser.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/LiBiOptimiser.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/LiBiOptimiser.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/LiBiOptimiser.o: LiBiNormSrc/Options.h
$(BUILD)/LiBiNormSrc/LiBiOptimiser.o: LiBiNormSrc/LiBiOptimiser.h
$(BUILD)/LiBiNormSrc/LiBiOptimiser.o: ../bioinformaticsLib/nelderMeadOptimiser.h
$(BUILD)/LiBiNormSrc/LiBiOptimiser.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/LiBiOptimiser.o: LiBiNormSrc/ModelData.h mcmcLib/mcmc.h
$(BUILD)/LiBiNormSrc/LiBiOptimiser.o: mcmcLib/params.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bamtools/src/api/BamReader.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bamtools/src/api/api_global.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bamtools/src/api/BamAux.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bioinformaticsLib/fastaFile.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: LiBiNormSrc/Options.h LiBiNormSrc/refSeqs.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: hisat2Lib/reference.h hisat2Lib/ref_read.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: hisat2Lib/alphabet.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: hisat2Lib/assert_helpers.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: hisat2Lib/word_io.h hisat2Lib/endian_swap.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: hisat2Lib/hisat2Lib.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: LiBiNormSrc/FeatureFileEx.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: LiBiNormSrc/Regions.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bioinformaticsLib/featureFile.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: LiBiNormSrc/GeneCountData.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bioinformaticsLib/libParser.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: mcmcLib/mcmc.h mcmcLib/params.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: ../bioinformaticsLib/nelderMeadOptimiser.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: LiBiNormSrc/LiBiTools.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: LiBiNormSrc/LiBiNormCore.h
$(BUILD)/LiBiNormSrc/LiBiTools.o: LiBiNormSrc/ModelData.h
$(BUILD)/LiBiNormSrc/LiBiVariation.o: ../bioinformaticsLib/libParser.h
$(BUILD)/LiBiNormSrc/LiBiVariation.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/LiBiVariation.o: LiBiNormSrc/LiBiVariation.h
$(BUILD)/LiBiNormSrc/LiBiVariation.o: LiBiNormSrc/LiBiNormCore.h
$(BUILD)/LiBiNormSrc/LiBiVariation.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/LiBiVariation.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/LiBiVariation.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/LiBiVariation.o: LiBiNormSrc/GeneCountData.h
$(BUILD)/LiBiNormSrc/LiBiVariation.o: mcmcLib/mcmc.h mcmcLib/params.h
$(BUILD)/LiBiNormSrc/LiBiVariation.o: ../bioinformaticsLib/nelderMeadOptimiser.h
$(BUILD)/LiBiNormSrc/LiBiVariation.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/LiBiNormSrc/LiBiVariation.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/LiBiVariation.o: LiBiNormSrc/Options.h
$(BUILD)/LiBiNormSrc/LiBiVariation.o: LiBiNormSrc/ModelData.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: LiBiNormSrc/MakeFastq.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/src/api/BamReader.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/src/api/api_global.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/src/api/BamAux.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/MakeFastq.o: ../bioinformaticsLib/fastaFile.h
$(BUILD)/LiBiNormSrc/ModelData.o: ../bioinformaticsLib/rand.h
$(BUILD)/LiBiNormSrc/ModelData.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/LiBiNormSrc/ModelData.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/ModelData.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/ModelData.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/ModelData.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/ModelData.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/ModelData.o: LiBiNormSrc/ModelData.h
$(BUILD)/LiBiNormSrc/ModelData.o: LiBiNormSrc/Options.h mcmcLib/mcmc.h
$(BUILD)/LiBiNormSrc/ModelData.o: mcmcLib/params.h
$(BUILD)/LiBiNormSrc/ModelData.o: ../bioinformaticsLib/nelderMeadOptimiser.h
$(BUILD)/LiBiNormSrc/refSeqs.o: LiBiNormSrc/refSeqs.h hisat2Lib/reference.h
$(BUILD)/LiBiNormSrc/refSeqs.o: hisat2Lib/ref_read.h hisat2Lib/alphabet.h
$(BUILD)/LiBiNormSrc/refSeqs.o: hisat2Lib/assert_helpers.h
$(BUILD)/LiBiNormSrc/refSeqs.o: hisat2Lib/word_io.h hisat2Lib/endian_swap.h
$(BUILD)/LiBiNormSrc/refSeqs.o: hisat2Lib/hisat2Lib.h
$(BUILD)/LiBiNormSrc/refSeqs.o: LiBiNormSrc/FeatureFileEx.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/LiBiNormSrc/refSeqs.o: LiBiNormSrc/Regions.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bamtools/src/api/BamReader.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bamtools/src/api/api_global.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bamtools/src/api/BamAux.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/LiBiNormSrc/refSeqs.o: LiBiNormSrc/Options.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bioinformaticsLib/featureFile.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/LiBiNormSrc/refSeqs.o: LiBiNormSrc/GeneCountData.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bioinformaticsLib/libParser.h
$(BUILD)/LiBiNormSrc/refSeqs.o: mcmcLib/mcmc.h mcmcLib/params.h
$(BUILD)/LiBiNormSrc/refSeqs.o: ../bioinformaticsLib/nelderMeadOptimiser.h
$(BUILD)/LiBiNormSrc/Regions.o: LiBiNormSrc/Regions.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bioinformaticsLib/printEx.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/src/api/BamReader.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/src/api/api_global.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/src/api/BamAux.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/LiBiNormSrc/Regions.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/LiBiNormSrc/Regions.o: LiBiNormSrc/Options.h

$(BUILD)/bioinformaticsLib/codFile.o: ../bioinformaticsLib/codFile.h
$(BUILD)/bioinformaticsLib/codFile.o: ../bioinformaticsLib/genomicPosition.h
$(BUILD)/bioinformaticsLib/codFile.o: ../bioinformaticsLib/libParser.h
$(BUILD)/bioinformaticsLib/codFile.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/dataVec.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/bioinformaticsLib/dataVec.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/dataVec.o: ../bioinformaticsLib/printEx.h
$(BUILD)/bioinformaticsLib/dataVec.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/bioinformaticsLib/dataVec.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/bioinformaticsLib/dataVec.o: ../bioinformaticsLib/libParser.h
$(BUILD)/bioinformaticsLib/dataVec.o: ../bioinformaticsLib/rand.h
$(BUILD)/bioinformaticsLib/fastaFile.o: ../bioinformaticsLib/printEx.h
$(BUILD)/bioinformaticsLib/fastaFile.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/bioinformaticsLib/fastaFile.o: ../bioinformaticsLib/libParser.h
$(BUILD)/bioinformaticsLib/fastaFile.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/fastaFile.o: ../bioinformaticsLib/fastaFile.h
$(BUILD)/bioinformaticsLib/featureFile.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/bioinformaticsLib/featureFile.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/bioinformaticsLib/featureFile.o: ../bioinformaticsLib/featureFile.h
$(BUILD)/bioinformaticsLib/featureFile.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/bioinformaticsLib/featureFile.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/bioinformaticsLib/featureFile.o: ../bioinformaticsLib/printEx.h
$(BUILD)/bioinformaticsLib/featureFile.o: ../bioinformaticsLib/libParser.h
$(BUILD)/bioinformaticsLib/featureFile.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/genbankFile.h
$(BUILD)/bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/containerEx.h
$(BUILD)/bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/printEx.h
$(BUILD)/bioinformaticsLib/genbankFile.o: ../bioinformaticsLib/libParser.h
$(BUILD)/bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/genomicPosition.h
$(BUILD)/bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/libParser.h
$(BUILD)/bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/genomicPosition.o: ../bioinformaticsLib/printEx.h
$(BUILD)/bioinformaticsLib/libCommon.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/libCommon.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/bioinformaticsLib/libCommon.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/bioinformaticsLib/libParser.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/libParser.o: ../bioinformaticsLib/libParser.h
$(BUILD)/bioinformaticsLib/nelderMeadOptimiser.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/nelderMeadOptimiser.o: ../bioinformaticsLib/nelderMeadOptimiser.h
$(BUILD)/bioinformaticsLib/nelderMeadOptimiser.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/bioinformaticsLib/nelderMeadOptimiser.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/bioinformaticsLib/nelderMeadOptimiser.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/bioinformaticsLib/nelderMeadOptimiser.o: ../bioinformaticsLib/printEx.h
$(BUILD)/bioinformaticsLib/printEx.o: ../bioinformaticsLib/printEx.h
$(BUILD)/bioinformaticsLib/printEx.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/bioinformaticsLib/printEx.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/bioinformaticsLib/rand.o: ../bioinformaticsLib/rand.h
$(BUILD)/bioinformaticsLib/rand.o: ../bioinformaticsLib/dataVec.h
$(BUILD)/bioinformaticsLib/rand.o: ../bioinformaticsLib/libCommon.h
$(BUILD)/bioinformaticsLib/rand.o: ../bioinformaticsLib/printEx.h
$(BUILD)/bioinformaticsLib/rand.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/bioinformaticsLib/smithWaterman.o: ../bioinformaticsLib/smithWaterman.h
$(BUILD)/bioinformaticsLib/stringEx.o: ../bioinformaticsLib/stringEx.h
$(BUILD)/bioinformaticsLib/stringEx.o: ../bioinformaticsLib/inQuotes.h
$(BUILD)/bamtools/src/api/BamAlignment.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/bamtools/src/api/BamAlignment.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/BamAlignment.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/BamAlignment.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/BamAlignment.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/BamMultiReader.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/BamReader.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/internal/bam/BamMultiReader_p.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/internal/bam/BamMultiMerger_p.h
$(BUILD)/bamtools/src/api/BamMultiReader.o: ../bamtools/src/api/algorithms/Sort.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/BamReader.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/internal/bam/BamReader_p.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/internal/bam/BamHeader_p.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/internal/bam/BamRandomAccessController_p.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/src/api/BamReader.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/api/BamWriter.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/api/internal/bam/BamWriter_p.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/src/api/BamWriter.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/bam/BamHeader_p.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/internal/bam/BamHeader_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/bam/BamHeader_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/bam/BamHeader_p.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/bamtools/src/api/internal/bam/BamHeader_p.o: ../bamtools/src/api/internal/bam/BamHeader_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamHeader_p.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/bamtools/src/api/internal/bam/BamHeader_p.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/bamtools/src/api/internal/bam/BamHeader_p.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/api/internal/bam/BamHeader_p.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/src/api/internal/bam/BamHeader_p.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/api/internal/bam/BamHeader_p.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/bamtools/src/api/internal/bam/BamHeader_p.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/api/internal/bam/BamHeader_p.o: ../bamtools/src/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamHeader_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/bam/BamHeader_p.o: ../bamtools/src/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/BamMultiReader.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/BamReader.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/SamConstants.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/algorithms/Sort.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/internal/bam/BamMultiReader_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamMultiReader_p.o: ../bamtools/src/api/internal/bam/BamMultiMerger_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/internal/bam/BamRandomAccessController_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/internal/bam/BamReader_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/BamReader.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/internal/bam/BamHeader_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/internal/index/BamIndexFactory_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamRandomAccessController_p.o: ../bamtools/src/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/BamReader.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/internal/bam/BamHeader_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/internal/bam/BamRandomAccessController_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/internal/bam/BamReader_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/internal/index/BamStandardIndex_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/internal/index/BamToolsIndex_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/internal/io/BamDeviceFactory_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamReader_p.o: ../bamtools/src/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamWriter_p.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/bamtools/src/api/internal/bam/BamWriter_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/bam/BamWriter_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/bam/BamWriter_p.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/internal/bam/BamWriter_p.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/bamtools/src/api/internal/bam/BamWriter_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/bam/BamWriter_p.o: ../bamtools/src/api/internal/bam/BamWriter_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamWriter_p.o: ../bamtools/src/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/src/api/internal/bam/BamWriter_p.o: ../bamtools/src/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/src/api/internal/index/BamIndexFactory_p.o: ../bamtools/src/api/internal/index/BamIndexFactory_p.h
$(BUILD)/bamtools/src/api/internal/index/BamIndexFactory_p.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/bamtools/src/api/internal/index/BamIndexFactory_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/index/BamIndexFactory_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/index/BamIndexFactory_p.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/internal/index/BamIndexFactory_p.o: ../bamtools/src/api/internal/index/BamStandardIndex_p.h
$(BUILD)/bamtools/src/api/internal/index/BamIndexFactory_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/index/BamIndexFactory_p.o: ../bamtools/src/api/internal/index/BamToolsIndex_p.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/internal/bam/BamReader_p.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/BamReader.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/internal/bam/BamHeader_p.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/internal/bam/BamRandomAccessController_p.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/internal/index/BamStandardIndex_p.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/internal/io/BamDeviceFactory_p.h
$(BUILD)/bamtools/src/api/internal/index/BamStandardIndex_p.o: ../bamtools/src/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/internal/bam/BamReader_p.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/BamReader.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/internal/bam/BamHeader_p.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/internal/bam/BamRandomAccessController_p.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/internal/index/BamToolsIndex_p.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/internal/io/BamDeviceFactory_p.h
$(BUILD)/bamtools/src/api/internal/index/BamToolsIndex_p.o: ../bamtools/src/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/src/api/internal/io/BamDeviceFactory_p.o: ../bamtools/src/api/internal/io/BamDeviceFactory_p.h
$(BUILD)/bamtools/src/api/internal/io/BamDeviceFactory_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/io/BamDeviceFactory_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/io/BamDeviceFactory_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/io/BamDeviceFactory_p.o: ../bamtools/src/api/internal/io/BamFile_p.h
$(BUILD)/bamtools/src/api/internal/io/BamDeviceFactory_p.o: ../bamtools/src/api/internal/io/ILocalIODevice_p.h
$(BUILD)/bamtools/src/api/internal/io/BamDeviceFactory_p.o: ../bamtools/src/api/internal/io/BamFtp_p.h
$(BUILD)/bamtools/src/api/internal/io/BamDeviceFactory_p.o: ../bamtools/src/api/internal/io/BamHttp_p.h
$(BUILD)/bamtools/src/api/internal/io/BamDeviceFactory_p.o: ../bamtools/src/api/internal/io/BamPipe_p.h
$(BUILD)/bamtools/src/api/internal/io/BamFile_p.o: ../bamtools/src/api/internal/io/BamFile_p.h
$(BUILD)/bamtools/src/api/internal/io/BamFile_p.o: ../bamtools/src/api/internal/io/ILocalIODevice_p.h
$(BUILD)/bamtools/src/api/internal/io/BamFile_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/io/BamFile_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/io/BamFile_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/io/BamFtp_p.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/internal/io/BamFtp_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/io/BamFtp_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/io/BamFtp_p.o: ../bamtools/src/api/internal/io/BamFtp_p.h
$(BUILD)/bamtools/src/api/internal/io/BamFtp_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/io/BamFtp_p.o: ../bamtools/src/api/internal/io/TcpSocket_p.h
$(BUILD)/bamtools/src/api/internal/io/BamFtp_p.o: ../bamtools/src/api/internal/io/HostInfo_p.h
$(BUILD)/bamtools/src/api/internal/io/BamFtp_p.o: ../bamtools/src/api/internal/io/HostAddress_p.h
$(BUILD)/bamtools/src/api/internal/io/BamFtp_p.o: ../bamtools/src/api/internal/io/RollingBuffer_p.h
$(BUILD)/bamtools/src/api/internal/io/BamFtp_p.o: ../bamtools/src/api/internal/io/ByteArray_p.h
$(BUILD)/bamtools/src/api/internal/io/BamHttp_p.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/internal/io/BamHttp_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/io/BamHttp_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/io/BamHttp_p.o: ../bamtools/src/api/internal/io/BamHttp_p.h
$(BUILD)/bamtools/src/api/internal/io/BamHttp_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/io/BamHttp_p.o: ../bamtools/src/api/internal/io/HttpHeader_p.h
$(BUILD)/bamtools/src/api/internal/io/BamHttp_p.o: ../bamtools/src/api/internal/io/TcpSocket_p.h
$(BUILD)/bamtools/src/api/internal/io/BamHttp_p.o: ../bamtools/src/api/internal/io/HostInfo_p.h
$(BUILD)/bamtools/src/api/internal/io/BamHttp_p.o: ../bamtools/src/api/internal/io/HostAddress_p.h
$(BUILD)/bamtools/src/api/internal/io/BamHttp_p.o: ../bamtools/src/api/internal/io/RollingBuffer_p.h
$(BUILD)/bamtools/src/api/internal/io/BamHttp_p.o: ../bamtools/src/api/internal/io/ByteArray_p.h
$(BUILD)/bamtools/src/api/internal/io/BamPipe_p.o: ../bamtools/src/api/internal/io/BamPipe_p.h
$(BUILD)/bamtools/src/api/internal/io/BamPipe_p.o: ../bamtools/src/api/internal/io/ILocalIODevice_p.h
$(BUILD)/bamtools/src/api/internal/io/BamPipe_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/io/BamPipe_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/io/BamPipe_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/io/BgzfStream_p.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/internal/io/BgzfStream_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/io/BgzfStream_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/io/BgzfStream_p.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/bamtools/src/api/internal/io/BgzfStream_p.o: ../bamtools/src/api/internal/io/BamDeviceFactory_p.h
$(BUILD)/bamtools/src/api/internal/io/BgzfStream_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/io/BgzfStream_p.o: ../bamtools/src/api/internal/io/BgzfStream_p.h
$(BUILD)/bamtools/src/api/internal/io/BgzfStream_p.o: ../bamtools/src/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/src/api/internal/io/ByteArray_p.o: ../bamtools/src/api/internal/io/ByteArray_p.h
$(BUILD)/bamtools/src/api/internal/io/ByteArray_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/io/ByteArray_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/io/HostAddress_p.o: ../bamtools/src/api/internal/io/HostAddress_p.h
$(BUILD)/bamtools/src/api/internal/io/HostAddress_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/io/HostAddress_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/io/HostInfo_p.o: ../bamtools/src/api/internal/io/HostInfo_p.h
$(BUILD)/bamtools/src/api/internal/io/HostInfo_p.o: ../bamtools/src/api/internal/io/HostAddress_p.h
$(BUILD)/bamtools/src/api/internal/io/HostInfo_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/io/HostInfo_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/io/HostInfo_p.o: ../bamtools/src/api/internal/io/NetUnix_p.h
$(BUILD)/bamtools/src/api/internal/io/HttpHeader_p.o: ../bamtools/src/api/internal/io/HttpHeader_p.h
$(BUILD)/bamtools/src/api/internal/io/HttpHeader_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/io/HttpHeader_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/io/ILocalIODevice_p.o: ../bamtools/src/api/internal/io/ILocalIODevice_p.h
$(BUILD)/bamtools/src/api/internal/io/ILocalIODevice_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/io/ILocalIODevice_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/io/ILocalIODevice_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/io/RollingBuffer_p.o: ../bamtools/src/api/internal/io/RollingBuffer_p.h
$(BUILD)/bamtools/src/api/internal/io/RollingBuffer_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/io/RollingBuffer_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/io/RollingBuffer_p.o: ../bamtools/src/api/internal/io/ByteArray_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_p.o: ../bamtools/src/api/internal/io/HostInfo_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_p.o: ../bamtools/src/api/internal/io/HostAddress_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_p.o: ../bamtools/src/api/internal/io/TcpSocketEngine_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_p.o: ../bamtools/src/api/internal/io/TcpSocket_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_p.o: ../bamtools/src/api/internal/io/RollingBuffer_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_p.o: ../bamtools/src/api/internal/io/ByteArray_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/src/api/internal/io/TcpSocketEngine_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/src/api/internal/io/HostAddress_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/src/api/internal/io/TcpSocket_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/src/api/internal/io/HostInfo_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/src/api/internal/io/RollingBuffer_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/src/api/internal/io/ByteArray_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocketEngine_unix_p.o: ../bamtools/src/api/internal/io/NetUnix_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocket_p.o: ../bamtools/src/api/internal/io/ByteArray_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocket_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocket_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocket_p.o: ../bamtools/src/api/internal/io/TcpSocket_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocket_p.o: ../bamtools/src/api/IBamIODevice.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocket_p.o: ../bamtools/src/api/internal/io/HostInfo_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocket_p.o: ../bamtools/src/api/internal/io/HostAddress_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocket_p.o: ../bamtools/src/api/internal/io/RollingBuffer_p.h
$(BUILD)/bamtools/src/api/internal/io/TcpSocket_p.o: ../bamtools/src/api/internal/io/TcpSocketEngine_p.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatParser_p.o: ../bamtools/src/api/SamConstants.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatParser_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatParser_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatParser_p.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatParser_p.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatParser_p.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatParser_p.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatParser_p.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatParser_p.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatParser_p.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatParser_p.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatParser_p.o: ../bamtools/src/api/internal/sam/SamFormatParser_p.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatParser_p.o: ../bamtools/src/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/src/api/SamConstants.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/api/internal/sam/SamFormatPrinter_p.o: ../bamtools/src/api/internal/sam/SamFormatPrinter_p.h
$(BUILD)/bamtools/src/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/src/api/SamConstants.h
$(BUILD)/bamtools/src/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/bamtools/src/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/bamtools/src/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/src/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/bamtools/src/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/src/api/internal/sam/SamHeaderValidator_p.h
$(BUILD)/bamtools/src/api/internal/sam/SamHeaderValidator_p.o: ../bamtools/src/api/internal/sam/SamHeaderVersion_p.h
$(BUILD)/bamtools/src/api/internal/utils/BamException_p.o: ../bamtools/src/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/src/api/SamHeader.o: ../bamtools/src/api/SamConstants.h
$(BUILD)/bamtools/src/api/SamHeader.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/SamHeader.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/SamHeader.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/bamtools/src/api/SamHeader.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/SamHeader.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/bamtools/src/api/SamHeader.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/api/SamHeader.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/src/api/SamHeader.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/api/SamHeader.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/bamtools/src/api/SamHeader.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/api/SamHeader.o: ../bamtools/src/api/internal/utils/BamException_p.h
$(BUILD)/bamtools/src/api/SamHeader.o: ../bamtools/src/api/internal/sam/SamFormatParser_p.h
$(BUILD)/bamtools/src/api/SamHeader.o: ../bamtools/src/api/internal/sam/SamFormatPrinter_p.h
$(BUILD)/bamtools/src/api/SamHeader.o: ../bamtools/src/api/internal/sam/SamHeaderValidator_p.h
$(BUILD)/bamtools/src/api/SamProgram.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/api/SamProgram.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/SamProgram.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/SamProgram.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/SamProgramChain.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/bamtools/src/api/SamProgramChain.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/SamProgramChain.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/SamProgramChain.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/api/SamProgramChain.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/SamReadGroup.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/api/SamReadGroup.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/SamReadGroup.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/SamReadGroup.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/SamReadGroupDictionary.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/src/api/SamReadGroupDictionary.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/SamReadGroupDictionary.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/SamReadGroupDictionary.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/api/SamReadGroupDictionary.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/SamSequence.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/api/SamSequence.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/SamSequence.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/SamSequence.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/api/SamSequenceDictionary.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/bamtools/src/api/SamSequenceDictionary.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/api/SamSequenceDictionary.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/api/SamSequenceDictionary.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/api/SamSequenceDictionary.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/toolkit/bamtools_sort.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/toolkit/bamtools_tool.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/SamConstants.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/api_global.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/shared/bamtools_global.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/BamMultiReader.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/BamReader.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/BamAlignment.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/BamAux.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/BamConstants.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/BamIndex.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/SamHeader.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/SamProgramChain.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/SamProgram.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/SamReadGroupDictionary.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/SamReadGroup.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/SamSequenceDictionary.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/SamSequence.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/BamWriter.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/api/algorithms/Sort.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/utils/bamtools_options.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/utils/bamtools_variant.h
$(BUILD)/bamtools/src/toolkit/bamtools_sort.o: ../bamtools/src/utils/utils_global.h
$(BUILD)/bamtools/src/utils/bamtools_options.o: ../bamtools/src/utils/bamtools_options.h
$(BUILD)/bamtools/src/utils/bamtools_options.o: ../bamtools/src/utils/bamtools_variant.h
$(BUILD)/bamtools/src/utils/bamtools_options.o: ../bamtools/src/utils/utils_global.h
$(BUILD)/bamtools/src/utils/bamtools_options.o: ../bamtools/src/shared/bamtools_global.h
