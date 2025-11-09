add_library(NPMugast SHARED MugastDetector.cxx MugastData.cxx MugastPhysics.cxx MugastSpectra.cxx)
target_link_libraries(NPMugast PUBLIC PUBLIC ROOT::RIO ROOT::Tree ROOT::Physics ROOT::Imt ROOT::Hist ROOT::Gpad NPRoot nptool )

if(MFM_FOUND)
  target_compile_definitions(NPMugast PUBLIC MFM_FOUND)
  target_link_libraries(NPMugast PUBLIC ${MFMCOMPLIB})
endif()

if(Geant4_FOUND)
  add_library(NPMugastGeant4 SHARED MugastGeant4.h MugastGeant4.cxx)
  target_link_libraries(NPMugastGeant4 PUBLIC NPMugast NPGeant4 nptool ${Geant4_LIBRARIES})
endif()

install(DIRECTORY ./ DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY} FILES_MATCHING PATTERN "*.h" )
ROOT_GENERATE_DICTIONARY(MugastDict MugastData.h MugastPhysics.h MugastSpectra.h LINKDEF Mugastlinkdef.hh MODULE NPMugast)
