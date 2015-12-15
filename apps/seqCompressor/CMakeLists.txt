# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/seqComrpessor/CMakeLists.txt
#
# CMakeLists.txt file for seqCompressor.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_apps_seqCompressor)
message (STATUS "Configuring apps/seqCompressor")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
#set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS_MAIN})
include_directories (SYSTEM ${SEQAN_INCLUDE_DIRS_DEPS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (seqCompressor seqCompressor.cpp 
                      seqCompressor_base.h
                      generateEvents.h
                      processEvents.h
)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (seqCompressor ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (seqCompressor)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install seqCompressor in ${PREFIX}/bin directory
install (TARGETS seqCompressor
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/seqCompressor for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (seqCompressor)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:seqComrpessor")
  set (CPACK_PACKAGE_NAME "seqCompressor")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "seqCompressor")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Rene Rahn <rene.rahn@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Rene Rahn <rene.rahn@fu-berlin.de>")

  seqan_configure_cpack_app (seqCompressor "seqCompressor")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:seqCompressor")

