cmake_minimum_required(VERSION 3.16.0 FATAL_ERROR)

set(CMAKE_SYSTEM_VERSION 10.0 CACHE STRING "" FORCE)

project(FEM++ CXX)

################################################################################
# Set target arch type if empty. Visual studio solution generator provides it.
################################################################################
if(NOT CMAKE_VS_PLATFORM_NAME)
    set(CMAKE_VS_PLATFORM_NAME "x64")
endif()
message("${CMAKE_VS_PLATFORM_NAME} architecture in use")

if(NOT ("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64"
     OR "${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86"))
    message(FATAL_ERROR "${CMAKE_VS_PLATFORM_NAME} arch is not supported!")
endif()

################################################################################
# Global configuration types
################################################################################
set(CMAKE_CONFIGURATION_TYPES
    "Debug"
    "Release"
    CACHE STRING "" FORCE
)

################################################################################
# Global compiler options
################################################################################
if(MSVC)
    # remove default flags provided with CMake for MSVC
    set(CMAKE_CXX_FLAGS "")
    set(CMAKE_CXX_FLAGS_DEBUG "")
    set(CMAKE_CXX_FLAGS_RELEASE "")
endif()

################################################################################
# Global linker options
################################################################################
if(MSVC)
    # remove default flags provided with CMake for MSVC
    set(CMAKE_EXE_LINKER_FLAGS "")
    set(CMAKE_MODULE_LINKER_FLAGS "")
    set(CMAKE_SHARED_LINKER_FLAGS "")
    set(CMAKE_STATIC_LINKER_FLAGS "")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_MODULE_LINKER_FLAGS_DEBUG "${CMAKE_MODULE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS}")
    set(CMAKE_STATIC_LINKER_FLAGS_DEBUG "${CMAKE_STATIC_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_MODULE_LINKER_FLAGS_RELEASE "${CMAKE_MODULE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS}")
    set(CMAKE_STATIC_LINKER_FLAGS_RELEASE "${CMAKE_STATIC_LINKER_FLAGS}")
endif()

################################################################################
# Nuget packages function stub.
################################################################################
function(use_package TARGET PACKAGE VERSION)
    message(WARNING "No implementation of use_package. Create yours. "
                    "Package \"${PACKAGE}\" with version \"${VERSION}\" "
                    "for target \"${TARGET}\" is ignored!")
endfunction()

################################################################################
# Common utils
################################################################################
include(CMake/Utils.cmake)

################################################################################
# Additional Global Settings(add specific info there)
################################################################################
include(CMake/GlobalSettingsInclude.cmake OPTIONAL)

################################################################################
# Use solution folders feature
################################################################################
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

################################################################################
# Sub-projects
################################################################################



























set(PROJECT_NAME FEM++)

################################################################################
# Source groups
################################################################################
set(no_group_source_files
    "../alglib-cpp/src/alglibinternal.cpp"
    "../alglib-cpp/src/alglibmisc.cpp"
    "../alglib-cpp/src/ap.cpp"
    "../alglib-cpp/src/dataanalysis.cpp"
    "../alglib-cpp/src/diffequations.cpp"
    "../alglib-cpp/src/fasttransforms.cpp"
    "../alglib-cpp/src/integration.cpp"
    "../alglib-cpp/src/interpolation.cpp"
    "../alglib-cpp/src/kernels_avx2.cpp"
    "../alglib-cpp/src/kernels_fma.cpp"
    "../alglib-cpp/src/kernels_sse2.cpp"
    "../alglib-cpp/src/linalg.cpp"
    "../alglib-cpp/src/optimization.cpp"
    "../alglib-cpp/src/solvers.cpp"
    "../alglib-cpp/src/specialfunctions.cpp"
    "../alglib-cpp/src/statistics.cpp"
    "../fileIO/fileIO2.cpp"
    "../fileIO/fileIO2.h"
    "../NLopt.v2.6.1.x86_64.tar/include/nlopt.h"
    "../NLopt.v2.6.1.x86_64.tar/include/nlopt.hpp"
    "../postprocessing/postprocessing.cpp"
    "../run_modes.cpp"
    "../utility/adp_methods.cpp"
    "../utility/refinement_tests.cpp"
)
source_group("" FILES ${no_group_source_files})

set(Header_Files
    "../Bayes/bayes.h"
    "../Kriging/kriging.h"
    "../run_modes.h"
)
source_group("Header Files" FILES ${Header_Files})

set(Header_Files__FILEIO
    "../utility/fileIO_functions.h"
    "../utility/fileIOCon.h"
)
source_group("Header Files\\FILEIO" FILES ${Header_Files__FILEIO})

set(Header_Files__HOPS
    "../HOPS/hops.h"
)
source_group("Header Files\\HOPS" FILES ${Header_Files__HOPS})

set(Header_Files__Postprocessing
    "../postprocessing/postprocessing.h"
)
source_group("Header Files\\Postprocessing" FILES ${Header_Files__Postprocessing})

set(Header_Files__QOI
    "../functions/qoi_error.h"
    "../functions/sensitivity.h"
    "../structure/adjacent.h"
    "../structure/element_error.h"
    "../structure/error_struct.h"
    "../utility/additional_basis.h"
    "../utility/qoi_helper.h"
)
source_group("Header Files\\QOI" FILES ${Header_Files__QOI})

set(Header_Files__Refinement
    "../utility/adp_methods.h"
    "../utility/h_refine_assist.h"
    "../utility/refinement.h"
    "../utility/refinement_tests.h"
)
source_group("Header Files\\Refinement" FILES ${Header_Files__Refinement})

set(Header_Files__Structure
    "../control/StructureControl.h"
    "../functions/basis.h"
    "../functions/Coordinates.h"
    "../functions/find_at_quad.h"
    "../functions/products.h"
    "../functions/unitVectorsM.h"
    "../integrals/BasisEval.h"
    "../structure/Domain.h"
    "../structure/Element.h"
    "../structure/Facet.h"
    "../structure/Material.h"
    "../structure/Point.h"
    "../structure/Scatter.h"
    "../structure/Tables.h"
    "../utility/constants.h"
    "../utility/functions.h"
    "../utility/memory_functions.h"
    "../utility/Plotter.h"
)
source_group("Header Files\\Structure" FILES ${Header_Files__Structure})

set(Header_Files__integration
    "../functions/elunknowns.h"
    "../functions/integralFunctions.h"
    "../functions/matrixFilling.h"
    "../integrals/Integral_c.h"
    "../integrals/Integral_d.h"
    "../integrals/Integral_g.h"
    "../structure/FrequencySweep.h"
    "../utility/sparseSolver.h"
)
source_group("Header Files\\integration" FILES ${Header_Files__integration})

set(Header_Files__matrices
    "../matricies/matrix2d.h"
    "../matricies/matrix3d.h"
    "../matricies/matrix4d.h"
)
source_group("Header Files\\matrices" FILES ${Header_Files__matrices})

set(Source_Files
    "../Bayes/bayes.cpp"
    "../control/StructureControl.cpp"
    "../functions/Coordinates.cpp"
    "../functions/integralFunctions.cpp"
    "../functions/matrixFilling.cpp"
    "../functions/unitVectorsM.cpp"
    "../integrals/BasisEval.cpp"
    "../integrals/Integral_c.cpp"
    "../integrals/Integral_d.cpp"
    "../integrals/Integral_g.cpp"
    "../Kriging/kriging.cpp"
    "../Main.cpp"
    "../structure/Domain.cpp"
    "../structure/element_error.cpp"
    "../structure/Facet.cpp"
    "../utility/constants.cpp"
    "../utility/functions.cpp"
    "../utility/h_refine_assist.cpp"
    "../utility/Plotter.cpp"
    "../utility/qoi_helper.cpp"
    "../utility/refinement.cpp"
    "../utility/sparseSolver.cpp"
)
source_group("Source Files" FILES ${Source_Files})

set(Source_Files__HOPS
    "../HOPS/hops.cpp"
)
source_group("Source Files\\HOPS" FILES ${Source_Files__HOPS})

set(ALL_FILES
    ${no_group_source_files}
    ${Header_Files}
    ${Header_Files__FILEIO}
    ${Header_Files__HOPS}
    ${Header_Files__Postprocessing}
    ${Header_Files__QOI}
    ${Header_Files__Refinement}
    ${Header_Files__Structure}
    ${Header_Files__integration}
    ${Header_Files__matrices}
    ${Source_Files}
    ${Source_Files__HOPS}
)

################################################################################
# Target
################################################################################
add_executable(${PROJECT_NAME} ${ALL_FILES})

use_props(${PROJECT_NAME} "${CMAKE_CONFIGURATION_TYPES}" "${DEFAULT_CXX_PROPS}")
set_target_properties(${PROJECT_NAME} PROPERTIES
    VS_GLOBAL_KEYWORD "Win32Proj"
)
################################################################################
# MSVC runtime library
################################################################################
get_property(MSVC_RUNTIME_LIBRARY_DEFAULT TARGET ${PROJECT_NAME} PROPERTY MSVC_RUNTIME_LIBRARY)
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    string(CONCAT "MSVC_RUNTIME_LIBRARY_STR"
        $<$<CONFIG:Release>:
            MultiThreaded
        >
        $<$<NOT:$<OR:$<CONFIG:Release>>>:${MSVC_RUNTIME_LIBRARY_DEFAULT}>
    )
elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
    string(CONCAT "MSVC_RUNTIME_LIBRARY_STR"
        $<$<CONFIG:Debug>:
            MultiThreadedDebugDLL
        >
        $<$<CONFIG:Release>:
            MultiThreadedDLL
        >
        $<$<NOT:$<OR:$<CONFIG:Debug>,$<CONFIG:Release>>>:${MSVC_RUNTIME_LIBRARY_DEFAULT}>
    )
endif()
set_target_properties(${PROJECT_NAME} PROPERTIES MSVC_RUNTIME_LIBRARY ${MSVC_RUNTIME_LIBRARY_STR})

################################################################################
# Include directories
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_include_directories(${PROJECT_NAME} PUBLIC
        "$<$<CONFIG:Release>:"
            "${CMAKE_CURRENT_SOURCE_DIR}/../nlopt/src/api;"
            "${CMAKE_CURRENT_SOURCE_DIR}/../alglib-cpp;"
            "${CMAKE_CURRENT_SOURCE_DIR}/../dlib"
        ">"
        "${CMAKE_CURRENT_SOURCE_DIR}/../eigen-eigen-5a0156e40feb"
    )
endif()

################################################################################
# Compile definitions
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        "$<$<CONFIG:Debug>:"
            "_DEBUG"
        ">"
        "$<$<CONFIG:Release>:"
            "NDEBUG"
        ">"
        "WIN32;"
        "_CONSOLE"
    )
endif()

################################################################################
# Compile and link options
################################################################################
if(MSVC)
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        target_compile_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Release>:
                /MP;
                /Ob2;
                /Oi;
                /GT;
                /Ot;
                /wd1224
            >
            ${DEFAULT_CXX_DEBUG_INFORMATION_FORMAT};
            ${DEFAULT_CXX_EXCEPTION_HANDLING}
        )
    elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
        target_compile_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Debug>:
                /Od
            >
            $<$<CONFIG:Release>:
                /MP;
                /Ob2;
                /Oi;
                /GT;
                /Ot;
                /wd1224
            >
            /W3;
            /Zi;
            ${DEFAULT_CXX_EXCEPTION_HANDLING}
        )
    endif()
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        target_link_options(${PROJECT_NAME} PRIVATE
            /DEBUG:FULL;
            /SUBSYSTEM:CONSOLE
        )
    elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
        target_link_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Release>:
                /OPT:REF;
                /OPT:ICF
            >
            /DEBUG;
            /MACHINE:X86;
            /SUBSYSTEM:CONSOLE;
            /INCREMENTAL
        )
    endif()
endif()

################################################################################
# Dependencies
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    set(ADDITIONAL_LIBRARY_DEPENDENCIES
        "$<$<CONFIG:Release>:"
            "nlopt"
        ">"
    )
endif()
target_link_libraries(${PROJECT_NAME} PRIVATE "${ADDITIONAL_LIBRARY_DEPENDENCIES}")

if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_link_directories(${PROJECT_NAME} PRIVATE
        "$<$<CONFIG:Release>:"
            "F:/Github/3DFEM-UQ-AR/cpp/nlopt/build/Release"
        ">"
    )
endif()

