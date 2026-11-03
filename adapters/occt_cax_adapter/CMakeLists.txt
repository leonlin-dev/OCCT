# ===========================================================================
# adapters/occt_cax_adapter/CMakeLists.txt
#
# Standalone and mountable CMake definition for occt_cax_adapter.
# Target: adapters::occt_cax (INTERFACE)
# Witness: smoke_occt_cax_adapter (Gate G0-C)
# ===========================================================================
cmake_minimum_required(VERSION 3.20)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(occt_cax_adapter LANGUAGES CXX)
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(OCCT_ADAPTER_STANDALONE ON)
else()
    set(OCCT_ADAPTER_STANDALONE OFF)
endif()

# Fallback path resolution for OCCT install
if(NOT DEFINED OCCT_INSTALL_DIR)
    if(DEFINED ENV{OCCT_INSTALL_DIR})
        set(OCCT_INSTALL_DIR "$ENV{OCCT_INSTALL_DIR}")
    else()
        set(OCCT_INSTALL_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../../codebase-baselines/aistudio/mathstudio/3rd-party/OCCT-7.8.1/install"
            CACHE PATH "Path to prebuilt OCCT 7.8.1 installation")
    endif()
endif()

# Fallback path resolution for caxlab
if(NOT DEFINED CAXLAB_ROOT)
    if(DEFINED ENV{CAXLAB_ROOT})
        set(CAXLAB_ROOT "$ENV{CAXLAB_ROOT}")
    else()
        set(CAXLAB_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../../../codebase-baselines/caxlab"
            CACHE PATH "Path to caxlab repository")
    endif()
endif()

set(CAXLAB_LOFT_INCLUDE "${CAXLAB_ROOT}/impls-loft/cpp14/include")
set(CAXLAB_LOFT_SUPPORT "${CAXLAB_ROOT}/impls-loft/cpp14/support")
set(EIGEN_INCLUDE_DIR "${CAXLAB_ROOT}/third_party/eigen-3.3.7")

# OCCT include and lib paths
set(OCCT_INCLUDE_DIR "${OCCT_INSTALL_DIR}/inc")
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(OCCT_LIB_DIR "${OCCT_INSTALL_DIR}/win64/vc14/lib")
    set(OCCT_BIN_DIR "${OCCT_INSTALL_DIR}/win64/vc14/bin")
else()
    set(OCCT_LIB_DIR "${OCCT_INSTALL_DIR}/win32/vc14/lib")
    set(OCCT_BIN_DIR "${OCCT_INSTALL_DIR}/win32/vc14/bin")
endif()

if(NOT TARGET adapters::occt_cax)
    add_library(adapter_occt_cax INTERFACE)
    add_library(adapters::occt_cax ALIAS adapter_occt_cax)

    target_include_directories(adapter_occt_cax INTERFACE
        "${CMAKE_CURRENT_LIST_DIR}"
        "${OCCT_INCLUDE_DIR}"
        "${CAXLAB_LOFT_INCLUDE}"
        "${CAXLAB_LOFT_SUPPORT}"
        "${EIGEN_INCLUDE_DIR}"
    )

    target_link_directories(adapter_occt_cax INTERFACE
        "${OCCT_LIB_DIR}"
    )

    target_link_libraries(adapter_occt_cax INTERFACE
        TKernel
        TKMath
        TKG2d
        TKG3d
        TKGeomBase
    )

    target_compile_definitions(adapter_occt_cax INTERFACE
        NOMINMAX
        _CRT_SECURE_NO_WARNINGS
    )
endif()

# Build smoke witness if standalone or requested
if(OCCT_ADAPTER_STANDALONE OR BUILD_TESTING)
    enable_testing()
    add_executable(smoke_occt_cax_adapter
        "${CMAKE_CURRENT_LIST_DIR}/smoke_occt_cax_adapter.cpp")
    target_link_libraries(smoke_occt_cax_adapter PRIVATE adapters::occt_cax)

    if(MSVC)
        target_compile_options(smoke_occt_cax_adapter PRIVATE /W4 /WX /permissive- /EHsc /bigobj)
    else()
        target_compile_options(smoke_occt_cax_adapter PRIVATE -Wall -Wextra -Werror -pedantic)
    endif()

    add_custom_command(TARGET smoke_occt_cax_adapter POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${OCCT_BIN_DIR}/TKernel.dll"
            "${OCCT_BIN_DIR}/TKMath.dll"
            "${OCCT_BIN_DIR}/TKG2d.dll"
            "${OCCT_BIN_DIR}/TKG3d.dll"
            "${OCCT_BIN_DIR}/TKGeomBase.dll"
            $<TARGET_FILE_DIR:smoke_occt_cax_adapter>
    )

    add_test(NAME smoke_occt_cax_adapter COMMAND smoke_occt_cax_adapter)
endif()
