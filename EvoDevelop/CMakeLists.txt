
set(CMAKE_CXX_STANDARD 20)

# OpenMP
if (APPLE)
    set(brew_libomp /opt/homebrew/Cellar/libomp/17.0.2)
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${brew_libomp}/include")
    set(OpenMP_LIBS ${brew_libomp}/lib/libomp.dylib)
    set(OpenMP_INCLUDE_DIRS ${brew_libomp}/include)
endif()

# OpenMesh
if (APPLE)
    #set(OpenMesh_ROOT ~/Documents/Code/OpenMesh-9.0.0)
    #set(OpenMesh_INCLUDE_DIRS ${OpenMesh_ROOT}/include)
    #set(OpenMesh_LIBS ${OpenMesh_ROOT}/lib/libOpenMeshCore.dylib)
    set(OpenMesh_LIBS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Library)
    set(OpenMesh_LIBS ${CMAKE_CURRENT_SOURCE_DIR}/Library/libOpenMeshCore.9.0.dylib)
endif()

if (WIN32)
    set(OpenMesh_LIBS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Library)
endif()


file(GLOB EVO_DEVELOP_SRCS
    "*.h"
    "*.cpp"
)

file(GLOB_RECURSE EVO_DEVELOP_COMPUTE_SRCS 
    "Compute/*.cpp"
    "Compute/*.h"
)

add_executable(EvoDevelop 
    ${EVO_DEVELOP_SRCS} 
    ${EVO_DEVELOP_COMPUTE_SRCS}
)

if (MSVC)
    target_compile_definitions(EvoDevelop PRIVATE -D_USE_MATH_DEFINES=1)
endif()

target_include_directories(EvoDevelop PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(EvoDevelop PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Compute)
target_include_directories(EvoDevelop PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Include/)
target_include_directories(EvoDevelop PUBLIC ${OpenMP_INCLUDE_DIRS})
target_include_directories(EvoDevelop PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Include/spdlog/include)

#target_include_directories(EvoDevelop PUBLIC ${OpenMesh_INCLUDE_DIRS})

target_link_libraries(EvoDevelop PUBLIC ${OpenMP_LIBS})

target_link_directories(EvoDevelop PUBLIC ${OpenMesh_LIBS_DIR})
target_link_libraries(EvoDevelop PUBLIC ${OpenMesh_LIBS})
