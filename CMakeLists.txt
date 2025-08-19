cmake_minimum_required(VERSION 3.21)
project(structure-formation)
set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "-ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE SOURCES_APP "${CMAKE_SOURCE_DIR}/*.cpp")

include(FetchContent)

# ImGui Bundle library
FetchContent_Declare(imgui_bundle GIT_REPOSITORY https://github.com/pthom/imgui_bundle.git GIT_TAG main)
FetchContent_MakeAvailable(imgui_bundle)

# muFFT library
add_subdirectory(muFFT)
list(APPEND LIBRARIES muFFT)
list(APPEND INCLUDES muFFT)

# glm library
FetchContent_Declare(glm GIT_REPOSITORY https://github.com/g-truc/glm.git GIT_TAG 1.0.1)
FetchContent_MakeAvailable(glm)
list(APPEND LIBRARIES glm::glm)
list(APPEND INCLUDES ${glm_SOURCE_DIR})

# Perlin noise library
FetchContent_Declare(perlin GIT_REPOSITORY https://github.com/Reputeless/PerlinNoise.git GIT_TAG v3.0.0)
FetchContent_MakeAvailable(perlin)
list(APPEND INCLUDES ${perlin_SOURCE_DIR})

set(CMAKE_SOURCE_DIR "src")
file(GLOB_RECURSE SOURCES "${CMAKE_SOURCE_DIR}/*.cpp")

imgui_bundle_add_app(${CMAKE_PROJECT_NAME} ${SOURCES})
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${LIBRARIES})
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${INCLUDES})
