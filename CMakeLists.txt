cmake_minimum_required (VERSION 3.18)
project (Ploto)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/3pty/glfw)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/3pty/imgui)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/3pty/gl3w)

set_target_properties(glfw    PROPERTIES FOLDER 3pty)
set_target_properties(imgui   PROPERTIES FOLDER 3pty)
set_target_properties(gl3w    PROPERTIES FOLDER 3pty)

add_executable(Ploto ${CMAKE_CURRENT_SOURCE_DIR}/ploto.cpp)
target_link_libraries(Ploto PRIVATE imgui glfw gl3w)
