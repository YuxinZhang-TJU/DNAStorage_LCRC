cmake_minimum_required(VERSION 3.10)

project(CompositeCode)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(${PROJECT_SOURCE_DIR}/src)

link_directories(lib)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

add_executable(Encode_OligoPool
    src/encode_oligoPool.cpp
)
target_link_libraries(Encode_OligoPool enc lcrc)

add_executable(Encode_LargeDNA
    src/encode_largeDNA.cpp
)
target_link_libraries(Encode_LargeDNA enc lcrc)
