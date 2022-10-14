cmake_minimum_required(VERSION 3.23)
project(realrootsofcubic)

set(CMAKE_CXX_STANDARD 23)

include_directories(${PROJECT_SOURCE_DIR}/excerpt/include)
add_executable(realrootsofcubic main.cpp)
target_link_libraries(realrootsofcubic PRIVATE ${PROJECT_SOURCE_DIR}/excerpt/libexcerpt.dll)
file(COPY ${PROJECT_SOURCE_DIR}/excerpt/libexcerpt.dll
        DESTINATION ${PROJECT_BINARY_DIR})