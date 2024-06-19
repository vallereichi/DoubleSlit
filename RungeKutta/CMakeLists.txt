cmake_minimum_required(VERSION 3.5)
project(doubleSlit2d)

add_executable(${PROJECT_NAME} doubleSlit.cpp)

add_subdirectory(src)
target_include_directories(${PROJECT_NAME} PRIVATE src/include)
target_link_directories(${PROJECT_NAME} PRIVATE src)
target_link_libraries(${PROJECT_NAME} Functions)