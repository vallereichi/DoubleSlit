cmake_minimum_required(VERSION 3.5)
project(doubleslit2d)

add_executable(${PROJECT_NAME} doubleSlit2d.cpp)

add_subdirectory(src)
target_include_directories(${PROJECT_NAME} PRIVATE src/include)
target_link_directories(${PROJECT_NAME} PRIVATE src)
target_link_libraries(${PROJECT_NAME} wavePacket)

