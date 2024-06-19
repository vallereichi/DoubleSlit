INCLUDE_DIRECTORIES(
    ./include
    ~/libs/eigen-3.4.0
)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
LINK_LIBRARIES(Eigen3::Eigen)

add_library(wavePacket wavePacket.cpp)