#This project implement the MCL algorithm
project(mcl)

aux_source_directory(. SRC_LIST)

include_directories("../utils/")
include_directories("./")

add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} ${OpenMP_LIBRARIES})


