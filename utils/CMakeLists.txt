#This project implement the LRW clustering algorithm
project(clstUtils)

aux_source_directory(csparse SRC_LIST)
aux_source_directory(. SRC_LIST)

include_directories("./")

add_library(${PROJECT_NAME} ${SRC_LIST})

#target_link_libraries(${PROJECT_NAME})


