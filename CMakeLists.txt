cmake_minimum_required(VERSION 2.8.3)
project(kinematics_calibration)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(RBDL REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(
 include
 ${EIGEN3_INCLUDE_DIRS}    
 ${RBDL_INCLUDE_DIR}
)


add_executable(${PROJECT_NAME}_node 
  src/main2.cpp
  src/robot_model/franka_panda_model.cpp
)

target_link_libraries(${PROJECT_NAME}_node
  ${catkin_LIBRARIES}  
  ${RBDL_LIBRARY}
  pthread
)

add_executable(${PROJECT_NAME}_node9
  src/main_KETI.cpp
  src/robot_model/franka_panda_model.cpp
)

target_link_libraries(${PROJECT_NAME}_node9
  ${catkin_LIBRARIES}  
  ${RBDL_LIBRARY}
  pthread
)

add_executable(hc_closed_chain_calibration_node
  src/hc_closed_chain.cpp
  src/robot_model/franka_panda_model.cpp
)

target_link_libraries(hc_closed_chain_calibration_node
  ${catkin_LIBRARIES}  
  ${RBDL_LIBRARY}
  pthread
)

add_executable(data_collector
  src/data_collector.cpp
  src/robot_model/franka_panda_model.cpp
)

target_link_libraries(data_collector
  ${catkin_LIBRARIES}  
  ${RBDL_LIBRARY}
  pthread
)
