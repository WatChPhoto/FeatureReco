cmake_minimum_required(VERSION 2.8)
project( UndistortImage )
find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )
add_executable( UndistortImage UndistortImage.cpp )
target_link_libraries( UndistortImage ${OpenCV_LIBS} )

