cmake_minimum_required(VERSION 2.8)
project( UndistortImage )
find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )
add_executable( UndistortImage UndistortImage.cpp )
target_link_libraries( UndistortImage ${OpenCV_LIBS} )


project( FindBoltLocations )
find_package( OpenCV REQUIRED )
find_package( ROOT 6.18 CONFIG REQUIRED)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
ADD_LIBRARY( featurereco_lib Configuration.cpp MedianTextReader.cpp featureFunctions.cpp bwlabel.cpp ellipse_detection.cpp to_text_functions.cpp)
include_directories( ${OpenCV_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})
add_executable( FindBoltLocations FindBoltLocations.cpp )
target_link_libraries( FindBoltLocations featurereco_lib ${OpenCV_LIBS} ${ROOT_LIBS} ROOT::Core ROOT::Hist)

