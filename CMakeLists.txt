cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(flagser)

enable_language(CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Threads REQUIRED)

find_package(HDF5)
IF(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIR})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_HDF5")
ENDIF()

IF(MSVC)
  set(CMAKE_CXX_RELEASE_FLAGS "${CMAKE_CXX_FLAGS} -Ox -D NDEBUG")
ELSE()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -Wall -pthread -pedantic-errors")
  set(CMAKE_CXX_RELEASE_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -D NDEBUG")
ENDIF()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_POSITION_INDEPENDANT}")
set(CMAKE_VS_PLATFORM_NAME "x64")
set(CMAKE_VS_PLATFORM_TOOLSET_HOST_ARCHITECTURE "x64")

# Also add the time of the file + output file into consideration
# Add the include dir to the watch list

add_custom_command(
    COMMAND ../build_custom_algorithms.sh
    DEPENDS parser/math_parser.py algorithms.math
    OUTPUT include/algorithms.h
    COMMENT "Check for custom algorithms"
)
add_custom_target(generate_custom_algorithms DEPENDS include/algorithms.h)


add_executable(flagser src/flagser.cpp include/algorithms.h)
add_dependencies(flagser generate_custom_algorithms)
set_target_properties(flagser PROPERTIES LINKER_LANGUAGE CXX)

add_executable(flagser_retrieve_persistence src/flagser.cpp include/algorithms.h)
add_dependencies(flagser_retrieve_persistence generate_custom_algorithms)
target_compile_definitions(flagser_retrieve_persistence PRIVATE RETRIEVE_PERSISTENCE)
set_target_properties(flagser_retrieve_persistence PROPERTIES LINKER_LANGUAGE CXX)


add_executable(flagser-memory src/flagser.cpp include/algorithms.h)
target_compile_definitions(flagser-memory PRIVATE KEEP_FLAG_COMPLEX_IN_MEMORY=1)
add_dependencies(flagser-memory generate_custom_algorithms)
set_target_properties(flagser-memory PROPERTIES LINKER_LANGUAGE CXX)

add_executable(ripser src/ripser.cpp)
set_target_properties(ripser PROPERTIES LINKER_LANGUAGE CXX)

add_executable(flagser-count src/flagser-count.cpp)
set_target_properties(flagser-count PROPERTIES LINKER_LANGUAGE CXX)

add_executable(er src/er.cpp)
set_target_properties(er PROPERTIES LINKER_LANGUAGE CXX)

IF(NOT MSVC)
add_executable(flagser-coefficients src/flagser.cpp)
add_dependencies(flagser-coefficients generate_custom_algorithms)
target_compile_definitions(flagser-coefficients PRIVATE USE_COEFFICIENTS=1)
set_target_properties(flagser-coefficients PROPERTIES LINKER_LANGUAGE CXX)

add_executable(flagser-coefficients-memory src/flagser.cpp include/algorithms.h)
target_compile_definitions(flagser-coefficients-memory PRIVATE USE_COEFFICIENTS=1 KEEP_FLAG_COMPLEX_IN_MEMORY=1)
add_dependencies(flagser-coefficients-memory generate_custom_algorithms)
set_target_properties(flagser-coefficients-memory PROPERTIES LINKER_LANGUAGE CXX)

add_executable(ripser-coefficients src/ripser.cpp)
target_compile_definitions(ripser-coefficients PRIVATE USE_COEFFICIENTS=1)
set_target_properties(ripser-coefficients PROPERTIES LINKER_LANGUAGE CXX)
ELSE()
  message("On Windows, compilation with coefficients is currently not supported.")
ENDIF()

message("")
message("** CMAKE_CXX_FLAGS: \t\t\t${CMAKE_CXX_FLAGS}")
message("** CMAKE_CXX_RELEASE_FLAGS: \t\t\t${CMAKE_CXX_RELEASE_FLAGS}")
message("** CMAKE_CXX_LINK_FLAGS: \t\t${CMAKE_CXX_LINK_FLAGS}")
message("")

IF(HDF5_FOUND)
TARGET_LINK_LIBRARIES(flagser ${HDF5_LIBRARIES} pthread)
TARGET_LINK_LIBRARIES(flagser_retrieve_persistence ${HDF5_LIBRARIES} pthread)
TARGET_LINK_LIBRARIES(flagser-memory ${HDF5_LIBRARIES} pthread)
TARGET_LINK_LIBRARIES(flagser-count ${HDF5_LIBRARIES} pthread)
TARGET_LINK_LIBRARIES(er ${HDF5_LIBRARIES} pthread)
TARGET_LINK_LIBRARIES(ripser ${HDF5_LIBRARIES} pthread)
IF(NOT MSVC)
TARGET_LINK_LIBRARIES(ripser-coefficients ${HDF5_LIBRARIES} pthread)
TARGET_LINK_LIBRARIES(flagser-coefficients ${HDF5_LIBRARIES} pthread)
TARGET_LINK_LIBRARIES(flagser-coefficients-memory ${HDF5_LIBRARIES} pthread)
ENDIF()
ENDIF()

IF(NOT MSVC)
install(TARGETS flagser flagser-coefficients flagser-memory flagser-coefficients-memory flagser-count RUNTIME DESTINATION bin)
ELSE()
install(TARGETS flagser flagser-memory flagser-count RUNTIME DESTINATION bin)
ENDIF()

add_custom_command(
    TARGET flagser POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:flagser> ${CMAKE_SOURCE_DIR})
add_custom_command(
    TARGET flagser_retrieve_persistence POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:flagser_retrieve_persistence> ${CMAKE_SOURCE_DIR})
add_custom_command(
    TARGET flagser-memory POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:flagser-memory> ${CMAKE_SOURCE_DIR})
add_custom_command(
    TARGET flagser-count POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:flagser-count> ${CMAKE_SOURCE_DIR})
add_custom_command(
    TARGET ripser POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:ripser> ${CMAKE_SOURCE_DIR})
add_custom_command(
    TARGET er POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:er> ${CMAKE_SOURCE_DIR}/tools)
IF(NOT MSVC)
add_custom_command(
    TARGET flagser-coefficients POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:flagser-coefficients> ${CMAKE_SOURCE_DIR})
add_custom_command(
    TARGET flagser-coefficients-memory POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:flagser-coefficients-memory> ${CMAKE_SOURCE_DIR})
add_custom_command(
    TARGET ripser-coefficients POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:ripser-coefficients> ${CMAKE_SOURCE_DIR})
ENDIF()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
    add_subdirectory(test)
    enable_testing()
endif()
