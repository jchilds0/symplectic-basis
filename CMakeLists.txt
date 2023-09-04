# To compile project 
# - run 'cmake -S . -B build/'
# - run 'cmake --build build/'

cmake_minimum_required(VERSION 3.12)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

project(SnapPea)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Set the location of the SnapPea kernel code
set(SNAPPEA_KERNEL ${PROJECT_SOURCE_DIR}/snappy_src/)

# Compiler options
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -pedantic -fanalyzer")

# Add include directories
include_directories(${SNAPPEA_KERNEL}/addl_code ${SNAPPEA_KERNEL}/headers ${SNAPPEA_KERNEL}/real_type ${SNAPPEA_KERNEL}/unix_kit ${SNAPPEA_KERNEL}/)
include_directories(${PROJECT_SOURCE_DIR}/symp_src)

# Add source files
file(GLOB KERNEL_SOURCES ${SNAPPEA_KERNEL}/*/*.c)
file(GLOB HEADER_FILES ${SNAPPEA_KERNEL}/*/*.h)
file(GLOB SYMP_SRC ${PROJECT_SOURCE_DIR}/symp_src/*.c ${PROJECT_SOURCE_DIR}/symp_src/*.h)
file(COPY ${PROJECT_SOURCE_DIR}/dev/CuspedCensusData DESTINATION ${PROJECT_SOURCE_DIR}/cmake-build-debug)      # Copy census data to be used by program

# Create executable target
add_executable(symplectic_basis_main ${PROJECT_SOURCE_DIR}/dev/symplectic_basis_main.c ${SYMP_SRC} ${KERNEL_SOURCES} ${HEADER_FILES} ${CUSPED_CENSUS})

# Link math library
target_link_libraries(symplectic_basis_main m)

# Define custom target to create BuildDate file
#add_custom_target(BuildDate ALL
#        COMMAND date > ${CMAKE_SOURCE_DIR}/BuildDate
#        DEPENDS ${KERNEL_SOURCES} ${HEADER_FILES}
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#        )

# testing binary
add_executable(symplectic_basis_test ${PROJECT_SOURCE_DIR}/dev/symplectic_basis_test.c ${SYMP_SRC} ${KERNEL_SOURCES} ${HEADER_FILES} ${CUSPED_CENSUS})
target_link_libraries(symplectic_basis_test m)

# enable testing functionality
enable_testing()

# define tests
add_test(
        NAME symplectic_test
        COMMAND $<TARGET_FILE:symplectic_basis_test>
)

# Clean target
add_custom_target(clean-all
        COMMAND ${CMAKE_COMMAND} -E remove symplectic_basis_main symplectic_basis_main.o
        COMMAND ${CMAKE_COMMAND} -E remove symplectic_basis_test symplectic_basis_test.o
        COMMAND ${CMAKE_COMMAND} -E remove_directory KernelObjects
        COMMAND ${CMAKE_COMMAND} -E remove *.pyc
        )

