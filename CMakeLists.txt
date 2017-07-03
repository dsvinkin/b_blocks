#for gaussian and poissonian BB
#ADD_DEFINITIONS(-Wall -g) # debugging mode
ADD_DEFINITIONS(-Wall -O3)

#list all source files here
ADD_EXECUTABLE(b_blocks
main_bb.cpp
bayesian_blocks_gaussian.h
bayesian_blocks_gaussian.cpp
bayesian_blocks.h
bayesian_blocks.cpp
light_curve.h
light_curve.cpp
thc_file.h
thc_file.cpp
)

TARGET_LINK_LIBRARIES(b_blocks gsl gslcblas)