cmake_minimum_required(VERSION 3.10)
project(SEIAR)

#set(CMAKE_CXX_STANDARD 11)
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

add_executable(mcmc_SEIAR SEIAR_model.h SEIAR_model.c particle_mcmc.h particle_mcmc.c mcmc_SEIAR.c)
add_executable(test_like test.c SEIAR_model.h SEIAR_model.c particle_mcmc.h particle_mcmc.c)

target_link_libraries(mcmc_SEIAR GSL::gsl GSL::gslcblas)
target_link_libraries(test_like GSL::gsl GSL::gslcblas)

