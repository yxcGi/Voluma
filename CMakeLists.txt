cmake_minimum_required(VERSION 3.16)
project(Voluma LANGUAGES CXX VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
add_compile_options(-Wall -Wextra -ffp-contract=off)

option(FVM_USE_MPI "Build with MPI" ON)
if(FVM_USE_MPI)
  find_package(MPI COMPONENTS CXX)
  if(NOT MPI_CXX_FOUND)
    message(WARNING "MPI not found, building serial version")
    set(FVM_USE_MPI OFF)
  endif()
endif()

file(GLOB_RECURSE FVM_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/fvm/*.cpp)
add_library(fvm STATIC ${FVM_SOURCES})
target_include_directories(fvm PUBLIC ${CMAKE_SOURCE_DIR}/src)
if(FVM_USE_MPI)
  target_compile_definitions(fvm PUBLIC FVM_USE_MPI)
  target_link_libraries(fvm PUBLIC MPI::MPI_CXX)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
file(GLOB_RECURSE FVM_EXAMPLES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/examples/*.cpp ${CMAKE_SOURCE_DIR}/apps/*.cpp)
foreach(src ${FVM_EXAMPLES})
  get_filename_component(name ${src} NAME_WE)
  add_executable(${name} ${src})
  target_link_libraries(${name} PRIVATE fvm)
endforeach()

enable_testing()
file(GLOB FVM_TESTS CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/tests/*.cpp)
foreach(src ${FVM_TESTS})
  get_filename_component(name ${src} NAME_WE)
  add_executable(${name} ${src})
  target_link_libraries(${name} PRIVATE fvm)
  add_test(NAME ${name} COMMAND ${name} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  if(FVM_USE_MPI)
    add_test(NAME ${name}_np3 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
             $<TARGET_FILE:${name}> WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
endforeach()

# 算例回归测试：运行小算例并用 tools/compare_*.py 按容差核对（需要 python3；圆柱、TGV 需要 numpy）
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  set(FVM_CHECK ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/check_case.py --exe $<TARGET_FILE:fvmFlow>)
  set(FVM_MPI_ARGS "")
  if(FVM_USE_MPI)
    set(FVM_MPI_ARGS --np 2 --mpiexec ${MPIEXEC_EXECUTABLE})
  endif()
  add_test(NAME case_pipe_laminar COMMAND ${FVM_CHECK} --case ${CMAKE_SOURCE_DIR}/cases/pipe/case.json
           --workdir ${CMAKE_BINARY_DIR}/cases/pipe ${FVM_MPI_ARGS}
           -- ${CMAKE_SOURCE_DIR}/tools/compare_pipe.py OUTPUT --tol 1.5)
  add_test(NAME case_cylinder_re40 COMMAND ${FVM_CHECK} --case ${CMAKE_SOURCE_DIR}/cases/cylinderRe40/case.json
           --workdir ${CMAKE_BINARY_DIR}/cases/cylinderRe40 ${FVM_MPI_ARGS}
           -- ${CMAKE_SOURCE_DIR}/tools/compare_cylinder_laminar.py OUTPUT --re 40 --check)
  add_test(NAME case_tgv3d_32 COMMAND ${FVM_CHECK} --case ${CMAKE_SOURCE_DIR}/cases/tgv3d/case32.json
           --workdir ${CMAKE_BINARY_DIR}/cases/tgv3d ${FVM_MPI_ARGS}
           -- ${CMAKE_SOURCE_DIR}/tools/compare_tgv.py OUTPUT/monitor.csv --check 30)
  # 边界条件算例（cases/channelInlet 等，见 docs/validation.md）
  set(FVM_CHECK_SCALAR ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/check_case.py --exe $<TARGET_FILE:fvmScalar>)
  add_test(NAME case_bc_parabolic COMMAND ${FVM_CHECK} --case ${CMAKE_SOURCE_DIR}/cases/channelInlet/parabolic.json
           --workdir ${CMAKE_BINARY_DIR}/cases/bc_parabolic ${FVM_MPI_ARGS}
           -- ${CMAKE_SOURCE_DIR}/tools/compare_channel_inlet.py OUTPUT --mode poiseuille --nu 0.05 --check)
  add_test(NAME case_bc_profile COMMAND ${FVM_CHECK} --case ${CMAKE_SOURCE_DIR}/cases/channelInlet/profile.json
           --workdir ${CMAKE_BINARY_DIR}/cases/bc_profile ${FVM_MPI_ARGS}
           -- ${CMAKE_SOURCE_DIR}/tools/compare_channel_inlet.py OUTPUT --mode poiseuille --nu 0.05 --check)
  add_test(NAME case_bc_fixedFluxPressure COMMAND ${FVM_CHECK} --case ${CMAKE_SOURCE_DIR}/cases/channelInlet/fixedFluxPressure.json
           --workdir ${CMAKE_BINARY_DIR}/cases/bc_fixedFluxPressure ${FVM_MPI_ARGS}
           -- ${CMAKE_SOURCE_DIR}/tools/compare_channel_inlet.py OUTPUT --mode poiseuille --nu 0.05 --check)
  add_test(NAME case_bc_flowRate COMMAND ${FVM_CHECK} --case ${CMAKE_SOURCE_DIR}/cases/channelInlet/flowRate.json
           --workdir ${CMAKE_BINARY_DIR}/cases/bc_flowRate ${FVM_MPI_ARGS}
           -- ${CMAKE_SOURCE_DIR}/tools/compare_channel_inlet.py OUTPUT --mode developed --nu 0.05 --check)
  add_test(NAME case_bc_totalPressure COMMAND ${FVM_CHECK} --case ${CMAKE_SOURCE_DIR}/cases/channelInlet/totalPressure.json
           --workdir ${CMAKE_BINARY_DIR}/cases/bc_totalPressure ${FVM_MPI_ARGS}
           -- ${CMAKE_SOURCE_DIR}/tools/compare_channel_inlet.py OUTPUT --mode bernoulli --check)
  add_test(NAME case_bc_timeTable COMMAND ${FVM_CHECK} --case ${CMAKE_SOURCE_DIR}/cases/channelInlet/timeTable.json
           --workdir ${CMAKE_BINARY_DIR}/cases/bc_timeTable ${FVM_MPI_ARGS}
           -- ${CMAKE_SOURCE_DIR}/tools/compare_channel_inlet.py OUTPUT --mode ramp --check)
  add_test(NAME case_bc_turbulentInlet COMMAND ${FVM_CHECK} --case ${CMAKE_SOURCE_DIR}/cases/channelInlet/rans.json
           --workdir ${CMAKE_BINARY_DIR}/cases/bc_turbulentInlet ${FVM_MPI_ARGS}
           -- ${CMAKE_SOURCE_DIR}/tools/check_fluxes.py OUTPUT --converged)
  add_test(NAME case_bc_rotatingWall COMMAND ${FVM_CHECK} --case ${CMAKE_SOURCE_DIR}/cases/taylorCouette/case.json
           --workdir ${CMAKE_BINARY_DIR}/cases/bc_rotatingWall ${FVM_MPI_ARGS}
           -- ${CMAKE_SOURCE_DIR}/tools/compare_taylor_couette.py OUTPUT --tol 0.01 --check)
  add_test(NAME case_bc_advective COMMAND ${FVM_CHECK_SCALAR} --case ${CMAKE_SOURCE_DIR}/cases/advectiveOutlet/case.json
           --workdir ${CMAKE_BINARY_DIR}/cases/bc_advective
           -- ${CMAKE_SOURCE_DIR}/tools/compare_advective.py OUTPUT --check)
  add_test(NAME case_bc_inletOutlet_backflow COMMAND ${FVM_CHECK} --case ${CMAKE_SOURCE_DIR}/cases/jetBox/case.json
           --workdir ${CMAKE_BINARY_DIR}/cases/bc_jetBox ${FVM_MPI_ARGS}
           -- ${CMAKE_SOURCE_DIR}/tools/check_fluxes.py OUTPUT --expect-backflow outlet)
  set_tests_properties(case_pipe_laminar case_cylinder_re40 case_tgv3d_32 case_bc_parabolic case_bc_profile
                       case_bc_fixedFluxPressure case_bc_flowRate case_bc_totalPressure case_bc_timeTable
                       case_bc_turbulentInlet case_bc_rotatingWall case_bc_advective case_bc_inletOutlet_backflow
                       PROPERTIES LABELS "case" TIMEOUT 1200)
endif()
