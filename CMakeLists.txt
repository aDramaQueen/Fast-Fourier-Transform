cmake_minimum_required(VERSION 3.5)
project(fft VERSION 0.8.0 LANGUAGES C)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

# Add files
add_executable(fft
   main.c
   fft.c
   testFFT.c
)

# Setup C compiler
set_property(TARGET fft PROPERTY C_STANDARD 11)
target_link_libraries(fft m)

include(CPack)
