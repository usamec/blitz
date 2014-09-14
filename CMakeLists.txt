cmake_minimum_required(VERSION 2.8)
project (Blitz)

list( APPEND CMAKE_CXX_FLAGS "-std=c++0x -g -O2 ${CMAKE_CXX_FLAGS}")

add_library(dna_seq dna_seq.cc)

add_executable(test_dna_seq test_dna_seq.cc)
target_link_libraries(test_dna_seq dna_seq)

add_executable(blitz blitz.cc)
target_link_libraries(blitz dna_seq)

add_executable(bpb debruijn-pacbio.cc)
target_link_libraries(bpb dna_seq)

add_executable(tc test-clustering4.cc)
target_link_libraries(tc dna_seq)

add_executable(tc3 test-clustering3.cc)
target_link_libraries(tc3 dna_seq)

add_executable(ti test-illum.cc)
target_link_libraries(ti dna_seq)
