if(BUILD_TESTS)
   #add_subdirectory("/Dependencies/Source/Catch_Project")
   include_directories(${CATCH_INCLUDE_DIR} ${headers})
   
   SET(TEST_FILES 
      test_problem.cpp
      assignment_problem.cpp
      dynamic_problem.cpp
      )

   enable_testing()

   add_executable(test_main test_main.cpp ${TEST_FILES})
   target_link_libraries(test_main CS2_CPP)

   #foreach( source_file ${TEST_FILES} )
   #   string( REPLACE ".cpp" "" object ${source_file} )
   #   add_library(${object} STATIC ${source_file}) # for some reason only shared library is linked to test_main
   #   target_link_libraries(test_main ${object})
   #endforeach( source_file ${TEST_FILES} )

   add_test(NAME CS2_CPP_TEST COMMAND test_main)
endif()

