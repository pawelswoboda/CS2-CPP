function(DOWNLOAD_FILE URL DOWNLOAD_TO)
   if(NOT EXISTS ${DOWNLOAD_TO})
      message("Download ${URL}")
      file(DOWNLOAD "${URL}" "${DOWNLOAD_TO}")
   endif()
endfunction(DOWNLOAD_FILE)

if(BUILD_TESTS)
   include_directories(${CATCH_INCLUDE_DIR} ${headers})
   
   SET(TEST_FILES 
      test_problem.cpp
      assignment_problem.cpp
      dynamic_problem.cpp
      )

   enable_testing()

   add_executable(test_main test_main.cpp ${TEST_FILES})

   if(BUILD_GTE_EVALUATION)
      target_sources(test_main PUBLIC evaluate_gte_instances.cpp)
      file(MAKE_DIRECTORY "instances")
      # Download gte instances
      SET(GTE_FILES 
         gte_bad.20
         gte_bad.40
         gte_bad.60
         gte_bad.160
         gte_bad.200
         gte_bad.460
         gte_bad.510
         gte_bad.1160
         gte_bad.1700
         gte_bad.6410
         gte_bad.6830
         gte_bad.15100
         gte_bad.15710
         gte_bad.35620
         gte_bad.49320
         gte_bad.60090
         gte_bad.65330
         gte_bad.176280
         gte_bad.298300
         gte_bad.451760
         gte_bad.469010
         gte_bad.508829
         )
      file(MAKE_DIRECTORY "instances/gte")
      foreach(INSTANCE ${GTE_FILES})
         DOWNLOAD_FILE("http://elib.zib.de/pub/Packages/mp-testdata/mincost/gte/${INSTANCE}" "test/instances/gte/${INSTANCE}")
      endforeach()
   endif()

   if(BUILD_NETGEN_EVALUATION)
      target_sources(test_main PUBLIC evaluate_netgen_instances.cpp)
      SET(NETGEN_FILES 
         big1.net   
         big2.net   
         big3.net   
         big4.net   
         big5.net   
         big6.net   
         big7.net   
         big8.net   
         cap1.net   
         cap10.net   
         cap11.net   
         cap12.net   
         cap13.net   
         cap14.net   
         cap15.net   
         cap16.net   
         cap17.net   
         cap18.net   
         cap19.net   
         cap2.net   
         cap20.net   
         cap21.net   
         cap22.net   
         cap23.net   
         cap24.net   
         cap25.net   
         cap26.net   
         cap27.net   
         cap28.net   
         cap29.net   
         cap3.net   
         cap30.net   
         cap31.net   
         cap32.net   
         cap33.net   
         cap34.net   
         cap35.net   
         cap36.net   
         cap37.net   
         cap38.net   
         cap39.net   
         cap4.net   
         cap40.net   
         cap41.net   
         cap5.net   
         cap6.net   
         cap7.net   
         cap8.net   
         cap9.net   
         stndrd1.net   
         stndrd10.net   
         stndrd16.net   
         stndrd17.net   
         stndrd18.net   
         stndrd19.net   
         stndrd2.net   
         stndrd20.net   
         stndrd21.net   
         stndrd22.net   
         stndrd23.net   
         stndrd24.net   
         stndrd25.net   
         stndrd26.net   
         stndrd27.net   
         stndrd28.net   
         stndrd29.net   
         stndrd3.net   
         stndrd30.net   
         stndrd31.net   
         stndrd32.net   
         stndrd33.net   
         stndrd34.net   
         stndrd35.net   
         stndrd36.net   
         stndrd37.net   
         stndrd38.net   
         stndrd39.net   
         stndrd4.net   
         stndrd40.net   
         stndrd41.net   
         stndrd42.net   
         stndrd43.net   
         stndrd44.net   
         stndrd45.net   
         stndrd46.net   
         stndrd47.net   
         stndrd48.net   
         stndrd5.net   
         stndrd50.net   
         stndrd51.net   
         stndrd52.net   
         stndrd53.net   
         stndrd54.net   
         stndrd6.net   
         stndrd7.net   
         stndrd8.net   
         stndrd9.net   
         transp1.net   
         transp10.net   
         transp11.net   
         transp12.net   
         transp13.net   
         transp14.net   
         transp2.net   
         transp3.net   
         transp4.net   
         transp5.net   
         transp6.net   
         transp7.net   
         transp8.net   
         transp9.net
         )
      file(MAKE_DIRECTORY "instances/netgen")
      foreach(INSTANCE ${NETGEN_FILES})
         DOWNLOAD_FILE("http://elib.zib.de/pub/Packages/mp-testdata/mincost/netg/${INSTANCE}" "test/instances/netgen/${INSTANCE}")
      endforeach()
   endif()

   add_test(NAME CS2_CPP_TEST COMMAND test_main)
endif()

