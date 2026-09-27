in /work/noaa/wrf-chem/pagowski/jedi/code/jedi-bundle/crtm/test/CMakeLists.txt

create parallel tests for AAOD to those for AOD ie.

  Test  #4: test_Unit_AerosolScatter_TL
  Test  #5: test_Unit_AerosolScatter_AD
  Test  #6: test_Unit_AerosolScatter_K

add

test_Unit_AerosolScatter_AAOD_TL
test_Unit_AerosolScatter_AAOD_AD
test_Unit_AerosolScatter_AAOD_K

also 
list( APPEND AAOD_Sensor_Ids
        cris-fsr_n21
        v.abi_g18
        cris399_npp
        v.abi_gr
        abi_g18
        airs_aqua
)

list (APPEND common_tests
     Simple
     AOD
     AAOD

and run to see if they pass