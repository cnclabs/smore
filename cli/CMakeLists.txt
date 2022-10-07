# file: CMakeLists.txt
# date: 2022-10-07


add_executable(app ${CMAKE_CURRENT_SOURCE_DIR}/app.cpp)
target_link_libraries(app smore)

add_executable(bpr ${CMAKE_CURRENT_SOURCE_DIR}/bpr.cpp)
target_link_libraries(bpr smore)

add_executable(deepwalk ${CMAKE_CURRENT_SOURCE_DIR}/deepwalk.cpp)
target_link_libraries(deepwalk smore)

add_executable(eco ${CMAKE_CURRENT_SOURCE_DIR}/eco.cpp)
target_link_libraries(eco smore)

add_executable(gcn ${CMAKE_CURRENT_SOURCE_DIR}/gcn.cpp)
target_link_libraries(gcn smore)

add_executable(hoprec ${CMAKE_CURRENT_SOURCE_DIR}/hoprec.cpp)
target_link_libraries(hoprec smore)

add_executable(hpe ${CMAKE_CURRENT_SOURCE_DIR}/hpe.cpp)
target_link_libraries(hpe smore)

add_executable(line ${CMAKE_CURRENT_SOURCE_DIR}/line.cpp)
target_link_libraries(line smore)

add_executable(mf ${CMAKE_CURRENT_SOURCE_DIR}/mf.cpp)
target_link_libraries(mf smore)

add_executable(nemf ${CMAKE_CURRENT_SOURCE_DIR}/nemf.cpp)
target_link_libraries(nemf smore)

add_executable(nerank ${CMAKE_CURRENT_SOURCE_DIR}/nerank.cpp)
target_link_libraries(nerank smore)

add_executable(skewopt ${CMAKE_CURRENT_SOURCE_DIR}/skewopt.cpp)
target_link_libraries(skewopt smore)

add_executable(textgcn ${CMAKE_CURRENT_SOURCE_DIR}/textgcn.cpp)
target_link_libraries(textgcn smore)

add_executable(textgcndev ${CMAKE_CURRENT_SOURCE_DIR}/textgcndev.cpp)
target_link_libraries(textgcndev smore)

add_executable(walklets ${CMAKE_CURRENT_SOURCE_DIR}/walklets.cpp)
target_link_libraries(walklets smore)

add_executable(warp ${CMAKE_CURRENT_SOURCE_DIR}/warp.cpp)
target_link_libraries(warp smore)
