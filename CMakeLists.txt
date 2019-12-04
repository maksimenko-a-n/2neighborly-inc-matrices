cmake_minimum_required(VERSION 2.8)

#--- имя проекта совпадает с именем каталога проекта
get_filename_component(project_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)

#--- имя главной цели. CXX значает, что проект на с++
project( ${project_name} C CXX )

set(CMAKE_CXX_STANDARD 14)

set(source_code src)
set(headers ${source_code}/headers)
set(implementation ${source_code}/implementation)

message(
"
config:
    source_code: ${source_code}
    headers: ${headers}
    implementation: ${implementation}
"
)

########################################################################

#--- обработка хедеров
if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${headers})
    list(APPEND path_headers ${CMAKE_CURRENT_SOURCE_DIR}/${headers})
    FILE (GLOB_RECURSE files_headers ${path_headers}/*)
else()
    message(SEND_ERROR "headers not found")
endif()

if (files_headers)
    message("files_headers not empty")
else()
    message(WARNING "files_headers empty")
endif()

########################################################################

#--- обработка cpp-файлов
if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${implementation})
    list(APPEND path_implementation ${CMAKE_CURRENT_SOURCE_DIR}/${implementation})
    FILE (GLOB_RECURSE files_implementation ${CMAKE_CURRENT_SOURCE_DIR}/${implementation}/*)
else()
    message(SEND_ERROR "implementation not found")
endif()

if (files_implementation)
    message("files_implementation not empty")
else()
    message(WARNING "files_implementation empty")
endif()

########################################################################

#add_compile_options(-O2)

add_executable(
    ${project_name}
    ${files_headers}
    ${files_implementation}
    ${source_code}/main.cpp
#    src/test.h
#    src/generate_table.h
)

target_compile_options(${project_name} PUBLIC -O2)

target_include_directories( ${project_name} PUBLIC
    ${path_headers}
)

message("")

#add_subdirectory(libs/bliss)

#target_include_directories( bliss PUBLIC
#    ${source_code}
#    )

#target_link_libraries(${project_name} bliss)

#target_include_directories(
#    bliss PRIVATE
#    libs/bliss
#)

