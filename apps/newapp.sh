#!/bin/sh

# This script adds a new application skeleton to the DiSk++
# application tree

if [ -z "$1" ]; then
    echo "Please specify a name for the new application"
    exit 1
fi

APPNAME=$1
APPNAME_UPPER=$(echo ${APPNAME} | tr '[a-z]' '[A-Z]')

if [ -e "$APPNAME" ]; then
    echo "Application \"${APPNAME}\" already exists"
    exit 1
fi


mkdir -p ${APPNAME}/src
touch ${APPNAME}/src/${APPNAME}.cpp

cat << HEREDOC >> ${APPNAME}/CMakeLists.txt
set(LINK_LIBS diskpp)

add_executable(${APPNAME} src/${APPNAME}.cpp)
target_link_libraries(${APPNAME} \${LINK_LIBS})
install(TARGETS ${APPNAME} RUNTIME DESTINATION bin)
HEREDOC

cat << HEREDOC >> CMakeLists.txt

option(BUILD_APP_${APPNAME_UPPER} "Build ${APPNAME} application" ON)
if (BUILD_APP_${APPNAME_UPPER})
    add_subdirectory(${APPNAME})
endif()
HEREDOC
