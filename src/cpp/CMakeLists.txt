# Project
add_library(${PROJECT_NAME} SHARED "${PROJECT_NAME}.cpp" "${PROJECT_NAME}.h")
set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES PUBLIC_HEADER ${PROJECT_NAME}.h
)
target_link_libraries(
    ${PROJECT_NAME}
    tardigrade_error_tools
    tardigrade_constitutive_tools
)

target_compile_options(${PROJECT_NAME} PUBLIC)

if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(
        ${PROJECT_NAME}
        PUBLIC
            "${tardigrade_vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
            "${tardigrade_error_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
            "${tardigrade_constitutive_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
    )
endif()

install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
