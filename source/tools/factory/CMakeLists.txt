
#set(INCLUDES_VAR TOOL_INCLUDES)

add_includes(${CMAKE_CURRENT_SOURCE_DIR})
add_sources(
        CSFactory.h
        CSRegistrar.h
        CSPlugin.h
        Parameter.h
        Parameter.cpp
        Types.h
        Types.cpp
    )
