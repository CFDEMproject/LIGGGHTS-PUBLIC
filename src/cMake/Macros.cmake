# based on: https://github.com/schrummy14/LIGGGHTS_Flexible_Fibers/blob/master/src/WINDOWS/CMake_patch.zip

MACRO(GET_SUBDIRS retval filename)
	FILE(GLOB file-list RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} */${filename})

	SET(dir-list "")
	FOREACH(item ${file-list})
		GET_FILENAME_COMPONENT(file-path ${item} PATH)
		SET(dir-list ${dir-list} ${file-path})
	ENDFOREACH()
	SET(${retval} ${dir-list})
ENDMACRO()

MACRO(GETDATETIME result format)
	STRING(TIMESTAMP ${result} ${format})
ENDMACRO()
