# based on: https://github.com/schrummy14/LIGGGHTS_Flexible_Fibers/blob/master/src/WINDOWS/CMake_patch.zip

INCLUDE(Macros)

# Style
MACRO(ADD_STYLE styleName fileName)
    # Get global style list
    GET_PROPERTY(style_${styleName}_local GLOBAL PROPERTY style_${styleName})
    # Append filename
    SET(style_${styleName}_local ${style_${styleName}_local} ${fileName})
	# Save to global style list
    SET_PROPERTY(GLOBAL PROPERTY style_${styleName} ${style_${styleName}_local})
ENDMACRO()

# Scan for Styles
MACRO(SCAN_STYLE searchExpression fileFilter styleName)
	# Get all files that match filter expression
	FILE(GLOB file_list RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${fileFilter})

	# Go through each file
	FOREACH(fileName ${file_list})
		# add full path
		SET(fileName ${CMAKE_CURRENT_SOURCE_DIR}/${fileName})

		# Read each file
		FILE(READ ${fileName} fileContent)
		# Does the file contain the searchExpression
		SET(pos -1)
		STRING(FIND "${fileContent}" ${searchExpression} pos)
		# If searchExpression is in fileContent
		IF(pos GREATER -1)
			ADD_STYLE(${styleName} ${fileName})
		ENDIF()
	ENDFOREACH()
ENDMACRO()

# Scan for all known styles
MACRO(SCAN_STYLES)
    SCAN_STYLE( ANGLE_CLASS     			angle_*.h      			angle)
    SCAN_STYLE( ATOM_CLASS      			atom_vec_*.h   			atom)
    SCAN_STYLE( BOND_CLASS      			bond_*.h       			bond)
    SCAN_STYLE( COMMAND_CLASS   			*.h          			command)
    SCAN_STYLE( COMPUTE_CLASS   			compute_*.h    			compute)
    SCAN_STYLE( DIHEDRAL_CLASS  			dihedral_*.h   			dihedral)
    SCAN_STYLE( DUMP_CLASS      			dump_*.h       			dump)
    SCAN_STYLE( FIX_CLASS       			fix_*.h        			fix)
    SCAN_STYLE( IMPROPER_CLASS  			improper_*.h   			improper)
    SCAN_STYLE( INTEGRATE_CLASS 			*.h          			integrate)
    SCAN_STYLE( KSPACE_CLASS    			*.h          			kspace)
    SCAN_STYLE( MINIMIZE_CLASS  			min_*.h        			minimize)
    SCAN_STYLE( PAIR_CLASS      			pair_*.h       			pair)
    SCAN_STYLE( REGION_CLASS    			region_*.h     			region)
    SCAN_STYLE( CFD_DATACOUPLING_CLASS 	    cfd_datacoupling_*.h	cfd_datacoupling)
    SCAN_STYLE( CFD_REGIONMODEL_CLASS  	    cfd_regionmodel_*.h  	cfd_regionmodel)
    SCAN_STYLE( LB_CLASS        			*.h          			lb)
    SCAN_STYLE( SPH_KERNEL_CLASS  		    sph_kernel_*.h  		sph_kernel)
    SCAN_STYLE( SURFACE_MODEL  		        surface_model_*.h  		surface_model)
    SCAN_STYLE( NORMAL_MODEL  		        normal_model_*.h  		normal_model)
    SCAN_STYLE( TANGENTIAL_MODEL  		    tangential_model_*.h  	tangential_model)
    SCAN_STYLE( ROLLING_MODEL  			    rolling_model_*.h  	    rolling_model)
    SCAN_STYLE( COHESION_MODEL  			cohesion_model_*.h  	cohesion_model)
    SCAN_STYLE( MESH_MOVER  			    mesh_mover_*.h  	    mesh_mover)
    SCAN_STYLE( MESH_MODULE  			    mesh_module_*.h  	    mesh_module)
    SCAN_STYLE( READER  			    	reader_*.h  	        reader)
ENDMACRO()

MACRO(GEN_STYLE styleName)
	# Get global style list
	GET_PROPERTY(styleFileList GLOBAL PROPERTY style_${styleName})

	# Filename
	SET(fileName ${CMAKE_CURRENT_SOURCE_DIR}/style_${styleName}.h)
	# Delete file if exists
	IF(EXISTS ${fileName})
	    FILE(REMOVE ${fileName})
	ENDIF()

	# Create File
	GETDATETIME(NOW "%Y-%m-%d %H:%M:%S")
	FILE(WRITE ${fileName} "/* created on ${NOW} */\n")

	# Append includes
	FOREACH(styleFile ${styleFileList})
		FILE(APPEND ${fileName} "\#include \"${styleFile}\"\n")
	ENDFOREACH()
ENDMACRO()

MACRO(GENERATE_STYLES)
    # Generate all style files
    GEN_STYLE(angle)
    GEN_STYLE(atom)
    GEN_STYLE(body)
    GEN_STYLE(bond)
    GEN_STYLE(cfd_datacoupling)
    GEN_STYLE(cfd_regionmodel)
    GEN_STYLE(cohesion_model)
    GEN_STYLE(command)
    GEN_STYLE(compute)
    GEN_STYLE(dihedral)
    GEN_STYLE(dump)
    GEN_STYLE(fix)
    GEN_STYLE(improper)
    GEN_STYLE(integrate)
    GEN_STYLE(kspace)
    GEN_STYLE(lb)
    GEN_STYLE(mesh_module)
    GEN_STYLE(mesh_mover)
    GEN_STYLE(minimize)
    GEN_STYLE(normal_model)
    GEN_STYLE(pair)
    GEN_STYLE(reader)
    GEN_STYLE(region)
    GEN_STYLE(rolling_model)
    GEN_STYLE(sph_kernel)
    GEN_STYLE(surface_model)
    GEN_STYLE(tangential_model)
ENDMACRO()