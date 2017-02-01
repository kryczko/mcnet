# Compile everything under 'dependencies'
add_subdirectory(deps)

# Add include directories
include_directories(
	"." 
	deps
	deps/yaml-cpp-release-0.5.3/include
)

# Find source files
aux_source_directory(
	"." 
	mcnet-src
)

aux_source_directory(
	"deps" 
	deps_src
)

# Add compilation target
add_executable( 
	mcnet
	${mcnet-src} 
	${deps_src}
)

add_library(mcnet-lib
        SHARED
	${mcnet-src} 

	${deps_src}
)

target_link_libraries (
	mcnet 
	yaml-cpp
)

target_link_libraries (
	mcnet-lib
	yaml-cpp
)