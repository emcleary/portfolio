DIR="build"
BUILD_TYPE="Release"
BOOL_VERBOSE="OFF"

all:
	mkdir -p $(DIR)
	cmake -S . -B $(DIR) -D CMAKE_BUILD_TYPE=$(BUILD_TYPE) -D CMAKE_VERBOSE_MAKEFILE:BOOL=$(BOOL_VERBOSE)
	cmake --build $(DIR)

clean:
	rm -vrf $(DIR)

verbose:
	$(eval BOOL_VERBOSE="ON")

debug:
	$(eval BUILD_TYPE="Debug")
