# Kokkos

## Install

> WARN: build and install directories for Kokkos should really be different.

1. Create a `~/Kokkos` dir.
2. `cd` into `~/Kokkos`. Download kokkos via  `git clone git@github.com:kokkos/kokkos.git`
3. Make sure to have modified Kokkos environment variables in `~/.bashrc`:

```
### Kokkos
export KOKKOS_SRC_DIR=${HOME}/Kokkos/kokkos
export KOKKOS_BUILD_DIR=${KOKKOS_SRC_DIR}/build
export KOKKOS_INSTALL_DIR=${HOME}/Kokkos/kokkos_install
```

Run `build_install_kokkos.sh`


## Cmake

`CMAKE_PREFIX_PATH` works as a build directive, rather than as an environment variable. Moreover, you may perform the build into a dedicated temporary directory (it's cleaner, because when done, you can remove that temporary directory and you get back a clean pristine source tree). [StackOverflow](https://stackoverflow.com/questions/8019505/how-to-set-the-cmake-prefix-path)

[Semicolon-separated list](https://cmake.org/cmake/help/latest/manual/cmake-language.7.html#cmake-language-lists) of directories specifying installation *prefixes* to be searched by the `find_package()`, `find_program()`, `find_library()`, `find_file()`, and `find_path()` commands. Each command will add appropriate subdirectories (like `bin`, `lib`, or `include`) as specified in its own documentation. [cmake.org](https://cmake.org/cmake/help/latest/variable/CMAKE_PREFIX_PATH.html)

By default this is empty. It is intended to be set by the project.

The system directories that are contained in `CMAKE_PREFIX_PATH` are locations that typically include installed software.
