"""
This is a simple build script that uses the clang++ (or g++) compiler to compile the C++ program. Call from the command line as follows: ´python ./dev/build.py´
You may have to use python3 instead. To see available options, use the -h flag: ´python ./dev/build.py -h´
Note, that this script is multithreaded. The /bin directory is created if it does not exist and all files in it are deleted.

Features:
- Syntax highlighting for C++ when logging errors to terminal
- Multithreaded compilation
- Logging (./bin/build.log)
- Automatically running the binary after compilation
"""

# Define the compilers and debugger
cpp_compiler = 'clang++'
cuda_compiler = 'nvcc'  # TODO: Add support for CUDA
debugger = 'lldb'
linker = cpp_compiler

# Define the source directory and output binary
src_dirs = ['./src', './test']
include_dirs = ['./include', './test', './mechanism']
build_dir = './bin'
output_binary_name = 'main'
submodules = dict(
    sundials = dict(
        include_dir = './submodules/sundials/install/include',
        lib_dir = './submodules/sundials/install/lib',
        core_lib = 'libsundials_core',
    ),
    hdf5 = dict(
        include_dir = './submodules/hdf5_build/include',
        lib_dir = './submodules/hdf5_build/bin',
        core_lib = 'libhdf5',
    ),
)

# Other settings
compiler_flags = [
    '-std=c++20',              # Use the C++20 standard
    '-march=native',           # Enable all instruction set extensions supported by the CPU
    '-ftree-vectorize',        # Enable tree vectorization
    #'-fvectorize',             # Enable auto-vectorization (only for clang)
    '-pedantic',               # Enforce strict standard compliance
    #'-funroll-loops',          # Enable loop unrolling
    #'-mlong-double-80',        # Use 80-bit long double. You may also try with -128, however it is less likely to work. Fails with clang++ on win11 with c++20
# Error handling flags: slight performance boost, harder error identification, might cause unexpected behavior
    #'-fno-exceptions',         # Disable exception handling
    #'-fno-math-errno',         # Do not set errno after math library functions
    #'-ffinite-math-only',      # Assume that floating-point arithmetic arguments and results are finite numbers   
]
linker_flags = [
    '-flto',                  # Enable link-time optimization
    '-Wl', '-O2',             # Linker optimization
]
common_flags = [
# Sanitizers: help identify issues in runtime, performance overhead
    '-fsanitize=address',     # Protect against memory errors (a.g. use after free)
    '-fsanitize=undefined',   # Protect against undefined behavior (e.g. integer overflow, invalid type casts)
    '-fno-omit-frame-pointer',# Keep frame pointer for better stack traces
    #'-fomit-frame-pointer',    # Omit frame pointer for functions that don't need one
]


import argparse
import json
from build_utility import Builder, Logger

def main():
    args = parse_args()

    global compiler_flags, linker_flags, common_flags
    if args.test:
        compiler_flags.append('-DTEST')
    if args.benchmark:
        compiler_flags.append('-DBENCHMARK')
    if args.warning:
        common_flags.extend(['-Wall', '-Wextra', '-Werror'])
    if args.debug:
        common_flags.append('-g')
        common_flags.append('-Og')
    if args.optimize2 and not args.debug:
        compiler_flags.append('-O2')
    if args.optimize3 and not args.optimize2 and not args.debug:
        compiler_flags.append('-O3')
    if not args.optimize2 and not args.optimize3:
        compiler_flags.append('-O0')

    builder = Builder(build_dir)
    log_global_settings(builder.logger)

    for submodule in submodules.values():
        include_dirs.append(submodule['include_dir'])
    cpp_files, cuda_files, include_flags = builder.gather_files(src_dirs, include_dirs)
    compiler_flags += common_flags + include_flags
    linker_flags += common_flags
    
    for cpp_file in cpp_files:
        builder.add_source_file(cpp_file, cpp_compiler, compiler_flags)
    #for cuda_file in cuda_files:
    #    builder.add_source_file(cuda_file, cuda_compiler, compiler_flags)

    object_files = builder.compile_all()
    builder.link(
        linker=linker,
        obj_files=object_files,
        submodules=submodules,
        output_binary_name=output_binary_name,
        linker_flags=linker_flags,
        shared=args.shared
    )

    if args.run and not args.shared:
        if args.debug:
            ret = builder.run_binary(debugger)
        else:
            ret = builder.run_binary()
        exit(ret)
    exit(0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build script")
    parser.add_argument('-t', '--test', action='store_true', help='build with test flag (TEST is defined)')
    parser.add_argument('-b', '--benchmark', action='store_true', help='build with benchmark flag (BENCHMARK is defined)')
    parser.add_argument('-r', '--run', action='store_true', help='run the binary after building')
    parser.add_argument('-w', '--warning', action='store_true', help='enable all warnings')
    parser.add_argument('-d', '--debug', action='store_true', help='compile in debug mode (also runs gdb if called with -r)')
    parser.add_argument('-s', '--shared', action='store_true', help='build as shared library (.dll/.so)')
    parser.add_argument('-o2', '-O2', '--optimize2', action='store_true', help='enable o2 optimization')
    parser.add_argument('-o3', '-O3', '--optimize3', action='store_true', help='enable o3 optimization')
    args = parser.parse_args()
    return args


def log_global_settings(logger: Logger):
    submodules_str = json.dumps(submodules, indent=4)
    logger.log(f'cpp_compiler = \'{cpp_compiler}\'')
    logger.log(f'cuda_compiler = \'{cuda_compiler}\'')
    logger.log(f'linker = \'{linker}\'')
    logger.log(f'src_dirs = {src_dirs}')
    logger.log(f'include_dirs = {include_dirs}')
    logger.log(f'build_dir = \'{build_dir}\'')
    logger.log(f'output_binary_name = \'{output_binary_name}\'')
    logger.log(f'submodules = {submodules_str}')
    logger.log(f'compiler_flags = {compiler_flags}')
    logger.log(f'linker_flags = {linker_flags}')
    logger.log(f'common_flags = {common_flags}\n')

if __name__ == "__main__":
    main()