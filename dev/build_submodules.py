"""
This is a simple script to check, download, build dependencies. Call from the command line as follows: ´python ./dev/build_submodules.py´. 
You may have to use python3 instead. Make sure, all dependencies are ok, and sundials is built correctly.

The following things should happen:
- Check if all necessary tools are installed (git, make, cmake, gcc, g++, gdb, clang, clang++, lldb)
- Check if all necessary python packages are installed (numpy, matplotlib, pygments)
- Clone all submodules (sundials)
- Build and install SUNDIALS from source
"""

submodules = ['submodules/sundials']

# Only for Windows
vs_path = 'C:/Program Files/Microsoft Visual Studio/2022/'  # Path to Visual Studio 2022

submodules = dict(
    # for building SUNDIALS:
    sundials = dict(
        dir = "./submodules/sundials/",
        build_dir = "./submodules/sundials/build/",
        install_dir = "./submodules/sundials/install/",
        src_dir = "./submodules/sundials/",
        cmake_command_list = [
            "-DCMAKE_BUILD_TYPE=RelWithDebInfo",
            "-DBUILD_ARKODE=ON",
            "-DBUILD_CVODE=ON",
            "-DBUILD_CVODES=OFF",
            "-DBUILD_IDA=OFF",
            "-DBUILD_IDAS=OFF",
            "-DBUILD_KINSOL=OFF",
            "-DBUILD_SUNLINSOL=ON",
            "-DBUILD_SUNMATRIX=ON",
            "-DBUILD_SHARED_LIBS=OFF",
            "-DBUILD_STATIC_LIBS=ON",
            "-DEXAMPLES_ENABLE_C=OFF",
            "-DEXAMPLES_ENABLE_CXX=OFF",
            "-DEXAMPLES_ENABLE_F2003=OFF",
            "-DEXAMPLES_INSTALL=OFF",   # requires '-DEXAMPLES_INSTALL_PATH=' + os.path.join(install_dir, "examples").replace("\\", "/")
            "-DENABLE_MPI=OFF",
            "-DENABLE_LAPACK=OFF",
            "-DENABLE_KLU=OFF",
            "-DENABLE_HYPRE=OFF",
            "-DENABLE_OPENMP=OFF",
            "-DENABLE_PTHREAD=OFF",
            "-DENABLE_SUPERLUMT=OFF",
            "-DENABLE_PETSC=OFF",
            #'-DKLU_INCLUDE_DIR=/usr/casc/sundials/apps/rh6/suitesparse/4.5.3/include',
            #'-DKLU_LIBRARY_DIR=/usr/casc/sundials/apps/rh6/suitesparse/4.5.3/lib',
            #'-DHYPRE_INCLUDE_DIR=/usr/casc/sundials/apps/rh6/hypre/2.11.1/include',
            #'-DHYPRE_LIBRARY=/usr/casc/sundials/apps/rh6/hypre/2.11.1/lib/libHYPRE.a',
            #'-DSUPERLUMT_INCLUDE_DIR=/usr/casc/sundials/apps/rh6/superlu_mt/SuperLU_MT_3.1/SRC',
            #'-DSUPERLUMT_LIBRARY_DIR=/usr/casc/sundials/apps/rh6/superlu_mt/SuperLU_MT_3.1/lib',
            #'-DSUPERLUMT_THREAD_TYPE=Pthread',
            #'-DPETSC_INCLUDE_DIR=/usr/casc/sundials/apps/rh6/petsc/3.7.2/include',
            #'-DPETSC_LIBRARY_DIR=/usr/casc/sundials/apps/rh6/petsc/3.7.2/lib'
        ]
    ),

    
)

import shutil
import os
import sys
import subprocess
from build_utility import red, light_red, green, yellow, blue, magenta, cyan, grey, reset, white, bold, italic


def main():
# Check dependencies
    check_tool_exists('git',        'sudo apt install git',         'https://git-scm.com/downloads/win')
    check_tool_exists('make',       'sudo apt install make',        'https://www.mingw-w64.org/downloads/')
    check_tool_exists('cmake',      'sudo apt install cmake',       'https://cmake.org/download/')
    check_tool_exists('gcc',        'sudo apt install gcc',         'https://www.mingw-w64.org/downloads/')
    check_tool_exists('g++',        'sudo apt install g++',         'https://www.mingw-w64.org/downloads/',         optinal=True)
    check_tool_exists('gdb',        'sudo apt install gdb',         'https://www.mingw-w64.org/downloads/',         optinal=True)
    check_tool_exists('clang',      'sudo apt install clang',       'https://releases.llvm.org/download.html',      optinal=True)
    check_tool_exists('clang++',    'sudo apt install clang',       'https://releases.llvm.org/download.html',      optinal=True)
    check_tool_exists('lldb',       'sudo apt install lldb',        'https://releases.llvm.org/download.html',      optinal=True)
    check_python_package('numpy')
    check_python_package('matplotlib')
    check_python_package('pygments')

# Check if we are in a Git repository
    if not os.path.exists('.git'):
        print(f'{bold}{red}Error:{reset} This is not a Git repository. Please clone the repository using Git instead of downloading it as a ZIP file.')
        return -1

# Clone submodules
    run_command('git submodule update --init --recursive')
    global submodules
    for name, settings in submodules.items():
        if not os.path.exists(settings['dir']):
            print(f'{bold}{red}Error:{reset} Submodule {name} not found at {settings["dir"]}. Please run "git submodule update --init --recursive"')
            return -1

# Build submodules
    ask_to_build('SUNDIALS', submodules['sundials'])

# end of main()


def run_command(command: str) -> bool:
    """Run a command in the shell. Returns True if the command was successful, False otherwise. 
    The command and its output are printed to the console neatly. """

    print(f'{bold}{blue}Running command:{reset} {command}')
    process = subprocess.Popen(command, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    try:
        process.stdin.write('y\n')
        process.stdin.flush()
    except:
        pass
    for line in process.stdout:
        print(f'    {line}', end='')
    for line in process.stderr:
        print(f'    {line}', end='')
    process.wait()
    if process.returncode != 0:
        print(f'    {bold}{red}Error:{reset} {italic}{command}{reset} failed with return code {process.returncode}')
        return False
    else:
        return True


apt_updated = False
def check_tool_exists(tool: str, install_linux: str, install_win: str, optinal: bool = False) -> bool:
    """Check if a tool is installed on the system. (e.g.: git, cmake, gcc...)"""

    if shutil.which(tool) is None:
        if optinal:
            print(f'{bold}{yellow}Warning:{reset} {tool} not found. Please install {tool}:')
        else:
            print(f'{bold}{red}Error:{reset} {tool} not found. Please install {tool}:')

        if os.name == 'posix' and install_linux is not None:
            print(f'    {install_linux}')
            print(f'    Would you like to install {tool}? [y/n]', end=' ')
            if input().lower() == 'y':
                global apt_updated
                if not apt_updated:
                    run_command('sudo apt update')
                    apt_updated = True
                run_command(install_linux)
        elif install_win is not None:
            print(f"    {install_win}")

        return False
    else:
        print(f'{tool}: {bold}{green}ok{reset}')
        return True


pip_updated = False
python_path = os.path.basename(sys.executable).split('.')[0]
def check_python_package(package: str) -> bool:
    """Check if a python package is installed on the system. (e.g.: numpy, matplotlib...)"""

    try:
        __import__(package)
        print(f'python.{package}: {bold}{green}ok{reset}')
        return True
    except ImportError:
        print(f'{bold}{red}Error:{reset} {package} not found. Please install {package}:')
        print(f'    pip install {package}')
        print(f'    Would you like to install {package}? [y/n]', end=' ')
        if input().lower() == 'y':
            global pip_updated
            if not pip_updated:
                check_tool_exists(python_path, '', '')
                check_tool_exists('pip', f'sudo apt install {python_path}-pip', f'Try command "{python_path} -m ensurepip"')
                run_command('pip install --upgrade pip')
                pip_updated = True
            run_command(f'pip install {package}')
        return False


def vs_build_command(vs_path: str, command: str) -> bool:
    """Initialize the Visual Studio build environment, to run commands like 'msbuild'. """

    if os.name == 'nt':
        if not os.path.exists(vs_path):
            print(f'{bold}{red}Error:{reset} Visual Studio 2022 not found at {vs_path}. Please install Visual Studio 2022, and update vs_path in this script.')
            return  False
        vs_dev_cmd_path = os.path.join(vs_path, 'Community/Common7/Tools/VsDevCmd.bat')
        if not os.path.exists(vs_dev_cmd_path):
            print(f'{bold}{red}Error:{reset} VsDevCmd.bat not found at {vs_dev_cmd_path}. You may run the later commands manually in the developer command prompt for VS 2022.')
            return False
        if not run_command(f'"{vs_dev_cmd_path}" && {command}'):
            return False
        return True


def build_with_cmake(name: str, src_dir: str, build_dir: str, install_dir: str, cmake_command_list: list, vs_path: str = vs_path) -> bool:
    """Build submodule (sundials, hdf5) from source with cmake. """

    # Check directories
    if not os.path.exists(src_dir):
        print(f'{bold}{red}Error:{reset} {name} source directory not found at {src_dir}.')
        return False
    if src_dir:
        src_dir = os.path.abspath(src_dir).replace('\\', '/')
    if build_dir:
        build_dir = os.path.abspath(build_dir).replace('\\', '/')
    if install_dir:
        install_dir = os.path.abspath(install_dir).replace('\\', '/')
    else:
        install_dir = build_dir
    print(f'    Source directory: {src_dir}')
    print(f'    Build directory: {build_dir}')
    print(f'    Install directory: {install_dir}')
    if os.path.exists(build_dir):
        shutil.rmtree(build_dir)
    os.makedirs(build_dir)
    if os.path.exists(install_dir):
        shutil.rmtree(install_dir)
    os.makedirs(install_dir)
    os.chdir(build_dir)

    # Run CMake
    generator = 'Unix Makefiles' if os.name == 'posix' else 'Visual Studio 17 2022'
    cmake_command_list = [
        'cmake',
        '-DCMAKE_INSTALL_PREFIX=' + install_dir,
        f'-G "{generator}"',
        '-A "x64"' if os.name == 'nt' else '',
    ] + cmake_command_list + [src_dir]
    command = ' '.join(cmake_command_list)
    if not run_command(command):
        return False
    
    # Build and install
    if os.name == 'posix':  # Linux
        if not run_command('make'):
            return False
        if not run_command('make install'):
            return False
    elif os.name == 'nt':   # Windows
        vs_build = os.path.join(build_dir, 'ALL_BUILD.vcxproj').replace('\\', '/')
        vs_install = os.path.join(build_dir, 'INSTALL.vcxproj').replace('\\', '/')
        if not os.path.exists(vs_build):
            print(f'{bold}{red}Error:{reset} ALL_BUILD.vcxproj not found at {vs_build}.')
            return False
        if not os.path.exists(vs_install):
            print(f'{bold}{red}Error:{reset} INSTALL.vcxproj not found at {vs_install}.')
            return False
        if not vs_build_command(vs_path, f'msbuild {vs_build}'):
            return False
        if not vs_build_command(vs_path, f'msbuild {vs_install}'):
            return False
        
    return True


def ask_to_build(name: str, settings: dict):
    print(f'{bold}Do you want to build {name} from source? [y/n]{reset}', end=' ')
    if input().strip().lower() == 'y':
        if not build_with_cmake(
            name,
            settings.get('src_dir', ''),
            settings.get('build_dir', ''),
            settings.get('install_dir', ''),
            settings.get('cmake_command_list', []),
            vs_path
        ):
            print(f'{bold}{red}Error:{reset} Building {name} failed.')
        else:
            print(f'{bold}{green}Success:{reset} {name} built successfully.')
    else:
        print(f'skipping building {name}')

if __name__ == "__main__":
    main()