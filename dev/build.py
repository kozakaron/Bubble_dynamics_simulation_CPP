"""
This is a simple build script that uses the clang++ (or g++) compiler to compile the C++ program. Call from the command line as follows: ´python ./dev/build.py´
You may have to use python3 instead. To see available options, use the -h flag: ´python ./dev/build.py -h´
Note, that this script is multithreaded. THe /bin directory is created if it does not exist and all files in it are deleted.

Features:
- Highlighted output for errors, warnings, and notes
- Syntax highlighting for C++ code in error messages
- Multithreaded compilation
"""

import argparse
import os
import subprocess
import glob
import time
import re
from concurrent.futures import ThreadPoolExecutor


# Define the source directory and output binary
src_dirs = ['./src', './test']
build_dir = './bin'
include_dirs = ['./include', './test', './src', './mechanism']
output_binary = 'main'


# Other settings
compiler = 'clang++'
compiler_flags = [
    '-std=c++17',              # Use the C++17 standard
    '-march=native',           # Enable all instruction set extensions supported by the CPU
    '-ftree-vectorize',        # Enable tree vectorization
    #'-fvectorize',             # Enable auto-vectorization (only for clang)
    '-pedantic',               # Enforce strict standard compliance
    #'-funroll-loops',          # Enable loop unrolling
    '-mlong-double-80',        # Use 80-bit long double. You may also try with -128, however it is less likely to work.
# Error handling flags: slight performance boost, harder error identification, might cause unexpected behavior
    #'-fno-exceptions',         # Disable exception handling
    #'-fno-math-errno',         # Do not set errno after math library functions
    #'-ffinite-math-only',      # Assume that floating-point arithmetic arguments and results are finite numbers   
]
linker_flags = [
    #'-flto',                  # Enable link-time optimization
    #'-Wl,-O2',                # Linker optimization
]
common_flags = [
    #'-v',                     # Verbose output
# Sanitizers: help identify issues in runtime, performance overhead
    '-fsanitize=address',     # Protect against memory errors (a.g. use after free)
    '-fsanitize=undefined',   # Protect against undefined behavior (e.g. integer overflow, invalid type casts)
    '-fno-omit-frame-pointer',# Keep frame pointer for better stack traces
    #'-fomit-frame-pointer',    # Omit frame pointer for functions that don't need one
]



# colors for syntax highlighting
red = '\033[91m'
light_red = '\033[31m'
green = '\033[92m'
yellow = '\033[93m'
blue = '\033[94m'
magenta = '\033[95m'
cyan = '\033[96m'
grey = '\033[90m'
reset = '\033[0m'
white = '\033[97m'
bold = '\033[1m'
italic = '\033[3m'
# dict to save compile time
compile_times = {}


def syntax_highlight_cpp(code_line):
    """Apply syntax highlighting to a single line of C++ code."""

    # Define regex patterns for different C++ components
    keywords = r'\b(break|case|catc|continue|default|delete|do|else|explicit|extern|false|for|friend|goto|if|inline|namespace|new|operator|private|protected|public|register|return|sizeof|static|switch|template|this|throw|true|try|typedef|typeid|typename|using|virtual|volatile|while)\b'
    types = r'\b(int|float|double|class|enum|struct|union|char|void|bool|string|short|long|signed|unsigned|const|constexpr|auto|long)\b'
    numbers = r'\b\d+\b'
    strings = r'\".*?\"|\'.*?\''
    comments = r'\/\/.*?$|\/\*.*?\*\/'
    symbols = r'[\[\]\{\}\(\)\.\,\;\:\+\-\*\/\&\!\=\%\~\^\>\<\?\#]'
    macros = r'\#\w+'

     # Create a new string of equal length with different characters for each component
    highlighted_line = list(' ' * len(code_line))

    def replace_with_char(match, char):
        start, end = match.span()
        for i in range(start, end):
            highlighted_line[i] = char

    # Apply syntax highlighting
    for pattern, char in [(symbols, 's'), (keywords, 'k'), (types, 't'), (numbers, 'n'), (strings, 'q'), (comments, 'c'), (macros, 'm')]:
        for match in re.finditer(pattern, code_line):
            replace_with_char(match, char)

    # Define ANSI escape codes for different components
    colors = {
        'k': bold + cyan,
        't': cyan,
        'n': magenta,
        'q': light_red,
        'c': green,
        's': grey,
        'm': bold + magenta,
        ' ': reset
    }

    highlighted_code = ''
    for char, highlight in zip(code_line, highlighted_line):
        highlighted_code += colors.get(highlight, reset) + char + reset
    highlighted_code += reset  # Reset color at the end

    return highlighted_code


def syntax_coloring(text, indent=True, linker=False):
    """Color the text with syntax highlighting."""
    
    ret = ''
    skip_next_line = False

    for line in text.split('\n'):
        if skip_next_line:
            skip_next_line = False
        elif 'errors generated.' in line or 'error generated.' in line:
            line = bold + line + reset
        elif '~~~~~^~~~~' in line or '^' in line:
            if line[:line.find('^')].isspace() or line[:line.find('^')].strip() == '':
                line = bold + line + reset
            if line[:line.find('~~~~~^~~~~')].isspace() or line[:line.find('~~~~~^~~~~')].strip() == '':
                line = bold + line + reset
                skip_next_line = True
        elif any(keyword in line for keyword in [' error: ', ' warning: ', ' note: ']):
            error_type = next(keyword for keyword in [' error: ', ' warning: ', ' note: '] if keyword in line)
            color = light_red if ' error: ' in line else yellow
            sections = line.split(error_type)
            line = re.sub(r"'(.*?)'", lambda m: ('\x1B[3m' + m.group() + '\x1B[0m'), ''.join(sections[1:])) # italic
            line = bold + sections[0] + reset + color + error_type + reset + line
        elif 'In file ' in line:
            pass
        else:
            if not linker:
                line = syntax_highlight_cpp(line)


        if indent:
            line = '|-> ' + line if any(keyword in line for keyword in ['note: ', ' error: ', ' warning: ']) else '|   ' + line
        ret += line + '\n'

    return ret[:-1]


def compile_source_file(args):
    """Compile a single source file into an object file."""

    global compiler
    global compile_times
    start = time.time()
    cpp_file, include_dirs, compiler_flags = args
    obj_file = os.path.join(build_dir, os.path.basename(cpp_file).replace('.cpp', '.o'))
    command_list = [compiler, '-c'] + include_dirs + compiler_flags + [cpp_file, '-o', obj_file]
    result = subprocess.run(command_list, capture_output=True, text=True)
    if result.returncode != 0:
        print(light_red + f'Error compiling ' + bold + f'{cpp_file}:\n' + reset + syntax_coloring(result.stderr))
    end = time.time()
    compile_times[cpp_file] = end - start
    return result.returncode



def link_object_files(obj_files, output_binary):
    """Link object files into a final binary."""

    global linker_flags
    global compiler
    global compile_times
    start = time.time()
    command_list = [compiler] + obj_files + ['-o', output_binary] + linker_flags
    result = subprocess.run(command_list, capture_output=True, text=True)
    if result.returncode != 0:
        print(red + f'Error linking:\n' + reset + syntax_coloring(result.stderr, linker=True))
    end = time.time()
    compile_times["linking"] = end - start
    return result.returncode



def main():
    parser = argparse.ArgumentParser(description="Build script")
    parser.add_argument('-t', '--test', action='store_true', help='build with test flag (TEST is defined)')
    parser.add_argument('-b', '--benchmark', action='store_true', help='build with benchmark flag (BENCHMARK is defined)')
    parser.add_argument('-r', '--run', action='store_true', help='run the binary after building')
    parser.add_argument('-w', '--warning', action='store_true', help='enable all warnings')
    parser.add_argument('-d', '--debug', action='store_true', help='compile in debug mode (also runs gdb if called with -r)')
    parser.add_argument('-s', '--shared', action='store_true', help='build as shared library (.dll/.so)')
    parser.add_argument('-o2', '--optimize2', action='store_true', help='enable o2 optimization')
    parser.add_argument('-o3', '--optimize3', action='store_true', help='enable o3 optimization')
    args = parser.parse_args()

    # manage directories
    global build_dir
    global output_binary
    global include_dirs
    global compiler_flags
    global linker_flags
    global common_flags
    global compiler
    compiler_flags += common_flags
    linker_flags += common_flags
    os.makedirs(build_dir, exist_ok=True)
    for file in os.listdir(build_dir):
        if file != '.gitkeep':
            os.remove(os.path.join(build_dir, file))
    cpp_files = []
    for src_dir in src_dirs:
        cpp_files += glob.glob(os.path.join(src_dir, '*.cpp'))  # Get all cpp files in the src directory
    include_dirs = [f'-I{d}' for d in include_dirs]
    output_binary = os.path.join(build_dir, output_binary)
    if 'clang++' not in compiler and 'g++' not in compiler:
        print(red + bold + 'Unsupported compiler: ' + reset + compiler)
        exit(1)

    # Compiler flags
    if args.test:
        compiler_flags.append('-DTEST')
        if '-fno-exceptions' in compiler_flags:
            compiler_flags.remove('-fno-exceptions')
    if args.benchmark:
        compiler_flags.append('-DBENCHMARK')
    if args.warning:
        compiler_flags.extend(['-Wall', '-Wextra', '-Werror'])
    if args.debug:
        compiler_flags.append('-g')
    if args.shared:
        if os.name == 'nt':
            output_binary = output_binary + '.dll'
        else:
            output_binary = output_binary + '.so'
    if args.optimize2:
        compiler_flags.append('-O2')
    if args.optimize3 and not args.optimize2:
        compiler_flags.append('-O3')

    # Compile each source file in parallel using ThreadPoolExecutor
    start = time.time()
    pool_args = [(cpp_file, include_dirs, compiler_flags) for cpp_file in cpp_files]
    with ThreadPoolExecutor() as executor:
        results = list(executor.map(compile_source_file, pool_args))

    if any(r != 0 for r in results):
        print(red + bold + 'Compilation failed' + reset)
        exit(1)

    # Link object files into the final binary
    obj_files = [os.path.join('bin', os.path.basename(cpp_file).replace('.cpp', '.o')) for cpp_file in cpp_files]
    if link_object_files(obj_files, output_binary) != 0:
        print(red + bold + 'Linking failed' + reset)
        exit(1)

    end = time.time()
    compile_time = end - start
    print(green + 'Build succeeded in' + bold + f'{compile_time: .2f}' + reset + green + ' seconds' + reset)
    output_file = os.path.join(build_dir, 'compile_times.txt')
    with open(output_file, 'w') as f:
        f.write(f'Build time: {compile_time:5.2f} s\n')
        f.write('Compile times:\n')
        for cpp_file, compile_time in compile_times.items():
            cpp_file = cpp_file.replace('\\\\', '/').replace('\\', '/')
            f.write(f'{cpp_file: >40}: {compile_time:5.2f} s\n')

    if args.run and not args.shared:
        if args.debug:
            print(f'{bold}GDB (gnu debugger) help:{reset}')
            print(f'    {bold}run{reset}       start/restart the program (r)')
            print(f'    {bold}continue{reset}  continue execution (c)')
            print(f'    {bold}quit{reset}      quit gdb (q)')
            print(f'  {bold}Breakpoints:{reset}')
            print(f'    {bold}break{reset}     set a breakpoint (b): {italic}break function_name{reset} or {italic}break file.h:22{reset}')
            print(f'    {bold}info{reset}      show information (i): {italic}info locals{reset} or {italic}info breakpoints{reset}')
            print(f'    {bold}delete{reset}    delete breakpoint (d): {italic}delete 1{reset}')
            print(f'  {bold}Stepping:{reset}')
            print(f'    {bold}next{reset}      step over (n)')
            print(f'    {bold}step{reset}      step into (s)')
            print(f'    {bold}finish{reset}    finish current function (f)')
            print(f'  {bold}Code, variables:{reset}')
            print(f'    {bold}print{reset}     print variable (p): {italic}print variable_name{reset}')
            print(f'    {bold}backtrace{reset} show stack trace (bt)')
            print(f'    {bold}list{reset}      show source code (l): {italic}list{reset} or {italic}list 1, 10{reset}')
            result = subprocess.run(['gdb', output_binary])
            exit(0)

        result = subprocess.run([output_binary])
        if (result.returncode != 0):
            if result.stderr != '' and result.stderr != None: print(result.stderr)
            print(red + bold + output_binary.replace('\\', '/') + reset  + red + f' returned with code {result.returncode}' + reset)
            exit(1)
        else:
            print(green + bold + output_binary.replace('\\', '/') + reset  + green + f' returned with code {result.returncode}' + reset)
            exit(0)
                    

if __name__ == "__main__":
    main()
