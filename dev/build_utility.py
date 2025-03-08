# colors for terminal output
red =       '\033[91m'
light_red = '\033[31m'
green =     '\033[92m'
yellow =    '\033[93m'
blue =      '\033[94m'
magenta =   '\033[95m'
cyan =      '\033[96m'
grey =      '\033[90m'
reset =     '\033[0m'
white =     '\033[97m'
bold =      '\033[1m'
italic =    '\033[3m'

debugger_help = f'''{bold}Debugger help:{reset}
    {bold}run{reset}       start/restart the program (r)
    {bold}continue{reset}  continue execution (c)
    {bold}quit{reset}      quit debugger (q)
  {bold}Breakpoints:{reset}
    {bold}break{reset}     set a breakpoint (b): {italic}break function_name{reset} or {italic}break file.h:22{reset}
    {bold}info{reset}      show information (i): {italic}info locals{reset} or {italic}info breakpoints{reset}
    {bold}delete{reset}    delete breakpoint (d): {italic}delete 1{reset}
  {bold}Stepping:{reset}
    {bold}next{reset}      step over (n)
    {bold}step{reset}      step into (s)
    {bold}finish{reset}    finish current function (f)
  {bold}Code, variables:{reset}
    {bold}print{reset}     print variable (p): {italic}print variable_name{reset}
    {bold}backtrace{reset} show stack trace (bt)
    {bold}list{reset}      show source code (l): {italic}list{reset} or {italic}list 1, 10{reset}'''



from pygments import highlight
from pygments.style import Style
from pygments.token import Token, Error, Other, Keyword, Literal, String, Text, Number, Operator, Punctuation, Comment, Generic
from pygments.lexers import CppLexer
from pygments.formatters import Terminal256Formatter
import re
import os

class Logger:
    def __init__(self, log_file: str):
        self.log_file = log_file
        self.log_str = ''
        _white =     '#fff'
        _black =     '#000'
        _grey =      '#aaa'
        _red =       '#f00'
        _green =     '#0a0'
        _blue =      '#00b3ff'
        _yellow =    '#ff0'
        _cyan =      '#0ff'
        _magenta =   '#f0f'
        _orange =    '#ffa600'
        class CppColorizer(Style):
            styles = {
                Token:                  '',
                Error:                  _red,
                Other:                  '',
                Keyword:                _magenta,
                Literal:                _green,
                String:                 _orange,
                Text:                   '',
                Number:                 _cyan,
                Operator:               _grey,
                Punctuation:            _grey,
                Comment:                _green,
                Generic:                '',
            }
        self.formatter = Terminal256Formatter(style=CppColorizer)
        self.lexer = CppLexer()


    def log_error(self, error: str, description: str = ''):
        global red, bold, reset
        if description == '' or description == None:
            print(red + bold + error + reset)
            self.log_str += f'ERROR: {error}\n'
        else:
            print(red + bold + error + reset + ': ' + description)
            self.log_str += f'ERROR: {error}: {description}\n'
    

    def log(self, message: str, formatted_message: str = ''):
        if formatted_message != '':
            print(formatted_message)
        self.log_str += message + '\n'


    def silent_log(self, message: str):
        self.log_str += message + '\n\n'


    def log_compiler_output(self, text: str, indent: bool = True, linker: bool = False):
        print(self._compiler_output_formatter(text, indent, linker))
        self.log_str += text + '\n'


    def write_to_file(self):
        if self.log_file is not None and self.log_file != '':
            with open(self.log_file, 'w') as file:
                file.write(self.log_str)


    def _syntax_highlight_cpp(self, code_line: str) -> str:
        return highlight(code_line, self.lexer, self.formatter)
    

    def _compiler_output_formatter(self, text: str, indent: bool = True, linker: bool = False):
        """Color the text output of compilers with syntax highlighting. Works for GCC and Clang formats."""
        
        colored_text = ''
        skip_next_line = False
        keywords = ['error: ', 'warning: ', 'note: ']

        for line in text.split('\n'):
            is_last_line = 'errors generated.' in line or 'error generated.' in line
            is_faulty_code_line = not linker and len(line.strip(' '))>0 and line.strip(' ')[0].isdigit()
            is_mistake_pointed_out = '~~~~~^~~~~' in line or '^' in line
            is_new_error = any(keyword in line for keyword in keywords)
            global bold, reset, light_red, yellow, italic

            if skip_next_line:
                skip_next_line = False
            elif is_last_line:
                line = bold + line + reset
            elif is_faulty_code_line:
                line = self._syntax_highlight_cpp(line).rstrip('\n')
            elif is_mistake_pointed_out:
                is_one_line = line[:line.find('^')].isspace() or line[:line.find('^')].strip() == ''
                is_two_line = line[:line.find('~~~~~^~~~~')].isspace() or line[:line.find('~~~~~^~~~~')].strip() == ''
                if is_one_line:
                    line = bold + line + reset
                if is_two_line:
                    line = bold + line + reset
                    skip_next_line = True
            elif is_new_error:
                error_type = next(keyword for keyword in keywords if keyword in line)
                color = light_red if ' error: ' in line else yellow
                sections = line.split(error_type)
                line = re.sub(r"'(.*?)'", lambda m: (italic + m.group() + reset), ''.join(sections[1:]))
                line = bold + sections[0].replace('\\', '/') + reset + color + error_type + reset + line
            elif 'In file ' in line:
                line = line.replace('\\', '/')
            else:
                pass

            if indent:
                line = '|-> ' + line if any(keyword in line for keyword in keywords) else '|   ' + line
            colored_text += line + '\n'

        return colored_text[:-1]



# Builder class
import subprocess
import shutil
import time
import os
import glob
from concurrent.futures import ThreadPoolExecutor

class Builder:
    def __init__(self, build_dir: str):
        self.build_dir = build_dir.replace('\\', '/')
        self.compile_times = dict()
        self.thread_pool_args = []
        self.output_binary = ''
        self.logger = Logger(os.path.join(build_dir, '.build.log'))
        if not os.path.exists(build_dir):
            os.mkdir(build_dir)
        else:
            self._clean_build_dir(build_dir)


    def _log_command(self, command_list: list[str]):
        command = os.getcwd().replace('\\', '/') + '> ' + ' '.join(command_list)
        self.logger.silent_log(command)


    def _clean_build_dir(self, build_dir: str):
        for root, dirs, files in os.walk(build_dir, topdown=False):
            for name in files:
                if name != '.gitkeep':
                    os.remove(os.path.join(root, name))
            for name in dirs:
                shutil.rmtree(os.path.join(root, name))


    def _check_tool_exists(self, tool: str) -> bool:
        if shutil.which(tool) is None:
            self.logger.log_error('tool not found', tool)
            return False
        else:
            return True
        

    def _compile_source_file(self, compiler: str, source_file: str, compiler_flags: list[str]) -> int:
        start = time.time()
        
        # checks
        if not os.path.exists(source_file):
            self.logger.log_error('source file not found', source_file)
            return -1
        if not self._check_tool_exists(compiler):
            return -1
        
        # compile
        base_name = os.path.basename(source_file)
        object_file = os.path.splitext(base_name)[0] + '.o'
        object_file = os.path.join(self.build_dir, object_file)
        command_list = [compiler, '-c'] + compiler_flags + [source_file, '-o', object_file]
        self._log_command(command_list)
        result = subprocess.run(command_list, capture_output=True, text=True)
        if result.returncode != 0:
            self.logger.log(message = 'error compiling ' + f'{source_file}:',
                            formatted_message = light_red + 'error compiling ' + bold + f'{source_file}:' + reset)
            self.logger.log_compiler_output(result.stderr, indent=True)

        # log compile time
        end = time.time()
        self.compile_times[source_file.replace('\\', '/')] = end - start
        return result.returncode
    

    def _compile_source_file_args(self, args: tuple[str, list[str], list[str]]) -> int:
        return self._compile_source_file(*args)


    def gather_files(self, src_dirs: list[str], include_dirs: list[str]) -> tuple[list[str], list[str], list[str]]:
        cpp_files = []
        cuda_files = []
        include_flags = []
        for src_dir in src_dirs:
            if not os.path.exists(src_dir):
                self.logger.log_error('source directory not found', src_dir)
            else:
                cpp_files += glob.glob(os.path.join(src_dir, '*.cpp').replace('\\', '/'))
                cuda_files += glob.glob(os.path.join(src_dir, '*.cu').replace('\\', '/'))
        for include_dir in include_dirs:
            if not os.path.exists(include_dir):
                self.logger.log_error('include directory not found', include_dir)
            else:
                include_flags.append('-I' + include_dir.replace("\\", "/"))

        return cpp_files, cuda_files, include_flags


    def add_source_file(self, source_file: str, compiler: str, compiler_flags: list[str]):
        self.thread_pool_args.append((compiler, source_file.replace('\\', '/'), compiler_flags))

    
    def compile_all(self) -> list[str]:
        self.start = time.time()
        object_files = []

        with ThreadPoolExecutor() as executor:
            results = list(executor.map(self._compile_source_file_args, self.thread_pool_args))
        if any(r != 0 for r in results):
            self.logger.log_error('compilation failed', '')
        else:
            for _, source_file, _ in self.thread_pool_args:
                base_name = os.path.basename(source_file)
                object_file = os.path.splitext(base_name)[0] + '.o'
                object_file = os.path.join(self.build_dir, object_file).replace('\\', '/')
                object_files.append(object_file)

        return object_files
    

    def link(self,
             linker: str,
             obj_files: list[str],
             lib_dirs: list[str],
             output_binary_name: str,
             linker_flags: list[str],
             shared: bool = False
            ) -> bool:
        start = time.time()

        #checks
        for obj_file in obj_files:
            if not os.path.exists(obj_file):
                self.logger.log_error('object file not found', obj_file)
                self.logger.write_to_file()
                return -1
        if not self._check_tool_exists(linker):
            self.logger.write_to_file()
            return -1
        if len(obj_files) == 0:
            self.logger.write_to_file()
            return -1
        
        # gather library directories
        for lib_dir in lib_dirs:
            if not os.path.exists(lib_dir):
                self.logger.log_error('library directory not found', lib_dir)
                self.logger.write_to_file()
                return -1
            linker_flags.append('-L' + os.path.abspath(lib_dir).replace('\\', '/'))
            
            core_lib = None
            for lib_file in os.listdir(lib_dir):
                if not lib_file.endswith('.a') and not lib_file.endswith('.so') and not lib_file.endswith('.dll'):
                    continue

                lib_file = os.path.join(lib_dir, lib_file)
                if 'core' in lib_file:
                    if core_lib is not None:
                        self.logger.log_error(f'multiple core libraries found: {core_lib}, {lib_file}', '')
                        self.logger.write_to_file()
                        return -1
                    core_lib = lib_file
                    continue
                
                linker_flags.append(lib_file)
            
            if core_lib is None:
                self.logger.log_error('core library not found', '')
                self.logger.write_to_file()
                return -1
            linker_flags.append(core_lib)            
            linker_flags.append('-lm')

        
        # determine output binary name
        output_binary_name = os.path.basename(output_binary_name).split('.')[0]
        if shared:
            output_binary = output_binary_name + ('.dll' if os.name == 'nt' else '.so')
        else:
            output_binary = output_binary_name + ('.exe' if os.name == 'nt' else '')
        output_binary = os.path.join(self.build_dir, output_binary).replace('\\', '/')
        self.output_binary = output_binary

        # link
        command_list = [linker] + obj_files + ['-o', output_binary] + linker_flags
        self._log_command(command_list)
        result = subprocess.run(command_list, capture_output=True, text=True)
        if result.returncode != 0:
            self.logger.log(message = 'Error linking:',
                            formatted_message = light_red + 'Error linking:' + reset)
            self.logger.log_compiler_output(result.stderr, indent=True, linker=True)
        end = time.time()

        # log compile times
        self.compile_times["linking"] = end - start
        self.compile_times["total"] = end - self.start
        compile_times_str = ''
        for cpp_file, compile_time in self.compile_times.items():
            compile_times_str += f'{cpp_file: >40}: {compile_time:5.2f} s\n'
        self.logger.log(message = 'bbject file compilation times: \n' + compile_times_str)
        if result.returncode != 0:
            self.logger.log_error('linking failed', '')
            self.logger.write_to_file()
            return False
        else:
            self.logger.log(message=f'build succeeded in {self.compile_times["total"]:.2f} seconds',
                            formatted_message=green + 'Build succeeded in ' + bold + f'{self.compile_times["total"]:.2f}' + reset + green + ' seconds' + reset)
            self.logger.write_to_file()
            return True
    

    def run_binary(self, debugger: str = '') -> int:
        if not os.path.exists(self.output_binary):
            self.logger.log_error('binary not found', self.output_binary)
            return -1
        
        if debugger == '' or debugger == None:
            result = subprocess.run([self.output_binary])
            if (result.returncode != 0):
                if result.stderr != '' and result.stderr != None: print(result.stderr)
                color = red
            else:
                color = green
            self.logger.log(message=f'program returned with code {result.returncode}',
                            formatted_message=color + bold + self.output_binary + ' program returned with code ' + str(result.returncode) + reset)
            return result.returncode
        else:
            if not self._check_tool_exists(debugger):
                return -1
            global debugger_help
            print(debugger_help)
            result = subprocess.run([debugger, self.output_binary])
            return result.returncode