from cantera import ck2yaml
import os

inp_dir = './mechanism/inp_files'
yaml_dir = './mechanism/yaml_files'
reset =     '\033[0m'
bold =      '\033[1m'

if __name__ == '__main__':
    inp_files = [file for file in os.listdir(inp_dir) if file.lower().endswith('.inp')]
    for inp_file in inp_files:
        inp_path = os.path.join(inp_dir, inp_file).replace('\\', '/').replace('//', '/')
        yaml_file = os.path.splitext(inp_file)[0] + '.yaml'
        yaml_path = os.path.join(yaml_dir, yaml_file).replace('\\', '/').replace('//', '/')
        print('_' * len(f'Converting {inp_file}'))
        print(f'{bold}Converting {inp_file}{reset}')

        ck2yaml.convert(
            input_file=inp_path,
            thermo_file=None,       # incuded in input_file
            transport_file=None,
            out_name=yaml_path,
            permissive=True,
            verbose=True
        )

        print('\n')