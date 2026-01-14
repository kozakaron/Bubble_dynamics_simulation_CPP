# How to add a new reaction mechanism?

It is possible to quickly adopt a new reaction mechanism without building the C++ code again. Reaction mechanisms can be provided in two ways:
 * **INP file**: the Chemkin II input file format, denoted as `.inp`. This is a legacy format from the Fortran era, using the CGS system often with erg or cal as energy units.
* **YAML file**: the Cantera YAML file format, denoted as `.yaml` or `.yml`. This is the modern format used by Cantera, using the SI system with J as energy unit.

This code base uses a JSON based format, extended with comments. You need to convert your mechanism to this format using the provided Python scripts. Why a third format? Cantera's YAML parser is actually available in C++, and I might make the change in the future. I did not want to suffer with another C++ dependency, and Cantera's Python interface is super easy to install and use. Furthermore, the reaction mechanisms are heavily preprocessed, and these steps are visible in the intermediate JSON files.

## From INP to YAML

Skip these steps if you already have your mechanism in YAML format.

 1. Copy your INP file to `Bubble_dynamics_simulation/mechanism/inp_files`. You may need to combine your thermodynamics and reaction mechanism files into a single INP file first.
 2. Make sure, that your current working directory is `Bubble_dynamics_simulation_CPP`. The terminal should show something like `.../Bubble_dynamics_simulation_CPP>`.
 3. Run the following command in your terminal: `python ./mechanism/convert_inp_to_yaml.py`. You might need to use python3.
 4. If the cantera module is not found, install Cantera using `pip install cantera`.
 5. All INP files in the `inp_files` folder will be converted to YAML files and stored in the `yaml_files` folder. If any error arises during converstaion, check the detailed cantera error messages, and modify your INP file accordingly.

## From YAML to JSON

 1. Copy your YAML file to `Bubble_dynamics_simulation/mechanism/yaml_files`.
 2. Make sure, that your current working directory is `Bubble_dynamics_simulation_CPP`. The terminal should show something like `.../Bubble_dynamics_simulation_CPP>`.
 3. Run the following command in your terminal: `python ./mechanism/convert_yaml_to_json.py`. You might need to use python3.
 4. If the cantera module is not found, install Cantera using `pip install cantera`.
 5. You will reference your mechanisms by the formatless filename of the resulting JSON. You may rename it to something simple and meaningful.