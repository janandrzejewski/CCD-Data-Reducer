# astro_reduction_tool

astro_reduction_tool is a command-line progam for astro data reduction.

## Requirements

- Python 3.10.8+
- Miniconda 22.11.1+
- Ccdproc 2.3.2+
- Astropy 5.2+
## Installation
install ccdproc using 
```commandline
conda install -c conda-forge ccdproc
```

install astropy using 
```commandline
conda install astropy
```

## Usage

### Data format

This program needs a specific data format.
(year_month_day/object_objectid)
and there must be subfolder with flats named 'flat'

### Reduction
To start use:
```commandline
python astro_reduction_tool.py -d main_folder_dir -md master_dark_dir -mb master_bias_dir
```

if you already made master_flats and you don't want to make new ones just add -mf dir to flat subfolder
## Development

Imports are organized with `isort`.

Code is formatted with `black`.

Master_flats are saved in flat folder.
