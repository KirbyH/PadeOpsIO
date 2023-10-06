# PadeOpsIO Overview

PadeOpsIO is a tool for loading, visualizing, and analyzing data from [PadeOps](https://github.com/FPAL-Stanford-University/PadeOps), an LES and DNS box code. 

## Dependencies

Aside from standard python library packages, PadeOpsIO requires: <br>
* `numpy` <br>
* `scipy` <br>
* `matplotlib` <br>

The `pyproject.toml` file lists the requirement Python >3.7, but if this causes issues it can be changed. `numpy` and `scipy` don't play nicely for some reason (to be resolved...) so they are left out of the `pyproject.toml` file. 

## Installation

Clone the repo from GitHub. This can be done with: 
```
pip install git+https://github.com/kirbyh/PadeOpsIO.git
```

Alternatively, for developers, clone the git repository locally and install it with editing enabled: 
```
# navigate to the installation directory
cd $install_dir

# clone the repo
git clone -b main https://github.com/kirbyh/PadeOpsIO.git

# then add with pip
pip install . -e
```

The module can be imported into a Python script with: 
```
import padeopsio as pio
```

## Usage

PadeOpsIO is used for visualization of and analysis for output data from PadeOps. For more, see the [quick start](https://github.com/kirbyh/PadeOpsIO/blob/main/padeopsio/padeopsio_quickstart.ipynb). 

Data can be instanteous data: 

![image](https://user-images.githubusercontent.com/8905274/197601106-86fd32e4-52dc-4cf5-bcc3-d1bd664cdc08.png)

or time-averaged: 

![image](https://user-images.githubusercontent.com/8905274/197600994-47325c6d-89f3-4d09-9a44-1a0822fe81b5.png)
