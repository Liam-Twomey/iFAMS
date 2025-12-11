# iFAMS

To download a stable version of iFAMS, go to https://github.com/prellgroup/iFAMS/releases/

- Documentation, including tutorials, can be found in the `documentation` directory.
- The source code can be found in the `source/` directory.

## Installing source code
There are three ways to install iFAMS source code as a python tool, for instance on a non-Windows machine or to edit the program. pipx or uv are recommended, because they add the tool to $PATH, allowing it to be called as a command.

Note: The files in this fork have been uncompressed and reorganized relative to the parent repo at [https://github.com//prellgroup/iFAMS], so no information in this readme is relevant to the upstream. 

### pipx installation
`pipx install git+https://github.com/Liam-Twomey/iFAMS` 
If the installation suceeds, the `ifams` command should now be available in the command line.

### uv installation
Alternately, you can use the faster, more modern, but less standard `uv` package manager ([homepage](https://docs.astral.sh/uv/)), and just `uv tool install git+https://github.com/Liam-Twomey/iFAMS`
After installation of the tool, the `ifams   command should be available in the command line.

### pip installation
```
git clone https://github.com/Liam-Twomey/iFAMS
cd iFAMS
python3 -m venv venv
source venv/bin/activate
pip install build
python3 -m build .
pip install -e .
```

iFAMS can then be run by `python3 source/iFAMS_GUI.py`
