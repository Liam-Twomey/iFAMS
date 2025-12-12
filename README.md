# iFAMS
iFAMS (Interactive Fourier-Transform Analysis for Mass Spectrometry) is a tool for white noise filtering and charge state deconvolution for polydisperse ESI-MS data.

## About this fork
This repo is a fork of iFAMS which has been slighly reorganized to allow easy installation with a package manager on any platform. The main/original repo of the software can be found at [https://github.com/prellgroup/iFAMS/releases/]

This was forked from version 6.3; the changes relative to that version are as follows:
- Decompress zip files of source code
- Reorganize source code and documentation files to play nice with the python packaging standard.
- Make import statements in code use relative paths, again to play nice with the python packaging standard.
- Add `pyproject.toml` and `source/__init__.py` to define dependencies, package build, and program initialization.
- Rename GUI file to work with packaging (periods before the file extension break the package build)

Other than import statements, the source code is unmodified.

## Organization of this repo

- Documentation, including tutorials, can be found in the `documentation/` directory.
- The source code can be found in the `source/` directory.
- Package definition and dependencies are in `pyproject.toml`.

## Installing source code
There are three ways to install iFAMS source code as a python tool, for instance on a non-Windows machine or to edit the program. `pipx` or `uv` are *strongly* recommended for general use, because they automate installation of the package in its own environment and addition to the system path. This makes the package *substantially* more stable and easy to use. The pip installation should only be used for development or modification of the code.

***Note:*** The files in this fork have been uncompressed and reorganized relative to the parent repo at [https://github.com//prellgroup/iFAMS], so no information in this readme is relevant to the upstream. 

### pipx installation
`pipx` installs the python app so it can be interacted with like a system program.
After [installing pipx](https://pipx.pypa.io/latest/installation/), iFAMS can be installed with:
```
pipx install git+https://github.com/Liam-Twomey/iFAMS
```
If the installation succeeds, the `ifams` command should now be available in the command line.

### uv installation
Alternately, you can use the faster, more modern, but less standard [uv](https://docs.astral.sh/uv/) package manager.
After [installing uv](https://docs.astral.sh/uv/getting-started/installation/) iFAMS can be installed with `uv tool install git+https://github.com/Liam-Twomey/iFAMS`. The `ifams` command should be now available in the command line. Alternately, if you want to try the program without installing it, `uvx git+https://github.com/Liam-Twomey/iFAMS`.

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
