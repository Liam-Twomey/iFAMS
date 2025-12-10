# iFAMS

To download a stable version of iFAMS, go to https://github.com/prellgroup/iFAMS/releases/

- Documentation, including tutorials, can be found in the `documentation` directory.
- The source code can be found in the `source/` directory.

## Installing source code
There are three ways to install iFAMS source code as a python tool, for instance on a non-Windows machine or to edit the program. pipx or uv are recommended, because they add the tool to $PATH, allowing it to be called as a command.

Note: these changes have not been merged to the parent repo at [https://github.com//prellgroup/iFAMS], so this will not work there yet. The source code here is unmodified from that repo, just uncompressed and reorganized for package management.

### pip installation

```
git clone https://github.com/Liam-Twomey/iFAMS
cd iFAMS
python3 -m venv venv
source venv/bin/activate
pip install build
python3 -m build .
pip install -e.
```

iFAMS can then be run by `python3 source/iFAMS_v6.3_Quant_GUI.py`

### pipx installation
`pipx install git+https://github.com/Liam-Twomey/iFAMS` 

### uv installation
Alternately, you can use the `uv` package manager ([homepage](https://docs.astral.sh/uv/)), and just `uv tool install git+https://github.com/Liam-Twomey/iFAMS`
