[← Documentation home](README.md) · [Main README](../README.md)

# Installation

## Fastest start: no terminal (Google Colab)

Colab runs Python in your web browser for free. You don't need to install anything.

1. Go to <https://colab.research.google.com>.
2. Click **File → Open notebook → GitHub**, paste `https://github.com/nothikemu/CFD`, and pick
   `examples/quickstart.ipynb` on branch `claude/sharp-carson-1l4xvs`.
   *(Alternatively, click **New notebook** and paste the cells from the box below.)*
3. Press **▶** on each grey box from top to bottom (or press **Shift + Enter**).
4. In the **SETTINGS** box, change the length, height, material or force, then press **▶** again.

> **Private repository?** Colab will ask you to sign in to GitHub. If the repo is public you won't be asked.

Manual Colab cells:

```python
# Cell 1: download + install
!git clone -b claude/sharp-carson-1l4xvs https://github.com/nothikemu/CFD.git
%cd CFD
!pip -q install -r requirements.txt
```

```python
# Cell 2: check everything works
!python -m pytest test_validation.py -q
```

```python
# Cell 3: run the benchmarks and show the pictures
!python app.py
from IPython.display import Image, display
display(Image("output/cantilever_von_mises.png"))
```

## Installing on your own computer

### You need Python 3.9 or newer

Check by opening a terminal (**Mac:** Terminal app · **Windows:** PowerShell · **Linux:** any terminal):

```bash
python3 --version        # Mac / Linux
python --version         # Windows
```

If it's missing or older than 3.9, install it from <https://www.python.org/downloads/>. On Windows, tick **"Add python.exe to PATH"** in the installer.

> Throughout this guide, **Windows users type `python` wherever it says `python3`**.

### Download the code

**With git:**
```bash
git clone -b claude/sharp-carson-1l4xvs https://github.com/nothikemu/CFD.git
cd CFD
```

**Without git:** on the GitHub page, switch the branch dropdown to `claude/sharp-carson-1l4xvs` and click **Code → Download ZIP**. Unzip it, then open a terminal inside the unzipped folder.

### (Recommended) make a private Python environment

This keeps the solver's packages separate from everything else on your computer.

```bash
# Mac / Linux
python3 -m venv .venv
source .venv/bin/activate
```
```powershell
# Windows PowerShell
python -m venv .venv
.venv\Scripts\Activate.ps1
```

If Windows refuses with "running scripts is disabled", run this once, then activate again:
`Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`

Run the activate line again every time you open a new terminal. To leave the environment, type `deactivate`.

### Install the packages

```bash
python3 -m pip install -r requirements.txt
```

This installs `numpy`, `scipy`, `matplotlib` and `pytest`. Nothing else is needed.

### Check it works

```bash
python3 -m pytest test_validation.py -v
```

You should see **`24 passed`**.

### Using an editor instead of a terminal (VS Code)

1. Install [VS Code](https://code.visualstudio.com/) and its **Python** extension.
2. **File → Open Folder…** and choose the `CFD` folder.
3. Pick your interpreter: press `Ctrl/Cmd + Shift + P`, run **Python: Select Interpreter**, and choose `.venv` if you made one.
4. Open any `.py` file (for example `examples/01_my_beam.py`) and click **▶ Run Python File** in the top-right corner.
5. Pictures appear in the `output/` folder in the sidebar. Click one to view it.

PyCharm, Spyder and Thonny work the same way: open the folder, open a file, press Run.

---

**Next:** run your first analysis with the [Tutorial](tutorial.md).
