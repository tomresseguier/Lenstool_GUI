# visualens

### A visual, fast interface to explore and manipulate astronomical images, catalogs, and lens models created using Lenstool.

`visualens` uses `PyQt` for fast handling of large astronomical images. Main lensing functionalities include:
- Fast, easy-to-run lensing calculations using `Lenstool`'s engine,
- Visual forward modeling based on `lenstronomy`.

---

### Installation:

Many `visualens` functions rely on `Lenstool`, which is easy to install through conda (highly recommended). We also recommend installing `lenstronomy` through conda rather than pip, before installing `visualens`.

Some IDEs like Spyder may cause conflicts with PyQt. We suggest installing your favorite IDE in your conda virtual environment before visualens.

**Installation steps:**

> `conda create -n visualens_env -c conda-forge lenstool lenstronomy python==3.12.2`

> `conda activate visualens_env`

Install your favorite python GUI

> `pip install visualens`

---

### Starting guide:

https://github.com/tomresseguier/Lenstool_GUI/tree/main/examples
