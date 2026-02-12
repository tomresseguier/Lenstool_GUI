# Lenstool GUI

A visual interface to explore astronomical images, catalogs, and lens models created using Lenstool.



Installation:

Many visualens functions rely on Lenstool, which is easy to install through conda (highly recommended).
We also recommend installing lenstronomy through conda rather than pip, because of common issues with the llvmlite (required by lenstronomy) install through pip(1), although this should not be a problem if using an older version of python (<=3.11).
If you use a certain IDE like Spyder and want to install it in your conda virtual environment, install it before visualens, as the other way around might result in conflicts between their PyQt requirements. The same may apply to other GUI, so make sure to install these first.

>> conda create -n visualens_env -c conda-forge lenstool lenstronomy python==3.12.2
>> conda activate visualens_env
## Install your favorite python GUI
>> pip install visualens


Notes:
(1) Installing lenstronomy through pip, or not installing it prior to installing visualens, may result in building failure due to a lack of pip installable wheels for llvmlite
