# FLiTs: Fast Line Tracer for Protoplanetary Disks

A line ray-tracer build to compute the atomic and molecular lines emerging from protoplanetary disks. It is made to work together with [ProDiMo](https://prodimo.iwf.oeaw.ac.at/). It can compute quickly many molecular lines in the infrared taking into account the rotation of the disk and possible line-blending effects. FLiTs is available on a collaborative basis.

## Disclaimer

By using FLiTs you agree to the terms of use that can be found in the [User Guide](https://github.com/michielmin/FLiTs/blob/master/doc/UserGuide.pdf). It basically means you offer **Michiel Min** co-author rights on any paper that uses results computed with FLiTs.

## Citation

If you use FLiTs for you research, please cite [Woitke P., Min M. et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A&A...618A..74W/abstract).

## Installation

Below you find instructions on how to install ``FLiTs``. The easiest is to use an existing conda environment for ``ProDiMo``, but you can also create a separate conda environment. The installation should work for Linux (gfortran, ifx) and Mac OS X (gfortran). Please use recent versions of the compilers. The following instructions are tested with gfortran 15.2 and ifx 2025.3.2.

### Existing conda environment for ProDiMo

If you used `conda` to install and setup `ProDiMo`, all you need to do is to activate your ProDiMo conda environment, clone this repository (see the green Code button) into a directory of your choice and do the following steps:

```bash
cd TO_YOURFLITS_DIRECTORY
cp Makefile.conda Makefile
make
make install 
conda activate YOUR_PRODIMO_ENVIRONMENT
FLiTs
```

If you use `ifx` as a Fortran compiler, use the `Makefile.conda.ifx` instead of `Makefile.conda`. It is important to re-activate your conda environment after the `make install` step, otherwise the changes to the environment variables are not applied.

The last command (`FLiTs`) is just to check if the compilation works. You should see some welcome message and some information on the compiler etc., and an error message at the end. That means the compilation was successful, the error message is expected as no input file was provided.

Whenever you update FLiTs, you need to compile the code with `make` again. It is best to use `make clean` first. The `make install` step does not need to be repeated.

### Separate conda environment

If you prefer to have a separate conda environment for FLiTs, you can create one via

```bash
conda create -n flits gfortran cfitsio git
conda activate flits
```

for gfortran. For ifx you can use:

```bash
conda create -n flits.ifx python cfitsio git ifx_linux-64 mkl 
conda activate flits.ifx
```

Then follow the same steps as above (git clone, make etc.).

### Manual installation

You can also install the dependencies manually, you need to have a Fortran compiler, Git and the cfitsio library. After you cloned the repository, you can start with the example `Makefile.FLiTs`; i.e. copy it to `Makefile` and edit it according to your system (e.g. paths etc.).
