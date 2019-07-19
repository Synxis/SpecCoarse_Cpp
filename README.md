# Spectral Coarsening
This is the C++ implementation of "Spectral Coarsening of Geometric Operators" [Liu et al. 2019]. The dependencies are libigl ```https://libigl.github.io``` and Spectra ```https://spectralib.org```, where the Spectra has been added to the `./external/` folder. 

### Installation
Our code is developed on MacOS, we provide the commands for installing SpecCoarse in MacOS: 

Assuming your libigl is installed and compiled in `/usr/local/libigl/`, the first step is to clone the SpecCoarse_Cpp repository to `/usr/local/`
```
cd /usr/local
git clone https://github.com/HTDerekLiu/SpecCoarse_Cpp.git
```

Now assume you are in the directory `/usr/local/SpecCoarse_Cpp/`, then you can compile the code via cmake
```
mkdir build
cd build
cmake ..
make
```

Once the installation is completed, you can run the code by
```
./specCoarsen_bin [path/to/XXX.obj] [m] [k] [lr] [lrReduce]
```


### bibtex
```
@article{Liu:SpecCoarse:2019,
  title = {Spectral Coarsening of Geometric Operators},
  author = {Hsueh-Ti Derek Liu and Alec Jacobson and Maks Ovsjanikov},
  year = {2019},
  journal = {ACM Transactions on Graphics}, 
}
```
