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
@article{Liu:2019:SCG:3306346.3322953,
  author = {Liu, Hsueh-Ti Derek and Jacobson, Alec and Ovsjanikov, Maks},
  title = {Spectral Coarsening of Geometric Operators},
  journal = {ACM Trans. Graph.},
  issue_date = {July 2019},
  volume = {38},
  number = {4},
  month = jul,
  year = {2019},
  issn = {0730-0301},
  pages = {105:1--105:13},
  articleno = {105},
  numpages = {13},
  url = {http://doi.acm.org/10.1145/3306346.3322953},
  doi = {10.1145/3306346.3322953},
  acmid = {3322953},
  publisher = {ACM},
  address = {New York, NY, USA},
  keywords = {geometry processing, numerical coarsening, spectral geometry},
} 
```
