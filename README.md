![APM](https://badgen.net/github/license/micromatch/micromatch)

# sympiler_eigen
This repository provide an Eigen interface for 
[Sympiler](https://github.com/sympiler/sympiler) and 
[NASOQ](https://github.com/sympiler/nasoq) 
codes. 
We are adding interfaces one by one and this is a work in progress. 
The current status:
- [x] Sparse supernodal Cholesky factorization. 
- [ ] NASOQ QP Solver
- [ ] Sparse lower triangular solver 
- [ ] Sparse upper triangular solver
- [ ] Sparse incomplete Cholesky zero
- [ ] Sparse incomplete LU zero
- [ ] Sparse matrix-vector multiplication
- [ ] Preconditioned GMRES
- [ ] Sparse Gauss-Seidel
- [ ] Sparse non-supernodal Cholesky
- [ ] Sparse supernodal LDL factorization
- [ ] Sparse simplicial (non-supernodal) LDL factorization
- [ ] ...


## Installation

### pre-requisites
Sympiler needs METIS and MKL as dependecies. 
If METIS does not exist, CMake will handle it. 
If METIS exists and it is not handled with the CMake, 
you can add its libray and include paths to `-DCMAKE_PREFIX_PATH=`.

If MKL does not exist, some packages might not be installed. 

After installing dependencies, you should follow:
```bash
   git clone --recursive https://github.com/sympiler/sympiler-eigen.git
   cd sympiler-eigen
   mkdir build
   cd build
   cmake ..
   make 
```
You can follow the [sympiler installation instructions](https://github.com/sympiler/sympiler). 




