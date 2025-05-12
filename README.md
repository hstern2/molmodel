# molmodel

command-line programs for molecular modeling

## Build Instructions

### macOS (Darwin)

requires Xcode command‑line tools (`clang`, `cmake`, etc.)

```bash
brew update
brew install fftw --build-from-source
brew install gsl
mkdir build && cd build
cmake .. && make
```

### Linux (Ubuntu)

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake pkg-config libopenblas-dev libfftw3-dev libgsl-dev
mkdir build && cd build
cmake .. && make
```

## Executables

* **geo**: geometry calculations
* **msim**: molecular dynamics simulations
* **rmsd**: root-mean-square deviation between two structures
* **hist**: histograms
