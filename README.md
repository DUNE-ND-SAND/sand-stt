Requirements
- [ROOT](https://root.cern/)
- [edep-sim](https://github.com/ClarkMcGrew/edep-sim)

Build

```
$ mkdir build
$ cd build
$ cmake -DCMAKE_INSTALL_PREFIX=<install-dir> <source-dir>
$ make
$ make install
```

Before running application
- To have dictionaties of the structs loaded at run time:
```
$ source setup.sh
```

Using with ROOT
- To load dictionaries of structs:
```
root [0] gSystem->Load("libStruct.so")
```

Digitization
- Create digits of STT and cells of calorimeter

```
$ Digitize <MC file> <reco file>
```

Reconstruction
- Track find and fit of STT track
- Clustering of calorimeter cells

```
$ Reconstruct <MC file> <reco file>
```

Analysis
- Evaluate parameters of particles
- Evaluate neutrino energy

```
Analyze <MC file> <reco file>
```

The description of the data format can be found [here](https://github.com/DUNE-ND-SAND/sand-stt/wiki/Data-Model)

The code format can be find [here](https://github.com/DUNE-ND-SAND/sand-stt/wiki/Code-Formatting)
