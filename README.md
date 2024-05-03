# seqfold.go

This is an incomplete Go port of the python [seqfold](https://github.com/Lattice-Automation/seqfold) project.

## Install

You need Go version 1.18 or later installed in your system.

```
go install gitlab.com/folago/seqfold.go/cmd/seqfold@latest
```

## Usage

```
usage: seqfold [-h] [-t FLOAT] [-d] [-r] SEQ

Predict the minimum free energy (kcal/mol) of a nucleic acid sequence

positional arguments:
  SEQ                   nucleic acid sequence to fold

optional arguments:
  -h, --help            show this help message and exit
  -t FLOAT, --celsius FLOAT
                        temperature in Celsius
  -d, --dot-bracket     write a dot-bracket of the MFE folding to stdout
  -r, --sub-structures  write each substructure of the MFE folding to stdout
```

## Contibuting 

This project has been included in https://github.com/bebop/poly.
Please contribute to that project.

