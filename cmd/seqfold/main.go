package main

import (
	"flag"
	"fmt"
	"os"

	"gitlab.com/folago/seqfold.go"
)

const usage = `usage: seqfold [-h] [-t FLOAT] [-d] [-r] SEQ

Predict the minimum free energy (kcal/mol) of a nucleic acid sequence

positional arguments:
  SEQ                   nucleic acid sequence to fold

optional arguments:
  -h, --help            show this help message and exit
  -t FLOAT, --celsius FLOAT
                        temperature in Celsius
  -d, --dot-bracket     write a dot-bracket of the MFE folding to stdout
  -r, --sub-structures  write each substructure of the MFE folding to stdout
`

// this was the python format string
// fmt = "{:>4} {:>4} {:>6}  {:<15}"
const alignFormatStr = "%4v %4v %6v %-15s\n"

func main() {
	var (
		seq    string
		temp   float64
		dot    bool
		substr bool
	)
	flag.BoolVar(&dot, "dot-bracket", false, "write a dot-bracket of the MFE folding to stdout")
	flag.BoolVar(&dot, "d", false, "write a dot-bracket of the MFE folding to stdout")

	flag.BoolVar(&substr, "sub-structure", false, "write each substructure of the MFE folding to stdout")
	flag.BoolVar(&substr, "r", false, "write each substructure of the MFE folding to stdout")

	flag.Float64Var(&temp, "celsius", 37.0, "temperature in Celsius")
	flag.Float64Var(&temp, "t", 37.0, "temperature in Celsius")

	flag.Usage = func() { fmt.Print(usage) }

	flag.Parse()

	seq = flag.Arg(0)

	if seq == "" {
		fmt.Print(usage)
		os.Exit(1)
	}
	structs := seqfold.Fold(seq, temp)
	if dot {
		fmt.Println(seq)
		fmt.Println(seqfold.DotBracket(structs))
	}
	if substr {
		fmt.Printf(alignFormatStr, "i", "j", "ddg", "description")
		for _, s := range structs {
			i, j, e, d := s.MultiString()
			fmt.Printf(alignFormatStr, i, j, e, d)
		}
	}

	energy := 0.0
	for _, s := range structs {
		energy += s.E
	}
	fmt.Println(seqfold.RoundFloat(energy, 2))
}
