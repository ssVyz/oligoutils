package oligoutils

import "bufio"
import "os"
import "io"
import "fmt"
import "strings"


type Seqr struct {
	Header string
	Seq string
}


const Iupac = "ATGCRYVBMHDSKN"

var matches = map[string][]string{
	"A": {"T"},
	"T": {"A"},
	"G": {"C"},
	"C": {"G"},
	"U": {"A"},

	"R": {"T", "C"},
	"Y": {"G", "A"},
	"S": {"C", "G"},
	"W": {"T", "A"},
	"K": {"C", "A"},
	"M": {"T", "G"},

	"B": {"G", "C", "A"},
	"D": {"T", "C", "A"},
	"H": {"T", "G", "A"},
	"V": {"T", "G", "C"},

	"N": {"A", "T", "G", "C"},
	"I": {"A", "T", "G", "C"},
}

var identical = map[string][]string{
	"A": {"A"},
	"T": {"T"},
	"G": {"G"},
	"C": {"C"},
	"U": {"T"},

	"R": {"A", "G"},
	"Y": {"C", "T"},
	"S": {"G", "C"},
	"W": {"A", "T"},
	"K": {"G", "T"},
	"M": {"A", "C"},

	"B": {"C", "G", "T"},
	"D": {"A", "G", "T"},
	"H": {"A", "C", "T"},
	"V": {"A", "C", "G"},

	"N": {"A", "T", "G", "C"},
	"I": {"A", "T", "G", "C"},
}


func isIdentical(bas byte, template byte) bool {
	if isValidBase(bas) != true || isValidBase(template) != true {
		return false
	}
	basConv := strings.ToUpper(string(bas))
	templConv := strings.ToUpper(string(template))
	matching := false

	lst, ok := identical[basConv]
	if !ok {return false;}

	for _, item := range lst {
		if item == templConv {
			matching = true
		}
	}
	return matching
}


func isComplementMatch(bas byte, template byte) bool {
	if isValidBase(bas) != true || isValidBase(template) != true {
		return false
	}
	basConv := strings.ToUpper(string(bas))
	templConv := strings.ToUpper(string(template))
	matching := false

	lst, ok := matches[basConv]
	if !ok {return false;}

	for _, item := range lst {
		if item == templConv {
			matching = true
		}
	}
	return matching
}


func isValidBase(b byte) bool {
	isValid := false
	bString := strings.ToUpper(string(b))
	for i := 0; i < len(Iupac); i++ {
		if bString == string(Iupac[i]) {
			isValid = true
		}
	}
	return isValid
}


func ParseFasta(path string) ([]Seqr, error) {
	// initialize the result slice
	var res []Seqr

	//file opening and assigning reader
	f, err := os.Open(path)
	if err != nil {
		return res, fmt.Errorf("Error while trying to open a file")
	}
	defer f.Close()

	r := bufio.NewReader(f)

	// ------------------------------
	// Go into the Seqr building loop

	firstPass := true
	buildHeader := true
	var currentRec = Seqr{}

	for true {
		// read a byte and handle EOF or errors
		c, err := r.ReadByte()
		if err == io.EOF {
			break
		} else if err != nil {
			return []Seqr{}, fmt.Errorf("Error while reading input: %w", err)
		}
		
		// Handle first pass > symbol OR a new sequence with reset
		if c == '>' && firstPass == true {
			firstPass = false
			continue
		} else if c != '>' && firstPass == true {
			return []Seqr{}, fmt.Errorf("fasta formatting, initial position is not >")
		} else if c == '>' && firstPass == false {
			res = append(res, currentRec)
			currentRec = Seqr{}
			buildHeader = true
			continue
		}

		// Handle backslash n and or backslash r delimiting the header
		if (c == '\n' || c == '\r') && buildHeader == true {
			buildHeader = false
			continue
		}

		// handle any letter or character
		if buildHeader == true {
			currentRec.Header = currentRec.Header + string(c)
			continue
		} else if buildHeader == false {
			if isValidBase(c) == true {
				currentRec.Seq = currentRec.Seq + strings.ToUpper(string(c))
			}
			continue
		}
	}
	return res, nil
}
