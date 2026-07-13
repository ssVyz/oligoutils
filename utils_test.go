package oligoutils

import "testing"
import "fmt"

func TestParseError(t *testing.T) {

	_, err := ParseFasta("not_exist.txt")
	if err != nil {
		fmt.Println(err)
	}
}

func TestIsValidBase(t *testing.T) {

	if isValidBase('A') != true {
		t.Fatalf("Valid for A is incorrect")
	}
	if isValidBase('Z') != false {
		t.Fatalf("Valid for A is incorrect")
	}
	if isValidBase('_') != false {
		t.Fatalf("Valid for A is incorrect")
	}

}

func TestParseFasta(t *testing.T) {

	content, _ := ParseFasta("example/exa1.fa")

	fmt.Printf("Reading first entry name: %v, and first entry sequence: %v\nLength of the slice: %v\n", content[0].Header, content[0].Seq, len(content))
}

func TestParseFastaMalformed(t *testing.T) {
	_, err := ParseFasta("example/malformed.txt")
	if err != nil {
		fmt.Println(err)
	} else {
		t.Fatalf("malformed did not result in an error")
	}
}

func TestBaseMatching(t *testing.T) {

	bases := "ART"
	template := "TCT"

	fmt.Println("Matching three bases:")
	fmt.Printf("input base %v, template base %v, matching: %v\n", string(bases[0]), string(template[0]), isComplementMatch(bases[0], template[0]))
	fmt.Printf("input base %v, template base %v, matching: %v\n", string(bases[1]), string(template[1]), isComplementMatch(bases[1], template[1]))
	fmt.Printf("input base %v, template base %v, matching: %v\n", string(bases[2]), string(template[2]), isComplementMatch(bases[2], template[2]))

}