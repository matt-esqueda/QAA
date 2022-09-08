# Author: Matt Esqueda mesqueda@uoregon.edu

__version__ = "0.6"


DNAbases = set('ATGCNatcgn')
RNAbases = set('AUGCNaucgn')


def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return ord(letter) - 33


def qual_score(phred_score: str) -> float:
    """Calcuates the average quality score of a phred string"""
    total_score = 0
    for letter in phred_score:
        total_score += convert_phred(letter)
    return total_score/len(phred_score)


def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests")


def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_DNA_seq(DNA), 'String contains invalid characters'
    
    DNA = DNA.upper()
    Cs = DNA.count('C')
    Gs = DNA.count('G')
    return (Cs + Gs)/len(DNA)


def oneline_fasta(file):
    ''' Reads a fasta file and writes to a new fasta file with the sequence line from each record contained in one line'''
    with open(output_file, 'w') as out:                                    
        with open(file, 'r') as fa:
            for i, line in enumerate(fa):
                line = line.strip('\n')
                if i == 0:
                    print(line, file=out)
                elif line[0] == '>' and i != 0:
                    print(file=out)
                    print(line, file=out)
                else:
                    print(line, end='', file=out)
            print(file=out)



    