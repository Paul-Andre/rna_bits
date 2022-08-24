from typing import List, Tuple, Union, DefaultDict

from collections import defaultdict

from .representation import *

class ParseError(Exception):
    def __init__(self, message=None, text=None, column=None, line=None):
        self.text = text
        self.column = column
        self.message = message
        self.line = line
        super().__init__(message)

    def print(self):
        #print(dir(self))
        print(self.message)
        if (self.text is not None and self.column is not None):
            # TODO: handle tabs and stuff
            print(self.text)
            print(" "*self.column + "^")

class ParseContext:
    def __init__(self, text=None, column=None, line=None):
        self.text = text
        self.column = column
        self.line = line
    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_value, exc_tb):
        if isinstance(exc_value, ParseError):
            if (self.column): exc_value.column = self.column
            if (self.line): exc_value.line = self.line
            if (self.text): exc_value.text = self.text
            raise exc_value
        # if isinstance(exc_value, Exception):
        #     raise ParseError(
        #             message = "Parse error: \n" + repr(exc_value),
        #             text = self.text,
        #             column = self.column,
        #             line = self.line,
        #             )

def isParen(a: str) -> bool:
    assert(len(a) == 1)
    return a.isalpha() or a in "{}()[]<>"


# All valid:
# <<<<[[[[....>>>>]]]]
# ((((AAAA....))))aaaa
# AAAA{{{{....aaaa}}}} 

def isOpenParen(c):
    assert(len(c) == 1)
    return (c in "({<[") or (c.isalpha() and c.isupper())

def isClosingParen(c):
    assert(len(c) == 1)
    return (c in ")}>]") or (c.isalpha() and c.islower())

def getOpenParen(c):
    assert(len(c) == 1)
    if c.isalpha():
        return c.upper()
    else:
        return {
            ")":"(",
            "}":"{",
            ">":"<",
            "]":"[",
        }[c]

def isDot(c):
    assert(len(c) == 1)
    return c == "."

def isSeparator(c):
    assert(len(c) == 1)
    return c == "&" or c == " "


def parseParens(s: Union[List[str], str]) -> List[Tuple[int,int]]:
    stacks: DefaultDict[str, List[int]] = defaultdict(list)
    links = []
    pos = 1
    for c in s:
        if isSeparator(c):
            pass
        else:
            if isDot(c):
                pass
            elif isOpenParen(c):
                stacks[c].append(pos)
            elif isClosingParen(c):
                d = getOpenParen(c)
                if len(stacks[d]) > 0:
                    open_pos = stacks[d].pop()
                    links.append((open_pos,  pos))
                else:
                    raise ParseError("Closing bracket " + repr(c) + " at position " + repr(i) + " has no openning bracket.",
                                     text=s,
                                     column=i)
            else:
                assert False, "Unrecognized character"
            pos+=1

    return links

def letterToNumber(a:str) -> int:
    #TODO: make more complicated strand names like AA, AB, AC work
    assert(len(a) == 1)
    if (a.isupper() and a.isalpha()):
        return ord(a) - ord("A") + 1
    if (a.islower() and a.isalpha()):
        return ord(a) - ord("a") + 1

def splitOffLetters(a:str) -> int:
    for (i,c) in enumerate(a):
        if c.isnumeric:
            break
    if (i == len(a)):
        raise ParseError("Expected nucleotide number")
    return (a[:i], a[i:])

# splits a string in a way where you can recover the positions
def splitAndPositions(string, splittingChar, startPos=0):
    last = 0
    i = 0
    while i<len(string):
        if string[i] == splittingChar:
            yield (last + startPos, string[last:i])
            last = i+1
        i+=1
    yield (last + startPos, string[last:i])


def parseMotifDefinition(line, strandLengths):
    """
    ex:
    directory: 1 2 3-6 12-14 A32-45 #23-#25
    The nucleotides in the motif pdb are mapped to the nucleotide numbers on
    the right-hand side. If there's more than 1 strand, the numbers need to
    have the strand name (A,B,C...) prepended, or # prepended to mean that we
    count the nucleotides in all strands together
    """
    # TODO make it parse in a nice way, with good error messages
    cumulStrandLengths = [0]
    for l in strandLengths:
        cumulStrandLengths.append(l+cumulStrandLengths[-1])

    filename, l = line.split(":")
    l = l.strip().split()

    nucIndices = []

    for p in l:
        if "-" in p:
            a,b = p.split("-")
        else:
            a = p
            b = p

        al, an = splitOffLetters(a)
        bl, bn = splitOffLetters(b)

        an = int(an)
        bn = int(bn)

        if al == "":
            if len(strandLengths) > 1:
                raise ParseError(
                "When more than 1 strand, need to either prepend a strand letter to the nucleotide "
                "numbers, or prepend them with '#' to indicate we number the nucleotides starting "
                "from the beginning. Ex: A2-10, B23-40, #34-30"
                )
            al = "A"

        if bl == "":
            bl = al
        

        def asdf(a,al,an):
            if al != "#":
                aln = letterToNumber(al)
                if an > strandLengths[aln-1] or an==0:
                    raise ParseError("Invalid nucleotide number " + a)
                apos = cumulStrandLengths[aln-1]+an
            else: 
                apos = an
            return apos

        apos = asdf(a,al,an)
        bpos = asdf(b,bl,bn)
        
        if bpos < apos:
            raise ParseError("Nucleotide {} comes before {}".format(b,a))

        nucIndices.extend(range(apos,bpos+1))

    return filename, nucIndices


def parse(infile):
    """
    Expected file structure:
GCUGGGAUGUUGGCUUAGAAGCAGCCAUCAUUUAAAGAGUGCGUAACAGCUCACCAGC
(((((.(((.((.............)).))).........((......))...)))))
1NJI.A.92: 42-49
1YL3.A.78: 12-26
2OM7.I.2: 5-7, 31-41, 50-54
    """
    first_line = infile.__next__().strip()
    second_line = infile.__next__().strip()

    print(first_line)
    print(second_line)

    # Check that the strands have the same length
    # TODO: use exception
    assert len(first_line) == len(second_line)
    for a,b in zip(first_line, second_line):
        assert isSeparator(a) == isSeparator(b)
        if not isSeparator(a):
            assert a in "ACGUacgu"
            assert isParen(a) or isDot(a)

    # Create Nucleotide objects
    nucsByIndex = [None] # None added to offset the array by 1
    index = 1
    strandId = 1
    posInStrand = 1
    prev = None
    for c in first_line:
        if isSeparator(c):
            strandId+=1
            posInStrand=0
            prev = None
            pass
        else:
            nuc = Nucleotide(index, (strandId, posInStrand), c)
            nuc.prev = prev
            if prev is not None:
                prev.next = nuc
            nucsByIndex.append(nuc)
            assert(nucsByIndex[index] is nuc)
            index+=1
            posInStrand+=1


    # Add basepairs from secondary structure
    with ParseContext(line=2):
        pairings = parseParens(second_line)

    for (i,j) in pairings:
        a = nucsByIndex[i]
        b = nucsByIndex[j]
        pair = Pair(a,b)
        a.pair = pair
        b.pair = pair
        a.paired = b
        b.paired = a


    # Detect stacks
    for nuc in nucsByIndex[1:]:
        if nuc.paired is None:
            continue
        if nuc.pair.a is not nuc:
            # Avoid visiting a pair from both sides
            continue
        if nuc.next is None:
            continue


    
        

    isMultistrand = (first_line.find("&") != -1)
    strandLengths = [ len(a) for a in first_line.split("&") ]
    cumulStrandLengths = [strandLengths[0]]
    for l in strandLengths[1:]:
        cumulStrandLengths.append(l+cumulStrandLengths[-1])

    # Parse the motif insertions
    # the file excluding the first two lines
    for i,line in enumerate(infile):
        with ParseContext(line=2+i):
            print(parseMotifDefinition(line, strandLengths))









