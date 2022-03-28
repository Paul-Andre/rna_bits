import sys

import parser

def main():
    if len(sys.argv) <= 1:
        print("Need input file")
        exit(0)

    indir = sys.argv[1]
    infile = open(indir)

    try:
        parsed = parser.parse(infile)
    except parser.ParseError as error:
        error.print()
        exit(1)

if __name__ == "__main__":
    main()

