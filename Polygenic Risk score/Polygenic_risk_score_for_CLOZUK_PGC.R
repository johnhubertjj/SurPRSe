### PRS with plink###

plink --file data --extract myrange.txt --range
All SNPs within that range will then be excluded or extracted. The format of myrange.txt should be, one range per line, whitespace-separated:
  CHR     Chromosome code (1-22, X, Y, XY, MT, 0)
BP1     Start of range, physical position in base units
BP2     End of range, as above
LABEL   Name of range/gene
For example,
2 30000000 35000000  R1
2 60000000 62000000  R2
X 10000000 20000000  R3

