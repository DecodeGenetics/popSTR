## BinSTR

Bins variant alleles together. Reads and writes VCFs.

### Usage

```sh
$ ./binstr -h
usage: binstr [-h] [--id2bins ID2BINS] [--output OUTPUT] input

positional arguments:
  input              Input VCF files. The VCF file may be uncompressed,
                     bgzipped, or gzipped.

optional arguments:
  -h, --help         show this help message and exit
  --id2bins ID2BINS  File containing Variant ID and comma-separated allele
                     indexes to bin with bins separated by semicolons. E.g.
                     chr11:12345:SG 0,1;2;3,6;4,5,7
  --output OUTPUT    Output VCF file.
```
