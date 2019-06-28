import gzip

def open_vcf_file(filename):
    _vcf_f = None

    if filename.endswith(".gz"):
        _vcf_f = gzip.open(filename, "r")
    else:
        _vcf_f = open(filename, "r")

    return _vcf_f


def read_vcf_header(vcf_f):
    header = ""
    samples_line = ""

    while True:
        line = vcf_f.readline().decode("utf-8")

        if line.startswith("##"):
            header += line
            continue

        samples_line = line
        break

    return header, samples_line
