from Bio import Entrez
from Bio import SeqIO
Entrez.email = "sefa1@umbc.edu"

def read_genbank(genome_accession_no, genbank_file=None):
    """Read genbank file. If the file is not given, based on the
    genome_accession_no, grab it from NCBI and parse it. Return Sequence Record
    object."""
    
    if genbank_file:
        print "reading genbank file %s" % genbank_file
        seq_record = SeqIO.read(genbank_file, "genbank")
    else:
        print "downloading and parsing genbank file for %s" % genome_accession_no
        handle = Entrez.efetch(db="nucleotide", rettype="gb",
                               retmode="text", id=genome_accession_no)
        seq_record = SeqIO.read(handle, "gb")
        handle.close()
    return seq_record

def extract_genes(seq_record):
    """Given BioPython SeqRecord object as argument, return the list of all
    genes where each gene is a SeqFeature object)"""
    return [f for f in seq_record.features if f.type == "gene"]

def extract_cds(seq_record):
    """Given BioPython SeqRecord object as argument, return the list of all
    coding sequences where each one is a SeqFeature object"""
    return [f for f in seq_record.features if f.type == "CDS"]

def reverse_complement(seq):
    return Seq(seq).reverse_complement().tostring()

def split_len(seq, length):
    """Given a string, returns a list containing _length_ sized pieces of the seq. For
    example, split_len('abcdefgh', 3) = ['abc', 'def', 'gh']"""
    return [seq[i:i+length] for i in range(0, len(seq), length)]
