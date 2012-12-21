from bioutils import *
from collections import namedtuple
import sys

# define promoter region as a named tuple:
# pos_start, pos_end - start and end position of the promoter (start < end)
# strand: strand of the transcribed gene
# gene, gene_locus_tag, gene_pos_start, gene_pos_end - name, locus tag and position
# of the transcribed gene
PromoterRegion = namedtuple("PromoterRegion",
                            "pos_start pos_end strand gene " +
                            "gene_locus_tag gene_pos_start gene_pos_end")

# use CodingRegion named tuple to represent sequneces
# The namedtuple CodingRegion contains three fields
# start, end -- start and end positions of the coding interval
# locus_tags -- locus tags of the genes in that region. It is a list as a
# contiguous region may contain multiple following genes. Instead of gene names,
# locus tags are used since some genes are not named in the genbank
# file. However, all coding sequences have the locus_tag information (hopefully).
CodingRegion = namedtuple("CodingRegion", "start end locus_tags")

# Intergenic region tuple contains follwing fieds
# start, end -- positions
# left_cr, right_cr -- left and right coding regions (CodingRegion objects)
IntergenicRegion = namedtuple("IntergenicRegion", "start end left_cr right_cr")


def extract_promoters(genome_accession_no, genbank_file=None,
                      extract_mode="fixed",
                      adaptive_mode_threshold=50,
                      promoter_start=-300, promoter_end=50,
                      output_file="promoter_regions.fasta"):
    """Extract promoter regions and report them in FASTA format.

    genome_accession_no -- genome accession number. If no genbank file is given,
    the genbank file associated with it is retrieved from NCBI.
    
    genbank_file -- genbank file to be parsed. If not given, genome_accession_no
    is used.
    
    extract_mode -- Three different beahviour is implemented.
    (1) fixed mode, the default one. Retrieve [-300, +50] region of every gene.
    (2) whole_chunk mode retrieves [-300, +50] region, just like fixed
    mode. The difference is that it is in the output set if and only if it does
    not overlap with another gene. For example let gene A be at [600, 1000], the
    promoter region for gene A be at [300, 650]. If there is another gene B at
    position [100, 350], since gene B and promoter for gene A overlap, it is not
    presented in the final list of promoter regions.
    (3) adaptive mode returns overlapped promoters as well, but only the
    non-overlapped sections of them. In the example given above, it returns
    promoter region of gene A as [351, 650], since the subsequence [300, 350]
    overlaps with gene B. Also, adaptive_mode_threshold parameter defines the
    minimum allowed length for promoter sequences.

    adaptive_mode_threshold -- In adaptive mode, the minimum length of allowed
    promoter region.

    promoter_start, promoter_end -- The promoter region relative to the gene
    start position. By default they are -300 and +50.

    output_file -- the output file name for the returned promoter sequences.
    The results are printed in fasta format.
    """

    # check arguments
    assert promoter_start < 0
    assert extract_mode in ["fixed", "whole_chunk", "adaptive"]
    
    # read genbank
    seq_record = read_genbank(genome_accession_no, genbank_file)
    # extract all genes (list of SeqFeature objects)
    genes = extract_genes(seq_record)
    # all genes should be sorted based on their positions on the genome
    assert all(ga.location.start <= gb.location.start
               for (ga, gb) in zip(genes, genes[1:]))
    
    # based on the extraction strategy, extract promoter regions from the
    # genome. Each promoter region is defined as a dictionary.
    promoter_regions = []
    
    for gindex, g in enumerate(genes):
        gene_start = g.location.nofuzzy_start
        gene_end = g.location.nofuzzy_end

        pr = None
        if (g.strand == 1 and
            gene_start + promoter_start >= 0 and
            gene_start + promoter_end < len(seq_record.seq)):
            pr_start = gene_start + promoter_start
            pr_end = gene_start + promoter_end
            assert pr_start < gene_start and pr_end > gene_start
            pr = PromoterRegion(pos_start = pr_start,
                                pos_end = pr_end,
                                strand = g.strand,
                                gene = ', '.join(g.qualifiers['gene']),
                                gene_locus_tag = ', '.join(g.qualifiers['locus_tag']),
                                gene_pos_start = g.location.nofuzzy_start,
                                gene_pos_end = g.location.nofuzzy_end)
            
        elif (g.strand == -1 and
              gene_end - promoter_start < len(seq_record.seq) and
              gene_end - promoter_end >= 0):
            pr_start = gene_end - promoter_end
            pr_end = gene_end - promoter_start
            assert pr_start < gene_end and pr_end > gene_end
            pr = PromoterRegion(pos_start = pr_start,
                                pos_end = pr_end,
                                strand = g.strand,
                                gene = g.qualifiers['gene'],
                                gene_locus_tag = g.qualifiers['locus_tag'],
                                gene_pos_start = g.location.nofuzzy_start,
                                gene_pos_end = g.location.nofuzzy_end)
            
        if not pr:
            continue

        if extract_mode == "fixed":
            promoter_regions.append(pr)

        elif extract_mode == "whole_chunks":
            # if pr overlaps with previous gene, discard it
            if (gindex > 0 and pr.pos_start <= genes[gindex-1].location.nofuzzy_end):
                pass
            elif (gindex < len(genes)-1 and
                  pr.pos_end >= genes[gindex+1].location.nofuzzy_start):
                pass
            else:
                # yey, it's valid
                promoter_regions.append(pr)

        elif extract_mode == "adaptive":
            pr_start = pr.pos_start
            pr_end = pr.pos_end

            if gindex > 0:
                pr_start = max(pr_start, genes[gindex-1].location.nofuzzy_end + 1)
            if gindex < len(genes)-1:
                pr_end = min(pr_end, genes[gindex+1].location.nofuzzy_start)
            # check the length of trimmed promoter region
            if pr_end - pr_start + 1 >= adaptive_mode_threshold:
                # valid
                # since tuple is immutable, create a new one
                new_pr = PromoterRegion(pos_start = pr_start,
                                        pos_end = pr_end,
                                        strand = pr.strand,
                                        gene = pr.gene,
                                        gene_locus_tag = pr.gene_locus_tag,
                                        gene_pos_start = pr.gene_pos_start,
                                        gene_pos_end = pr.gene_pos_end)
                promoter_regions.append(new_pr)

        # write output in FASTA format
        print "writing output to '%s' in fasta format." % output_file
        with open(output_file, 'w') as f:
            for pr in promoter_regions:
                f.write(">%s\n" % str(pr))
                # get promoter sequence
                pseq = seq_record.seq[pr.pos_start:pr.pos_end]
                if pr.strand == -1:
                    # reverse complement
                    pseq = pseq.reverse_complement()
                # print 100 chars each line
                f.write("%s\n" % '\n'.join(split_len(pseq.tostring(), 100)))

    return promoter_regions

def extract_intergenics(genome_accession_no, genbank_file=None,
                        coding_output_file="coding_regions.fasta",
                        intergenic_output_file="intergenic_regions.fasta"):
    """Split coding regions and non-coding regions of a genome and report them
    in FASTA format.

    genome_accession_no -- genome accession number, if the genbank file is not
    given, the genbank file with this accession number is retrieved from NCBI.

    genbank_file -- The file to be used. If given, the previous argument
    _genome_accession_no_ is obsolete.

    coding_output_file -- the output file name for coding region sequences. The
    output is written in FASTA format. Each FASTA entry corresponds to a
    contiguous region in genome. Its header line has the start/end positions and
    locus tags of genes in that region.
    
    intergenic_output_file -- the output file name for intergenic regions. The
    output is in FASTA format. Each entry corresponds to a contiguous region in
    genome. 
    """

    # no surprises, read genome file first
    seq_record = read_genbank(genome_accession_no, genbank_file)
    # extract all coding sequences (list of SeqFeature objects)
    cds = extract_cds(seq_record)
    # make sure all coding regions are sorted by position
    cds.sort(key=lambda self: self.location.nofuzzy_start)
    # find all coding region intervals
    assert len(cds) > 0

    # define first coding region interval
    cr = CodingRegion(start=cds[0].location.nofuzzy_start,
                      end=cds[0].location.nofuzzy_end,
                      locus_tags=[cds[0].qualifiers['locus_tag'][0]])
    coding_regions = [cr]
    
    for i in cds[1:]:
        i_start = i.location.nofuzzy_start
        i_end = i.location.nofuzzy_end
        i_locus_tags = i.qualifiers['locus_tag'][0]
        if coding_regions[-1].end < i_start:
            # new interval
            cr = CodingRegion(start=i_start,
                              end=i_end,
                              locus_tags=[i_locus_tags])
            coding_regions.append(cr)
        else:
            # merge this into the last coding interval
            # implementation detail: since namedtuples are immutable, create new
            # one and replace it with the last element in the list
            cr = CodingRegion(coding_regions[-1].start,
                              i_end,
                              coding_regions[-1].locus_tags + [i_locus_tags])
            coding_regions[-1] = cr

    # calculate intergenic regions
    assert len(coding_regions) > 0
    intergenic_regions = []
    # first intergenic region
    if (coding_regions[0].start > 0):
        intergenic_regions.append(IntergenicRegion(start=0,
                                                   end=coding_regions[0].start,
                                                   left_cr=None,
                                                   right_cr=coding_regions[0]))
    # include all others
    for left_cr, right_cr in zip(coding_regions, coding_regions[1:]):
        intergenic_regions.append(IntergenicRegion(start=left_cr.end,
                                                   end=right_cr.start,
                                                   left_cr=left_cr,
                                                   right_cr=right_cr))
    # aaaand the last one
    if (coding_regions[-1].end < len(seq_record.seq)):
        intergenic_regions.append(IntergenicRegion(start=coding_regions[-1].end,
                                                   end=len(seq_record.seq),
                                                   left_cr=coding_regions[-1],
                                                   right_cr=None))
    
    # writing time!
    print "Writing coding regions into the file %s" % coding_output_file
    with open(coding_output_file, 'w') as f:
        for cr in coding_regions:
            cseq = seq_record.seq[cr.start:cr.end]
            f.write(">%s\n" % str(cr))
            f.write("%s\n" % '\n'.join(split_len(cseq.tostring(), 100)))

    # write intergenic regions too.
    print "Writing intergenic regions into the file %s" % intergenic_output_file
    with open(intergenic_output_file, 'w') as f:
        for ir in intergenic_regions:
            iseq = seq_record.seq[ir.start:ir.end]
            f.write(">%s\n" % str(ir))
            f.write("%s\n" % '\n'.join(split_len(iseq.tostring(), 100)))
        
            
print "loaded"

