from edge.blast import blast_genome
from Bio.Seq import Seq


def compute_pcr_product(primer_a_sequence, primer_a_blastres,
                        primer_b_sequence, primer_b_blastres):
    """
    Computes a PCR product based on two blast results.
    """

    MIN_IDENTITIES = 0.90
    MIN_BINDING_LENGTH = 10

    # primers must be on the same fragment
    if primer_a_blastres.fragment_id != primer_b_blastres.fragment_id:
        return None

    # primers must align to different senses
    primer_a_strand = primer_a_blastres.strand()
    primer_b_strand = primer_b_blastres.strand()
    if primer_a_strand == primer_b_strand:
        return None

    # forward primer aligns to the sense strand, binds to the antisense strand,
    # and elongates along the sense strand and copies sense strand bases; it
    # must align at the 3' end of the primer against the sense strand.
    #
    # reverse primer aligns to the antisense strand, binds to the sense strand,
    # and elongates along the antisense strand and copies antisense strand
    # bases; it must align at the 3' end of the primer against the antisense
    # strand.

    if primer_a_strand == 1:
        fwd_primer = primer_a_sequence
        fwd_primer_res = primer_a_blastres
        rev_primer = primer_b_sequence
        rev_primer_res = primer_b_blastres
    else:
        fwd_primer = primer_b_sequence
        fwd_primer_res = primer_b_blastres
        rev_primer = primer_a_sequence
        rev_primer_res = primer_a_blastres

    if fwd_primer_res.query_end != len(fwd_primer) or\
       rev_primer_res.query_end != len(rev_primer):
        return None

    # primers must align well
    if primer_a_blastres.identity_ratio() < MIN_IDENTITIES or\
       primer_b_blastres.identity_ratio() < MIN_IDENTITIES or\
       primer_a_blastres.alignment_length() < MIN_BINDING_LENGTH or\
       primer_b_blastres.alignment_length() < MIN_BINDING_LENGTH:
        return None

    # cannot produce product if elongated regions do not overlap
    if fwd_primer_res.subject_end >= rev_primer_res.subject_end:
        return None

    # get sequence between primers
    fragment = primer_a_blastres.fragment.indexed_fragment()
    product_mid = fragment.get_sequence(bp_lo=fwd_primer_res.subject_end+1,
                                        bp_hi=rev_primer_res.subject_end-1)
    product = '%s%s%s' % (fwd_primer, product_mid, str(Seq(rev_primer).reverse_complement()))
    return product


def pcr_from_genome(genome, primer_a_sequence, primer_b_sequence):
    """
    PCR on genome, using two primer sequences. Returns tuple of PCR product
    sequence, primer a blast results, and primer b blast results. Produuct
    sequence may be None if primers do not bind unique or in a way that can
    result in a PCR product.

    A PCR can create a product if each primer uniquely binds to the genome,
    primers bind to the same fragment, on sense and antisense strands
    respectively, and non-overlapping in their sense strand binding positions
    """

    primer_a_results = blast_genome(genome, 'blastn', primer_a_sequence)
    primer_b_results = blast_genome(genome, 'blastn', primer_b_sequence)

    pcr_products = []
    for a_res in primer_a_results:
        for b_res in primer_b_results:
            product = compute_pcr_product(primer_a_sequence, a_res,
                                          primer_b_sequence, b_res)
            if product is not None:
                pcr_products.append(product)

    if len(pcr_products) == 1:
        return (pcr_products[0], primer_a_results, primer_b_results)
    else:
        return (None, primer_a_results, primer_b_results)
