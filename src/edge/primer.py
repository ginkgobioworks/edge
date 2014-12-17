import primer3.bindings as p3
from edge.models import Fragment


PRIMER3_MIN_PRIMER = 18
PRIMER3_MAX_PRIMER = 25

PRIMER3_DEFAULTS = {
    'PRIMER_OPT_SIZE': 20,
    'PRIMER_PICK_INTERNAL_OLIGO': 1,
    'PRIMER_INTERNAL_MAX_SELF_END': 8,
    'PRIMER_MIN_SIZE': PRIMER3_MIN_PRIMER,
    'PRIMER_MAX_SIZE': PRIMER3_MAX_PRIMER,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 57.0,
    'PRIMER_MAX_TM': 63.0,
    'PRIMER_MIN_GC': 20.0,
    'PRIMER_MAX_GC': 80.0,
    'PRIMER_MAX_POLY_X': 100,
    'PRIMER_INTERNAL_MAX_POLY_X': 100,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_DNA_CONC': 5000.0,  # nM
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_SELF_ANY': 12,
    'PRIMER_MAX_SELF_END': 8,
    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    'PRIMER_PAIR_MAX_COMPL_END': 8,
}


def primer3_design(fragment, product_start_bp, product_len,
                   front_distance, back_distance, primer3_opts):
    """
    Design primer using primer3 pipeline. Returns list of designed primers.
    """

    fragment = fragment.indexed_fragment()

    print 'desired product: %s, %s' % (product_start_bp, product_len)
    print 'fragment: %s' % fragment.sequence

    if product_start_bp+product_len-1 <= fragment.length:
        product = fragment.get_sequence(bp_lo=product_start_bp, bp_hi=product_start_bp+product_len-1)
    else:
        end_bp = fragment.circ_bp(product_start_bp+product_len-1)
        product = fragment.get_sequence(bp_lo=product_start_bp)+fragment.get_sequence(bp_hi=end_bp)

    min_primer_product_size = PRIMER3_MIN_PRIMER*2+len(product)
    max_primer_product_size = PRIMER3_MAX_PRIMER*2+len(product)+front_distance+back_distance

    bp_lo = fragment.circ_bp(product_start_bp-front_distance-PRIMER3_MAX_PRIMER)
    bp_hi = fragment.circ_bp(product_start_bp+product_len-1+back_distance+PRIMER3_MAX_PRIMER)

    template = None
    if bp_hi < bp_lo:
        assert fragment.circular is True
        template_p1 = fragment.get_sequence(bp_lo=bp_lo)
        template_p2 = fragment.get_sequence(bp_hi=bp_hi)
        template = template_p1+template_p2
    else:
        template = fragment.get_sequence(bp_lo=bp_lo, bp_hi=bp_hi)

    product_start = template.index(product)+1

    opts = {}
    opts.update(PRIMER3_DEFAULTS)
    opts.update(primer3_opts)
    opts['PRIMER_PRODUCT_SIZE_RANGE'] = [[min_primer_product_size, max_primer_product_size]]
    opts['PRIMER_TASK'] = 'pick_cloning_primers'
    opts['PRIMER_PICK_ANYWAY'] = 1

    print 'template %s' % template
    print 'product start %s, len %s' % (product_start, len(product))
    print 'product size range: %s' % opts['PRIMER_PRODUCT_SIZE_RANGE']

    r = p3.designPrimers(
        dict(SEQUENCE_ID='_', SEQUENCE_TEMPLATE=str(template),
             SEQUENCE_TARGET=[product_start, len(product)]), opts)

    return r
