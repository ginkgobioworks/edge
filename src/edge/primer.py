import re
import os
import tempfile
import subprocess
from django.conf import settings
from Bio.Seq import Seq


PRIMER3_MIN_PRIMER = 18
PRIMER3_MAX_PRIMER = 25

PRIMER3_INT_DEFAULTS = {
    'PRIMER_OPT_SIZE': 20,
    'PRIMER_PICK_INTERNAL_OLIGO': 1,
    'PRIMER_INTERNAL_MAX_SELF_END': 8,
    'PRIMER_MIN_SIZE': PRIMER3_MIN_PRIMER,
    'PRIMER_MAX_SIZE': PRIMER3_MAX_PRIMER,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_SELF_ANY': 12,
    'PRIMER_MAX_SELF_END': 8,
    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    'PRIMER_PAIR_MAX_COMPL_END': 8,
    'PRIMER_MAX_POLY_X': 100,
    'PRIMER_INTERNAL_MAX_POLY_X': 100,
}

PRIMER3_FLOAT_DEFAULTS = {
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 57.0,
    'PRIMER_MAX_TM': 63.0,
    'PRIMER_MIN_GC': 20.0,
    'PRIMER_MAX_GC': 80.0,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_DNA_CONC': 5000.0,  # nM
}


def parse_primer3_output(d):
    primers = {}

    for k in d:
        m = re.search(br'\d+', k)
        if m:
            i = int(m.group(0))
            if i not in primers:
                primers[i] = {}
            k2 = re.sub(b'_%d_' % i, b'_', k)
            primers[i][k2.decode()] = d[k].decode()

    return [primers[k] for k in sorted(primers.keys())]


def primer3_run(opts):
    fn = None
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        fn = f.name
        for k in opts:
            f.write('%s=%s\n' % (k, opts[k]))
        f.write('PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s\n' % settings.PRIMER3_CONFIG_DIR)
        f.write('=')

    cmd = "%s %s" % (settings.PRIMER3_BIN, fn)
    out = subprocess.check_output(cmd.split(' '))
    os.unlink(fn)
    r = {}
    for ll in out.split(b'\n'):
        if b'=' in ll:
            ll = ll.split(b'=')
            r[ll[0]] = b'='.join(ll[1:])

    return parse_primer3_output(r)


def design_primers_from_template(template, roi_start, roi_len, junctions, primer3_opts):
    primer3_opts = {} if primer3_opts is None else primer3_opts
    min_primer_product_size = min(PRIMER3_MIN_PRIMER * 2 + roi_len, len(template))
    max_primer_product_size = len(template)

    opts = {}
    opts.update(PRIMER3_INT_DEFAULTS)
    opts.update(PRIMER3_FLOAT_DEFAULTS)
    for k, v in primer3_opts.items():
        if k in PRIMER3_INT_DEFAULTS:
            opts[k] = int(v)
        elif k in PRIMER3_FLOAT_DEFAULTS:
            opts[k] = float(v)

    opts.update(dict(PRIMER_PRODUCT_SIZE_RANGE='%s-%s' % (min_primer_product_size,
                                                          max_primer_product_size),
                     SEQUENCE_ID='_',
                     SEQUENCE_TEMPLATE=template,
                     SEQUENCE_TARGET='%s,%s' % (roi_start, roi_len)))

    primers = primer3_run(opts)
    if junctions and len(junctions) > 0:
        for primer in primers:
            p1 = primer['PRIMER_LEFT_SEQUENCE']
            p2 = primer['PRIMER_RIGHT_SEQUENCE']
            p2_rc = str(Seq(p2).reverse_complement())
            p1_i = template.lower().index(p1.lower())
            p2_i = template.lower().index(p2_rc.lower())
            primer['PRIMER_LEFT_SEQUENCE_DISTANCE_TO_JUNCTION'] = junctions[0] - (p1_i + len(p1))
            primer['PRIMER_RIGHT_SEQUENCE_DISTANCE_TO_JUNCTION'] = p2_i - junctions[-1]

    return primers


def design_primers(fragment, roi_start_bp, roi_len, upstream_window, downstream_window, opts):
    """
    Design primer using primer3 pipeline. Returns list of designed primers.
    """

    opts = {} if opts is None else opts
    fragment = fragment.indexed_fragment()

    roi = fragment.get_sequence(bp_lo=roi_start_bp, bp_hi=roi_start_bp + roi_len - 1)
    bp_lo = roi_start_bp - upstream_window
    bp_hi = roi_start_bp + roi_len - 1 + downstream_window
    template = fragment.get_sequence(bp_lo=bp_lo, bp_hi=bp_hi)

    roi_start = template.index(roi) + 1
    return design_primers_from_template(template, roi_start, len(roi),
                                        [upstream_window, upstream_window + len(roi) - 1], opts)
