import re
import os
import tempfile
import subprocess
from django.conf import settings


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


def parse_primer3_output(d):
    primers = {}

    for k in d:
        m = re.search('\d+', k)
        if m:
            i = int(m.group(0))
            if i not in primers:
                primers[i] = {}
            k2 = re.sub('_%s_' % i, '_', k)
            primers[i][k2] = d[k]

    return [primers[k] for k in sorted(primers.keys())]


def primer3_run(opts):
    fn = None
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        fn = f.name
        for k in opts:
            f.write('%s=%s\n' % (k, opts[k]))
        f.write('PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s/primer3_config/\n' % settings.PRIMER3_DIR)
        f.write('=')

    cmd = "%s/primer3_core %s" % (settings.PRIMER3_DIR, fn)
    out = subprocess.check_output(cmd.split(' '))
    os.unlink(fn)
    r = {}
    for l in out.split('\n'):
        if '=' in l:
            l = l.split('=')
            r[l[0]]='='.join(l[1:])

    return parse_primer3_output(r)


def primer3_design(template, roi_start, roi_len, opts):

    min_primer_product_size = min(PRIMER3_MIN_PRIMER*2+roi_len, len(template))
    max_primer_product_size = len(template)

    opts = {}
    opts.update(PRIMER3_DEFAULTS)
    opts.update(opts)
    opts.update(dict(PRIMER_PRODUCT_SIZE_RANGE='%s-%s' % (min_primer_product_size,
                                                          max_primer_product_size),
                     SEQUENCE_ID='_',
                     SEQUENCE_TEMPLATE=template,
                     SEQUENCE_TARGET='%s,%s' % (roi_start, roi_len)))

    return primer3_run(opts)
    

def design_primers(fragment, roi_start_bp, roi_len, upstream_window, downstream_window, opts,
                   replacement_sequence=None):
    """
    Design primer using primer3 pipeline. Returns list of designed primers.
    """

    fragment = fragment.indexed_fragment()

    if roi_start_bp+roi_len-1 <= fragment.length:
        roi = fragment.get_sequence(bp_lo=roi_start_bp, bp_hi=roi_start_bp+roi_len-1)
    else:
        end_bp = fragment.circ_bp(roi_start_bp+roi_len-1)
        roi = fragment.get_sequence(bp_lo=roi_start_bp)+fragment.get_sequence(bp_hi=end_bp)

    bp_lo = roi_start_bp-upstream_window
    bp_hi = roi_start_bp+roi_len-1+downstream_window
    if fragment.circular is True:
        bp_lo = fragment.circ_bp(bp_lo)
        bp_hi = fragment.circ_bp(bp_hi)
    else:
        bp_lo = 1 if bp_lo < 1 else bp_lo
        bp_hi = fragment.length if bp_hi > fragment.length else bp_hi

    template = None
    if bp_hi < bp_lo:
        assert fragment.circular is True
        template_p1 = fragment.get_sequence(bp_lo=bp_lo)
        template_p2 = fragment.get_sequence(bp_hi=bp_hi)
        template = template_p1+template_p2
    else:
        template = fragment.get_sequence(bp_lo=bp_lo, bp_hi=bp_hi)

    roi_start = template.index(roi)+1
    if replacement_sequence is not None:
        re.sub(roi, replacement_sequence, template)
        roi = replacement_sequence
    return primer3_design(template, roi_start, len(roi), opts)
