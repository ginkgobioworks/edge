from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


class IO(object):
    """
    Class for exporting genome sequences and features.
    """

    def __init__(self, genome):
        self.__genome = genome.indexed_genome()

    def to_fasta(self, filename):
        """
        Export to FASTA format, saving to the specified filename.
        """
        outf = open(filename, 'w')
        for fragment in self.__genome.fragments.all():
            fragment = fragment.indexed_fragment()
            outf.write('>%s\n' % (fragment.name,))
            outf.write(fragment.sequence)
            outf.write('\n')
        outf.close()

    def to_gff(self, filename):
        """
        Export to GFF format, saving to the specified filename.
        """
        records = []

        for fragment in self.__genome.fragments.all():
            fragment = fragment.indexed_fragment()
            seq = Seq(fragment.sequence)
            rec = SeqRecord(seq, "%s" % (fragment.name,))
            features = []

            for annotation in fragment.annotations():
                # FeatureLocation first bp is AfterPosition, so -1
                loc = FeatureLocation(annotation.base_first-1, annotation.base_last)
                qualifiers = {'name': annotation.feature.name}
                feature = SeqFeature(loc,
                                     type=annotation.feature.type,
                                     strand=1,
                                     qualifiers=qualifiers)
                features.append(feature)

            rec.features = features
            records.append(rec)

        with open(filename, "w") as out_handle:
            GFF.write(records, out_handle, include_fasta=True)
