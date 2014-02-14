from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


class IO(object):
    """
    Class for exporting fragment sequences and features.
    """

    def __init__(self, fragment):
        """
        Constructor. fragment should be the fragment object to export.
        """

        self.__fragment = fragment.indexed_fragment()

    def to_fasta(self, filename):
        """
        Export to FASTA format, saving to the specified filename.
        """
        outf = open(filename, 'w')
        outf.write('>fragment %s|%s\n' % (self.__fragment.id, self.__fragment.name))
        outf.write(self.__fragment.sequence)
        outf.write('\n')
        outf.close()

    def to_gff(self, filename):
        """
        Export to GFF format, saving to the specified filename.
        """
        seq = Seq(self.__fragment.sequence)
        rec = SeqRecord(seq, "fragment %s: %s" % (self.__fragment.id, self.__fragment.name))
        features = []

        for annotation in self.__fragment.annotations():
            # FeatureLocation first bp is AfterPosition, so -1
            loc = FeatureLocation(annotation.base_first-1, annotation.base_last)
            qualifiers = {'name': annotation.feature.name}
            feature = SeqFeature(loc, type=annotation.feature.type, strand=1, qualifiers=qualifiers)
            features.append(feature)

        rec.features = features

        with open(filename, "w") as out_handle:
            GFF.write([rec], out_handle)
