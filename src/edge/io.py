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

    def to_gff_file(self, file):
        """
        Export to GFF format, saving to provided file like object.
        """
        records = []

        for fragment in self.__genome.fragments.all():
            fragment = fragment.indexed_fragment()
            seq = Seq(fragment.sequence)
            rec = SeqRecord(seq, "%s" % (fragment.name,))
            features = []

            for annotation in fragment.annotations():
                # FeatureLocation first bp is AfterPosition, so -1
                loc = FeatureLocation(annotation.base_first - 1, annotation.base_last)
                qualifiers = annotation.feature.qualifiers

                # "phase" is how GFF annotations frame shifts. CDS with phase 1
                # means annotation is 1 bp frameshift. GFF parser we use parses
                # GFF features into BioPython SeqFeature and stores the "phase"
                # in the qualifiers. Here we guard against manual annotation
                # that set phase to None and change it to 0 instead, because it
                # will be output as None otherwise.
                if "phase" in qualifiers and qualifiers["phase"] is None:
                    qualifiers["phase"] = 0
                elif "Phase" in qualifiers and qualifiers["Phase"] is None:
                    qualifiers["Phase"] = 0

                qualifiers.update({'name': annotation.feature_name})
                strand = annotation.feature.strand
                feature = SeqFeature(loc,
                                     type=annotation.feature.type,
                                     strand=0 if strand is None else strand,
                                     qualifiers=qualifiers)
                features.append(feature)

            rec.features = features
            records.append(rec)

        GFF.write(records, file, include_fasta=True)

    def to_gff(self, filename):
        """
        Export to GFF format, saving to the specified filename.
        """
        with open(filename, "w") as out_handle:
            self.to_gff_file(out_handle)
