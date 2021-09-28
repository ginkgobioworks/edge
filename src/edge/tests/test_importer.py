import os
import tempfile
import unittest.mock as mock

from django.test import TestCase

from edge import import_gff
from edge.models import Genome


class ImporterTest(TestCase):
    def test_import_gff_procedure_creates_genome_and_annotations(self):

        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tcds\t30\t80\t.\t-\t.\tID=i2;Name=f2
chrI\tTest\trbs\t20\t28\t.\t+\t.\tID=i3
chrII\tTest\tgene\t40\t60\t.\t-\t.\tID=f4;gene=g4
chrII\tTest\tgene\t20\t80\t.\t+\t.\tID=i5;Name=f5
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
>chrII
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
            f.write(data)
            f.close()

            self.assertEquals(Genome.objects.filter(name="TestGenome").count(), 0)
            import_gff("TestGenome", f.name)
            self.assertEquals(Genome.objects.filter(name="TestGenome").count(), 1)

            # import again with same name does not work
            self.assertRaises(Exception, import_gff, "TestGenome", f.name)

            # can import again with different name
            self.assertEquals(Genome.objects.filter(name="TestGenome2").count(), 0)
            import_gff("TestGenome2", f.name)
            self.assertEquals(Genome.objects.filter(name="TestGenome2").count(), 1)

            genome = Genome.objects.get(name="TestGenome")

            os.unlink(f.name)

        # created one fragment for each sequence in GFF file
        self.assertCountEqual(
            [fr.name for fr in genome.fragments.all()], ["chrI", "chrII"]
        )
        chrI = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrI"
        ][0]
        self.assertEquals(len(chrI.sequence), 160)
        self.assertEquals(len(chrI.annotations()), 2)
        chrII = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrII"
        ][0]
        self.assertEquals(len(chrII.sequence), 160)
        self.assertEquals(len(chrII.annotations()), 2)

    def test_import_gff_creates_fragments_and_annotate_features(self):

        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tcds\t30\t80\t.\t-\t.\tID=i2;Name=f2
chrI\tTest\trbs\t20\t28\t.\t+\t.\tID=i3
chrII\tTest\tgene\t40\t60\t.\t-\t.\tID=f4;gene=g4
chrII\tTest\tgene\t20\t80\t.\t+\t.\tID=i5;Name=f5
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
>chrII
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
            f.write(data)
            f.close()
            genome = Genome.import_gff("Foo", f.name)
            os.unlink(f.name)

        # created one fragment for each sequence in GFF file
        self.assertCountEqual(
            [fr.name for fr in genome.fragments.all()], ["chrI", "chrII"]
        )

        # verify chrI fragment
        chrI = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrI"
        ][0]
        self.assertEquals(len(chrI.sequence), 160)
        # verify skips annotation on entire sequence
        self.assertEquals(len(chrI.annotations()), 2)
        self.assertEquals(chrI.annotations()[0].base_first, 20)
        self.assertEquals(chrI.annotations()[0].base_last, 28)
        self.assertEquals(
            chrI.annotations()[0].feature.name, "i3"
        )  # no name, loaded ID
        self.assertEquals(chrI.annotations()[0].feature.strand, 1)
        self.assertEquals(chrI.annotations()[1].base_first, 30)
        self.assertEquals(chrI.annotations()[1].base_last, 80)
        self.assertEquals(chrI.annotations()[1].feature.name, "f2")
        self.assertEquals(chrI.annotations()[1].feature.strand, -1)

        # verify chrII fragment
        chrII = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrII"
        ][0]
        self.assertEquals(len(chrII.sequence), 160)
        # consecutive annotations merged even though they span multiple chunks
        self.assertEquals(len(chrII.annotations()), 2)
        self.assertEquals(chrII.annotations()[0].base_first, 20)
        self.assertEquals(chrII.annotations()[0].base_last, 80)
        self.assertEquals(chrII.annotations()[0].feature.name, "f5")
        self.assertEquals(chrII.annotations()[0].feature.strand, 1)
        self.assertEquals(chrII.annotations()[1].base_first, 40)
        self.assertEquals(chrII.annotations()[1].base_last, 60)
        self.assertEquals(
            chrII.annotations()[1].feature.name, "g4"
        )  # has gene, use gene name
        self.assertEquals(chrII.annotations()[1].feature.strand, -1)

    def test_import_feature_starting_at_first_base(self):

        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tcds\t1\t80\t.\t-\t.\tID=i2;Name=f2
chrI\tTest\trbs\t20\t28\t.\t+\t.\tID=i3
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
            f.write(data)
            f.close()
            genome = Genome.import_gff("Foo", f.name)
            os.unlink(f.name)

        # verify chrI fragment
        chrI = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrI"
        ][0]
        self.assertEquals(len(chrI.sequence), 160)
        # verify skips annotation on entire sequence
        self.assertEquals(len(chrI.annotations()), 2)
        self.assertEquals(chrI.annotations()[0].base_first, 1)
        self.assertEquals(chrI.annotations()[0].base_last, 80)
        self.assertEquals(chrI.annotations()[0].feature.name, "f2")
        self.assertEquals(chrI.annotations()[0].feature.strand, -1)
        self.assertEquals(chrI.annotations()[1].base_first, 20)
        self.assertEquals(chrI.annotations()[1].base_last, 28)
        self.assertEquals(
            chrI.annotations()[1].feature.name, "i3"
        )  # no name, loaded ID
        self.assertEquals(chrI.annotations()[1].feature.strand, 1)

    def test_import_partially_overlapping_features(self):

        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tcds\t19\t21\t.\t-\t.\tID=i2;Name=f2
chrI\tTest\trbs\t20\t28\t.\t+\t.\tID=i3
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
            f.write(data)
            f.close()
            genome = Genome.import_gff("Foo", f.name)
            os.unlink(f.name)

        # verify chrI fragment
        chrI = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrI"
        ][0]
        self.assertEquals(len(chrI.sequence), 160)
        # verify skips annotation on entire sequence
        self.assertEquals(len(chrI.annotations()), 2)
        self.assertEquals(chrI.annotations()[1].base_first, 20)
        self.assertEquals(chrI.annotations()[1].base_last, 28)
        self.assertEquals(
            chrI.annotations()[1].feature.name, "i3"
        )  # no name, loaded ID
        self.assertEquals(chrI.annotations()[1].feature.strand, 1)
        self.assertEquals(chrI.annotations()[0].base_first, 19)
        self.assertEquals(chrI.annotations()[0].base_last, 21)
        self.assertEquals(chrI.annotations()[0].feature.name, "f2")
        self.assertEquals(chrI.annotations()[0].feature.strand, -1)

    def test_import_overlapping_features(self):

        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tcds\t20\t28\t.\t-\t.\tID=i2;Name=f2
chrI\tTest\trbs\t20\t28\t.\t+\t.\tID=i3
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
            f.write(data)
            f.close()
            genome = Genome.import_gff("Foo", f.name)
            os.unlink(f.name)

        # verify chrI fragment
        chrI = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrI"
        ][0]
        self.assertEquals(len(chrI.sequence), 160)
        # verify skips annotation on entire sequence
        self.assertEquals(len(chrI.annotations()), 2)
        self.assertEquals(chrI.annotations()[1].base_first, 20)
        self.assertEquals(chrI.annotations()[1].base_last, 28)
        self.assertEquals(
            chrI.annotations()[1].feature.name, "i3"
        )  # no name, loaded ID
        self.assertEquals(chrI.annotations()[1].feature.strand, 1)
        self.assertEquals(chrI.annotations()[0].base_first, 20)
        self.assertEquals(chrI.annotations()[0].base_last, 28)
        self.assertEquals(chrI.annotations()[0].feature.name, "f2")
        self.assertEquals(chrI.annotations()[0].feature.strand, -1)

    def test_import_feature_ending_at_last_base(self):

        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tcds\t20\t28\t.\t-\t.\tID=i2;Name=f2
chrI\tTest\trbs\t20\t160\t.\t+\t.\tID=i3
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
            f.write(data)
            f.close()
            genome = Genome.import_gff("Foo", f.name)
            os.unlink(f.name)

        # verify chrI fragment
        chrI = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrI"
        ][0]
        self.assertEquals(len(chrI.sequence), 160)
        # verify skips annotation on entire sequence
        self.assertEquals(len(chrI.annotations()), 2)
        self.assertEquals(chrI.annotations()[1].base_first, 20)
        self.assertEquals(chrI.annotations()[1].base_last, 160)
        self.assertEquals(
            chrI.annotations()[1].feature.name, "i3"
        )  # no name, loaded ID
        self.assertEquals(chrI.annotations()[1].feature.strand, 1)
        self.assertEquals(chrI.annotations()[0].base_first, 20)
        self.assertEquals(chrI.annotations()[0].base_last, 28)
        self.assertEquals(chrI.annotations()[0].feature.name, "f2")
        self.assertEquals(chrI.annotations()[0].feature.strand, -1)


class QualifierTest(TestCase):
    def import_with_qualifiers(self, qualifiers, phase="."):
        data = """##gff-version 3
chrI\tTest\tcds\t30\t80\t.\t-\t%s\t%s
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
""" % (
            phase,
            qualifiers,
        )
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
            f.write(data)
            f.close()
            self.genome = Genome.import_gff("Foo", f.name)
            os.unlink(f.name)

    def test_stores_phase(self):
        qualifiers = "name=g2"
        self.import_with_qualifiers(qualifiers)
        chrI = [
            f.indexed_fragment()
            for f in self.genome.fragments.all()
            if f.name == "chrI"
        ][0]
        self.assertEquals(len(chrI.annotations()), 1)
        self.assertEquals(chrI.annotations()[0].feature.name, "g2")
        self.assertEquals(
            chrI.annotations()[0].feature.qualifiers, dict(name=["g2"], source=["Test"])
        )

        qualifiers = "name=g2"
        self.import_with_qualifiers(qualifiers, phase=2)
        chrI = [
            f.indexed_fragment()
            for f in self.genome.fragments.all()
            if f.name == "chrI"
        ][0]
        self.assertEquals(len(chrI.annotations()), 1)
        self.assertEquals(chrI.annotations()[0].feature.name, "g2")
        self.assertEquals(
            chrI.annotations()[0].feature.qualifiers,
            dict(name=["g2"], phase=["2"], source=["Test"]),
        )

    def test_keeps_qualifiers(self):
        qualifiers = "name=g2;locus_tag=b0002;aliases=a,b,c"
        self.import_with_qualifiers(qualifiers)
        chrI = [
            f.indexed_fragment()
            for f in self.genome.fragments.all()
            if f.name == "chrI"
        ][0]
        self.assertEquals(len(chrI.annotations()), 1)
        self.assertEquals(chrI.annotations()[0].feature.name, "g2")
        self.assertEquals(
            chrI.annotations()[0].feature.qualifiers,
            dict(
                name=["g2"],
                locus_tag=["b0002"],
                aliases=["a", "b", "c"],
                source=["Test"],
            ),
        )

    def test_uses_name_qualifier_as_name_over_gene_qualifier(self):
        qualifiers = "ID=i2;gene=g2;Name=f2"
        self.import_with_qualifiers(qualifiers)
        chrI = [
            f.indexed_fragment()
            for f in self.genome.fragments.all()
            if f.name == "chrI"
        ][0]
        self.assertEquals(len(chrI.annotations()), 1)
        self.assertEquals(chrI.annotations()[0].feature.name, "f2")

    def test_uses_name_qualifier_as_name_over_locus_tag_qualifier(self):
        qualifiers = "ID=i2;Name=f2;locus_tag=l2"
        self.import_with_qualifiers(qualifiers)
        chrI = [
            f.indexed_fragment()
            for f in self.genome.fragments.all()
            if f.name == "chrI"
        ][0]
        self.assertEquals(len(chrI.annotations()), 1)
        self.assertEquals(chrI.annotations()[0].feature.name, "f2")

    def test_uses_locus_qualifier_as_name_as_name_over_id(self):
        qualifiers = "ID=i2;locus_tag=l2"
        self.import_with_qualifiers(qualifiers)
        chrI = [
            f.indexed_fragment()
            for f in self.genome.fragments.all()
            if f.name == "chrI"
        ][0]
        self.assertEquals(len(chrI.annotations()), 1)
        self.assertEquals(chrI.annotations()[0].feature.name, "l2")

    def test_uses_id_as_name_if_nothing_else_available(self):
        qualifiers = "ID=i2"
        self.import_with_qualifiers(qualifiers)
        chrI = [
            f.indexed_fragment()
            for f in self.genome.fragments.all()
            if f.name == "chrI"
        ][0]
        self.assertEquals(len(chrI.annotations()), 1)
        self.assertEquals(chrI.annotations()[0].feature.name, "i2")

    def test_uses_feature_type_as_name_if_no_id(self):
        qualifiers = ""
        self.import_with_qualifiers(qualifiers)
        chrI = [
            f.indexed_fragment()
            for f in self.genome.fragments.all()
            if f.name == "chrI"
        ][0]
        self.assertEquals(len(chrI.annotations()), 1)
        self.assertEquals(chrI.annotations()[0].feature.name, "cds")


class JoinImporterTest(TestCase):

    def test_import_gff_CDS_subfragments(self):

        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tgene\t30\t80\t.\t+\t.\tID=i2g;Name=f2g
chrI\tTest\tCDS\t30\t80\t.\t+\t.\tID=i2;Name=f2
chrI\tTest\tCDS\t30\t41\t.\t+\t.\tParent=i2
chrI\tTest\tCDS\t50\t55\t.\t+\t.\tParent=i2
chrI\tTest\tCDS\t60\t80\t.\t+\t.\tParent=i2
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
            f.write(data)
            f.close()
            genome = Genome.import_gff("Foo", f.name)
            os.unlink(f.name)

        # created one fragment for each sequence in GFF file
        self.assertCountEqual(
            [fr.name for fr in genome.fragments.all()], ["chrI"]
        )
        chrI = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrI"
        ][0]
        self.assertEqual(len(chrI.sequence), 160)
        self.assertEqual(len(chrI.annotations()), 4)
        self.assertEqual(chrI.annotations()[0].base_first, 30)
        self.assertEqual(chrI.annotations()[0].base_last, 80)
        self.assertEqual(chrI.annotations()[0].feature.name, 'f2g')
        self.assertEqual(chrI.annotations()[0].feature.type, 'gene')
        self.assertEqual(chrI.annotations()[1].base_first, 30)
        self.assertEqual(chrI.annotations()[1].base_last, 41)
        self.assertEqual(chrI.annotations()[1].feature.name, 'f2')
        self.assertEqual(chrI.annotations()[1].feature.type, 'CDS')
        self.assertEqual(chrI.annotations()[2].base_first, 50)
        self.assertEqual(chrI.annotations()[2].base_last, 55)
        self.assertEqual(chrI.annotations()[2].feature.name, 'f2')
        self.assertEqual(chrI.annotations()[2].feature.type, 'CDS')
        self.assertEqual(chrI.annotations()[3].base_first, 60)
        self.assertEqual(chrI.annotations()[3].base_last, 80)
        self.assertEqual(chrI.annotations()[3].feature.name, 'f2')
        self.assertEqual(chrI.annotations()[3].feature.type, 'CDS')

    def test_import_gff_CDS_subfragments_overlap(self):
        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tgene\t30\t80\t.\t+\t.\tID=i2g;Name=f2g
chrI\tTest\tCDS\t30\t80\t.\t+\t.\tID=i2;Name=f2
chrI\tTest\tCDS\t30\t41\t.\t+\t.\tParent=i2
chrI\tTest\tCDS\t41\t50\t.\t+\t.\tParent=i2
chrI\tTest\tCDS\t56\t61\t.\t+\t.\tParent=i2
chrI\tTest\tCDS\t60\t80\t.\t+\t.\tParent=i2
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
            f.write(data)
            f.close()
            genome = Genome.import_gff("Foo", f.name)
            os.unlink(f.name)

        # created one fragment for each sequence in GFF file
        self.assertCountEqual(
            [fr.name for fr in genome.fragments.all()], ["chrI"]
        )
        chrI = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrI"
        ][0]
        self.assertEqual(len(chrI.sequence), 160)
        self.assertEqual(len(chrI.annotations()), 5)
        self.assertEqual(chrI.annotations()[0].base_first, 30)
        self.assertEqual(chrI.annotations()[0].base_last, 80)
        self.assertEqual(chrI.annotations()[0].feature.name, 'f2g')
        self.assertEqual(chrI.annotations()[0].feature.type, 'gene')
        self.assertNotEqual(chrI.annotations()[0].feature.id, chrI.annotations()[1].feature.id)
        self.assertEqual(chrI.annotations()[1].base_first, 30)
        self.assertEqual(chrI.annotations()[1].base_last, 41)
        self.assertEqual(chrI.annotations()[1].feature.name, 'f2')
        self.assertEqual(chrI.annotations()[1].feature.type, 'CDS')
        self.assertEqual(chrI.annotations()[2].base_first, 41)
        self.assertEqual(chrI.annotations()[2].base_last, 50)
        self.assertEqual(chrI.annotations()[2].feature.name, 'f2')
        self.assertEqual(chrI.annotations()[2].feature.type, 'CDS')
        self.assertEqual(chrI.annotations()[3].base_first, 56)
        self.assertEqual(chrI.annotations()[3].base_last, 61)
        self.assertEqual(chrI.annotations()[3].feature.name, 'f2')
        self.assertEqual(chrI.annotations()[3].feature.type, 'CDS')
        self.assertEqual(chrI.annotations()[4].base_first, 60)
        self.assertEqual(chrI.annotations()[4].base_last, 80)
        self.assertEqual(chrI.annotations()[4].feature.name, 'f2')
        self.assertEqual(chrI.annotations()[4].feature.type, 'CDS')
        self.assertEqual(chrI.annotations()[1].feature.id, chrI.annotations()[2].feature.id)
        self.assertEqual(chrI.annotations()[2].feature.id, chrI.annotations()[3].feature.id)
        self.assertEqual(chrI.annotations()[3].feature.id, chrI.annotations()[4].feature.id)

    @mock.patch("edge.models.fragment.Annotation.from_chunk_feature_and_location_array")
    def test_import_gff_CDS_subfragments_overlap_check_chunks(self, cf_fcl_mock):
        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tgene\t30\t80\t.\t+\t.\tID=i2g;Name=f2g
chrI\tTest\tCDS\t30\t80\t.\t+\t.\tID=i2;Name=f2
chrI\tTest\tCDS\t30\t41\t.\t+\t.\tParent=i2
chrI\tTest\tCDS\t41\t50\t.\t+\t.\tParent=i2
chrI\tTest\tCDS\t56\t61\t.\t+\t.\tParent=i2
chrI\tTest\tCDS\t60\t80\t.\t+\t.\tParent=i2
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
            f.write(data)
            f.close()
            genome = Genome.import_gff("Foo", f.name)
            os.unlink(f.name)

        # created one fragment for each sequence in GFF file
        self.assertCountEqual(
            [fr.name for fr in genome.fragments.all()], ["chrI"]
        )
        chrI = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrI"
        ][0]
        chrI.annotations()
        gene_cfs = []
        cds_cfs = []
        args, _ = cf_fcl_mock.call_args
        for cf, fcl in args[0]:
            if cf.feature.type == 'gene':
                gene_cfs.append(cf)
            elif cf.feature.type == 'CDS':
                cds_cfs.append(cf)
        gene_chunk_starts = sorted([cf.fcl_base_first for cf in gene_cfs])
        gene_chunk_ends = sorted([cf.fcl_base_last for cf in gene_cfs])
        cds_chunk_starts = sorted([cf.fcl_base_first for cf in cds_cfs])
        cds_chunk_ends = sorted([cf.fcl_base_last for cf in cds_cfs])

        self.assertEqual(len(chrI.sequence), 160)
        self.assertEqual(gene_chunk_starts, [30, 41, 42, 51, 56, 60, 62])
        self.assertEqual(gene_chunk_ends, [40, 41, 50, 55, 59, 61, 80])
        self.assertEqual(cds_chunk_starts, [30, 41, 41, 42, 56, 60, 60, 62])
        self.assertEqual(cds_chunk_ends, [40, 41, 41, 50, 59, 61, 61, 80])

    def test_import_gff_CDS_subfragments_SGD_CDS(self):
        data = """##gff-version 3
chrI\tSGD\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tSGD\tgene\t20\t60\t.\t+\t.\tID=A1;Name=A1;gene=gene_A1
chrI\tSGD\tCDS\t20\t37\t.\t+\t0\tParent=A1_mRNA;Name=A1_CDS;orf_classification=Verified
chrI\tSGD\tintron\t38\t39\t.\t+\t.\tParent=A1_mRNA;Name=A1_intron;orf_classification=Verified
chrI\tSGD\tCDS\t40\t60\t.\t+\t0\tParent=A1_mRNA;Name=A1_CDS;orf_classification=Verified
chrI\tSGD\tmRNA\t20\t60\t.\t+\t.\tID=A1_mRNA;Name=A1_mRNA;Parent=A1
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
            f.write(data)
            f.close()
            genome = Genome.import_gff("Foo", f.name)
            os.unlink(f.name)

        # created one fragment for each sequence in GFF file
        self.assertCountEqual(
            [fr.name for fr in genome.fragments.all()], ["chrI"]
        )
        chrI = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrI"
        ][0]
        self.assertEqual(len(chrI.sequence), 160)
        self.assertEqual(len(chrI.annotations()), 5)
        self.assertEqual(chrI.annotations()[0].base_first, 20)
        self.assertEqual(chrI.annotations()[0].base_last, 60)
        self.assertEqual(chrI.annotations()[0].feature.name, 'A1')
        self.assertEqual(chrI.annotations()[0].feature.type, 'gene')
        self.assertNotEqual(chrI.annotations()[0].feature.id, chrI.annotations()[1].feature.id)
        self.assertNotEqual(chrI.annotations()[0].feature.id, chrI.annotations()[2].feature.id)
        self.assertNotEqual(chrI.annotations()[0].feature.id, chrI.annotations()[3].feature.id)
        self.assertEqual(chrI.annotations()[1].base_first, 20)
        self.assertEqual(chrI.annotations()[1].base_last, 37)
        self.assertEqual(chrI.annotations()[1].feature.name, 'A1_CDS')
        self.assertEqual(chrI.annotations()[1].feature.type, 'CDS')
        self.assertNotEqual(chrI.annotations()[1].feature.id, chrI.annotations()[2].feature.id)
        self.assertNotEqual(chrI.annotations()[1].feature.id, chrI.annotations()[3].feature.id)
        self.assertEqual(chrI.annotations()[2].base_first, 20)
        self.assertEqual(chrI.annotations()[2].base_last, 60)
        self.assertEqual(chrI.annotations()[2].feature.name, 'A1_mRNA')
        self.assertEqual(chrI.annotations()[2].feature.type, 'mRNA')
        self.assertNotEqual(chrI.annotations()[2].feature.id, chrI.annotations()[3].feature.id)
        self.assertEqual(chrI.annotations()[3].base_first, 38)
        self.assertEqual(chrI.annotations()[3].base_last, 39)
        self.assertEqual(chrI.annotations()[3].feature.name, 'A1_intron')
        self.assertEqual(chrI.annotations()[3].feature.type, 'intron')
        self.assertEqual(chrI.annotations()[4].base_first, 40)
        self.assertEqual(chrI.annotations()[4].base_last, 60)
        self.assertEqual(chrI.annotations()[4].feature.name, 'A1_CDS')
        self.assertEqual(chrI.annotations()[4].feature.type, 'CDS')
        self.assertEqual(chrI.annotations()[1].feature.id, chrI.annotations()[4].feature.id)

    def test_import_subfeatures_simple_reverse_coordinates(self):

        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tgene\t20\t65\t.\t-\t.\tID=i2;Name=f2
chrI\tTest\tcds\t20\t28\t.\t-\t.\tParent=i2;Name=f2_cds
chrI\tTest\tcds\t58\t65\t.\t-\t.\tParent=i2;Name=f2_cds
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
            f.write(data)
            f.close()
            genome = Genome.import_gff("Foo", f.name)
            os.unlink(f.name)

        # verify chrI fragment
        chrI = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrI"
        ][0]
        self.assertEqual(len(chrI.sequence), 160)
        self.assertEqual(len(chrI.annotations()), 3)
        self.assertEqual(chrI.annotations()[0].base_first, 20)
        self.assertEqual(chrI.annotations()[0].base_last, 65)
        self.assertEqual(chrI.annotations()[0].feature_base_first, 1)
        self.assertEqual(chrI.annotations()[0].feature_base_last, 46)
        self.assertEqual(chrI.annotations()[0].feature.name, "f2")
        self.assertEquals(chrI.annotations()[0].feature.strand, -1)
        self.assertEqual(chrI.annotations()[1].base_first, 20)
        self.assertEqual(chrI.annotations()[1].base_last, 28)
        self.assertEqual(chrI.annotations()[1].feature_base_first, 9)
        self.assertEqual(chrI.annotations()[1].feature_base_last, 17)
        self.assertEqual(chrI.annotations()[1].feature.name, "f2_cds")
        self.assertEquals(chrI.annotations()[1].feature.strand, -1)
        self.assertEqual(chrI.annotations()[2].base_first, 58)
        self.assertEqual(chrI.annotations()[2].base_last, 65)
        self.assertEqual(chrI.annotations()[2].feature_base_first, 1)
        self.assertEqual(chrI.annotations()[2].feature_base_last, 8)
        self.assertEqual(chrI.annotations()[2].feature.name, "f2_cds")
        self.assertEquals(chrI.annotations()[2].feature.strand, -1)

    def test_import_subfeatures_overlap_reverse_coordinates(self):

        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
chrI\tTest\tgene\t20\t65\t.\t-\t.\tID=i2;Name=f2
chrI\tTest\tcds\t20\t28\t.\t-\t.\tParent=i2;Name=f2_cds
chrI\tTest\tintron\t29\t40\t.\t-\t.\tParent=i2;Name=f2_intron
chrI\tTest\tcds\t41\t60\t.\t-\t.\tParent=i2;Name=f2_cds
chrI\tTest\tcds\t58\t65\t.\t-\t.\tParent=i2;Name=f2_cds
###
##FASTA
>chrI
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTAAC
ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
"""

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
            f.write(data)
            f.close()
            genome = Genome.import_gff("Foo", f.name)
            os.unlink(f.name)

        # verify chrI fragment
        chrI = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrI"
        ][0]
        self.assertEqual(len(chrI.sequence), 160)
        self.assertEqual(len(chrI.annotations()), 5)
        self.assertEqual(chrI.annotations()[0].base_first, 20)
        self.assertEqual(chrI.annotations()[0].base_last, 65)
        self.assertEqual(chrI.annotations()[0].feature_base_first, 1)
        self.assertEqual(chrI.annotations()[0].feature_base_last, 46)
        self.assertEqual(chrI.annotations()[0].feature.name, "f2")
        self.assertEquals(chrI.annotations()[0].feature.strand, -1)
        self.assertEqual(chrI.annotations()[1].base_first, 20)
        self.assertEqual(chrI.annotations()[1].base_last, 28)
        self.assertEqual(chrI.annotations()[1].feature_base_first, 29)
        self.assertEqual(chrI.annotations()[1].feature_base_last, 37)
        self.assertEqual(chrI.annotations()[1].feature.name, "f2_cds")
        self.assertEquals(chrI.annotations()[1].feature.strand, -1)
        self.assertEqual(chrI.annotations()[2].base_first, 29)
        self.assertEqual(chrI.annotations()[2].base_last, 40)
        self.assertEqual(chrI.annotations()[2].feature_base_first, 1)
        self.assertEqual(chrI.annotations()[2].feature_base_last, 12)
        self.assertEqual(chrI.annotations()[2].feature.name, "f2_intron")
        self.assertEquals(chrI.annotations()[2].feature.strand, -1)
        self.assertEqual(chrI.annotations()[3].base_first, 41)
        self.assertEqual(chrI.annotations()[3].base_last, 60)
        self.assertEqual(chrI.annotations()[3].feature_base_first, 9)
        self.assertEqual(chrI.annotations()[3].feature_base_last, 28)
        self.assertEqual(chrI.annotations()[3].feature.name, "f2_cds")
        self.assertEquals(chrI.annotations()[3].feature.strand, -1)
        self.assertEqual(chrI.annotations()[4].base_first, 58)
        self.assertEqual(chrI.annotations()[4].base_last, 65)
        self.assertEqual(chrI.annotations()[4].feature_base_first, 1)
        self.assertEqual(chrI.annotations()[4].feature_base_last, 8)
        self.assertEqual(chrI.annotations()[4].feature.name, "f2_cds")
        self.assertEquals(chrI.annotations()[4].feature.strand, -1)
