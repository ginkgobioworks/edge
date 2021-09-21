import os, tempfile
import unittest.mock as mock

from django.test import TestCase

from edge import import_gff
from edge.models import Genome

class JoinImporterTest(TestCase):
    """
    Example
    -------
    U00096.3	feature	gene	57364	58179	.	+	.	db_xref=EcoGene:EG11570;gene=djlA;gene_synonym=ECK0056,JW0054,yabH;locus_tag=b0055
    U00096.3	feature	CDS	57364	58179	.	+	0	codon_start=1;db_xref=GI:1786241,ASAP:ABE-0000187,UniProtKB/Swiss-Prot:P31680,EcoGene:EG11570;function=phenotype%3B Not classified;gene=djlA;gene_synonym=ECK0056,JW0054,yabH;locus_tag=b0055;note=GO_process: GO:0006457 - protein folding;product=DnaJ-like protein%2C membrane anchored;protein_id=AAC73166.1;transl_table=11;translation=MQYWGKIIGVAVALLMGGGFWGVVLGLLIGHMFDKARSRKMAWFANQRERQALFFATTFEVMGHLTKSKGRVTEADIHIASQLMDRMNLHGASRTAAQNAFRVGKSDNYPLREKMRQFRSVCFGRFDLIRMFLEIQIQAAFADGSLHPNERAVLYVIAEELGISRAQFDQFLRMMQGGAQFGGGYQQQTGGGNWQQAQRGPTLEDACNVLGVKPTDDATTIKRAYRKLMSEHHPDKLVAKGLPPEMMEMAKQKAQEIQQAYELIKQQKGFK
    U00096.3	feature	gene	58474	59269	.	+	.	db_xref=EcoGene:EG12610;gene=yabP;gene_synonym=ECK0057,JW0055,yabQ;locus_tag=b4659;pseudo=
    U00096.3	feature	CDS	58474	59269	.	+	0	ID=biopygen1;codon_start=1;db_xref=ASAP:ABE-0000192,ASAP:ABE-0000194,UniProtKB/Swiss-Prot:P39220,EcoGene:EG12610;gene=yabP;gene_synonym=ECK0057,JW0055,yabQ;locus_tag=b4659;note=pseudogene;pseudo=;transl_table=11
    U00096.3	feature	CDS	58474	59052	.	+	0	Parent=biopygen1
    U00096.3	feature	CDS	59228	59269	.	+	0	Parent=biopygen1
    """

    def test_import_gff_CDS_subfragments(self):

        data = """##gff-version 3
chrI\tTest\tchromosome\t1\t160\t.\t.\t.\tID=i1;Name=f1
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
        self.assertEquals(len(chrI.sequence), 160)
        self.assertEquals(len(chrI.annotations()), 1)
        self.assertEquals(chrI.annotations()[0].base_first, 30)
        self.assertEquals(chrI.annotations()[0].base_last, 80)
        self.assertEquals(chrI.annotations()[0].feature.name, 'f2')

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
        self.assertEquals(len(chrI.sequence), 160)
        self.assertEquals(len(chrI.annotations()), 2)
        self.assertEquals(chrI.annotations()[0].base_first, 30)
        self.assertEquals(chrI.annotations()[0].base_last, 80)
        self.assertEquals(chrI.annotations()[1].base_first, 30)
        self.assertEquals(chrI.annotations()[1].base_last, 80)
        self.assertEquals(chrI.annotations()[0].feature.name, 'f2g')
        self.assertEquals(chrI.annotations()[1].feature.name, 'f2')
        self.assertEquals(chrI.annotations()[0].feature.id, 2)
        self.assertEquals(chrI.annotations()[1].feature.id, 3)


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
chrII\tTest\tchromosome\t1\t100\t.\t.\t.\tID=i3;Name=f3
chrII\tTest\trbs\t50\t70\t.\t+\t.\tID=i4r;Name=f4r
chrII\tTest\tCDS\t20\t90\t.\t+\t.\tID=i4;Name=f4
chrII\tTest\tCDS\t20\t45\t.\t+\t.\tParent=i4
chrII\tTest\tCDS\t80\t90\t.\t+\t.\tParent=i4
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
            [fr.name for fr in genome.fragments.all()], ["chrI", "chrII"]
        )
        chrI = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrI"
        ][0]
        
        chrI.annotations()
        gene_cfs = []
        cds_cfs = []
        args, _ = cf_fcl_mock.call_args
        for cf, _ in args[0]:
            if cf.feature_id == 4:
                gene_cfs.append(cf)
            elif cf.feature_id == 5:
                cds_cfs.append(cf)
        gene_starts = sorted([cf.fcl_base_first for cf in gene_cfs])
        gene_ends = sorted([cf.fcl_base_last for cf in gene_cfs])
        cds_starts = sorted([cf.fcl_base_first for cf in cds_cfs])
        cds_ends = sorted([cf.fcl_base_last for cf in cds_cfs])
    
        self.assertEquals(len(chrI.sequence), 160)
        self.assertEquals(gene_starts, [30, 41, 51, 56, 60])
        self.assertEquals(gene_ends, [41, 50, 55, 61, 80])
        self.assertEquals(cds_starts, [30, 41, 56, 60])
        self.assertEquals(cds_ends, [41, 50, 61, 80])

        chrII = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "chrII"
        ][0]
        chrII.annotations()
        rbs_cfs = []
        cds_cfs = []
        args, _ = cf_fcl_mock.call_args
        for cf, fcl in args[0]:
            if cf.feature_id == 6:
                rbs_cfs.append(cf)
            elif cf.feature_id == 7:
                cds_cfs.append(cf)
        rbs_starts = sorted([cf.fcl_base_first for cf in rbs_cfs])
        rbs_ends = sorted([cf.fcl_base_last for cf in rbs_cfs])
        cds_starts = sorted([cf.fcl_base_first for cf in cds_cfs])
        cds_ends = sorted([cf.fcl_base_last for cf in cds_cfs])
    
        self.assertEquals(len(chrII.sequence), 100)
        self.assertEquals(rbs_starts, [50])
        self.assertEquals(rbs_ends, [70])
        self.assertEquals(cds_starts, [20, 80])
        self.assertEquals(cds_ends, [45, 90])