from django.test import TestCase

from edge import import_gff
from edge.models import Genome

class JoinImporterTest(TestCase):
    def test_import_join_cds(self):
        """
        U00096.3	feature	gene	57364	58179	.	+	.	db_xref=EcoGene:EG11570;gene=djlA;gene_synonym=ECK0056,JW0054,yabH;locus_tag=b0055
        U00096.3	feature	CDS	57364	58179	.	+	0	codon_start=1;db_xref=GI:1786241,ASAP:ABE-0000187,UniProtKB/Swiss-Prot:P31680,EcoGene:EG11570;function=phenotype%3B Not classified;gene=djlA;gene_synonym=ECK0056,JW0054,yabH;locus_tag=b0055;note=GO_process: GO:0006457 - protein folding;product=DnaJ-like protein%2C membrane anchored;protein_id=AAC73166.1;transl_table=11;translation=MQYWGKIIGVAVALLMGGGFWGVVLGLLIGHMFDKARSRKMAWFANQRERQALFFATTFEVMGHLTKSKGRVTEADIHIASQLMDRMNLHGASRTAAQNAFRVGKSDNYPLREKMRQFRSVCFGRFDLIRMFLEIQIQAAFADGSLHPNERAVLYVIAEELGISRAQFDQFLRMMQGGAQFGGGYQQQTGGGNWQQAQRGPTLEDACNVLGVKPTDDATTIKRAYRKLMSEHHPDKLVAKGLPPEMMEMAKQKAQEIQQAYELIKQQKGFK
        U00096.3	feature	gene	58474	59269	.	+	.	db_xref=EcoGene:EG12610;gene=yabP;gene_synonym=ECK0057,JW0055,yabQ;locus_tag=b4659;pseudo=
        U00096.3	feature	CDS	58474	59269	.	+	0	ID=biopygen1;codon_start=1;db_xref=ASAP:ABE-0000192,ASAP:ABE-0000194,UniProtKB/Swiss-Prot:P39220,EcoGene:EG12610;gene=yabP;gene_synonym=ECK0057,JW0055,yabQ;locus_tag=b4659;note=pseudogene;pseudo=;transl_table=11
        U00096.3	feature	CDS	58474	59052	.	+	0	Parent=biopygen1
        U00096.3	feature	CDS	59052	59228	.	+	0	Parent=biopygen1
        U00096.3	feature	CDS	59228	59269	.	+	0	Parent=biopygen1

        If the row in GFF is of type "CDS" and is a joint of 2 fragments, you would have one Feature and two Chunk_Feature rows.
        --> CDS 58474-59269 should have a Feature of (58474-59269) and 3 Chunk_Features of (58474-59052, 59052-59228, 59228-59269)
        """

        import_gff("TestGenome", "example/test_join.gff")
        genome = Genome.objects.get(name="TestGenome")
        U00096_3 = [
            fr.indexed_fragment() for fr in genome.fragments.all() if fr.name == "U00096.3"
        ][0]
        self.assertEqual(1, 2)