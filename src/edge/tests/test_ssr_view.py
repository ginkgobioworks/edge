import os
import json

from django.test import TestCase
from edge.ssr.crelox import Sites
from edge.models import Genome, Fragment, Genome_Fragment


class SSRViewTest(TestCase):

    def test_check_api_without_insert(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            Sites.loxP + "t"*100 + Sites.loxP + "c"*1000
        )
        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        res = self.client.post(
            "/edge/genomes/%s/ssr/" % parent_genome.id,
            data=json.dumps(
                dict(donor=None, is_donor_circular=None, reaction="crelox", create=False)
            ),
            content_type="application/json",
        )
        self.assertEquals(res.status_code, 200)
        r = json.loads(res.content)

        self.assertEquals(len(r), 1)
        self.assertEquals(r[0]["recombination"]["type"], "Excision")
        self.assertEquals(len(r[0]["genomic_locations"]), 2)
        self.assertEquals(r[0]["genomic_locations"][0]["fragment_id"], f.id)
        self.assertEquals(r[0]["genomic_locations"][0]["start"], 1)
        self.assertEquals(r[0]["genomic_locations"][1]["fragment_id"], f.id)
        self.assertEquals(r[0]["genomic_locations"][1]["start"], len(Sites.loxP)+100+1)

    def test_check_api_with_insert(self):
        assert False

    def test_create_genome_api_without_insert(self):
        assert False

    def test_create_genome_api_with_insert(self):
        assert False

    def test_unsupported_reaction(self):
        assert False
