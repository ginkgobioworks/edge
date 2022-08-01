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
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "t"*100 + Sites.lox71 + "c"*1000
        )
        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        donor = "a"*10 + Sites.lox66 + "g"*100

        res = self.client.post(
            "/edge/genomes/%s/ssr/" % parent_genome.id,
            data=json.dumps(
                dict(donor=donor, is_donor_circular=True, reaction="crelox", create=False)
            ),
            content_type="application/json",
        )
        self.assertEquals(res.status_code, 200)
        r = json.loads(res.content)

        self.assertEquals(len(r), 1)
        self.assertEquals(r[0]["recombination"]["type"], "Integration")
        self.assertEquals(len(r[0]["genomic_locations"]), 1)
        self.assertEquals(r[0]["genomic_locations"][0]["fragment_id"], f.id)
        self.assertEquals(r[0]["genomic_locations"][0]["start"], 101)

    def test_create_genome_api_without_insert(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            Sites.loxP + "t"*100 + Sites.loxP + "c"*1000
        )
        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        res = self.client.post(
            "/edge/genomes/%s/ssr/" % parent_genome.id,
            data=json.dumps(
                dict(donor=None, is_donor_circular=None, reaction="crelox",
                     genome_name="foo foo bar bar", create=True)
            ),
            content_type="application/json",
        )
        self.assertEquals(res.status_code, 201)
        r = json.loads(res.content)

        child = Genome.objects.get(pk=r["id"])
        self.assertEquals(child.name, "foo foo bar bar")
        modified_sequence = child.fragments.all()[0].indexed_fragment().sequence
        self.assertEquals("t"*100 not in modified_sequence, True)

    def test_create_genome_api_with_insert(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "t"*100 + Sites.lox71 + "c"*1000
        )
        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        donor = "a"*10 + Sites.lox66 + "g"*100

        res = self.client.post(
            "/edge/genomes/%s/ssr/" % parent_genome.id,
            data=json.dumps(
                dict(donor=donor, is_donor_circular=True, reaction="crelox",
                     genome_name="foo foo bar bar", create=True)
            ),
            content_type="application/json",
        )
        self.assertEquals(res.status_code, 201)
        r = json.loads(res.content)

        child = Genome.objects.get(pk=r["id"])
        self.assertEquals(child.name, "foo foo bar bar")
        modified_sequence = child.fragments.all()[0].indexed_fragment().sequence
        self.assertEquals("g"*100 in modified_sequence, True)
