import os
import json
from Bio.Seq import Seq
from django.test import TestCase
from edge.recombine import find_swap_region, recombine
from edge.models import Genome, Fragment, Genome_Fragment, Operation
from edge.blastdb import build_all_genome_dbs, fragment_fasta_fn


class GenomeRecombinationTest(TestCase):

    def build_genome(self, *templates):
        g = Genome(name='Foo')
        g.save()
        for seq in templates:
            f = Fragment.create_with_sequence('Bar', seq)
            Genome_Fragment(genome=g, fragment=f, inherited=False).save()
            try:
                os.unlink(fragment_fasta_fn(f))
            except:
                pass
        build_all_genome_dbs(refresh=True)
        return g

    def test_finds_correct_region_for_swapping(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = ''.join([front_bs, replaced, back_bs])

        arm_len = min(len(front_bs), len(back_bs))
        g = self.build_genome(template)
        r = find_swap_region(g, cassette, arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(r[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(r[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(r[0].start, len(upstream)+1)
        self.assertEquals(r[0].end, len(template)-len(downstream))
        self.assertEquals(r[0].sequence, ''.join([front_bs, middle, back_bs]))
        self.assertEquals(r[0].cassette_reversed, False)
        self.assertEquals(r[0].front_arm, front_bs[0:arm_len])
        self.assertEquals(r[0].back_arm, back_bs[-arm_len:])

    def test_finds_correct_region_for_swapping_with_reverse_complement_cassette(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = str(Seq(''.join([front_bs, replaced, back_bs])).reverse_complement())

        arm_len = min(len(front_bs), len(back_bs))
        g = self.build_genome(template)
        r = find_swap_region(g, cassette, arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(r[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(r[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(r[0].start, len(upstream)+1)
        self.assertEquals(r[0].end, len(template)-len(downstream))
        self.assertEquals(r[0].sequence, ''.join([front_bs, middle, back_bs]))
        self.assertEquals(r[0].cassette_reversed, True)
        self.assertEquals(r[0].front_arm, str(Seq(back_bs[-arm_len:]).reverse_complement()))
        self.assertEquals(r[0].back_arm, str(Seq(front_bs[0:arm_len]).reverse_complement()))

    def test_recombines_correctly(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = ''.join([front_bs, replaced, back_bs])

        arm_len = min(len(front_bs), len(back_bs))
        g = self.build_genome(template)
        c = recombine(g, cassette, arm_len)

        self.assertNotEqual(g.id, c.id)
        self.assertEquals(c.fragments.all()[0].indexed_fragment().sequence,
                          ''.join([upstream, cassette, downstream]))

    def test_creates_operation(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = ''.join([front_bs, replaced, back_bs])

        arm_len = min(len(front_bs), len(back_bs))
        g = self.build_genome(template)
        self.assertEquals(Operation.objects.count(), 0)
        c = recombine(g, cassette, arm_len)
        self.assertEquals(Operation.objects.count(), 1)
        self.assertEquals(c.operations.all()[0].type, Operation.RECOMBINATION[0])
        self.assertEquals(c.operations.all()[0].params,
                          json.dumps(dict(cassette=cassette, homology_arm_length=arm_len)))

    def test_recombines_with_reverse_complement_cassette_correctly(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = str(Seq(''.join([front_bs, replaced, back_bs])).reverse_complement())

        arm_len = min(len(front_bs), len(back_bs))
        g = self.build_genome(template)
        c = recombine(g, cassette, arm_len)

        self.assertNotEqual(g.id, c.id)
        self.assertEquals(c.fragments.all()[0].indexed_fragment().sequence,
                          ''.join([upstream, front_bs, replaced, back_bs, downstream]))

    def test_find_swap_region_api(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = ''.join([front_bs, replaced, back_bs])

        arm_len = min(len(front_bs), len(back_bs))
        g = self.build_genome(template)

        res = self.client.post('/edge/genomes/%s/recombination/' % g.id,
                               data=json.dumps(dict(cassette=cassette,
                                                    homology_arm_length=arm_len,
                                                    create=False)),
                               content_type='application/json')
        self.assertEquals(res.status_code, 200)
        r = json.loads(res.content)
        self.assertEquals(len(r), 1)
        self.assertEquals(r[0]['fragment_id'], g.fragments.all()[0].id)
        self.assertEquals(r[0]['fragment_name'], g.fragments.all()[0].name)
        self.assertEquals(r[0]['start'], len(upstream)+1)
        self.assertEquals(r[0]['end'], len(template)-len(downstream))
        self.assertEquals(r[0]['sequence'], ''.join([front_bs, middle, back_bs]))
        self.assertEquals(r[0]['cassette_reversed'], False)
        self.assertEquals(r[0]['front_arm'], front_bs[0:arm_len])
        self.assertEquals(r[0]['back_arm'], back_bs[-arm_len:])

    def test_recombination_api(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = ''.join([front_bs, replaced, back_bs])

        arm_len = min(len(front_bs), len(back_bs))
        g = self.build_genome(template)

        res = self.client.post('/edge/genomes/%s/recombination/' % g.id,
                               data=json.dumps(dict(cassette=cassette,
                                                    homology_arm_length=arm_len,
                                                    create=True,
                                                    genome_name='FooBar')),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        r = json.loads(res.content)
        self.assertEquals(r['name'], 'FooBar')

        c = Genome.objects.get(pk=r['id'])
        self.assertEquals(c.fragments.all()[0].indexed_fragment().sequence,
                          ''.join([upstream, cassette, downstream]))
