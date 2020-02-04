import os
import json

from Bio.Seq import Seq
from django.test import TestCase

import edge.recombine
from edge.recombine import find_swap_region, find_root_genome, recombine, remove_overhangs
from edge.blastdb import build_all_genome_dbs, fragment_fasta_fn
from edge.models import Genome, Fragment, Genome_Fragment, Operation


class RemoveOverhangsTest(TestCase):
    def test_removes_front_overhang(self):
        self.assertEquals(remove_overhangs('(atg/)aa'), 'aa')

    def test_removes_back_overhang(self):
        self.assertEquals(remove_overhangs('aa(atg/)'), 'aa')

    def test_removes_front_and_back(self):
        self.assertEquals(remove_overhangs('(atg/)aa(atg/)'), 'aa')

    def test_does_not_remove_internal_overhang(self):
        self.assertEquals(remove_overhangs('(atg/)a(atg/)a(atg/)'), 'a(atg/)a')

    def test_does_not_remove_unclosed_overhang(self):
        self.assertEquals(remove_overhangs('(atg/aa'), '(atg/aa')
        self.assertEquals(remove_overhangs('atg/aa)'), 'atg/aa)')

    def test_works_with_single_char_input(self):
        self.assertEquals(remove_overhangs(')'), ')')
        self.assertEquals(remove_overhangs('('), '(')


class GenomeRecombinationTest(TestCase):
    def setUp(self):
        self.old_check_junction_lu = edge.recombine.CHECK_JUNCTION_LEFT_UP
        self.old_check_junction_ld = edge.recombine.CHECK_JUNCTION_LEFT_DN
        self.old_check_junction_ru = edge.recombine.CHECK_JUNCTION_RIGHT_UP
        self.old_check_junction_rd = edge.recombine.CHECK_JUNCTION_RIGHT_DN

        edge.recombine.CHECK_JUNCTION_LEFT_UP = 10
        edge.recombine.CHECK_JUNCTION_LEFT_DN = 40
        edge.recombine.CHECK_JUNCTION_RIGHT_UP = 40
        edge.recombine.CHECK_JUNCTION_RIGHT_DN = 10

        self.old_single_cross_over_gap_max = edge.recombine.SINGLE_CROSSOVER_MAX_GAP
        edge.recombine.SINGLE_CROSSOVER_MAX_GAP = 10

    def tearDown(self):
        edge.recombine.CHECK_JUNCTION_LEFT_UP = self.old_check_junction_lu
        edge.recombine.CHECK_JUNCTION_LEFT_DN = self.old_check_junction_ld
        edge.recombine.CHECK_JUNCTION_RIGHT_UP = self.old_check_junction_ru
        edge.recombine.CHECK_JUNCTION_RIGHT_DN = self.old_check_junction_rd
        edge.recombine.SINGLE_CROSSOVER_MAX_GAP = self.old_single_cross_over_gap_max

    def build_genome(self, circular, *templates):
        g = Genome(name='Foo')
        g.save()
        for seq in templates:
            f = Fragment.create_with_sequence('Bar', seq, circular=circular)
            Genome_Fragment(genome=g, fragment=f, inherited=False).save()
            try:
                os.unlink(fragment_fasta_fn(f))
            except OSError:
                pass
        build_all_genome_dbs(refresh=True)
        return Genome.objects.get(pk=g.id)

    def build_genome_and_ancestors(self, names):
        # build a genome of name in names with each as its subsquent genome's parent
        g_parent = Genome(name=names[0])
        g_parent.save()
        g_child = g_parent
        for name in names[:-1]:
            g_child = Genome(name=name)
            g_child.parent = g_parent
            g_child.save()
            g_parent = g_child
        return g_child

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
        g = self.build_genome(False, template)
        r = find_swap_region(g, cassette, arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(r[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(r[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(r[0].start, len(upstream) + 1)
        self.assertEquals(r[0].end, len(template) - len(downstream))
        self.assertEquals(r[0].sequence, ''.join([front_bs, middle, back_bs]))
        self.assertEquals(r[0].cassette_reversed, False)
        self.assertEquals(r[0].front_arm, front_bs[0:arm_len])
        self.assertEquals(r[0].back_arm, back_bs[-arm_len:])
        self.assertEquals(r[0].is_double_crossover, True)

    def test_finding_swap_region_across_circular_boundary(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        cassette = ''.join([front_bs, replaced, back_bs])
        arm_len = min(len(front_bs), len(back_bs))

        template = ''.join([middle[8:], back_bs, downstream, upstream, front_bs, middle[0:8]])
        g = self.build_genome(True, template)

        r = find_swap_region(g, cassette, arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(r[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(r[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(r[0].start, len(template) - 8 - len(front_bs) + 1)
        self.assertEquals(r[0].end, len(middle) - 8 + len(back_bs))
        self.assertEquals(r[0].sequence, ''.join([front_bs, middle, back_bs]))
        self.assertEquals(r[0].cassette_reversed, False)
        self.assertEquals(r[0].front_arm, front_bs[0:arm_len])
        self.assertEquals(r[0].back_arm, back_bs[-arm_len:])

    def test_finding_swap_region_when_front_arm_is_across_circular_boundary(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        cassette = ''.join([front_bs, replaced, back_bs])
        arm_len = min(len(front_bs), len(back_bs))

        template = ''.join([front_bs[8:], middle, back_bs, downstream, upstream, front_bs[0:8]])
        g = self.build_genome(True, template)

        r = find_swap_region(g, cassette, arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(r[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(r[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(r[0].start, len(template) - 8 + 1)
        self.assertEquals(r[0].end, len(front_bs + middle + back_bs) - 8)
        self.assertEquals(r[0].sequence, ''.join([front_bs, middle, back_bs]))
        self.assertEquals(r[0].cassette_reversed, False)
        self.assertEquals(r[0].front_arm, front_bs[0:arm_len])
        self.assertEquals(r[0].back_arm, back_bs[-arm_len:])

    def test_finding_swap_region_when_back_arm_is_across_circular_boundary(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        cassette = ''.join([front_bs, replaced, back_bs])
        arm_len = min(len(front_bs), len(back_bs))

        template = ''.join([back_bs[8:], downstream, upstream, front_bs, middle, back_bs[0:8]])
        g = self.build_genome(True, template)

        r = find_swap_region(g, cassette, arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(r[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(r[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(r[0].start, len(back_bs) - 8 + len(downstream + upstream) + 1)
        self.assertEquals(r[0].end, len(back_bs) - 8)
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

        cassette = str(Seq(''.join([front_bs, replaced, back_bs])).reverse_complement())
        arm_len = min(len(front_bs), len(back_bs))

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        g = self.build_genome(False, template)

        r = find_swap_region(g, cassette, arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(r[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(r[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(r[0].start, len(upstream) + 1)
        self.assertEquals(r[0].end, len(template) - len(downstream))
        self.assertEquals(r[0].sequence, ''.join([front_bs, middle, back_bs]))
        self.assertEquals(r[0].cassette_reversed, True)
        self.assertEquals(r[0].front_arm, str(Seq(back_bs[-arm_len:]).reverse_complement()))
        self.assertEquals(r[0].back_arm, str(Seq(front_bs[0:arm_len]).reverse_complement()))

    def test_finding_reverse_complement_region_across_circular_boundary(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        cassette = str(Seq(''.join([front_bs, replaced, back_bs])).reverse_complement())
        arm_len = min(len(front_bs), len(back_bs))

        template = ''.join([middle[8:], back_bs, downstream, upstream, front_bs, middle[0:8]])
        g = self.build_genome(True, template)

        r = find_swap_region(g, cassette, arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(r[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(r[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(r[0].start, len(template) - 8 - len(front_bs) + 1)
        self.assertEquals(r[0].end, len(middle) - 8 + len(back_bs))
        self.assertEquals(r[0].sequence, ''.join([front_bs, middle, back_bs]))
        self.assertEquals(r[0].cassette_reversed, True)
        self.assertEquals(r[0].front_arm, str(Seq(back_bs[-arm_len:]).reverse_complement()))
        self.assertEquals(r[0].back_arm, str(Seq(front_bs[0:arm_len]).reverse_complement()))

    def test_find_root_genome(self):
        genome = self.build_genome_and_ancestors(['Foo', 'Bar', 'Foobar'])
        root_genome = find_root_genome(genome)
        self.assertEqual(root_genome.name, 'Foo')

    def test_finding_reverse_complement_region_when_front_arm_is_across_circular_boundary(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        cassette = str(Seq(''.join([front_bs, replaced, back_bs])).reverse_complement())
        arm_len = min(len(front_bs), len(back_bs))

        template = ''.join([front_bs[8:], middle, back_bs, downstream, upstream, front_bs[0:8]])
        g = self.build_genome(True, template)

        r = find_swap_region(g, cassette, arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(r[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(r[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(r[0].start, len(template) - 8 + 1)
        self.assertEquals(r[0].end, len(front_bs + middle + back_bs) - 8)
        self.assertEquals(r[0].sequence, ''.join([front_bs, middle, back_bs]))
        self.assertEquals(r[0].cassette_reversed, True)
        self.assertEquals(r[0].front_arm, str(Seq(back_bs[-arm_len:]).reverse_complement()))
        self.assertEquals(r[0].back_arm, str(Seq(front_bs[0:arm_len]).reverse_complement()))

    def test_finding_reverse_complement_region_when_back_arm_is_across_circular_boundary(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        cassette = str(Seq(''.join([front_bs, replaced, back_bs])).reverse_complement())
        arm_len = min(len(front_bs), len(back_bs))

        template = ''.join([back_bs[8:], downstream, upstream, front_bs, middle, back_bs[0:8]])
        g = self.build_genome(True, template)

        r = find_swap_region(g, cassette, arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(r[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(r[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(r[0].start, len(back_bs) - 8 + len(downstream + upstream) + 1)
        self.assertEquals(r[0].end, len(back_bs) - 8)
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
        replaced = "a" * 100

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = ''.join([front_bs, replaced, back_bs])

        arm_len = min(len(front_bs), len(back_bs))
        g = self.build_genome(False, template)
        c = recombine(g, cassette, arm_len)

        self.assertNotEqual(g.id, c.id)
        self.assertEquals(c.fragments.all()[0].indexed_fragment().sequence,
                          ''.join([upstream, cassette, downstream]))

    def test_recombines_ignoring_extra_bases_upstream_and_downstream_of_cassette(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggagcgacgtagtctgcatctgatgcatgcactac"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagggcatcgtactactgatgcatgcacactgacgta"
        downstream = "gttaaggcgcgaacat"
        replaced = "a" * 100

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = ''.join(['c' * 6 + front_bs, replaced, back_bs + 'c' * 6])

        arm_len = int(min(len(front_bs), len(back_bs)) / 2)
        g = self.build_genome(False, template)
        c = recombine(g, cassette, arm_len)

        self.assertNotEqual(g.id, c.id)
        self.assertEquals(c.fragments.all()[0].indexed_fragment().sequence,
                          ''.join([upstream, front_bs, replaced, back_bs, downstream]))

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
        g = self.build_genome(False, template)
        self.assertEquals(Operation.objects.count(), 0)
        c = recombine(g, cassette, arm_len)
        self.assertEquals(Operation.objects.count(), 1)
        self.assertEquals(c.operation_set.all()[0].type, Operation.RECOMBINATION[0])
        self.assertEquals(c.operation_set.all()[0].params,
                          json.dumps(dict(cassette=cassette, homology_arm_length=arm_len)))

    def test_annotates_cassette(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = ''.join([front_bs, replaced, back_bs])
        arm_len = min(len(front_bs), len(back_bs))
        g = self.build_genome(False, template)

        a = g.fragments.all()[0].indexed_fragment().annotations()
        self.assertEquals(len(a), 0)

        c = recombine(g, cassette, arm_len)

        a = c.fragments.all()[0].indexed_fragment().annotations()
        self.assertEquals(len(a), 1)
        self.assertEquals(a[0].base_first, len(upstream) + 1)
        self.assertEquals(a[0].base_last, len(upstream + cassette))
        self.assertEquals(a[0].feature_base_first, 1)
        self.assertEquals(a[0].feature_base_last, len(cassette))
        self.assertEquals(a[0].feature.strand, 1)
        self.assertEquals(a[0].feature.operation.type, Operation.RECOMBINATION[0])
        self.assertEquals(a[0].feature.operation.genome, c)

    def test_annotates_reversed_cassette(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = str(Seq(''.join([front_bs, replaced, back_bs])).reverse_complement())
        arm_len = min(len(front_bs), len(back_bs))
        g = self.build_genome(False, template)

        a = g.fragments.all()[0].indexed_fragment().annotations()
        self.assertEquals(len(a), 0)

        c = recombine(g, cassette, arm_len)

        a = c.fragments.all()[0].indexed_fragment().annotations()
        self.assertEquals(len(a), 1)
        self.assertEquals(a[0].base_first, len(upstream) + 1)
        self.assertEquals(a[0].base_last, len(upstream + cassette))
        self.assertEquals(a[0].feature_base_first, 1)
        self.assertEquals(a[0].feature_base_last, len(cassette))
        # on reverse strand
        self.assertEquals(a[0].feature.strand, -1)
        self.assertEquals(a[0].feature.operation.type, Operation.RECOMBINATION[0])
        self.assertEquals(a[0].feature.operation.genome, c)

    def test_integrates_and_annotates_cassette_across_circular_boundary(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        cassette = ''.join([front_bs, replaced, back_bs])
        arm_len = min(len(front_bs), len(back_bs))

        template = ''.join([middle[8:], back_bs, downstream, upstream, front_bs, middle[0:8]])
        g = self.build_genome(True, template)
        c = recombine(g, cassette, arm_len)

        self.assertNotEqual(g.id, c.id)
        self.assertEquals(c.fragments.all()[0].indexed_fragment().sequence,
                          ''.join([downstream, upstream, cassette]))

        a = c.fragments.all()[0].indexed_fragment().annotations()
        self.assertEquals(len(a), 1)
        self.assertEquals(a[0].base_first, len(downstream + upstream) + 1)
        self.assertEquals(a[0].base_last, len(downstream + upstream + cassette))
        self.assertEquals(a[0].feature_base_first, 1)
        self.assertEquals(a[0].feature_base_last, len(cassette))
        self.assertEquals(a[0].feature.strand, 1)
        self.assertEquals(a[0].feature.operation.type, Operation.RECOMBINATION[0])
        self.assertEquals(a[0].feature.operation.genome, c)

    def test_recombine_when_front_arm_is_across_circular_boundary(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        cassette = ''.join([front_bs, replaced, back_bs])
        arm_len = min(len(front_bs), len(back_bs))

        template = ''.join([front_bs[8:], middle, back_bs, downstream, upstream, front_bs[0:8]])
        g = self.build_genome(True, template)
        c = recombine(g, cassette, arm_len)

        self.assertNotEqual(g.id, c.id)
        self.assertEquals(c.fragments.all()[0].indexed_fragment().sequence,
                          ''.join([downstream, upstream, cassette]))

    def test_recombine_when_back_arm_is_across_circular_boundary(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        cassette = ''.join([front_bs, replaced, back_bs])
        arm_len = min(len(front_bs), len(back_bs))

        template = ''.join([back_bs[8:], downstream, upstream, front_bs, middle, back_bs[0:8]])
        g = self.build_genome(True, template)
        c = recombine(g, cassette, arm_len)

        self.assertNotEqual(g.id, c.id)
        self.assertEquals(c.fragments.all()[0].indexed_fragment().sequence,
                          ''.join([downstream, upstream, cassette]))

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
        g = self.build_genome(False, template)
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
        g = self.build_genome(False, template)

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
        self.assertEquals(r[0]['start'], len(upstream) + 1)
        self.assertEquals(r[0]['end'], len(template) - len(downstream))
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
        g = self.build_genome(False, template)

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

    def test_recombines_multiple_times_on_different_fragments(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaaaaa"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = ''.join([front_bs, replaced, back_bs])

        f1 = 't' * 20 + template + 'c' * 20 + template + 'c' * 30
        f2 = 't' * 40 + template + 'c' * 15 + template + 'c' * 20

        arm_len = min(len(front_bs), len(back_bs))
        g = self.build_genome(False, f1, f2)
        self.assertEquals(g.fragments.count(), 2)

        c = recombine(g, cassette, arm_len)
        self.assertEquals(c.fragments.count(), 2)

        sequences = [f.indexed_fragment().sequence for f in c.fragments.all()]
        sequences = sorted(sequences, key=lambda s: len(s))

        self.assertEquals(sequences[0],
                          't' * 20 + upstream + cassette + downstream
                          + 'c' * 20 + upstream + cassette + downstream + 'c' * 30)
        self.assertEquals(sequences[1],
                          't' * 40 + upstream + cassette + downstream
                          + 'c' * 15 + upstream + cassette + downstream + 'c' * 20)

    def test_multiple_recombinations_on_multiple_fragments(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaaaaa"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = ''.join([front_bs, replaced, back_bs])

        f1 = 't' * 20 + template + 'c' * 20 + template + 'c' * 30
        f2 = 't' * 40 + template + 'c' * 15 + template + 'c' * 20

        arm_len = min(len(front_bs), len(back_bs))
        g = self.build_genome(False, f1, f2)
        self.assertEquals(g.fragments.count(), 2)

        # first recombination on g, producing first child c
        c1 = recombine(g, cassette, arm_len)
        self.assertEquals(c1.fragments.count(), 2)
        sequences = [f.indexed_fragment().sequence for f in c1.fragments.all()]
        sequences = sorted(sequences, key=lambda s: len(s))
        self.assertEquals(sequences[0],
                          't' * 20 + upstream + cassette + downstream
                          + 'c' * 20 + upstream + cassette + downstream + 'c' * 30)
        self.assertEquals(sequences[1],
                          't' * 40 + upstream + cassette + downstream
                          + 'c' * 15 + upstream + cassette + downstream + 'c' * 20)

        # check original genome sequences are correct
        g = Genome.objects.get(pk=g.id)
        self.assertEquals(g.fragments.count(), 2)
        sequences = [f.indexed_fragment().sequence for f in g.fragments.all()]
        sequences = sorted(sequences, key=lambda s: len(s))
        self.assertEquals(sequences[0], f1)
        self.assertEquals(sequences[1], f2)

        # second recombination on g, producing second child c2
        c2 = recombine(g, cassette, arm_len)
        self.assertEquals(c2.fragments.count(), 2)
        sequences = [f.indexed_fragment().sequence for f in c2.fragments.all()]
        sequences = sorted(sequences, key=lambda s: len(s))
        self.assertEquals(sequences[0],
                          't' * 20 + upstream + cassette + downstream
                          + 'c' * 20 + upstream + cassette + downstream + 'c' * 30)
        self.assertEquals(sequences[1],
                          't' * 40 + upstream + cassette + downstream
                          + 'c' * 15 + upstream + cassette + downstream + 'c' * 20)

        # check second recombination has not affected the results of the first recombination
        c1 = Genome.objects.get(pk=c1.id)
        self.assertEquals(c1.fragments.count(), 2)
        sequences = [f.indexed_fragment().sequence for f in c1.fragments.all()]
        sequences = sorted(sequences, key=lambda s: len(s))
        self.assertEquals(sequences[0],
                          't' * 20 + upstream + cassette + downstream
                          + 'c' * 20 + upstream + cassette + downstream + 'c' * 30)
        self.assertEquals(sequences[1],
                          't' * 40 + upstream + cassette + downstream
                          + 'c' * 15 + upstream + cassette + downstream + 'c' * 20)

    def test_recombines_multiple_times_on_circular_fragment(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaaaaa"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = ''.join([front_bs, replaced, back_bs])

        f = (
            middle[0:8]
            + back_bs
            + downstream
            + 't' * 20
            + template
            + 'c' * 20
            + upstream
            + front_bs
            + middle[8:]
        )
        arm_len = min(len(front_bs), len(back_bs))
        g = self.build_genome(True, f)
        c = recombine(g, cassette, arm_len)
        self.assertEquals(c.fragments.all()[0].indexed_fragment().sequence,
                          downstream + 't' * 20 + upstream + cassette + downstream
                          + 'c' * 20 + upstream + cassette)

    def test_multiple_recombines_return_same_child(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = ''.join([front_bs, replaced, back_bs])

        arm_len = min(len(front_bs), len(back_bs))
        g = self.build_genome(False, template)

        res = self.client.post('/edge/genomes/%s/recombination/' % g.id,
                               data=json.dumps(dict(cassette=cassette,
                                                    homology_arm_length=arm_len,
                                                    create=True,
                                                    genome_name='FooBar')),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        r = json.loads(res.content)
        c1 = r['id']

        res = self.client.post('/edge/genomes/%s/recombination/' % g.id,
                               data=json.dumps(dict(cassette=cassette,
                                                    homology_arm_length=arm_len,
                                                    create=True,
                                                    genome_name='FooBar')),
                               content_type='application/json')
        # returns 200 not 201
        self.assertEquals(res.status_code, 200)
        r = json.loads(res.content)
        c2 = r['id']
        self.assertEquals(c1, c2)

    def test_multiple_recombines_return_active_child(self):
        upstream = "gagattgtccgcgtttt"
        front_bs = "catagcgcacaggacgcggag"
        middle = "cggcacctgtgagccg"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacat"
        replaced = "aaaaaaaaaaaaaaaaaaa"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = ''.join([front_bs, replaced, back_bs])

        arm_len = min(len(front_bs), len(back_bs))
        g = self.build_genome(False, template)

        res = self.client.post('/edge/genomes/%s/recombination/' % g.id,
                               data=json.dumps(dict(cassette=cassette,
                                                    homology_arm_length=arm_len,
                                                    create=True,
                                                    genome_name='FooBar')),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        r = json.loads(res.content)
        c = Genome.objects.get(pk=r['id'])
        c.active = False
        c.save()

        res = self.client.post('/edge/genomes/%s/recombination/' % g.id,
                               data=json.dumps(dict(cassette=cassette,
                                                    homology_arm_length=arm_len,
                                                    create=True,
                                                    genome_name='FooBar')),
                               content_type='application/json')
        self.assertEquals(res.status_code, 200)
        r = json.loads(res.content)
        c = Genome.objects.get(pk=r['id'])
        self.assertEquals(c.active, True)

    def __test_verification_primers(self, template, middle, cassette, arm_len, is_reversed):
        from edge.pcr import pcr_from_genome

        g = self.build_genome(False, template)
        r = find_swap_region(g, cassette, arm_len, design_primers=True)

        self.assertEquals(len(r), 1)
        self.assertEquals(len(r[0].verification_cassette), 5)
        self.assertEquals(len(r[0].verification_front), 5)
        self.assertEquals(len(r[0].verification_back), 5)

        # cassette verification primers should work on unmodified genome
        for primer in r[0].verification_cassette:
            p = pcr_from_genome(g, primer['PRIMER_LEFT_SEQUENCE'], primer['PRIMER_RIGHT_SEQUENCE'])
            self.assertNotEqual(p[0], None)
            self.assertEquals(p[0].index(middle) >= 0, True)

        # front verification primers should NOT produce product
        for primer in r[0].verification_front:
            p = pcr_from_genome(g, primer['PRIMER_LEFT_SEQUENCE'], primer['PRIMER_RIGHT_SEQUENCE'])
            self.assertEqual(p[0], None)

        # back verification primers should NOT produce product
        for primer in r[0].verification_back:
            p = pcr_from_genome(g, primer['PRIMER_LEFT_SEQUENCE'], primer['PRIMER_RIGHT_SEQUENCE'])
            self.assertEqual(p[0], None)

        # do recombination, then try primers again on modified genome

        c = recombine(g, cassette, arm_len)

        for f in c.fragments.all():
            try:
                os.unlink(fragment_fasta_fn(f))
            except OSError:
                pass
        build_all_genome_dbs(refresh=True)
        # reload to get blastdb
        c = Genome.objects.get(pk=c.id)

        if is_reversed:
            cassette = str(Seq(cassette).reverse_complement())

        # cassette verification primers should work on modified genome, finding cassette
        for primer in r[0].verification_cassette:
            p = pcr_from_genome(c, primer['PRIMER_LEFT_SEQUENCE'], primer['PRIMER_RIGHT_SEQUENCE'])
            self.assertNotEqual(p[0], None)
            self.assertEquals(p[0].index(cassette) >= 0, True)

        # front verification primers should find a product including front of cassette
        for primer in r[0].verification_front:
            p = pcr_from_genome(c, primer['PRIMER_LEFT_SEQUENCE'], primer['PRIMER_RIGHT_SEQUENCE'])
            self.assertNotEqual(p[0], None)
            self.assertEquals(p[0].index(cassette[0:edge.recombine.CHECK_JUNCTION_LEFT_DN]) >= 0,
                              True)

        # back verification primers should find a product including back of cassette
        for primer in r[0].verification_back:
            p = pcr_from_genome(c, primer['PRIMER_LEFT_SEQUENCE'], primer['PRIMER_RIGHT_SEQUENCE'])
            self.assertNotEqual(p[0], None)
            self.assertEquals(p[0].index(cassette[-edge.recombine.CHECK_JUNCTION_RIGHT_UP:]) >= 0,
                              True)

    def test_finds_verification_primers_for_swap_region(self):
        upstream = "gagattgtccgcgttttagctgatacgtacgtgtcgatcgacttgcgtatctgatcatctgacgtagat"
        front_bs = "catagcgcacaggacgcggag"
        middle = "ccagtcgctgaggcagtcgatgcaggcatcgatcaggctggcacctgtgagccgagctcacgtatgcatcatcattga"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacatagctagtactagtcacgtagtcatttgtcgtacgtacgtattgagtcatca"
        replaced = "gactacgatcagtcgtagtaacgcgtagcgtagtcagcgtacacgtacgtagacgacgtacatgcatcgtactgtatc"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = ''.join([front_bs, replaced, back_bs])
        arm_len = min(len(front_bs), len(back_bs))

        self.__test_verification_primers(template, middle, cassette, arm_len, False)

    def test_finds_verification_primers_with_reverse_complement_cassette(self):
        upstream = "gagattgtccgcgttttagctgatacgtacgtgtcgatcgacttgcgtatctgatcatctgacgtagat"
        front_bs = "catagcgcacaggacgcggag"
        middle = "ccagtcgctgaggcagtcgatgcaggcatcgatcaggctggcacctgtgagccgagctcacgtatgcatcatcattga"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacatagctagtactagtcacgtagtcatttgtcgtacgtacgtattgagtcatca"
        replaced = "gactacgatcagtcgtagtaacgcgtagcgtagtcagcgtacacgtacgtagacgacgtacatgcatcgtactgtatc"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = str(Seq(''.join([front_bs, replaced, back_bs])).reverse_complement())
        arm_len = min(len(front_bs), len(back_bs))

        self.__test_verification_primers(template, middle, cassette, arm_len, True)

    def test_does_not_return_front_back_junction_verification_primers_if_they_are_not_useful(self):
        upstream = "gagattgtccgcgttttagctgatacgtacgtgtcgatcgacttgcgtatctgatcatctgacgtagat"
        front_bs = "catagcgcacaggacgcggag"
        middle = "ccagtcgctgaggcagtcgatgcaggcatcgatcaggctggcacctgtgagccgagctcacgtatgcatcatcattga"
        back_bs = "taatgaccccgaagcagg"
        downstream = "gttaaggcgcgaacatagctagtactagtcacgtagtcatttgtcgtacgtacgtattgagtcatca"
        replaced = "gactacgatcagtcgtagtaacgcgtagcgtagtcagcgtacacgtacgtagacgacgtacatgcatcgtactgtatc"

        template = ''.join([upstream, front_bs, middle, back_bs, downstream])
        cassette = ''.join([front_bs, replaced, back_bs])
        arm_len = min(len(front_bs), len(back_bs))

        g = self.build_genome(False, template)
        r = find_swap_region(g, cassette, arm_len, design_primers=True)

        self.assertEquals(len(r), 1)
        self.assertEquals(len(r[0].verification_cassette), 5)
        self.assertEquals(len(r[0].verification_front), 5)
        self.assertEquals(len(r[0].verification_back), 5)

        # cassette same as region to be replaced, front and back verification
        # primers are not useful
        cassette = ''.join([front_bs, middle, back_bs])
        r = find_swap_region(g, cassette, arm_len, design_primers=True)

        self.assertEquals(len(r), 1)
        self.assertEquals(len(r[0].verification_cassette), 5)
        self.assertEquals(len(r[0].verification_front), 0)
        self.assertEquals(len(r[0].verification_back), 0)


class SingleCrossoverTest(TestCase):
    def setUp(self):
        self.old_check_junction_lu = edge.recombine.CHECK_JUNCTION_LEFT_UP
        self.old_check_junction_ld = edge.recombine.CHECK_JUNCTION_LEFT_DN
        self.old_check_junction_ru = edge.recombine.CHECK_JUNCTION_RIGHT_UP
        self.old_check_junction_rd = edge.recombine.CHECK_JUNCTION_RIGHT_DN

        edge.recombine.CHECK_JUNCTION_LEFT_UP = 10
        edge.recombine.CHECK_JUNCTION_LEFT_DN = 40
        edge.recombine.CHECK_JUNCTION_RIGHT_UP = 40
        edge.recombine.CHECK_JUNCTION_RIGHT_DN = 10

        self.old_single_cross_over_gap_max = edge.recombine.SINGLE_CROSSOVER_MAX_GAP
        self.new_max_gap = 20
        edge.recombine.SINGLE_CROSSOVER_MAX_GAP = self.new_max_gap

    def tearDown(self):
        edge.recombine.CHECK_JUNCTION_LEFT_UP = self.old_check_junction_lu
        edge.recombine.CHECK_JUNCTION_LEFT_DN = self.old_check_junction_ld
        edge.recombine.CHECK_JUNCTION_RIGHT_UP = self.old_check_junction_ru
        edge.recombine.CHECK_JUNCTION_RIGHT_DN = self.old_check_junction_rd
        edge.recombine.SINGLE_CROSSOVER_MAX_GAP = self.old_single_cross_over_gap_max

    def build_genome(self, circular, *templates):
        g = Genome(name='Foo')
        g.save()
        for seq in templates:
            f = Fragment.create_with_sequence('Bar', seq, circular=circular)
            Genome_Fragment(genome=g, fragment=f, inherited=False).save()
            try:
                os.unlink(fragment_fasta_fn(f))
            except OSError:
                pass
        build_all_genome_dbs(refresh=True)
        return Genome.objects.get(pk=g.id)

    def test_finds_single_crossover_region(self):
        upstream = 'gagattgtccgcgtttt'
        locus = 'catagcgcacaggacgcggagtaggcgtagtcggttgatctgatgtc'
        downstream = 'gttaaggcgcgaacat'
        insertion = 'aaaaaaaaaaaaaaaaaaa'
        locus_len = len(locus)
        bs_len = int(locus_len / 2)

        template = ''.join([upstream, locus, downstream])
        cassette = ''.join([locus[locus_len - bs_len:], insertion, locus[0:bs_len]])

        g = self.build_genome(False, template)
        arm_len = bs_len - 2
        r = find_swap_region(g, cassette, arm_len)
        self.assertEquals(len(r), 1)
        self.assertEquals(r[0].fragment_id, g.fragments.all()[0].id)
        self.assertEquals(r[0].fragment_name, g.fragments.all()[0].name)
        self.assertEquals(r[0].is_double_crossover, False)
        self.assertEquals(r[0].start, 20)
        self.assertEquals(r[0].end, 62)

    def test_single_crossover_integrates_correctly(self):
        upstream = 'gagattgtccgcgtttt'
        locus = 'catagcgcacaggacgcggagtaggcgtagtcggttgatctgatgtc'
        downstream = 'gttaaggcgcgaacat'
        insertion = 'aaaaaaaaaaaaaaaaaaa'
        locus_len = len(locus)
        bs_len = int(locus_len / 2)

        template = ''.join([upstream, locus, downstream])
        cassette = ''.join([locus[locus_len - bs_len:], insertion, locus[0:bs_len]])

        g = self.build_genome(False, template)
        c = recombine(g, cassette, bs_len - 2)

        self.assertNotEqual(g.id, c.id)
        self.assertEquals(c.fragments.all()[0].indexed_fragment().sequence,
                          ''.join([upstream, locus, insertion, locus, downstream]))

    def test_single_crossover_integrates_correctly_with_gap_in_homology(self):
        upstream = 'gagattgtccgcgtttt'
        locus = 'catagcgcacaggacgcggagtaggcgtagtcggttgatctgatgtc'
        downstream = 'gttaaggcgcgaacat'
        insertion = 'aaaaaaaaaaaaaaaaaaa'
        locus_len = len(locus)
        gap = int(self.new_max_gap / 2)
        arm_short = 2
        bs_len = int(locus_len / 2) - (gap - arm_short)

        template = ''.join([upstream, locus, downstream])
        cassette = ''.join([locus[locus_len - bs_len:], insertion, locus[0:bs_len]])

        g = self.build_genome(False, template)
        c = recombine(g, cassette, bs_len - arm_short)

        self.assertNotEqual(g.id, c.id)
        self.assertEquals(c.fragments.all()[0].indexed_fragment().sequence,
                          ''.join([upstream, locus, insertion, locus, downstream]))

    def test_single_crossover_integrates_correctly_with_reverse_complement_of_locus(self):
        upstream = 'gagattgtccgcgtttt'
        locus = 'catagcgcacaggacgcggagtaggcgtagtcggttgatctgatgtc'
        downstream = 'gttaaggcgcgaacat'
        insertion = 'aaaaaaaaaaaaaaaaaaa'
        locus_len = len(locus)
        bs_len = int(locus_len / 2)

        template = ''.join([upstream, locus, downstream])
        cassette = ''.join([locus[locus_len - bs_len:], insertion, locus[0:bs_len]])
        cassette = str(Seq(cassette).reverse_complement())

        g = self.build_genome(False, template)
        c = recombine(g, cassette, bs_len - 2)

        self.assertNotEqual(g.id, c.id)
        self.assertEquals(c.fragments.all()[0].indexed_fragment().sequence,
                          ''.join([upstream, locus, insertion, locus, downstream]))

    def test_single_crossover_integrates_correctly_with_reverse_complement_and_gap_in_locus(self):
        upstream = 'gagattgtccgcgtttt'
        locus = 'catagcgcacaggacgcggagtaggcgtagtcggttgatctgatgtc'
        downstream = 'gttaaggcgcgaacat'
        insertion = 'aaaaaaaaaaaaaaaaaaa'
        locus_len = len(locus)
        gap = int(self.new_max_gap / 2)
        arm_short = 2
        bs_len = int(locus_len / 2) - (gap - arm_short)

        template = ''.join([upstream, locus, downstream])
        cassette = ''.join([locus[locus_len - bs_len:], insertion, locus[0:bs_len]])
        cassette = str(Seq(cassette).reverse_complement())

        g = self.build_genome(False, template)
        c = recombine(g, cassette, bs_len - arm_short)

        self.assertNotEqual(g.id, c.id)
        self.assertEquals(c.fragments.all()[0].indexed_fragment().sequence,
                          ''.join([upstream, locus, insertion, locus, downstream]))
