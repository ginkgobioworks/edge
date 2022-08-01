import os
import json

from Bio.Seq import Seq
from django.test import TestCase

from edge.models import Genome, Fragment, Genome_Fragment
from edge.ssr import (
  rc,
  find_query_locations_on_duplicated_template,
  find_site_locations_on_sequence,
  Event,
  Excision,
  Inversion,
  Integration,
  RMCE,
  SiteLocation,
  Recombination,
  Reaction,
)


class FindQueryOnSequenceTest(TestCase):
    def test_finds_locations_on_sequence(self):
        site = "aggc"
        template = "a"*100+site+"g"*100

        locs = find_query_locations_on_duplicated_template(template, len(template), site)
        self.assertEquals(len(locs), 1)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site)
        self.assertEquals(locs[0].start_0based, 100)

    def test_mark_location_as_on_fragment(self):
        site = "aggc"
        template = "a"*100+site+"g"*100

        locs = find_query_locations_on_duplicated_template(template, len(template), site, fragment_obj=object())
        self.assertEquals(len(locs), 1)
        self.assertEquals(locs[0].on_insert, False)

    def test_finds_multiple_locations(self):
        site = "aggc"
        template = "a"*100+site+"g"*100+site+"c"*100

        locs = find_query_locations_on_duplicated_template(template, len(template), site)
        self.assertEquals(len(locs), 2)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site)
        self.assertEquals(locs[0].start_0based, 100)
        self.assertEquals(locs[1].on_insert, True)
        self.assertEquals(locs[1].site, site)
        self.assertEquals(locs[1].start_0based, 204)

    def test_ignores_starting_position_after_sequence_length_argument(self):
        site = "aggc"
        template = "a"*100+site+"g"*100+site+"c"*100

        locs = find_query_locations_on_duplicated_template(template, len(template), site)
        self.assertEquals(len(locs), 2)

        locs = find_query_locations_on_duplicated_template(template, 102, site)
        self.assertEquals(len(locs), 1)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site)
        self.assertEquals(locs[0].start_0based, 100)


class FindLocationsTest(TestCase):
    def test_finds_locations_on_sequence(self):
        site = "aggc"
        template = "a"*100+site+"g"*100

        locs = find_site_locations_on_sequence(template, False, [site])
        self.assertEquals(len(locs), 1)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site)
        self.assertEquals(locs[0].start_0based, 100)

    def test_finds_multiple_sites_on_sequence(self):
        site1 = "aggc"
        site2 = "tggc"
        template = "a"*100+site1+"g"*100+site2+"t"*100

        locs = find_site_locations_on_sequence(template, False, [site1, site2])
        self.assertEquals(len(locs), 2)
        locs = sorted(locs, key=lambda l: l.site)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site1)
        self.assertEquals(locs[0].start_0based, 100)
        self.assertEquals(locs[1].on_insert, True)
        self.assertEquals(locs[1].site, site2)
        self.assertEquals(locs[1].start_0based, 204)

    def test_finds_reverse_complement_locations_on_sequence(self):
        site1 = "aggc"
        site2 = "tggc"
        self.assertEquals(rc(site2), "gcca")
        template = "a"*100+site1+"g"*100+rc(site2)+"t"*100

        locs = find_site_locations_on_sequence(template, False, [site1, site2])
        self.assertEquals(len(locs), 2)
        locs = sorted(locs, key=lambda l: l.site)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site1)
        self.assertEquals(locs[0].start_0based, 100)
        self.assertEquals(locs[1].on_insert, True)
        self.assertEquals(locs[1].site, rc(site2))
        self.assertEquals(locs[1].start_0based, 204)

    def test_finds_forward_location_across_circular_boundary(self):
        site1 = "aggc"
        template = site1[2:]+"g"*100+site1[:2]

        # cannot find anything if sequence is not circular
        locs = find_site_locations_on_sequence(template, False, [site1])
        self.assertEquals(len(locs), 0)

        locs = find_site_locations_on_sequence(template, True, [site1])
        self.assertEquals(len(locs), 1)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site1)
        self.assertEquals(locs[0].start_0based, 102)

    def test_finds_reverse_location_across_circular_boundary(self):
        site1 = "aggc"
        template = rc(site1)[2:]+"g"*100+rc(site1)[:2]

        # cannot find anything if sequence is not circular
        locs = find_site_locations_on_sequence(template, False, [site1])
        self.assertEquals(len(locs), 0)

        locs = find_site_locations_on_sequence(template, True, [site1])
        self.assertEquals(len(locs), 1)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, rc(site1))
        self.assertEquals(locs[0].start_0based, 102)


class ReactionTest(TestCase):

    def test_gathers_locations_on_genome(self):
        class FakeFragment(object):
            def __init__(self, n, s, c):
                self.id = n
                self.sequence = s
                self.circular = c

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc"]

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [FakeRecombination()]

        f = FakeFragment(14, "c"*100+"aggc"+"t"*100, False)
        r = FakeReaction(None, None, None, [f])
        r.determine_site_locations()

        self.assertEquals(len(r.locations), 1)
        self.assertEquals(r.locations[0].on_insert, False)
        self.assertEquals(r.locations[0].fragment_id, 14)
        self.assertEquals(r.locations[0].site, "aggc")
        self.assertEquals(r.locations[0].start_0based, 100)

    def test_gathers_locations_on_multiple_fragments(self):
        class FakeFragment(object):
            def __init__(self, n, s, c):
                self.id = n
                self.sequence = s
                self.circular = c

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc"]

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [FakeRecombination()]

        f1 = FakeFragment(14, "c"*100+"aggc"+"t"*100, False)
        f2 = FakeFragment(15, "c"*33+"aggc"+"t"*100, False)
        r = FakeReaction(None, None, None, [f1, f2])
        r.determine_site_locations()

        self.assertEquals(len(r.locations), 2)
        self.assertEquals(r.locations[0].on_insert, False)
        self.assertEquals(r.locations[0].fragment_id, 14)
        self.assertEquals(r.locations[0].site, "aggc")
        self.assertEquals(r.locations[0].start_0based, 100)
        self.assertEquals(r.locations[1].on_insert, False)
        self.assertEquals(r.locations[1].fragment_id, 15)
        self.assertEquals(r.locations[1].site, "aggc")
        self.assertEquals(r.locations[1].start_0based, 33)

    def test_gathers_locations_on_insert_and_genome(self):
        class FakeFragment(object):
            def __init__(self, n, s, c):
                self.id = n
                self.sequence = s
                self.circular = c

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc"]

            def required_insert_sites(self):
                return ["ttag"]

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [FakeRecombination()]

        f = FakeFragment(14, "c"*100+"aggc"+"t"*100, False)
        r = FakeReaction(None, "g"*33+"ttag"+"c"*100, True, [f])
        r.determine_site_locations()

        self.assertEquals(len(r.locations), 2)
        self.assertEquals(r.locations[0].on_insert, True)
        self.assertEquals(r.locations[0].fragment_id, None)
        self.assertEquals(r.locations[0].site, "ttag")
        self.assertEquals(r.locations[0].start_0based, 33)
        self.assertEquals(r.locations[1].on_insert, False)
        self.assertEquals(r.locations[1].fragment_id, 14)
        self.assertEquals(r.locations[1].site, "aggc")
        self.assertEquals(r.locations[1].start_0based, 100)

    def test_group_locations_into_multiple_events(self):
        class FakeFragment(object):
            def __init__(self, n, s, c):
                self.id = n
                self.sequence = s
                self.circular = c

        class FakeEvent(Event):
            pass

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc"]

            def required_insert_sites(self):
                return ["ttag"]

            def events(self, locations, errors):
                return self.generate_events(locations, FakeEvent, errors)

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [FakeRecombination()]

        f = FakeFragment(14, "c"*100+"aggc"+"t"*100+"aggc"+"g"*100, False)
        r = FakeReaction(None, "g"*33+"ttag"+"c"*100, True, [f])
        r.group_into_events()

        self.assertEquals(len(r.events), 2)

        self.assertEquals(type(r.events[0]), FakeEvent)
        locations = sorted(r.events[0].locations, key=lambda l: l.on_insert)
        self.assertEquals(len(locations), 2)
        self.assertEquals(locations[0].on_insert, False)
        self.assertEquals(locations[0].fragment_id, 14)
        self.assertEquals(locations[0].site, "aggc")
        self.assertEquals(locations[0].start_0based, 100)
        self.assertEquals(locations[1].on_insert, True)
        self.assertEquals(locations[1].fragment_id, None)
        self.assertEquals(locations[1].site, "ttag")
        self.assertEquals(locations[1].start_0based, 33)

        self.assertEquals(type(r.events[1]), FakeEvent)
        locations = sorted(r.events[1].locations, key=lambda l: l.on_insert)
        self.assertEquals(len(locations), 2)
        self.assertEquals(locations[0].on_insert, False)
        self.assertEquals(locations[0].fragment_id, 14)
        self.assertEquals(locations[0].site, "aggc")
        self.assertEquals(locations[0].start_0based, 204)
        self.assertEquals(locations[1].on_insert, True)
        self.assertEquals(locations[1].fragment_id, None)
        self.assertEquals(locations[1].site, "ttag")
        self.assertEquals(locations[1].start_0based, 33)


class RecombinationTest(TestCase):

    def test_recognizes_matched_locations_on_genome(self):

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc"]

        class FakeFragment(object):
            def __init__(self, n):
                self.id = n

        locations = [
            SiteLocation("atgc", FakeFragment(11), None, 8),
            SiteLocation("aggc", FakeFragment(11), None, 12),
        ]

        matched = FakeRecombination().possible_locations(locations, [])
        self.assertEquals(len(matched), 1)
        self.assertEquals(len(matched[0]), 1)
        self.assertEquals(matched[0][0].site, "aggc")

    def test_recognizes_matched_locations_in_reverse(self):

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["gcct"]

        class FakeFragment(object):
            def __init__(self, n):
                self.id = n

        locations = [
            SiteLocation("atgc", FakeFragment(11), None, 8),
            SiteLocation("aggc", FakeFragment(11), None, 12),
        ]

        matched = FakeRecombination().possible_locations(locations, [])
        self.assertEquals(len(matched), 1)
        self.assertEquals(len(matched[0]), 1)
        self.assertEquals(matched[0][0].site, "aggc")

    def test_recognizes_multiple_matched_locations(self):

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc"]

        class FakeFragment(object):
            def __init__(self, n):
                self.id = n

        locations = [
            SiteLocation("atgc", FakeFragment(11), None, 8),
            SiteLocation("aggc", FakeFragment(11), None, 12),
            SiteLocation("gcct", FakeFragment(11), None, 16),
        ]

        matched = FakeRecombination().possible_locations(locations, [])
        self.assertEquals(len(matched), 2)
        self.assertEquals(len(matched[0]), 1)
        self.assertEquals(matched[0][0].site, "aggc")
        self.assertEquals(len(matched[1]), 1)
        self.assertEquals(matched[1][0].site, "gcct")

    def test_recognizes_site_combinations(self):

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc", "gaat", "gggc"]

        class FakeFragment(object):
            def __init__(self, n):
                self.id = n

        locations = [
            SiteLocation("atgc", FakeFragment(11), None, 8),
            SiteLocation("aggc", FakeFragment(11), None, 12),
            SiteLocation("gaat", FakeFragment(11), None, 26),
            SiteLocation("gggc", FakeFragment(11), None, 88),
        ]

        matched = FakeRecombination().possible_locations(locations, [])
        self.assertEquals(len(matched), 1)
        self.assertEquals(len(matched[0]), 3)
        self.assertEquals(matched[0][0].site, "aggc")
        self.assertEquals(matched[0][0].start_0based, 12)
        self.assertEquals(matched[0][1].site, "gaat")
        self.assertEquals(matched[0][1].start_0based, 26)
        self.assertEquals(matched[0][2].site, "gggc")
        self.assertEquals(matched[0][2].start_0based, 88)

    def test_recognizes_site_combinations_in_reverse(self):

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc", "gaat", "gggc"]

        class FakeFragment(object):
            def __init__(self, n):
                self.id = n

        locations = [
            SiteLocation("atgc", FakeFragment(11), None, 8),
            SiteLocation("gccc", FakeFragment(11), None, 12),
            SiteLocation("attc", FakeFragment(11), None, 26),
            SiteLocation("gcct", FakeFragment(11), None, 88),
        ]

        matched = FakeRecombination().possible_locations(locations, [])
        self.assertEquals(len(matched), 1)
        self.assertEquals(len(matched[0]), 3)
        self.assertEquals(matched[0][0].site, "gccc")
        self.assertEquals(matched[0][0].start_0based, 12)
        self.assertEquals(matched[0][1].site, "attc")
        self.assertEquals(matched[0][1].start_0based, 26)
        self.assertEquals(matched[0][2].site, "gcct")
        self.assertEquals(matched[0][2].start_0based, 88)

    def test_does_not_recognize_combination_where_sites_are_reversed_but_in_forward_order(self):

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc", "gaat", "gggc"]

        class FakeFragment(object):
            def __init__(self, n):
                self.id = n

        locations = [
            SiteLocation("atgc", FakeFragment(11), None, 8),
            SiteLocation("gcct", FakeFragment(11), None, 12),
            SiteLocation("attc", FakeFragment(11), None, 26),
            SiteLocation("gccc", FakeFragment(11), None, 88),
        ]

        matched = FakeRecombination().possible_locations(locations, [])
        self.assertEquals(len(matched), 0)

    def test_recognizes_required_site_on_insert(self):

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc", "gaat"]

            def required_insert_sites(self):
                return ["gggc"]

        class FakeFragment(object):
            def __init__(self, n):
                self.id = n

        locations = [
            SiteLocation("gggc", None, object(), 12),
            SiteLocation("atgc", FakeFragment(11), None, 8),
            SiteLocation("attc", FakeFragment(11), None, 26),
            SiteLocation("gcct", FakeFragment(11), None, 88),
        ]

        matched = FakeRecombination().possible_locations(locations, [])
        self.assertEquals(len(matched), 1)
        self.assertEquals(len(matched[0]), 3)
        self.assertEquals(matched[0][0].site, "attc")
        self.assertEquals(matched[0][0].start_0based, 26)
        self.assertEquals(matched[0][0].on_insert, False)
        self.assertEquals(matched[0][1].site, "gcct")
        self.assertEquals(matched[0][1].start_0based, 88)
        self.assertEquals(matched[0][1].on_insert, False)
        self.assertEquals(matched[0][2].site, "gggc")
        self.assertEquals(matched[0][2].start_0based, 12)
        self.assertEquals(matched[0][2].on_insert, True)

    def test_recognizes_required_site_in_reverse_on_insert(self):

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc", "gaat"]

            def required_insert_sites(self):
                return ["gggc"]

        class FakeFragment(object):
            def __init__(self, n):
                self.id = n

        locations = [
            SiteLocation("gccc", None, object(), 12),
            SiteLocation("atgc", FakeFragment(11), None, 8),
            SiteLocation("attc", FakeFragment(11), None, 26),
            SiteLocation("gcct", FakeFragment(11), None, 88),
        ]

        matched = FakeRecombination().possible_locations(locations, [])
        self.assertEquals(len(matched), 1)
        self.assertEquals(len(matched[0]), 3)
        self.assertEquals(matched[0][0].site, "attc")
        self.assertEquals(matched[0][0].start_0based, 26)
        self.assertEquals(matched[0][0].on_insert, False)
        self.assertEquals(matched[0][1].site, "gcct")
        self.assertEquals(matched[0][1].start_0based, 88)
        self.assertEquals(matched[0][1].on_insert, False)
        self.assertEquals(matched[0][2].site, "gccc")
        self.assertEquals(matched[0][2].start_0based, 12)
        self.assertEquals(matched[0][2].on_insert, True)


class ReactionEventsTest(TestCase):

    def test_runs_multiple_events_on_same_fragment_with_adjusted_bps(self):

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc"]

            def events(self, locations, errors):
                return self.generate_events(locations, FakeEvent, errors)

        event_run_on_fragments = []
        bp_shift_per_event = 17

        class FakeEvent(Event):
            def run(self, fragment, insert, insert_circular):
                event_run_on_fragments.append(fragment)
                return bp_shift_per_event

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [FakeRecombination()]

        parent_genome = Genome(name="foo")
        parent_genome.save()
        parent_fragment = Fragment.create_with_sequence("bar", "a"*100+"aggc"+"t"*100+"aggc"+"c"*100+"aggc"+"g"*100)
        Genome_Fragment(genome=parent_genome, fragment=parent_fragment, inherited=False).save()

        r = FakeReaction(parent_genome, None, None)
        child_genome = r.run_reaction("far")

        self.assertNotEqual(child_genome.id, parent_genome.id)
        self.assertNotEqual(child_genome.fragments.all()[0].id, parent_fragment.id)
        self.assertEqual(child_genome.name, "far")

        self.assertEquals(len(r.events), 3)
        self.assertEquals(len(event_run_on_fragments), 3)

        self.assertEquals(r.events[0].fragment_id, parent_fragment.id)
        self.assertEquals(r.events[0].genomic_start_0based, 100)
        self.assertEquals(r.events[0].genomic_adjusted_start_0based, 100)
        self.assertEquals(event_run_on_fragments[0].id, child_genome.fragments.all()[0].id)

        self.assertEquals(r.events[1].fragment_id, parent_fragment.id)
        self.assertEquals(r.events[1].genomic_start_0based, 204)
        self.assertEquals(r.events[1].genomic_adjusted_start_0based, 204+bp_shift_per_event)
        self.assertEquals(event_run_on_fragments[1].id, child_genome.fragments.all()[0].id)

        self.assertEquals(r.events[2].fragment_id, parent_fragment.id)
        self.assertEquals(r.events[2].genomic_start_0based, 308)
        self.assertEquals(r.events[2].genomic_adjusted_start_0based, 308+2*bp_shift_per_event)
        self.assertEquals(event_run_on_fragments[2].id, child_genome.fragments.all()[0].id)


class ExcisionRecombinationTest(TestCase):
    def test_can_excise_multiple_places_on_fragment(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "attg"+"t"*100+"cttg" +\
            "c"*1000 +\
            "attg"+"g"*100+"cttg" +\
            "g"*1000 +\
            "attg"+"a"*100+"cttg"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Excision("attg", "cttg", "attc")]

        r = FakeReaction(parent_genome, None, None)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "attc"+"c"*1000+"attc"+"g"*1000+"attc")

    def test_can_excise_reverse_sites(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "cttg"+"t"+"attg" +\
            "c"*1000 +\
            "caag"+"g"*100+"caat"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Excision("attg", "cttg", "attc")]

        r = FakeReaction(parent_genome, None, None)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "cttg"+"t"+"attg"+"c"*1000+"gaat")

    def test_works_with_left_and_right_sites_of_different_length(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "cttga"+"t"+"att" +\
            "c"*1000 +\
            "acaag"+"g"*100+"caat"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Excision("attg", "cttgt", "att"),
                        Excision("cttga", "att", "ggg")]

        r = FakeReaction(parent_genome, None, None)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "ggg"+"c"*1000+"aat")


class InversionRecombinationTest(TestCase):
    def test_can_invert_multiple_places_on_fragment(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "attg"+"t"*100+"cttg" +\
            "c"*1000 +\
            "attg"+"g"*100+"cttg" +\
            "g"*1000 +\
            "attg"+"a"*100+"cttg"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Inversion("attg", "cttg", "aatt", "ccgg")]

        r = FakeReaction(parent_genome, None, None)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "aatt"+"a"*100+"ccgg" +\
                                      "c"*1000 +\
                                      "aatt"+"c"*100+"ccgg" +\
                                      "g"*1000 +\
                                      "aatt"+"t"*100+"ccgg")

    def test_can_invert_reverse_sites(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "cttg"+"t"+"attg" +\
            "c"*1000 +\
            "caag"+"g"*100+"caat"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Inversion("attg", "cttg", "atta", "cggc")]

        r = FakeReaction(parent_genome, None, None)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "cttg"+"t"+"attg"+"c"*1000+"gccg"+"c"*100+"taat")

    def test_works_with_left_and_right_sites_of_different_length(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "cttga"+"t"+"att" +\
            "c"*1000 +\
            "acaag"+"g"*100+"caat"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Inversion("attg", "cttgt", "att", "gcccgg"),
                        Inversion("cttga", "att", "gg", "cccccc")]

        r = FakeReaction(parent_genome, None, None)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "gg"+"a"+"cccccc"+"c"*1000+"ccgggc"+"c"*100+"aat")


class IntegrationRecombinationTest(TestCase):

    def test_can_integrate_multiple_places_on_fragment(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "attg"+"t"*1000+"attg"+"c"*100
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Integration("attg", "attg", "attc", "atcc")]

        r = FakeReaction(parent_genome, "ggg"+"attg"+"a"*100, True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "attc"+"a"*100+"ggg"+"atcc" +\
                                      "t"*1000 +\
                                      "attc"+"a"*100+"ggg"+"atcc" +\
                                      "c"*100)

    def test_can_integrate_reverse_sites_on_genome(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "attg"+"t"*1000+"caat"+"c"*100
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Integration("attg", "caat", "attc", "atcc")]

        r = FakeReaction(parent_genome, "ggg"+"attg"+"a"*100, True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "ggat"+"ccc"+"t"*100+"gaat" +\
                                      "t"*1000 +\
                                      "attc"+"a"*100+"ggg"+"atcc" +\
                                      "c"*100)

    def test_can_integrate_reverse_sites_on_insert(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "attg"+"t"*1000+"attg"+"c"*100
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Integration("attg", "attg", "attc", "atcc")]

        r = FakeReaction(parent_genome, "ggg"+"caat"+"a"*100, True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "attc"+"ccc"+"t"*100+"atcc" +\
                                      "t"*1000 +\
                                      "attc"+"ccc"+"t"*100+"atcc" +\
                                      "c"*100)

    def test_can_integrate_reverse_sites_on_genome_and_insert(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "t"*1000+"attg"+"c"*1000
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Integration("caat", "caat", "attc", "atcc")]

        r = FakeReaction(parent_genome, "ggg"+"attg"+"a"*100, True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "t"*1000 +\
                                      "ggat"+"a"*100+"ggg"+"gaat" +\
                                      "c"*1000)

    def test_works_with_left_and_right_sites_of_different_length(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "t"*1000+"attg"+"c"*1000
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Integration("cttg", "attg", "aattc", "atcccc")]

        r = FakeReaction(parent_genome, "ggg"+"cttg"+"a"*100, True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "t"*1000 +\
                                      "aattc"+"a"*100+"ggg"+"atcccc" +\
                                      "c"*1000)

    def test_can_integrate_insert_with_different_site(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "t"*1000+"attg"+"c"*1000
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Integration("cttg", "attg", "attc", "atcc")]

        r = FakeReaction(parent_genome, "ggg"+"cttg"+"a"*100, True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "t"*1000 +\
                                      "attc"+"a"*100+"ggg"+"atcc" +\
                                      "c"*1000)

    def test_can_integrate_insert_with_site_across_junction(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "t"*1000+"attg"+"c"*1000
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Integration("cttg", "attg", "attc", "atcc")]

        r = FakeReaction(parent_genome, "ttg"+"a"*100+"c", True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "t"*1000 +\
                                      "attc"+"a"*100+"atcc" +\
                                      "c"*1000)


class RMCERecombinationTest(TestCase):

    def test_can_integrate_multiple_places_on_fragment(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "ggaa"+"c"*100+"gcaa"+"t"*1000+"ggaa"+"c"*100+"gcaa"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [RMCE("attg", "attc", "ggaa", "gcaa", "cctt", "ccta")]

        r = FakeReaction(parent_genome, "ggg"+"attg"+"a"*100+"attc", True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "cctt"+"a"*100+"ccta" +\
                                      "t"*1000 +\
                                      "cctt"+"a"*100+"ccta")

    def test_can_integrate_reverse_sites_on_genome(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "ttgc"+"c"*100+"ttcc"+"t"*1000+"ggaa"+"c"*100+"gcaa"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [RMCE("attg", "attc", "ggaa", "gcaa", "ccta", "cctt")]

        r = FakeReaction(parent_genome, "ggg"+"attg"+"a"*100+"attc", True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "aagg"+"t"*100+"tagg" +\
                                      "t"*1000 +\
                                      "ccta"+"a"*100+"cctt")

    def test_can_integrate_reverse_sites_on_insert(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "ggaa"+"c"*100+"gcaa"+"t"*1000+"ggaa"+"c"*100+"gcaa"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [RMCE("attg", "attc", "ggaa", "gcaa", "ccta", "cctt")]

        r = FakeReaction(parent_genome, "ggg"+"gaat"+"a"*100+"caat", True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "ccta"+"t"*100+"cctt" +\
                                      "t"*1000 +\
                                      "ccta"+"t"*100+"cctt")

    def test_can_integrate_reverse_sites_on_genome_and_insert(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "a"+"ttgc"+"c"*100+"ttcc"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [RMCE("attg", "attc", "ggaa", "gcaa", "cctt", "ccta")]

        r = FakeReaction(parent_genome, "ggg"+"gaat"+"a"*100+"caat", True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "a"+"tagg"+"a"*100+"aagg")

    def test_works_with_left_and_right_sites_of_different_length(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence("bar",
            "gaa"+"c"*100+"gcaa"+"t"*1000
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [RMCE("attg", "attcc", "gaa", "gcaa", "cctta", "cctaaa")]

        r = FakeReaction(parent_genome, "ggg"+"attg"+"a"*100+"attcc", True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "cctta"+"a"*100+"cctaaa"+"t"*1000)