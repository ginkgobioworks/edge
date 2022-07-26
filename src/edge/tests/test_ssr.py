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

            def events(self, locations):
                return self.generate_events(locations, FakeEvent)

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

        matched = FakeRecombination().possible_locations(locations)
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

        matched = FakeRecombination().possible_locations(locations)
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

        matched = FakeRecombination().possible_locations(locations)
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

        matched = FakeRecombination().possible_locations(locations)
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

        matched = FakeRecombination().possible_locations(locations)
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

        matched = FakeRecombination().possible_locations(locations)
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

        matched = FakeRecombination().possible_locations(locations)
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

        matched = FakeRecombination().possible_locations(locations)
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

            def events(self, locations):
                return self.generate_events(locations, FakeEvent)

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
        self.assertEquals(r.events[0].start_0based, 100)
        self.assertEquals(r.events[0].adjusted_start_0based, 100)
        self.assertEquals(event_run_on_fragments[0].id, child_genome.fragments.all()[0].id)

        self.assertEquals(r.events[1].fragment_id, parent_fragment.id)
        self.assertEquals(r.events[1].start_0based, 204)
        self.assertEquals(r.events[1].adjusted_start_0based, 204+bp_shift_per_event)
        self.assertEquals(event_run_on_fragments[1].id, child_genome.fragments.all()[0].id)

        self.assertEquals(r.events[2].fragment_id, parent_fragment.id)
        self.assertEquals(r.events[2].start_0based, 308)
        self.assertEquals(r.events[2].adjusted_start_0based, 308+2*bp_shift_per_event)
        self.assertEquals(event_run_on_fragments[2].id, child_genome.fragments.all()[0].id)


class IntegrationRecombinationTest(TestCase):
    pass


class ExcisionRecombinationTest(TestCase):
    pass


class InversionRecombinationTest(TestCase):
    pass


class RMCERecombinationTest(TestCase):
    pass


class CreLoxTest(TestCase):
    pass


"""
        g = Genome(name="foo")
        g.save()
        f = Fragment.create_with_sequence("bar", "a"*100+"atgg"+"t"*100)
        Genome_Fragment(genome=g, fragment=f, inherited=False).save()

        l = SiteLocation("atgg", f, None, 100)
        self.assertEquals(l.prev_chunk, None)
        self.assertEquals(l.first_chunk, None)
        self.assertEquals(l.last_chunk, None)
        self.assertEquals(l.next_chunk, None)

        l.chunkify()
        self.assertEquals(l.prev_chunk.sequence, "a"*100)
        self.assertEquals(l.first_chunk.sequence, "atgg")
        self.assertEquals(l.last_chunk.sequence, "atgg")
        self.assertEquals(l.next_chunk.sequence, "t"*100)
"""
