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
  SITE_SEPARATION_MAX as SSMAX
)


class FindQueryOnSequenceTest(TestCase):
    def test_finds_locations_on_sequence(self):
        site = "aggc"
        template = "a" * 100 + site + "g" * 100

        locs = find_query_locations_on_duplicated_template(template, len(template), False, site)
        self.assertEquals(len(locs), 1)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site)
        self.assertEquals(locs[0].start_0based, 100)
        self.assertEquals(locs[0].is_fragment_or_insert_circular, False)

    def test_finds_locations_on_circular_sequence_and_remember_circularity_on_location_object(self):
        site = "aggc"
        template = "a" * 100 + site + "g" * 100

        locs = find_query_locations_on_duplicated_template(template, len(template), True, site)
        self.assertEquals(len(locs), 1)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site)
        self.assertEquals(locs[0].start_0based, 100)
        self.assertEquals(locs[0].is_fragment_or_insert_circular, True)

    def test_mark_location_as_on_fragment(self):
        site = "aggc"
        template = "a" * 100 + site + "g" * 100

        locs = find_query_locations_on_duplicated_template(
            template, len(template), False, site, fragment_obj=object()
        )
        self.assertEquals(len(locs), 1)
        self.assertEquals(locs[0].on_insert, False)

    def test_finds_multiple_locations(self):
        site = "aggc"
        template = "a" * 100 + site + "g" * 100 + site + "c" * 100

        locs = find_query_locations_on_duplicated_template(template, len(template), False, site)
        self.assertEquals(len(locs), 2)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site)
        self.assertEquals(locs[0].start_0based, 100)
        self.assertEquals(locs[1].on_insert, True)
        self.assertEquals(locs[1].site, site)
        self.assertEquals(locs[1].start_0based, 204)

    def test_ignores_starting_position_after_sequence_length_argument(self):
        site = "aggc"
        template = "a" * 100 + site + "g" * 100 + site + "c" * 100

        locs = find_query_locations_on_duplicated_template(template, len(template), False, site)
        self.assertEquals(len(locs), 2)

        locs = find_query_locations_on_duplicated_template(template, 102, False, site)
        self.assertEquals(len(locs), 1)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site)
        self.assertEquals(locs[0].start_0based, 100)


class FindLocationsTest(TestCase):
    def test_finds_locations_on_sequence(self):
        site = "aggc"
        template = "a" * 100 + site + "g" * 100

        locs = find_site_locations_on_sequence(template, False, [site])
        self.assertEquals(len(locs), 1)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site)
        self.assertEquals(locs[0].start_0based, 100)

    def test_finds_multiple_sites_on_sequence(self):
        site1 = "aggc"
        site2 = "tggc"
        template = "a" * 100 + site1 + "g" * 100 + site2 + "t" * 100

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
        template = "a" * 100 + site1 + "g" * 100 + rc(site2) + "t" * 100

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
        template = site1[2:] + "g" * 100 + site1[:2]

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
        template = rc(site1)[2:] + "g" * 100 + rc(site1)[:2]

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

        f = FakeFragment(14, "c" * 100 + "aggc" + "t" * 100, False)
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

        f1 = FakeFragment(14, "c" * 100 + "aggc" + "t" * 100, False)
        f2 = FakeFragment(15, "c" * 33 + "aggc" + "t" * 100, False)
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

        f = FakeFragment(14, "c" * 100 + "aggc" + "t" * 100, False)
        r = FakeReaction(None, "g" * 33 + "ttag" + "c" * 100, True, [f])
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

        f = FakeFragment(14, "c" * 100 + "aggc" + "t" * 100 + "aggc" + "g" * 100, False)
        r = FakeReaction(None, "g" * 33 + "ttag" + "c" * 100, True, [f])
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

    def test_prioritize_based_on_allowed_array_order(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "attg" + "t" * 1000 + "attc" + "c" * 100
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        recomb_rmce = RMCE("attg", "attc", "attg", "attc", "gggg", "gggg")
        recomb_inte = Integration("attg", "attg", "cccc", "cccc")

        class FakeReaction1(Reaction):
            @staticmethod
            def allowed():
                return [recomb_inte]

        class FakeReaction2(Reaction):
            @staticmethod
            def allowed():
                return [recomb_rmce, recomb_inte]

        donor = "ggg" + "attg" + "a" * 100 + "attc"

        # integration only
        r = FakeReaction1(parent_genome, donor, True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(
            f.sequence,
            "cccc" + "a" * 100 + "attc" + "ggg" + "cccc" + "t" * 1000 + "attc" + "c" * 100
        )

        # RMCE takes precedence because it's defined earlier
        r = FakeReaction2(parent_genome, donor, True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "gggg" + "a" * 100 + "gggg" + "c" * 100)


class RecombinationTest(TestCase):

    def test_recognizes_matched_locations_on_genome(self):

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc"]

        class FakeFragment(object):
            def __init__(self, n):
                self.id = n

        locations = [
            SiteLocation("atgc", FakeFragment(11), None, False, 8),
            SiteLocation("aggc", FakeFragment(11), None, False, 12),
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
            SiteLocation("atgc", FakeFragment(11), None, False, 8),
            SiteLocation("aggc", FakeFragment(11), None, False, 12),
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
            SiteLocation("atgc", FakeFragment(11), None, False, 8),
            SiteLocation("aggc", FakeFragment(11), None, False, 12),
            SiteLocation("gcct", FakeFragment(11), None, False, 16),
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
            SiteLocation("atgc", FakeFragment(11), None, False, 8),
            SiteLocation("aggc", FakeFragment(11), None, False, 12),
            SiteLocation("gaat", FakeFragment(11), None, False, 26),
            SiteLocation("gggc", FakeFragment(11), None, False, 88),
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
            SiteLocation("atgc", FakeFragment(11), None, False, 8),
            SiteLocation("gccc", FakeFragment(11), None, False, 12),
            SiteLocation("attc", FakeFragment(11), None, False, 26),
            SiteLocation("gcct", FakeFragment(11), None, False, 88),
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
            SiteLocation("atgc", FakeFragment(11), None, False, 8),
            SiteLocation("gcct", FakeFragment(11), None, False, 12),
            SiteLocation("attc", FakeFragment(11), None, False, 26),
            SiteLocation("gccc", FakeFragment(11), None, False, 88),
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
            SiteLocation("gggc", None, object(), False, 12),
            SiteLocation("atgc", FakeFragment(11), None, False, 8),
            SiteLocation("attc", FakeFragment(11), None, False, 26),
            SiteLocation("gcct", FakeFragment(11), None, False, 88),
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
            SiteLocation("gccc", None, object(), False, 12),
            SiteLocation("atgc", FakeFragment(11), None, False, 8),
            SiteLocation("attc", FakeFragment(11), None, False, 26),
            SiteLocation("gcct", FakeFragment(11), None, False, 88),
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

    def test_recognizes_required_site_on_insert_across_circular_boundary(self):

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc"]

            def required_insert_sites(self):
                return ["gggc", "ggga"]

        class FakeFragment(object):
            def __init__(self, n):
                self.id = n

        locations = [
            SiteLocation("gggc", None, object(), True, 12),
            SiteLocation("ggga", None, object(), True, 2),
            SiteLocation("aggc", FakeFragment(11), None, False, 33),
        ]

        matched = FakeRecombination().possible_locations(locations, [])
        self.assertEquals(len(matched), 1)
        self.assertEquals(len(matched[0]), 3)
        self.assertEquals(matched[0][0].site, "aggc")
        self.assertEquals(matched[0][0].start_0based, 33)
        self.assertEquals(matched[0][0].on_insert, False)
        self.assertEquals(matched[0][1].site, "gggc")
        self.assertEquals(matched[0][1].start_0based, 12)
        self.assertEquals(matched[0][1].on_insert, True)
        self.assertEquals(matched[0][2].site, "ggga")
        self.assertEquals(matched[0][2].start_0based, 2)
        self.assertEquals(matched[0][2].on_insert, True)

    def test_does_not_recognize_sites_across_chromosomes_as_one_set(self):

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc", "gaat", "gggc"]

        class FakeFragment(object):
            def __init__(self, n):
                self.id = n

        locations = [
            SiteLocation("atgc", FakeFragment(11), None, False, 8),
            SiteLocation("aggc", FakeFragment(11), None, False, 12),
            SiteLocation("gaat", FakeFragment(11), None, False, 26),
            SiteLocation("gggc", FakeFragment(12), None, False, 88),
        ]

        matched = FakeRecombination().possible_locations(locations, [])
        self.assertEquals(len(matched), 0)

    def test_does_recognize_different_combinations_of_sites_on_different_chromosomes(self):

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc", "gaat", "gggc"]

        class FakeFragment(object):
            def __init__(self, n):
                self.id = n

        locations = [
            SiteLocation("atgc", FakeFragment(11), None, False, 8),
            SiteLocation("aggc", FakeFragment(11), None, False, 12),
            SiteLocation("gaat", FakeFragment(11), None, False, 26),
            SiteLocation("gggc", FakeFragment(11), None, False, 88),
            SiteLocation("aggc", FakeFragment(12), None, False, 42),
            SiteLocation("gaat", FakeFragment(12), None, False, 56),
            SiteLocation("gggc", FakeFragment(12), None, False, 78),
        ]

        matched = FakeRecombination().possible_locations(locations, [])
        self.assertEquals(len(matched), 2)
        self.assertEquals(len(matched[0]), 3)
        self.assertEquals(matched[0][0].site, "aggc")
        self.assertEquals(matched[0][0].start_0based, 12)
        self.assertEquals(matched[0][0].fragment_id, 11)
        self.assertEquals(len(matched[1]), 3)
        self.assertEquals(matched[1][0].site, "aggc")
        self.assertEquals(matched[1][0].start_0based, 42)
        self.assertEquals(matched[1][0].fragment_id, 12)

    def test_ignores_sites_too_far_part_and_considers_those_that_are_closer_together(self):

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc", "aggc"]

        class FakeFragment(object):
            def __init__(self, n):
                self.id = n

        locations = [
            SiteLocation("aggc", FakeFragment(11), None, False, 12),
            SiteLocation("aggc", FakeFragment(11), None, False, 12 + SSMAX + 1),
            SiteLocation("aggc", FakeFragment(11), None, False, 12 + SSMAX + 50),
            SiteLocation("aggc", FakeFragment(11), None, False, 12 + SSMAX + 50 + SSMAX + 1),
            SiteLocation("aggc", FakeFragment(11), None, False,
                         12 + SSMAX + 50 + SSMAX + 1 + SSMAX + 1),
            SiteLocation("aggc", FakeFragment(12), None, False, 42),
            SiteLocation("aggc", FakeFragment(12), None, False, 78),
        ]

        matched = FakeRecombination().possible_locations(locations, [])
        self.assertEquals(len(matched), 2)
        self.assertEquals(len(matched[0]), 2)
        self.assertEquals(matched[0][0].site, "aggc")
        self.assertEquals(matched[0][0].start_0based, 12 + SSMAX + 1)
        self.assertEquals(matched[0][0].fragment_id, 11)
        self.assertEquals(matched[0][1].site, "aggc")
        self.assertEquals(matched[0][1].start_0based, 12 + SSMAX + 50)
        self.assertEquals(matched[0][1].fragment_id, 11)
        self.assertEquals(len(matched[1]), 2)
        self.assertEquals(matched[1][0].site, "aggc")
        self.assertEquals(matched[1][0].start_0based, 42)
        self.assertEquals(matched[1][0].fragment_id, 12)
        self.assertEquals(matched[1][1].site, "aggc")
        self.assertEquals(matched[1][1].start_0based, 78)
        self.assertEquals(matched[1][1].fragment_id, 12)

    def test_ignores_site_combinations_with_multiple_events(self):

        class FakeRecombination(Recombination):
            def required_genome_sites(self):
                return ["aggc", "aggc"]

        class FakeFragment(object):
            def __init__(self, n):
                self.id = n

        locations = [
            SiteLocation("aggc", FakeFragment(11), None, False, 12),
            SiteLocation("aggc", FakeFragment(11), None, False, 22),
            SiteLocation("aggc", FakeFragment(11), None, False, 33),
            SiteLocation("aggc", FakeFragment(12), None, False, 42),
            SiteLocation("aggc", FakeFragment(12), None, False, 78),
        ]

        matched = FakeRecombination().possible_locations(locations, [])
        # 22 and 33 combo is ignored
        self.assertEquals(len(matched), 2)
        self.assertEquals(len(matched[0]), 2)
        self.assertEquals(matched[0][0].site, "aggc")
        self.assertEquals(matched[0][0].start_0based, 12)
        self.assertEquals(matched[0][0].fragment_id, 11)
        self.assertEquals(matched[0][1].site, "aggc")
        self.assertEquals(matched[0][1].start_0based, 22)
        self.assertEquals(matched[0][1].fragment_id, 11)
        self.assertEquals(len(matched[1]), 2)
        self.assertEquals(matched[1][0].site, "aggc")
        self.assertEquals(matched[1][0].start_0based, 42)
        self.assertEquals(matched[1][0].fragment_id, 12)
        self.assertEquals(matched[1][1].site, "aggc")
        self.assertEquals(matched[1][1].start_0based, 78)
        self.assertEquals(matched[1][1].fragment_id, 12)


class ExcisionRecombinationTest(TestCase):
    def test_can_excise_multiple_places_on_fragment(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        seed_sequence = "attg" + "t" * 100 + "cttg" + "c" * 1000 + "attg" + "g" * 100 +\
            "cttg" + "g" * 1000 + "attg" + "a" * 100 + "cttg"
        f = Fragment.create_with_sequence("bar", seed_sequence)

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Excision("attg", "cttg", "attc")]

        r = FakeReaction(parent_genome, None, None)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "attc" + "c" * 1000 + "attc" + "g" * 1000 + "attc")

    def test_can_excise_reverse_sites(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "cttg" + "t" + "attg" + "c" * 1000 + "caag" + "g" * 100 + "caat"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Excision("attg", "cttg", "attc")]

        r = FakeReaction(parent_genome, None, None)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "cttg" + "t" + "attg" + "c" * 1000 + "gaat")

    def test_works_with_left_and_right_sites_of_different_length(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "cttga" + "t" + "att" + "c" * 1000 + "acaag" + "g" * 100 + "caat"
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
        self.assertEquals(f.sequence, "ggg" + "c" * 1000 + "aat")


class InversionRecombinationTest(TestCase):

    def test_can_invert_using_same_site(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "attg" + "atg" + "caat" + "c" * 1000
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Inversion("attg", "caat", "attg", "caat")]

        r = FakeReaction(parent_genome, None, None)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "attg" + "cat" + "caat" + "c" * 1000)

    def test_can_invert_multiple_places_on_fragment(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        seed_sequence = "attg" + "t" * 100 + "cttg" + "c" * 1000 + "attg" +\
            "g" * 100 + "cttg" + "g" * 1000 + "attg" + "a" * 100 + "cttg"
        f = Fragment.create_with_sequence("bar", seed_sequence)

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Inversion("attg", "cttg", "aatt", "ccgg")]

        r = FakeReaction(parent_genome, None, None)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()

        new_sequence = "aatt" + "a" * 100 + "ccgg" + "c" * 1000 + "aatt" +\
            "c" * 100 + "ccgg" + "g" * 1000 + "aatt" + "t" * 100 + "ccgg"
        self.assertEquals(f.sequence, new_sequence)

    def test_can_invert_reverse_sites(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "cttg" + "t" + "attg" + "c" * 1000 + "caag" + "g" * 100 + "caat"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Inversion("attg", "cttg", "atta", "cggc")]

        r = FakeReaction(parent_genome, None, None)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(
            f.sequence,
            "cttg" + "t" + "attg" + "c" * 1000 + "gccg" + "c" * 100 + "taat"
        )

    def test_works_with_left_and_right_sites_of_different_length(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "cttga" + "t" + "att" + "c" * 1000 + "acaag" + "g" * 100 + "caat"
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
        self.assertEquals(
            f.sequence,
            "gg" + "a" + "cccccc" + "c" * 1000 + "ccgggc" + "c" * 100 + "aat"
        )


class IntegrationRecombinationTest(TestCase):

    def test_can_integrate_multiple_places_on_fragment(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "attg" + "t" * 1000 + "attg" + "c" * 100
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Integration("attg", "attg", "attc", "atcc")]

        r = FakeReaction(parent_genome, "ggg" + "attg" + "a" * 100, True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()

        new_sequence = "attc" + "a" * 100 + "ggg" + "atcc" + "t" * 1000 + "attc" +\
            "a" * 100 + "ggg" + "atcc" + "c" * 100
        self.assertEquals(f.sequence, new_sequence)

    def test_can_integrate_reverse_sites_on_genome(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "attg" + "t" * 1000 + "caat" + "c" * 100
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Integration("attg", "caat", "attc", "atcc")]

        r = FakeReaction(parent_genome, "ggg" + "attg" + "a" * 100, True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()

        new_sequence = "ggat" + "ccc" + "t" * 100 + "gaat" + "t" * 1000 +\
            "attc" + "a" * 100 + "ggg" + "atcc" + "c" * 100
        self.assertEquals(f.sequence, new_sequence)

    def test_can_integrate_reverse_sites_on_insert(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "attg" + "t" * 1000 + "attg" + "c" * 100
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Integration("attg", "attg", "attc", "atcc")]

        r = FakeReaction(parent_genome, "ggg" + "caat" + "a" * 100, True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()

        new_sequence = "attc" + "ccc" + "t" * 100 + "atcc" + "t" * 1000 +\
            "attc" + "ccc" + "t" * 100 + "atcc" + "c" * 100
        self.assertEquals(f.sequence, new_sequence)

    def test_can_integrate_reverse_sites_on_genome_and_insert(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "t" * 1000 + "attg" + "c" * 1000
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Integration("caat", "caat", "attc", "atcc")]

        r = FakeReaction(parent_genome, "ggg" + "attg" + "a" * 100, True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(
            f.sequence,
            "t" * 1000 + "ggat" + "a" * 100 + "ggg" + "gaat" + "c" * 1000
        )

    def test_works_with_left_and_right_sites_of_different_length(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "t" * 1000 + "attg" + "c" * 1000
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Integration("cttg", "attg", "aattc", "atcccc")]

        r = FakeReaction(parent_genome, "ggg" + "cttg" + "a" * 100, True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(
            f.sequence,
            "t" * 1000 + "aattc" + "a" * 100 + "ggg" + "atcccc" + "c" * 1000
        )

    def test_can_integrate_insert_with_different_site(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "t" * 1000 + "attg" + "c" * 1000
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Integration("cttg", "attg", "attc", "atcc")]

        r = FakeReaction(parent_genome, "ggg" + "cttg" + "a" * 100, True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(
            f.sequence,
            "t" * 1000 + "attc" + "a" * 100 + "ggg" + "atcc" + "c" * 1000
        )

    def test_can_integrate_insert_with_site_across_circular_boundary(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "t" * 1000 + "attg" + "c" * 1000
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [Integration("cttg", "attg", "attc", "atcc")]

        r = FakeReaction(parent_genome, "ttg" + "a" * 100 + "c", True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(
            f.sequence,
            "t" * 1000 + "attc" + "a" * 100 + "atcc" + "c" * 1000
        )


class RMCERecombinationTest(TestCase):

    def test_can_integrate_multiple_places_on_fragment(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "ggaa" + "c" * 100 + "gcaa" + "t" * 1000 + "ggaa" + "c" * 100 + "gcaa"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [RMCE("attg", "attc", "ggaa", "gcaa", "cctt", "ccta")]

        r = FakeReaction(parent_genome, "ggg" + "attg" + "a" * 100 + "attc", True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(
            f.sequence,
            "cctt" + "a" * 100 + "ccta" + "t" * 1000 + "cctt" + "a" * 100 + "ccta"
        )

    def test_can_integrate_reverse_sites_on_genome(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "ttgc" + "c" * 100 + "ttcc" + "t" * 1000 + "ggaa" + "c" * 100 + "gcaa"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [RMCE("attg", "attc", "ggaa", "gcaa", "ccta", "cctt")]

        r = FakeReaction(parent_genome, "ggg" + "attg" + "a" * 100 + "attc", True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(
            f.sequence,
            "aagg" + "t" * 100 + "tagg" + "t" * 1000 + "ccta" + "a" * 100 + "cctt"
        )

    def test_can_integrate_reverse_sites_on_insert(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "ggaa" + "c" * 100 + "gcaa" + "t" * 1000 + "ggaa" + "c" * 100 + "gcaa"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [RMCE("attg", "attc", "ggaa", "gcaa", "ccta", "cctt")]

        r = FakeReaction(parent_genome, "ggg" + "gaat" + "a" * 100 + "caat", True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(
            f.sequence,
            "ccta" + "t" * 100 + "cctt" + "t" * 1000 + "ccta" + "t" * 100 + "cctt"
        )

    def test_can_integrate_reverse_sites_on_genome_and_insert(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "a" + "ttgc" + "c" * 100 + "ttcc"
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [RMCE("attg", "attc", "ggaa", "gcaa", "cctt", "ccta")]

        r = FakeReaction(parent_genome, "ggg" + "gaat" + "a" * 100 + "caat", True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "a" + "tagg" + "a" * 100 + "aagg")

    def test_works_with_left_and_right_sites_of_different_length(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "gaa" + "c" * 100 + "gcaa" + "t" * 1000
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [RMCE("attg", "attcc", "gaa", "gcaa", "cctta", "cctaaa")]

        r = FakeReaction(parent_genome, "ggg" + "attg" + "a" * 100 + "attcc", True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "cctta" + "a" * 100 + "cctaaa" + "t" * 1000)

    def test_can_integrate_with_sites_across_circular_boundary_on_insert(self):
        parent_genome = Genome(name="foo")
        parent_genome.save()
        f = Fragment.create_with_sequence(
            "bar",
            "ggaa" + "c" * 100 + "gcaa" + "t" * 1000
        )

        Genome_Fragment(genome=parent_genome, fragment=f, inherited=False).save()

        class FakeReaction(Reaction):
            @staticmethod
            def allowed():
                return [RMCE("attg", "attc", "ggaa", "gcaa", "cctt", "ccta")]

        r = FakeReaction(parent_genome, "g" + "a" * 10 + "attc" + "ggg" + "att", True)
        child_genome = r.run_reaction("far")
        f = child_genome.fragments.all()[0].indexed_fragment()
        self.assertEquals(f.sequence, "cctt" + "a" * 10 + "ccta" + "t" * 1000)
