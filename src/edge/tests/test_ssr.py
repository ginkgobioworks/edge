import os
import json

from Bio.Seq import Seq
from django.test import TestCase

from edge.ssr import (
  rc,
  find_query_locations_on_duplicated_template,
  find_site_locations_on_sequence,
  Recombination,
  Reaction,
)

from edge.models import Genome, Fragment, Genome_Fragment


class FindQueryOnSequenceTest(TestCase):
    def test_finds_locations_on_sequence(self):
        site = "aggc"
        template = "a"*100+site+"g"*100

        locs = find_query_locations_on_duplicated_template(template, len(template), site)
        self.assertEquals(len(locs), 1)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site)
        self.assertEquals(locs[0].start, 100)

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
        self.assertEquals(locs[0].start, 100)
        self.assertEquals(locs[1].on_insert, True)
        self.assertEquals(locs[1].site, site)
        self.assertEquals(locs[1].start, 204)

    def test_ignores_starting_position_after_sequence_length_argument(self):
        site = "aggc"
        template = "a"*100+site+"g"*100+site+"c"*100

        locs = find_query_locations_on_duplicated_template(template, len(template), site)
        self.assertEquals(len(locs), 2)

        locs = find_query_locations_on_duplicated_template(template, 102, site)
        self.assertEquals(len(locs), 1)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site)
        self.assertEquals(locs[0].start, 100)


class FindLocationsTest(TestCase):
    def test_finds_locations_on_sequence(self):
        site = "aggc"
        template = "a"*100+site+"g"*100

        locs = find_site_locations_on_sequence(template, False, [site])
        self.assertEquals(len(locs), 1)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site)
        self.assertEquals(locs[0].start, 100)

    def test_finds_multiple_sites_on_sequence(self):
        site1 = "aggc"
        site2 = "tggc"
        template = "a"*100+site1+"g"*100+site2+"t"*100

        locs = find_site_locations_on_sequence(template, False, [site1, site2])
        self.assertEquals(len(locs), 2)
        locs = sorted(locs, key=lambda l: l.site)
        self.assertEquals(locs[0].on_insert, True)
        self.assertEquals(locs[0].site, site1)
        self.assertEquals(locs[0].start, 100)
        self.assertEquals(locs[1].on_insert, True)
        self.assertEquals(locs[1].site, site2)
        self.assertEquals(locs[1].start, 204)

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
        self.assertEquals(locs[0].start, 100)
        self.assertEquals(locs[1].on_insert, True)
        self.assertEquals(locs[1].site, rc(site2))
        self.assertEquals(locs[1].start, 204)

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
        self.assertEquals(locs[0].start, 102)

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
        self.assertEquals(locs[0].start, 102)


class RecombinationTest(TestCase):
    def test_recognizes_matched_locations(self):
        pass

    def test_recognizes_matched_locations_in_reverse(self):
        pass


class ReactionTest(TestCase):
    def test_gathers_locations_on_insert_and_genome(self):
        pass

    def group_locations_into_events(self):
        pass


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
