import json
import re

from django import forms
from django.test import TestCase
from django.urls import reverse
import unittest.mock as mock
from edge.models import Genome, Fragment, Genome_Fragment


class GenomeListTest(TestCase):
    def test_empty_db(self):
        url = reverse("genome_list")
        res = self.client.get(url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [])

    def test_add_genome(self):
        url = reverse("genome_list")
        res = self.client.post(
            url,
            data=json.dumps(dict(name="foo", notes="bar")),
            content_type="application/json",
        )
        self.assertEquals(res.status_code, 201)
        uri = json.loads(res.content)["uri"]
        m = re.match(r"^/edge/genomes/(\d+)/$", uri)
        genome_id = int(m.group(1))
        self.assertNotEqual(re.match(r"^/edge/genomes/\d+/$", uri), None)
        self.assertEquals(
            json.loads(res.content),
            {
                "fragments": [],
                "id": genome_id,
                "name": "foo",
                "notes": "bar",
                "parent_id": None,
                "parent_name": "",
                "uri": uri,
            },
        )

    def test_derives_new_genome_after_adding_fragment(self):
        g1 = Genome(name="Foo")
        g1.save()
        url = reverse("derive-genome-with-new-fragments", kwargs={"genome_id": g1.id})
        res = self.client.post(
            url,
            data=json.dumps(
                [{"name": "test-fragment", "sequence": "AGCTAGCTTCGATCGA"}]
            ),
            content_type="application/json",
        )
        self.assertEquals(res.status_code, 201)

        child = Genome.objects.get(parent=g1.id)
        self.assertNotEquals(child.id, g1.id)

        fragment = child.fragments.first()
        self.assertEquals(
            json.loads(res.content),
            {
                "fragments": [
                    {
                        "id": fragment.id,
                        "uri": "/edge/fragments/{}/".format(fragment.id),
                        "name": fragment.name,
                        "circular": fragment.circular,
                        "parent_id": None,
                        "length": 16,
                    }
                ],
                "id": child.id,
                "name": child.name,
                "notes": None,
                "parent_id": g1.id,
                "parent_name": g1.name,
                "uri": "/edge/genomes/{}/".format(child.id),
            },
        )

    def test_doesnt_derives_new_genome_on_invalid_fragment(self):
        g1 = Genome(name="Foo")
        g1.save()
        url = reverse("derive-genome-with-new-fragments", kwargs={"genome_id": g1.id})
        with self.assertRaises(forms.ValidationError) as exception:
            self.client.post(
                url,
                data=json.dumps(
                    [
                        {"name": "valid-fragment", "sequence": "AGCTAGCTTCGATCGA"},
                        {"name": "invalid-fragment",},
                    ]
                ),
                content_type="application/json",
            )
            self.assertIn("sequence", exception.exception.error_dict)

        # Ensure that when an error is hit, no child genome was derived from the initially
        # valid fragments
        self.assertEqual(Genome.objects.filter(parent=g1.id).count(), 0)

    def test_derives_new_genome_with_multiple_fragments(self):
        g1 = Genome(name="Foo")
        g1.save()
        url = reverse("derive-genome-with-new-fragments", kwargs={"genome_id": g1.id})
        res = self.client.post(
            url,
            data=json.dumps(
                [
                    {"name": "test-fragment", "sequence": "AGCTAGCTTCGATCGA"},
                    {
                        "name": "circular-fragment",
                        "sequence": "AGCTAGCTTCGATCGAAGCTATTATATCGATA",
                        "circular": True,
                    },
                ]
            ),
            content_type="application/json",
        )
        self.assertEquals(res.status_code, 201)

        child = Genome.objects.get(parent=g1.id)
        self.assertNotEquals(child.id, g1.id)

        fragments = [
            {
                "id": fragment.id,
                "uri": "/edge/fragments/{}/".format(fragment.id),
                "name": fragment.name,
                "circular": fragment.circular,
                "parent_id": None,
                "length": fragment.indexed_fragment().length,
            }
            for fragment in child.fragments.all()
        ]
        self.assertEquals(
            json.loads(res.content),
            {
                "fragments": fragments,
                "id": child.id,
                "name": child.name,
                "notes": None,
                "parent_id": g1.id,
                "parent_name": g1.name,
                "uri": "/edge/genomes/{}/".format(child.id),
            },
        )

    def test_can_use_uri_from_add_genome_to_fetch_genome(self):
        url = reverse("genome_list")
        res = self.client.post(
            url,
            data=json.dumps(dict(name="foo", notes="bar")),
            content_type="application/json",
        )
        re2 = self.client.get(json.loads(res.content)["uri"])
        self.assertEquals(re2.status_code, 200)
        self.assertEquals(json.loads(res.content), json.loads(re2.content))

    def test_finds_genomes_with_specified_fragment_ids(self):
        g1 = Genome(name="Foo")
        g1.save()
        g2 = Genome(name="Bar")
        g2.save()
        f1 = Fragment(circular=True, name="FooF1")
        f1.save()
        f2 = Fragment(circular=True, name="FooF2")
        f2.save()
        f3 = Fragment(circular=True, name="FooF3", parent=f2)
        f3.save()
        Genome_Fragment(genome=g1, fragment=f1, inherited=False).save()
        Genome_Fragment(genome=g1, fragment=f2, inherited=False).save()
        Genome_Fragment(genome=g2, fragment=f1, inherited=True).save()
        Genome_Fragment(genome=g2, fragment=f3, inherited=False).save()

        # no filter, return both genomes
        url = reverse("genome_list")
        res = self.client.get(url)
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertCountEqual([g["id"] for g in d], [g1.id, g2.id])

        # looking for f1 and f2
        res = self.client.get("%s?f=%d&f=%d" % (url, f1.id, f2.id))
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertCountEqual([g["id"] for g in d], [g1.id])

        # looking for f1 and f3
        res = self.client.get("%s?f=%d&f=%d" % (url, f1.id, f3.id))
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertCountEqual([g["id"] for g in d], [g2.id])

        # looking for f2 and f3
        res = self.client.get("%s?f=%d&f=%d" % (url, f2.id, f3.id))
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEquals(d, [])

        # looking for f1
        res = self.client.get("%s?f=%d" % (url, f1.id,))
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEquals(d, [])

        # bad input, return []
        res = self.client.get("%s?f=[1,2,3]" % url)
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEquals(d, [])

    def test_finds_genomes_with_name(self):
        from edge.models import Genome

        a = Genome(name="Foo")
        a.save()
        Genome(name="Bar %s" % a.id).save()

        # no filter, return both genomes
        url = reverse("genome_list")
        res = self.client.get(url)
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertCountEqual([g["name"] for g in d], ["Foo", "Bar %s" % a.id])

        # finds genome by ID and query
        res = self.client.get("%s?q=%s" % (url, a.id))
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertCountEqual([g["name"] for g in d], ["Foo", "Bar %s" % a.id])

        # finds one
        res = self.client.get("%s?q=oo" % url)
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertCountEqual([g["name"] for g in d], ["Foo"])

        # finds none
        res = self.client.get("%s?q=ooo" % url)
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertCountEqual([g["name"] for g in d], [])

    def test_genome_list_paginates(self):
        from edge.models import Genome

        Genome(name="Foo").save()
        Genome(name="Bar").save()
        Genome(name="Far").save()
        Genome(name="Baz").save()

        url = reverse("genome_list")
        res = self.client.get(url)
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEqual(len(d), 4)

        res = self.client.get("%s?s=1" % url)
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEqual(len(d), 3)

        res = self.client.get("%s?s=1&p=2" % url)
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEqual(len(d), 2)

        res = self.client.get("%s?p=2" % url)
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEqual(len(d), 2)

    def test_genome_list_does_not_return_inactive_genomes(self):
        from edge.models import Genome

        Genome(name="Foo", active=False).save()
        Genome(name="Bar", active=False).save()
        Genome(name="Far").save()
        Genome(name="Baz").save()

        url = reverse("genome_list")
        res = self.client.get(url)
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEqual(len(d), 2)
        self.assertEqual(b"Foo" in res.content, False)
        self.assertEqual(b"Bar" in res.content, False)
        self.assertEqual(b"Far" in res.content, True)
        self.assertEqual(b"Baz" in res.content, True)


class GenomeTest(TestCase):
    def setUp(self):
        url = reverse("genome_list")
        res = self.client.post(
            url,
            data=json.dumps(dict(name="foo", notes="bar")),
            content_type="application/json",
        )
        uri = json.loads(res.content)["uri"]
        m = re.match(r"^/edge/genomes/(\d+)/$", uri)
        self.assertNotEqual(m, None)
        self.genome_uri = uri
        self.genome_id = int(m.group(1))

    def test_get_genome(self):
        res = self.client.get(self.genome_uri)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(
            json.loads(res.content),
            {
                "fragments": [],
                "id": self.genome_id,
                "name": "foo",
                "notes": "bar",
                "parent_id": None,
                "parent_name": "",
                "uri": self.genome_uri,
            },
        )

    def test_returns_404_if_genome_does_not_exist(self):
        url = reverse("genome", kwargs=dict(genome_id=98765))
        res = self.client.get(url)
        self.assertEquals(res.status_code, 404)

    def test_add_non_circular_fragment(self):
        data = dict(name="chrI", sequence="AGCTAGCTTCGATCGA")
        res = self.client.post(
            self.genome_uri + "fragments/",
            data=json.dumps(data),
            content_type="application/json",
        )
        self.assertEquals(res.status_code, 201)
        uri = json.loads(res.content)["uri"]
        m = re.match(r"^/edge/fragments/(\d+)/$", uri)
        fragment_id = int(m.group(1))
        self.assertEquals(
            json.loads(res.content),
            {
                "circular": False,
                "id": fragment_id,
                "length": len(data["sequence"]),
                "name": data["name"],
                "uri": uri,
                "parent_id": None,
            },
        )

    def test_add_circular_fragment(self):
        data = dict(name="chrI", sequence="AGCTAGCTTCGATCGA", circular=True)
        res = self.client.post(
            self.genome_uri + "fragments/",
            data=json.dumps(data),
            content_type="application/json",
        )
        self.assertEquals(res.status_code, 201)
        uri = json.loads(res.content)["uri"]
        m = re.match(r"^/edge/fragments/(\d+)/$", uri)
        fragment_id = int(m.group(1))
        self.assertEquals(
            json.loads(res.content),
            {
                "circular": True,
                "id": fragment_id,
                "length": len(data["sequence"]),
                "name": data["name"],
                "uri": uri,
                "parent_id": None,
            },
        )

    def test_get_genome_returns_fragments(self):
        data = dict(name="chrI", sequence="AGCTAGCTTCGATCGA", circular=True)
        res = self.client.post(
            self.genome_uri + "fragments/",
            data=json.dumps(data),
            content_type="application/json",
        )
        uri = json.loads(res.content)["uri"]
        m = re.match(r"^/edge/fragments/(\d+)/$", uri)
        fragment_id = int(m.group(1))
        res = self.client.get(self.genome_uri)
        self.assertEquals(
            json.loads(res.content)["fragments"],
            [
                {
                    "circular": True,
                    "id": fragment_id,
                    "length": len(data["sequence"]),
                    "name": data["name"],
                    "uri": uri,
                    "parent_id": None,
                }
            ],
        )

    def test_can_use_uri_from_add_fragment_to_fetch_fragment(self):
        data = dict(name="chrI", sequence="AGCTAGCTTCGATCGA")
        res = self.client.post(
            self.genome_uri + "fragments/",
            data=json.dumps(data),
            content_type="application/json",
        )
        re2 = self.client.get(json.loads(res.content)["uri"])
        self.assertEquals(re2.status_code, 200)
        self.assertEquals(json.loads(res.content), json.loads(re2.content))


class FragmentTest(TestCase):
    def setUp(self):
        url = reverse("fragment_list")
        self.sequence = "AGCTAGCTTCGATCGA"
        self.name = "foo"
        res = self.client.post(
            url,
            data=json.dumps(dict(name=self.name, sequence=self.sequence)),
            content_type="application/json",
        )
        uri = json.loads(res.content)["uri"]
        m = re.match(r"^/edge/fragments/(\d+)/$", uri)
        self.fragment_id = int(m.group(1))
        self.uri = uri

    def test_get_all_user_defined_fragments(self):
        url = reverse("fragment_list")
        res = self.client.get(url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(
            json.loads(res.content),
            [
                {
                    "id": self.fragment_id,
                    "name": self.name,
                    "length": len(self.sequence),
                    "uri": self.uri,
                    "circular": False,
                    "parent_id": None,
                }
            ],
        )

    def test_does_not_return_genomic_fragments(self):
        from edge.models import Genome, Genome_Fragment

        url = reverse("fragment_list")
        res = self.client.get(url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), 1)

        g = Genome(name="Foo")
        g.save()
        Genome_Fragment(genome=g, fragment_id=self.fragment_id, inherited=False).save()

        res = self.client.get(url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), 0)

    def test_does_not_return_inactive_fragments(self):
        from edge.models import Fragment

        url = reverse("fragment_list")
        res = self.client.get(url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), 1)

        fragment = Fragment.objects.get(pk=self.fragment_id)
        fragment.active = False
        fragment.save()

        res = self.client.get(url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), 0)

    def test_get_fragment(self):
        res = self.client.get(self.uri)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(
            json.loads(res.content),
            {
                "id": self.fragment_id,
                "name": self.name,
                "length": len(self.sequence),
                "uri": self.uri,
                "circular": False,
                "parent_id": None,
            },
        )

    def test_get_fragment_sequence(self):
        url = reverse("fragment_sequence", kwargs=dict(fragment_id=self.fragment_id))
        res = self.client.get(url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(
            json.loads(res.content),
            {
                "sequence": self.sequence,
                "base_first": 1,
                "base_last": len(self.sequence),
            },
        )

    def test_get_fragment_sequence_by_position(self):
        url = reverse("fragment_sequence", kwargs=dict(fragment_id=self.fragment_id))
        res = self.client.get("%s?f=3&l=10" % url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(
            json.loads(res.content),
            {"sequence": self.sequence[2:10], "base_first": 3, "base_last": 10},
        )

    def test_returns_404_if_fragment_does_not_exist(self):
        url = reverse("fragment", kwargs=dict(fragment_id=98765))
        res = self.client.get(url)
        self.assertEquals(res.status_code, 404)

    def test_add_annotation_on_forward_strand(self):
        data = dict(base_first=2, base_last=9, name="proC", type="promoter", strand=1)
        url = reverse("fragment_annotations", kwargs=dict(fragment_id=self.fragment_id))
        res = self.client.post(
            url, data=json.dumps(data), content_type="application/json"
        )
        self.assertEquals(res.status_code, 201)
        self.assertEquals(json.loads(res.content), {})
        res = self.client.get(self.uri + "annotations/")
        self.assertEquals(res.status_code, 200)
        self.assertEquals(
            json.loads(res.content),
            [
                {
                    "base_first": 2,
                    "base_last": 9,
                    "strand": 1,
                    "feature_base_first": 1,
                    "feature_base_last": 8,
                    "feature": {
                        "id": mock.ANY,
                        "length": 8,
                        "name": "proC",
                        "type": "promoter",
                        "qualifiers": {}
                    },
                    "name": "proC",
                    "type": "promoter",
                    "feature_full_length": 8,
                    "qualifiers": {},
                }
            ],
        )

    def test_add_annotation_with_qualifiers(self):
        data = dict(
            base_first=2,
            base_last=9,
            name="proC",
            type="promoter",
            strand=1,
            qualifiers=dict(gene="PROC"),
        )
        url = reverse("fragment_annotations", kwargs=dict(fragment_id=self.fragment_id))
        res = self.client.post(
            url, data=json.dumps(data), content_type="application/json"
        )
        self.assertEquals(res.status_code, 201)
        self.assertEquals(json.loads(res.content), {})
        res = self.client.get(self.uri + "annotations/")
        self.assertEquals(res.status_code, 200)
        self.assertEquals(
            json.loads(res.content),
            [
                {
                    "base_first": 2,
                    "base_last": 9,
                    "strand": 1,
                    "feature_base_first": 1,
                    "feature_base_last": 8,
                    "feature": {
                        "id": mock.ANY,
                        "length": 8,
                        "name": "proC",
                        "type": "promoter",
                        "qualifiers": {"gene": "PROC"}
                    },
                    "name": "proC",
                    "type": "promoter",
                    "qualifiers": {"gene": "PROC"},
                    "feature_full_length": 8,
                }
            ],
        )

    def test_add_annotation_on_reverse_strand(self):
        data = dict(base_first=3, base_last=10, name="proC", type="promoter", strand=-1)
        url = reverse("fragment_annotations", kwargs=dict(fragment_id=self.fragment_id))
        res = self.client.post(
            url, data=json.dumps(data), content_type="application/json"
        )
        self.assertEquals(res.status_code, 201)
        self.assertEquals(json.loads(res.content), {})
        res = self.client.get(self.uri + "annotations/")
        self.assertEquals(res.status_code, 200)
        self.assertEquals(
            json.loads(res.content),
            [
                {
                    "base_first": 3,
                    "base_last": 10,
                    "strand": -1,
                    "feature_base_first": 1,
                    "feature_base_last": 8,
                    "feature": {
                        "id": mock.ANY,
                        "length": 8,
                        "name": "proC",
                        "type": "promoter",
                        "qualifiers": {}
                    },
                    "name": "proC",
                    "type": "promoter",
                    "qualifiers": {},
                    "feature_full_length": 8,
                }
            ],
        )

    def test_annotate_multiple_chunks_with_qualifiers(self):
        data = dict(
            bases=((2, 4), (6, 9)),
            name="proC",
            type="promoter",
            strand=1,
            qualifiers=dict(gene="PROC"),
        )
        url = reverse(
            "fragment_annotate_chunks", kwargs=dict(fragment_id=self.fragment_id)
        )
        res = self.client.post(
            url, data=json.dumps(data), content_type="application/json"
        )
        self.assertEquals(res.status_code, 201)
        self.assertEquals(json.loads(res.content), {})
        res = self.client.get(self.uri + "annotations/")
        self.assertEquals(res.status_code, 200)
        self.assertEquals(
            json.loads(res.content),
            [
                {
                    "base_first": 2,
                    "base_last": 4,
                    "strand": 1,
                    "feature_base_first": 1,
                    "feature_base_last": 3,
                    "feature": {
                        "id": mock.ANY,
                        "length": 7,
                        "name": "proC",
                        "type": "promoter",
                        "qualifiers": {"gene": "PROC"}
                    },
                    "name": "proC",
                    "type": "promoter",
                    "qualifiers": {"gene": "PROC"},
                    "feature_full_length": 7,
                },
                {
                    "base_first": 6,
                    "base_last": 9,
                    "strand": 1,
                    "feature_base_first": 4,
                    "feature_base_last": 7,
                    "feature": {
                        "id": mock.ANY,
                        "length": 7,
                        "name": "proC",
                        "type": "promoter",
                        "qualifiers": {"gene": "PROC"}
                    },
                    "name": "proC",
                    "type": "promoter",
                    "qualifiers": {"gene": "PROC"},
                    "feature_full_length": 7,
                },
            ],
        )

    def test_get_multiple_annotations_that_are_overlapping(self):
        data = dict(base_first=2, base_last=9, name="proC", type="promoter", strand=1)
        url = reverse("fragment_annotations", kwargs=dict(fragment_id=self.fragment_id))
        res = self.client.post(
            url, data=json.dumps(data), content_type="application/json"
        )
        data = dict(base_first=3, base_last=10, name="proD", type="promoter", strand=-1)
        res = self.client.post(
            url, data=json.dumps(data), content_type="application/json"
        )
        res = self.client.get(self.uri + "annotations/")
        self.assertEquals(res.status_code, 200)
        self.assertEquals(
            json.loads(res.content),
            [
                {
                    "base_first": 2,
                    "base_last": 9,
                    "strand": 1,
                    "feature_base_first": 1,
                    "feature_base_last": 8,
                    "feature": {
                        "id": mock.ANY,
                        "length": 8,
                        "name": "proC",
                        "type": "promoter",
                        "qualifiers": {}
                    },
                    "name": "proC",
                    "type": "promoter",
                    "qualifiers": {},
                    "feature_full_length": 8,
                },
                {
                    "base_first": 3,
                    "base_last": 10,
                    "strand": -1,
                    "feature_base_first": 1,
                    "feature_base_last": 8,
                    "feature": {
                        "id": mock.ANY,
                        "length": 8,
                        "name": "proD",
                        "type": "promoter",
                        "qualifiers": {}
                    },
                    "name": "proD",
                    "type": "promoter",
                    "qualifiers": {},
                    "feature_full_length": 8,
                },
            ],
        )

    def test_get_annotations_by_region(self):
        data = dict(base_first=2, base_last=8, name="proC", type="promoter", strand=1)
        url = reverse("fragment_annotations", kwargs=dict(fragment_id=self.fragment_id))
        res = self.client.post(
            url, data=json.dumps(data), content_type="application/json"
        )
        data = dict(
            base_first=10, base_last=13, name="proD", type="promoter", strand=-1
        )
        res = self.client.post(
            url, data=json.dumps(data), content_type="application/json"
        )
        res = self.client.get(url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(
            json.loads(res.content),
            [
                {
                    "base_first": 2,
                    "base_last": 8,
                    "strand": 1,
                    "feature_base_first": 1,
                    "feature_base_last": 7,
                    "feature": {
                        "id": mock.ANY,
                        "length": 7,
                        "name": "proC",
                        "type": "promoter",
                        "qualifiers": {}
                    },
                    "name": "proC",
                    "type": "promoter",
                    "qualifiers": {},
                    "feature_full_length": 7,
                },
                {
                    "base_first": 10,
                    "base_last": 13,
                    "strand": -1,
                    "feature_base_first": 1,
                    "feature_base_last": 4,
                    "feature": {
                        "id": mock.ANY,
                        "length": 4,
                        "name": "proD",
                        "type": "promoter",
                        "qualifiers": {}
                    },
                    "name": "proD",
                    "type": "promoter",
                    "qualifiers": {},
                    "feature_full_length": 4,
                },
            ],
        )
        res = self.client.get("%s?f=8&l=14" % url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(
            json.loads(res.content),
            [
                {
                    "base_first": 2,
                    "base_last": 8,
                    "feature": {
                        "id": mock.ANY,
                        "length": 7,
                        "name": "proC",
                        "type": "promoter",
                        "qualifiers": {}
                    },
                    "strand": 1,
                    "feature_base_first": 1,
                    "feature_base_last": 7,
                    "name": "proC",
                    "type": "promoter",
                    "qualifiers": {},
                    "feature_full_length": 7,
                },
                {
                    "base_first": 10,
                    "base_last": 13,
                    "strand": -1,
                    "feature_base_first": 1,
                    "feature_base_last": 4,
                    "feature": {
                        "id": mock.ANY,
                        "length": 4,
                        "name": "proD",
                        "type": "promoter",
                        "qualifiers": {}
                    },
                    "name": "proD",
                    "type": "promoter",
                    "qualifiers": {},
                    "feature_full_length": 4,
                },
            ],
        )
        res = self.client.get("%s?f=9&l=14" % url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(
            json.loads(res.content),
            [
                {
                    "base_first": 10,
                    "base_last": 13,
                    "strand": -1,
                    "feature_base_first": 1,
                    "feature_base_last": 4,
                    "feature": {
                        "id": mock.ANY,
                        "length": 4,
                        "name": "proD",
                        "type": "promoter",
                        "qualifiers": {}
                    },
                    "name": "proD",
                    "type": "promoter",
                    "qualifiers": {},
                    "feature_full_length": 4,
                }
            ],
        )

    def test_limit_max_annotations_to_fetch(self):
        from edge.models import Fragment
        import random

        fragment = Fragment.objects.get(pk=self.fragment_id)
        fragment = fragment.indexed_fragment()
        flen = fragment.length
        nannotations = 10

        # create some annotations
        for n in range(0, nannotations):
            bf = random.randint(1, flen)
            bl = random.randint(bf, flen)
            fragment.annotate(bf, bl, "Feature %s" % (n,), "Feature", 1)

        url = reverse("fragment_annotations", kwargs=dict(fragment_id=self.fragment_id))
        res = self.client.get(url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), nannotations)

        # limit number of annotations
        res = self.client.get("%s?m=1" % url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), 1)

        res = self.client.get("%s?m=%s" % (url, nannotations,))
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), nannotations)

        res = self.client.get("%s?m=0" % url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), 0)


class GenomeAnnotationsTest(TestCase):
    def setUp(self):
        url = reverse("genome_list")
        res = self.client.post(
            url,
            data=json.dumps(dict(name="foo", notes="bar")),
            content_type="application/json",
        )
        self.genome_id = json.loads(res.content)["id"]
        self.genome_uri = json.loads(res.content)["uri"]

        data = dict(name="chrI", sequence="AGCTAGCTTCGATCGA")
        url = reverse("genome_fragments", kwargs=dict(genome_id=self.genome_id))
        res = self.client.post(
            url, data=json.dumps(data), content_type="application/json"
        )
        self.fragment_uri = json.loads(res.content)["uri"]
        self.fragment_data = json.loads(res.content)
        m = re.match(r"^/edge/fragments/(\d+)/$", self.fragment_uri)
        self.fragment_id = int(m.group(1))

    def test_returns_errors_if_genome_not_indexed_and_creates_indexed_genome(self):
        genome = Genome.objects.get(pk=self.genome_id)
        for fr in genome.fragments.all():
            fr.fragment_chunk_location_set.all().delete()

        # not indexed before call
        self.assertEquals(genome.has_location_index, False)

        # gets 200 with error
        url = reverse("genome_annotations", kwargs=dict(genome_id=self.genome_id))
        res = self.client.get("%s?q=proC" % url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(
            "Missing genome indices" in json.loads(res.content)["error"], True
        )

        # in test mode, the .delay call via celery is inlined
        genome = Genome.objects.get(pk=self.genome_id)
        self.assertEquals(genome.has_location_index, True)

    def test_find_annotation(self):
        data = dict(base_first=2, base_last=9, name="proC", type="promoter", strand=1)
        url = reverse("fragment_annotations", kwargs=dict(fragment_id=self.fragment_id))
        res = self.client.post(
            url, data=json.dumps(data), content_type="application/json"
        )

        url = reverse("genome_annotations", kwargs=dict(genome_id=self.genome_id))
        res = self.client.get("%s?q=proC" % url)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(
            json.loads(res.content),
            [
                [
                    self.fragment_data,
                    [
                        {
                            "base_first": 2,
                            "base_last": 9,
                            "feature": {
                                "id": mock.ANY,
                                "length": 8,
                                "name": "proC",
                                "type": "promoter",
                                "qualifiers": {}
                            },
                            "strand": 1,
                            "feature_base_first": 1,
                            "feature_base_last": 8,
                            "name": "proC",
                            "type": "promoter",
                            "qualifiers": {},
                            "feature_full_length": 8,
                        }
                    ],
                ]
            ],
        )

    def test_find_annotation_by_qualifier_field(self):
        data = dict(
            base_first=2,
            base_last=4,
            name="Some annotation",
            type="promoter",
            strand=1,
            qualifiers=dict(product=["Foobar"]),
        )
        url = reverse("fragment_annotations", kwargs=dict(fragment_id=self.fragment_id))
        res = self.client.post(
            url, data=json.dumps(data), content_type="application/json"
        )

        data = dict(
            base_first=5,
            base_last=7,
            name="Another annotation",
            type="promoter",
            strand=1,
            qualifiers=dict(product=["Atg20p"]),
        )
        res = self.client.post(
            url, data=json.dumps(data), content_type="application/json"
        )

        url = reverse("genome_annotations", kwargs=dict(genome_id=self.genome_id))
        res = self.client.get("%s?q=Atg20p&field=product" % url)
        self.assertEquals(res.status_code, 200)
        print(json.loads(res.content)[0])
        self.assertEquals(
            json.loads(res.content),
            [
                [
                    self.fragment_data,
                    [
                        {
                            "base_first": 5,
                            "base_last": 7,
                            "feature": {
                                "id": mock.ANY,
                                "length": 3,
                                "name": "Another annotation",
                                "type": "promoter",
                                "qualifiers": {"product": ["Atg20p"]}
                            },
                            "strand": 1,
                            "feature_base_first": 1,
                            "feature_base_last": 3,
                            "name": "Another annotation",
                            "type": "promoter",
                            "qualifiers": {"product": ["Atg20p"]},
                            "feature_full_length": 3,
                        }
                    ],
                ]
            ],
        )


class GenomeImportTest(TestCase):
    def test_import_works(self):
        with open("edge/tests/fixtures/ecoli-mg1655-simple.gff") as fp:
            url = reverse("import")
            res = self.client.post(url, {"name": "ecoli", "attachment": fp})
            self.assertEquals(len(json.loads(res.content)["imported_genomes"]), 1)
