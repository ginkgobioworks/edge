import json
import re
from django.test import TestCase

from edge.models import Genome, Fragment, Genome_Fragment


class GenomeListTest(TestCase):

    def test_empty_db(self):
        res = self.client.get('/edge/genomes/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [])

    def test_add_genome(self):
        res = self.client.post('/edge/genomes/', data=json.dumps(dict(name='foo', notes='bar')),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        uri = json.loads(res.content)['uri']
        m = re.match(r'^/edge/genomes/(\d+)/$', uri)
        genome_id = int(m.group(1))
        self.assertNotEqual(re.match(r'^/edge/genomes/\d+/$', uri), None)
        self.assertEquals(json.loads(res.content), {
            "fragments": [],
            "id": genome_id,
            "name": "foo",
            "notes": "bar",
            "parent_id": None,
            "parent_name": '',
            "uri": uri
        })

    def test_can_use_uri_from_add_genome_to_fetch_genome(self):
        res = self.client.post('/edge/genomes/', data=json.dumps(dict(name='foo', notes='bar')),
                               content_type='application/json')
        re2 = self.client.get(json.loads(res.content)['uri'])
        self.assertEquals(re2.status_code, 200)
        self.assertEquals(json.loads(res.content), json.loads(re2.content))

    def test_finds_genomes_with_specified_fragment_ids(self):
        g1 = Genome(name='Foo')
        g1.save()
        g2 = Genome(name='Bar')
        g2.save()
        f1 = Fragment(circular=True, name='FooF1')
        f1.save()
        f2 = Fragment(circular=True, name='FooF2')
        f2.save()
        f3 = Fragment(circular=True, name='FooF3', parent=f2)
        f3.save()
        Genome_Fragment(genome=g1, fragment=f1, inherited=False).save()
        Genome_Fragment(genome=g1, fragment=f2, inherited=False).save()
        Genome_Fragment(genome=g2, fragment=f1, inherited=True).save()
        Genome_Fragment(genome=g2, fragment=f3, inherited=False).save()

        # no filter, return both genomes
        res = self.client.get('/edge/genomes/')
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertCountEqual([g['id'] for g in d], [g1.id, g2.id])

        # looking for f1 and f2
        res = self.client.get('/edge/genomes/?f=%d&f=%d' % (f1.id, f2.id))
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertCountEqual([g['id'] for g in d], [g1.id])

        # looking for f1 and f3
        res = self.client.get('/edge/genomes/?f=%d&f=%d' % (f1.id, f3.id))
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertCountEqual([g['id'] for g in d], [g2.id])

        # looking for f2 and f3
        res = self.client.get('/edge/genomes/?f=%d&f=%d' % (f2.id, f3.id))
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEquals(d, [])

        # looking for f1
        res = self.client.get('/edge/genomes/?f=%d' % (f1.id,))
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEquals(d, [])

        # bad input, return []
        res = self.client.get('/edge/genomes/?f=[1,2,3]')
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEquals(d, [])

    def test_finds_genomes_with_name(self):
        from edge.models import Genome

        a = Genome(name='Foo')
        a.save()
        Genome(name='Bar %s' % a.id).save()

        # no filter, return both genomes
        res = self.client.get('/edge/genomes/')
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertCountEqual([g['name'] for g in d], ['Foo', 'Bar %s' % a.id])

        # finds genome by ID and query
        res = self.client.get('/edge/genomes/?q=%s' % a.id)
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertCountEqual([g['name'] for g in d], ['Foo', 'Bar %s' % a.id])

        # finds one
        res = self.client.get('/edge/genomes/?q=oo')
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertCountEqual([g['name'] for g in d], ['Foo'])

        # finds none
        res = self.client.get('/edge/genomes/?q=ooo')
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertCountEqual([g['name'] for g in d], [])

    def test_genome_list_paginates(self):
        from edge.models import Genome

        Genome(name='Foo').save()
        Genome(name='Bar').save()
        Genome(name='Far').save()
        Genome(name='Baz').save()

        res = self.client.get('/edge/genomes/')
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEqual(len(d), 4)

        res = self.client.get('/edge/genomes/?s=1')
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEqual(len(d), 3)

        res = self.client.get('/edge/genomes/?s=1&p=2')
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEqual(len(d), 2)

        res = self.client.get('/edge/genomes/?p=2')
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEqual(len(d), 2)

    def test_genome_list_does_not_return_inactive_genomes(self):
        from edge.models import Genome

        Genome(name='Foo', active=False).save()
        Genome(name='Bar', active=False).save()
        Genome(name='Far').save()
        Genome(name='Baz').save()

        res = self.client.get('/edge/genomes/')
        self.assertEquals(res.status_code, 200)
        d = json.loads(res.content)
        self.assertEqual(len(d), 2)
        self.assertEqual(b'Foo' in res.content, False)
        self.assertEqual(b'Bar' in res.content, False)
        self.assertEqual(b'Far' in res.content, True)
        self.assertEqual(b'Baz' in res.content, True)


class GenomeTest(TestCase):

    def setUp(self):
        res = self.client.post('/edge/genomes/', data=json.dumps(dict(name='foo', notes='bar')),
                               content_type='application/json')
        uri = json.loads(res.content)['uri']
        m = re.match(r'^/edge/genomes/(\d+)/$', uri)
        self.assertNotEqual(m, None)
        self.genome_uri = uri
        self.genome_id = int(m.group(1))

    def test_get_genome(self):
        res = self.client.get(self.genome_uri)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), {
            "fragments": [],
            "id": self.genome_id,
            "name": "foo",
            "notes": "bar",
            "parent_id": None,
            "parent_name": '',
            "uri": self.genome_uri
        })

    def test_returns_404_if_genome_does_not_exist(self):
        res = self.client.get('/edge/genomes/98765/')
        self.assertEquals(res.status_code, 404)

    def test_add_non_circular_fragment(self):
        data = dict(name='chrI', sequence='AGCTAGCTTCGATCGA')
        res = self.client.post(self.genome_uri + 'fragments/', data=json.dumps(data),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        uri = json.loads(res.content)['uri']
        m = re.match(r'^/edge/fragments/(\d+)/$', uri)
        fragment_id = int(m.group(1))
        self.assertEquals(json.loads(res.content), {
            "circular": False,
            "id": fragment_id,
            "length": len(data['sequence']),
            "name": data['name'],
            "uri": uri,
            "parent_id": None,
        })

    def test_add_circular_fragment(self):
        data = dict(name='chrI', sequence='AGCTAGCTTCGATCGA', circular=True)
        res = self.client.post(self.genome_uri + 'fragments/', data=json.dumps(data),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        uri = json.loads(res.content)['uri']
        m = re.match(r'^/edge/fragments/(\d+)/$', uri)
        fragment_id = int(m.group(1))
        self.assertEquals(json.loads(res.content), {
            "circular": True,
            "id": fragment_id,
            "length": len(data['sequence']),
            "name": data['name'],
            "uri": uri,
            "parent_id": None,
        })

    def test_get_genome_returns_fragments(self):
        data = dict(name='chrI', sequence='AGCTAGCTTCGATCGA', circular=True)
        res = self.client.post(self.genome_uri + 'fragments/',
                               data=json.dumps(data), content_type='application/json')
        uri = json.loads(res.content)['uri']
        m = re.match(r'^/edge/fragments/(\d+)/$', uri)
        fragment_id = int(m.group(1))
        res = self.client.get(self.genome_uri)
        self.assertEquals(json.loads(res.content)['fragments'], [{
            "circular": True,
            "id": fragment_id,
            "length": len(data['sequence']),
            "name": data['name'],
            "uri": uri,
            "parent_id": None,
        }])

    def test_can_use_uri_from_add_fragment_to_fetch_fragment(self):
        data = dict(name='chrI', sequence='AGCTAGCTTCGATCGA')
        res = self.client.post(self.genome_uri + 'fragments/', data=json.dumps(data),
                               content_type='application/json')
        re2 = self.client.get(json.loads(res.content)['uri'])
        self.assertEquals(re2.status_code, 200)
        self.assertEquals(json.loads(res.content), json.loads(re2.content))


class FragmentTest(TestCase):

    def setUp(self):
        self.sequence = 'AGCTAGCTTCGATCGA'
        self.name = 'foo'
        res = self.client.post('/edge/fragments/',
                               data=json.dumps(dict(name=self.name, sequence=self.sequence)),
                               content_type='application/json')
        uri = json.loads(res.content)['uri']
        m = re.match(r'^/edge/fragments/(\d+)/$', uri)
        self.fragment_id = int(m.group(1))
        self.uri = uri

    def test_get_all_user_defined_fragments(self):
        res = self.client.get('/edge/fragments/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [{
            "id": self.fragment_id,
            "name": self.name,
            "length": len(self.sequence),
            "uri": self.uri,
            "circular": False,
            "parent_id": None
        }])

    def test_does_not_return_genomic_fragments(self):
        from edge.models import Genome, Genome_Fragment

        res = self.client.get('/edge/fragments/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), 1)

        g = Genome(name='Foo')
        g.save()
        Genome_Fragment(genome=g, fragment_id=self.fragment_id, inherited=False).save()

        res = self.client.get('/edge/fragments/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), 0)

    def test_does_not_return_inactive_fragments(self):
        from edge.models import Fragment

        res = self.client.get('/edge/fragments/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), 1)

        fragment = Fragment.objects.get(pk=self.fragment_id)
        fragment.active = False
        fragment.save()

        res = self.client.get('/edge/fragments/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), 0)

    def test_get_fragment(self):
        res = self.client.get(self.uri)
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), {
            "id": self.fragment_id,
            "name": self.name,
            "length": len(self.sequence),
            "uri": self.uri,
            "circular": False,
            "parent_id": None
        })

    def test_get_fragment_sequence(self):
        res = self.client.get(self.uri + 'sequence/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), {
            'sequence': self.sequence,
            'base_first': 1,
            'base_last': len(self.sequence)
        })

    def test_get_fragment_sequence_by_position(self):
        res = self.client.get(self.uri + 'sequence/?f=3&l=10')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), {
            'sequence': self.sequence[2:10],
            'base_first': 3,
            'base_last': 10
        })

    def test_returns_404_if_fragment_does_not_exist(self):
        res = self.client.get('/edge/fragments/98765/')
        self.assertEquals(res.status_code, 404)

    def test_add_annotation_on_forward_strand(self):
        data = dict(base_first=2, base_last=9, name='proC', type='promoter', strand=1)
        res = self.client.post(self.uri + 'annotations/', data=json.dumps(data),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        self.assertEquals(json.loads(res.content), {})
        res = self.client.get(self.uri + 'annotations/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [{
            "base_first": 2,
            "base_last": 9,
            "name": 'proC',
            "type": 'promoter',
            "strand": 1,
            "qualifiers": {},
            "feature_full_length": 8,
            "feature_base_first": 1,
            "feature_base_last": 8,
        }])

    def test_add_annotation_with_qualifiers(self):
        data = dict(base_first=2, base_last=9, name='proC', type='promoter', strand=1,
                    qualifiers=dict(gene="PROC"))
        res = self.client.post(self.uri + 'annotations/', data=json.dumps(data),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        self.assertEquals(json.loads(res.content), {})
        res = self.client.get(self.uri + 'annotations/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [{
            "base_first": 2,
            "base_last": 9,
            "name": 'proC',
            "type": 'promoter',
            "strand": 1,
            "qualifiers": {"gene": "PROC"},
            "feature_full_length": 8,
            "feature_base_first": 1,
            "feature_base_last": 8,
        }])

    def test_add_annotation_on_reverse_strand(self):
        data = dict(base_first=3, base_last=10, name='proC', type='promoter', strand=-1)
        res = self.client.post(self.uri + 'annotations/', data=json.dumps(data),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        self.assertEquals(json.loads(res.content), {})
        res = self.client.get(self.uri + 'annotations/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [{
            "base_first": 3,
            "base_last": 10,
            "name": 'proC',
            "type": 'promoter',
            "strand": -1,
            "qualifiers": {},
            "feature_full_length": 8,
            "feature_base_first": 1,
            "feature_base_last": 8,
        }])

    def test_get_multiple_annotations(self):
        data = dict(base_first=2, base_last=9, name='proC', type='promoter', strand=1)
        res = self.client.post(self.uri + 'annotations/', data=json.dumps(data),
                               content_type='application/json')
        data = dict(base_first=3, base_last=10, name='proD', type='promoter', strand=-1)
        res = self.client.post(self.uri + 'annotations/', data=json.dumps(data),
                               content_type='application/json')
        res = self.client.get(self.uri + 'annotations/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [{
            "base_first": 2,
            "base_last": 9,
            "name": 'proC',
            "type": 'promoter',
            "strand": 1,
            "qualifiers": {},
            "feature_full_length": 8,
            "feature_base_first": 1,
            "feature_base_last": 8,
        }, {
            "base_first": 3,
            "base_last": 10,
            "name": 'proD',
            "type": 'promoter',
            "strand": -1,
            "qualifiers": {},
            "feature_full_length": 8,
            "feature_base_first": 1,
            "feature_base_last": 8,
        }])

    def test_get_annotations_by_region(self):
        data = dict(base_first=2, base_last=8, name='proC', type='promoter', strand=1)
        res = self.client.post(self.uri + 'annotations/', data=json.dumps(data),
                               content_type='application/json')
        data = dict(base_first=10, base_last=13, name='proD', type='promoter', strand=-1)
        res = self.client.post(self.uri + 'annotations/', data=json.dumps(data),
                               content_type='application/json')
        res = self.client.get(self.uri + 'annotations/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [{
            "base_first": 2,
            "base_last": 8,
            "name": 'proC',
            "type": 'promoter',
            "strand": 1,
            "qualifiers": {},
            "feature_full_length": 7,
            "feature_base_first": 1,
            "feature_base_last": 7,
        }, {
            "base_first": 10,
            "base_last": 13,
            "name": 'proD',
            "type": 'promoter',
            "strand": -1,
            "qualifiers": {},
            "feature_full_length": 4,
            "feature_base_first": 1,
            "feature_base_last": 4,
        }])
        res = self.client.get(self.uri + 'annotations/?f=8&l=14')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [{
            "base_first": 2,
            "base_last": 8,
            "name": 'proC',
            "type": 'promoter',
            "strand": 1,
            "qualifiers": {},
            "feature_full_length": 7,
            "feature_base_first": 1,
            "feature_base_last": 7,
        }, {
            "base_first": 10,
            "base_last": 13,
            "name": 'proD',
            "type": 'promoter',
            "strand": -1,
            "qualifiers": {},
            "feature_full_length": 4,
            "feature_base_first": 1,
            "feature_base_last": 4,
        }])
        res = self.client.get(self.uri + 'annotations/?f=9&l=14')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [{
            "base_first": 10,
            "base_last": 13,
            "name": 'proD',
            "type": 'promoter',
            "strand": -1,
            "qualifiers": {},
            "feature_full_length": 4,
            "feature_base_first": 1,
            "feature_base_last": 4,
        }])

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
            fragment.annotate(bf, bl, 'Feature %s' % (n,), 'Feature', 1)

        res = self.client.get(self.uri + 'annotations/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), nannotations)

        # limit number of annotations
        res = self.client.get(self.uri + 'annotations/?m=1')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), 1)

        res = self.client.get(self.uri + 'annotations/?m=%s' % (nannotations,))
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), nannotations)

        res = self.client.get(self.uri + 'annotations/?m=0')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(len(json.loads(res.content)), 0)


class GenomeAnnotationsTest(TestCase):

    def setUp(self):
        res = self.client.post('/edge/genomes/', data=json.dumps(dict(name='foo', notes='bar')),
                               content_type='application/json')
        self.genome_uri = json.loads(res.content)['uri']

        data = dict(name='chrI', sequence='AGCTAGCTTCGATCGA')
        res = self.client.post(self.genome_uri + 'fragments/', data=json.dumps(data),
                               content_type='application/json')
        self.fragment_uri = json.loads(res.content)['uri']
        self.fragment_data = json.loads(res.content)

    def test_find_annotation(self):
        data = dict(base_first=2, base_last=9, name='proC', type='promoter', strand=1)
        res = self.client.post(self.fragment_uri + 'annotations/', data=json.dumps(data),
                               content_type='application/json')

        res = self.client.get(self.genome_uri + 'annotations/?q=proC')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [[
            self.fragment_data, [
                {"base_first": 2,
                 "base_last": 9,
                 "name": 'proC',
                 "type": 'promoter',
                 "strand": 1,
                 "qualifiers": {},
                 "feature_full_length": 8,
                 "feature_base_first": 1,
                 "feature_base_last": 8}]
        ]])

    def test_find_annotation_by_qualifier_field(self):
        data = dict(base_first=2, base_last=4, name='Some annotation', type='promoter', strand=1,
                    qualifiers=dict(product=["Foobar"]))
        res = self.client.post(self.fragment_uri + 'annotations/', data=json.dumps(data),
                               content_type='application/json')

        data = dict(base_first=5, base_last=7, name='Another annotation', type='promoter', strand=1,
                    qualifiers=dict(product=["Atg20p"]))
        res = self.client.post(self.fragment_uri + 'annotations/', data=json.dumps(data),
                               content_type='application/json')

        res = self.client.get(self.genome_uri + 'annotations/?q=Atg20p&field=product')
        self.assertEquals(res.status_code, 200)
        print(json.loads(res.content))
        self.assertEquals(json.loads(res.content), [[
            self.fragment_data, [
                {"base_first": 5,
                 "base_last": 7,
                 "name": 'Another annotation',
                 "type": 'promoter',
                 "strand": 1,
                 "qualifiers": {'product': ['Atg20p']},
                 "feature_full_length": 3,
                 "feature_base_first": 1,
                 "feature_base_last": 3}]
        ]])


class GenomeImportTest(TestCase):

    def test_import_works(self):
        with open('edge/tests/fixtures/ecoli-mg1655-simple.gff') as fp:
            res = self.client.post('/edge/import_genome/', {'name': 'ecoli', 'attachment': fp})
            self.assertEquals(len(json.loads(res.content)['imported_genomes']), 1)
