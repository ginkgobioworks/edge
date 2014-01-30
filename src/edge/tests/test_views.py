import os
import json
import tempfile
from django.test import TestCase


class GenomeListTest(TestCase):

    def test_empty_db(self):
        res = self.client.get('/edge/genomes/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [])

    def test_add_genome(self):
        res = self.client.post('/edge/genomes/', data=json.dumps(dict(name='foo', notes='bar')),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        print json.loads(res.content)
        self.assertEquals(json.loads(res.content), {
            "fragments": [],
            "id": 1,
            "name": "foo",
            "notes": "bar",
            "parent_id": None,
            "uri": "/edge/genomes/1/",
        })

    def test_can_use_uri_from_add_genome_to_fetch_genome(self):
        res = self.client.post('/edge/genomes/', data=json.dumps(dict(name='foo', notes='bar')),
                               content_type='application/json')
        re2 = self.client.get(json.loads(res.content)['uri'])
        self.assertEquals(re2.status_code, 200)
        self.assertEquals(json.loads(res.content), json.loads(re2.content))


class GenomeTest(TestCase):

    def setUp(self):
        res = self.client.post('/edge/genomes/', data=json.dumps(dict(name='foo', notes='bar')),
                               content_type='application/json')

    def test_get_genome(self):
        res = self.client.get('/edge/genomes/1/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), {
            "fragments": [],
            "id": 1,
            "name": "foo",
            "notes": "bar",
            "parent_id": None,
            "uri": '/edge/genomes/1/',
        })

    def test_returns_404_if_genome_does_not_exist(self):
        res = self.client.get('/edge/genomes/100/')
        self.assertEquals(res.status_code, 404)

    def test_add_non_circular_fragment(self):
        data = dict(name='chrI', sequence='AGCTAGCTTCGATCGA')
        res = self.client.post('/edge/genomes/1/fragments/', data=json.dumps(data),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        self.assertEquals(json.loads(res.content), {
            "circular": False,
            "id": 1,
            "length": len(data['sequence']),
            "name": data['name'],
            "uri": "/edge/fragments/1/"
        })

    def test_add_circular_fragment(self):
        data = dict(name='chrI', sequence='AGCTAGCTTCGATCGA', circular=True)
        res = self.client.post('/edge/genomes/1/fragments/', data=json.dumps(data),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        self.assertEquals(json.loads(res.content), {
            "circular": True,
            "id": 1,
            "length": len(data['sequence']),
            "name": data['name'],
            "uri": "/edge/fragments/1/"
        })

    def test_get_genome_returns_fragments(self):
        data = dict(name='chrI', sequence='AGCTAGCTTCGATCGA', circular=True)
        res = self.client.post('/edge/genomes/1/fragments/',
                               data=json.dumps(data), content_type='application/json')
        print json.loads(res.content)
        res = self.client.get('/edge/genomes/1/')
        print json.loads(res.content)
        self.assertEquals(json.loads(res.content)['fragments'], [{
            "circular": True,
            "id": 1,
            "length": len(data['sequence']),
            "name": data['name'],
            "uri": "/edge/fragments/1/"
        }])

    def test_can_use_uri_from_add_fragment_to_fetch_fragment(self):
        data = dict(name='chrI', sequence='AGCTAGCTTCGATCGA')
        res = self.client.post('/edge/genomes/1/fragments/', data=json.dumps(data),
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

    def test_get_all_non_genomic_fragments(self):
        res = self.client.get('/edge/fragments/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [{
            "id": 1,
            "name": self.name,
            "length": len(self.sequence),
            "uri": '/edge/fragments/1/',
            "circular": False,
        }])

    def test_get_fragment(self):
        res = self.client.get('/edge/fragments/1/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), {
            "id": 1,
            "name": self.name,
            "length": len(self.sequence),
            "uri": '/edge/fragments/1/',
            "circular": False,
        })

    def test_get_fragment_sequence(self):
        res = self.client.get('/edge/fragments/1/sequence/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), {
            'sequence': self.sequence,
            'base_first': 1,
            'base_last': len(self.sequence)
        })

    def test_get_fragment_sequence_by_position(self):
        res = self.client.get('/edge/fragments/1/sequence/?f=3&l=10')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), {
            'sequence': self.sequence[2:10],
            'base_first': 3,
            'base_last': 10
        })

    def test_returns_404_if_fragment_does_not_exist(self):
        res = self.client.get('/edge/fragments/100/')
        self.assertEquals(res.status_code, 404)

    def test_add_annotation_on_forward_strand(self):
        data = dict(base_first=2, base_last=9, name='proC', type='promoter', strand=1)
        res = self.client.post('/edge/fragments/1/annotations/', data=json.dumps(data),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        self.assertEquals(json.loads(res.content), {})
        res = self.client.get('/edge/fragments/1/annotations/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [{
            "base_first": 2,
            "base_last": 9,
            "name": 'proC',
            "type": 'promoter',
            "strand": 1,
            "feature_full_length": 8,
            "feature_base_first": 1,
            "feature_base_last": 8,
        }])

    def test_add_annotation_on_reverse_strand(self):
        data = dict(base_first=3, base_last=10, name='proC', type='promoter', strand=-1)
        res = self.client.post('/edge/fragments/1/annotations/', data=json.dumps(data),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        self.assertEquals(json.loads(res.content), {})
        res = self.client.get('/edge/fragments/1/annotations/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [{
            "base_first": 3,
            "base_last": 10,
            "name": 'proC',
            "type": 'promoter',
            "strand": -1,
            "feature_full_length": 8,
            "feature_base_first": 1,
            "feature_base_last": 8,
        }])

    def test_get_multiple_annotations(self):
        data = dict(base_first=2, base_last=9, name='proC', type='promoter', strand=1)
        res = self.client.post('/edge/fragments/1/annotations/', data=json.dumps(data),
                               content_type='application/json')
        data = dict(base_first=3, base_last=10, name='proD', type='promoter', strand=-1)
        res = self.client.post('/edge/fragments/1/annotations/', data=json.dumps(data),
                               content_type='application/json')
        res = self.client.get('/edge/fragments/1/annotations/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [{
            "base_first": 2,
            "base_last": 9,
            "name": 'proC',
            "type": 'promoter',
            "strand": 1,
            "feature_full_length": 8,
            "feature_base_first": 1,
            "feature_base_last": 8,
        }, {
            "base_first": 3,
            "base_last": 10,
            "name": 'proD',
            "type": 'promoter',
            "strand": -1,
            "feature_full_length": 8,
            "feature_base_first": 1,
            "feature_base_last": 8,
        }])

    def test_get_annotations_by_region(self):
        data = dict(base_first=2, base_last=8, name='proC', type='promoter', strand=1)
        res = self.client.post('/edge/fragments/1/annotations/', data=json.dumps(data),
                               content_type='application/json')
        data = dict(base_first=10, base_last=13, name='proD', type='promoter', strand=-1)
        res = self.client.post('/edge/fragments/1/annotations/', data=json.dumps(data),
                               content_type='application/json')
        res = self.client.get('/edge/fragments/1/annotations/')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [{
            "base_first": 2,
            "base_last": 8,
            "name": 'proC',
            "type": 'promoter',
            "strand": 1,
            "feature_full_length": 7,
            "feature_base_first": 1,
            "feature_base_last": 7,
        }, {
            "base_first": 10,
            "base_last": 13,
            "name": 'proD',
            "type": 'promoter',
            "strand": -1,
            "feature_full_length": 4,
            "feature_base_first": 1,
            "feature_base_last": 4,
        }])
        res = self.client.get('/edge/fragments/1/annotations/?f=8&l=14')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [{
            "base_first": 2,
            "base_last": 8,
            "name": 'proC',
            "type": 'promoter',
            "strand": 1,
            "feature_full_length": 7,
            "feature_base_first": 1,
            "feature_base_last": 7,
        }, {
            "base_first": 10,
            "base_last": 13,
            "name": 'proD',
            "type": 'promoter',
            "strand": -1,
            "feature_full_length": 4,
            "feature_base_first": 1,
            "feature_base_last": 4,
        }])
        res = self.client.get('/edge/fragments/1/annotations/?f=9&l=14')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [{
            "base_first": 10,
            "base_last": 13,
            "name": 'proD',
            "type": 'promoter',
            "strand": -1,
            "feature_full_length": 4,
            "feature_base_first": 1,
            "feature_base_last": 4,
        }])


class GenomeAnnotationsTest(TestCase):

    def setUp(self):
        res = self.client.post('/edge/genomes/', data=json.dumps(dict(name='foo', notes='bar')),
                               content_type='application/json')
        data = dict(name='chrI', sequence='AGCTAGCTTCGATCGA')
        res = self.client.post('/edge/genomes/1/fragments/', data=json.dumps(data),
                               content_type='application/json')
        self.fragment_data = json.loads(res.content)
        data = dict(base_first=2, base_last=9, name='proC', type='promoter', strand=1)
        res = self.client.post('/edge/fragments/1/annotations/', data=json.dumps(data),
                               content_type='application/json')

    def test_find_annotation(self):
        res = self.client.get('/edge/genomes/1/annotations/?q=proC')
        self.assertEquals(res.status_code, 200)
        self.assertEquals(json.loads(res.content), [[
            self.fragment_data, [
                {"base_first": 2,
                 "base_last": 9,
                 "name": 'proC',
                 "type": 'promoter',
                 "strand": 1,
                 "feature_full_length": 8,
                 "feature_base_first": 1,
                 "feature_base_last": 8}]
        ]])


class GenomeFragmentTest(TestCase):

    def setUp(self):
        res = self.client.post('/edge/genomes/', data=json.dumps(dict(name='foo', notes='bar')),
                               content_type='application/json')
        self.sequence = 'AGCTAGCTTCGATCGA'
        self.fragment_name = 'chrI'
        data = dict(name=self.fragment_name, sequence=self.sequence)
        res = self.client.post('/edge/genomes/1/fragments/', data=json.dumps(data),
                               content_type='application/json')
        self.uri = json.loads(res.content)['uri'].replace('/edge', '/edge/genomes/1')

    def test_insert_bases(self):
        data = dict(name='chrI_m', op='insert_bases', before_bp=3, sequence='GATACA')
        res = self.client.put(self.uri, data=json.dumps(data),
                              content_type='application/json')
        self.assertEquals(res.status_code, 201)
        expected = {
            "fragments": [{
                'circular': False,
                'id': 2,
                'length': len(self.sequence)+6,
                'name': self.fragment_name,
                'uri': '/edge/fragments/2/',
                'changes': [[1, len(self.sequence)+6]],
            }],
            "id": 2,
            "name": "chrI_m",
            "notes": None,
            "parent_id": 1,
            "uri": "/edge/genomes/2/"
        }
        self.assertEquals(json.loads(res.content), expected)
        res = self.client.get('/edge/fragments/2/sequence/')
        self.assertEquals(json.loads(res.content)['sequence'],
                          self.sequence[0:2]+'GATACA'+self.sequence[2:])

    def test_remove_bases(self):
        data = dict(name='chrI_m', op='remove_bases', before_bp=3, length=4)
        res = self.client.put(self.uri, data=json.dumps(data),
                              content_type='application/json')
        self.assertEquals(res.status_code, 201)
        self.assertEquals(json.loads(res.content), {
            "fragments": [{
                'circular': False,
                'id': 2,
                'length': len(self.sequence)-4,
                'name': self.fragment_name,
                'uri': '/edge/fragments/2/',
                'changes': [[1, len(self.sequence)-4]],
            }],
            "id": 2,
            "name": "chrI_m",
            "notes": None,
            "parent_id": 1,
            "uri": "/edge/genomes/2/"
        })
        res = self.client.get('/edge/fragments/2/sequence/')
        self.assertEquals(json.loads(res.content)['sequence'],
                          self.sequence[0:2]+self.sequence[6:])

    def test_replace_bases(self):
        data = dict(name='chrI_m', op='replace_bases', before_bp=3, length=4, sequence='GATACA')
        res = self.client.put(self.uri, data=json.dumps(data),
                              content_type='application/json')
        self.assertEquals(res.status_code, 201)
        self.assertEquals(json.loads(res.content), {
            "fragments": [{
                'circular': False,
                'id': 2,
                'length': len(self.sequence)+2,
                'name': self.fragment_name,
                'uri': '/edge/fragments/2/',
                'changes': [[1, len(self.sequence)+2]],
            }],
            "id": 2,
            "name": "chrI_m",
            "notes": None,
            "parent_id": 1,
            "uri": "/edge/genomes/2/"
        })
        res = self.client.get('/edge/fragments/2/sequence/')
        self.assertEquals(json.loads(res.content)['sequence'],
                          self.sequence[0:2]+'GATACA'+self.sequence[6:])

    def test_insert_fragment(self):
        sequence = 'GATACA'
        name = 'bar'
        res = self.client.post('/edge/fragments/', data=json.dumps(dict(name=name, sequence=sequence)),
                               content_type='application/json')
        fragment_id = json.loads(res.content)['id']

        data = dict(name='chrI_m', op='insert_fragment', before_bp=3, fragment_id=fragment_id)
        res = self.client.put(self.uri, data=json.dumps(data),
                              content_type='application/json')
        self.assertEquals(res.status_code, 201)
        self.assertEquals(json.loads(res.content), {
            "fragments": [{
                'circular': False,
                'id': 3,
                'length': len(self.sequence)+6,
                'name': self.fragment_name,
                'uri': '/edge/fragments/3/',
                'changes': [[1, len(self.sequence)+6]],
            }],
            "id": 2,
            "name": "chrI_m",
            "notes": None,
            "parent_id": 1,
            "uri": "/edge/genomes/2/"
        })
        res = self.client.get('/edge/fragments/3/sequence/')
        self.assertEquals(json.loads(res.content)['sequence'],
                          self.sequence[0:2]+'GATACA'+self.sequence[2:])

    def test_replace_with_fragment(self):
        sequence = 'GATACA'
        name = 'bar'
        res = self.client.post('/edge/fragments/', data=json.dumps(dict(name=name, sequence=sequence)),
                               content_type='application/json')
        fragment_id = json.loads(res.content)['id']

        data = dict(name='chrI_m', op='replace_with_fragment', before_bp=3, length=4,
                    fragment_id=fragment_id)
        res = self.client.put(self.uri, data=json.dumps(data),
                              content_type='application/json')
        self.assertEquals(res.status_code, 201)
        self.assertEquals(json.loads(res.content), {
            "fragments": [{
                'circular': False,
                'id': 3,
                'length': len(self.sequence)+2,
                'name': self.fragment_name,
                'uri': '/edge/fragments/3/',
                'changes': [[1, len(self.sequence)+2]],
            }],
            "id": 2,
            "name": "chrI_m",
            "notes": None,
            "parent_id": 1,
            "uri": "/edge/genomes/2/"
        })
        res = self.client.get('/edge/fragments/3/sequence/')
        self.assertEquals(json.loads(res.content)['sequence'],
                          self.sequence[0:2]+'GATACA'+self.sequence[6:])
