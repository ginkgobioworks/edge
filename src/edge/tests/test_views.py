import os
import json
import unittest
import edge_srv
import tempfile


class GenomeListTest(unittest.TestCase):

    def setUp(self):
        self.db_fd, edge_srv.app.config['DATABASE'] = tempfile.mkstemp()
        edge_srv.app.config['TESTING'] = True
        self.app = edge_srv.app.test_client()
        edge_srv.init_db()

    def tearDown(self):
        os.close(self.db_fd)
        os.unlink(edge_srv.app.config['DATABASE'])

    def test_empty_db(self):
        res = self.app.get('/genomes')
        self.assertEquals(res.status, '200 OK')
        self.assertEquals(json.loads(res.data), [])

    def test_add_genome(self):
        res = self.app.post('/genomes', data=json.dumps(dict(name='foo', notes='bar')),
                            content_type='application/json')
        self.assertEquals(res.status, '201 CREATED')
        self.assertEquals(json.loads(res.data), {
            "fragments": [],
            "id": 1,
            "name": "foo",
            "notes": "bar",
            "parent_id": None,
            "uri": "/genomes/1",
        })

    def test_can_use_uri_from_add_genome_to_fetch_genome(self):
        res = self.app.post('/genomes', data=json.dumps(dict(name='foo', notes='bar')),
                            content_type='application/json')
        re2 = self.app.get(json.loads(res.data)['uri'])
        self.assertEquals(re2.status, '200 OK')
        self.assertEquals(json.loads(res.data), json.loads(re2.data))


class GenomeTest(unittest.TestCase):

    def setUp(self):
        self.db_fd, edge_srv.app.config['DATABASE'] = tempfile.mkstemp()
        edge_srv.app.config['TESTING'] = True
        self.app = edge_srv.app.test_client()
        edge_srv.init_db()
        res = self.app.post('/genomes', data=json.dumps(dict(name='foo', notes='bar')),
                            content_type='application/json')

    def tearDown(self):
        os.close(self.db_fd)
        os.unlink(edge_srv.app.config['DATABASE'])

    def test_get_genome(self):
        res = self.app.get('/genomes/1')
        self.assertEquals(res.status, '200 OK')
        self.assertEquals(json.loads(res.data), {
            "fragments": [],
            "id": 1,
            "name": "foo",
            "notes": "bar",
            "parent_id": None,
            "uri": '/genomes/1',
        })

    def test_returns_404_if_genome_does_not_exist(self):
        res = self.app.get('/genomes/100')
        self.assertEquals(res.status, '404 NOT FOUND')

    def test_add_non_circular_fragment(self):
        data = dict(name='chrI', sequence='AGCTAGCTTCGATCGA')
        res = self.app.post('/genomes/1/fragments', data=json.dumps(data),
                            content_type='application/json')
        self.assertEquals(res.status, '201 CREATED')
        self.assertEquals(json.loads(res.data), {
            "circular": False,
            "id": 1,
            "length": len(data['sequence']),
            "name": data['name'],
            "uri": "/fragments/1"
        })

    def test_add_circular_fragment(self):
        data = dict(name='chrI', sequence='AGCTAGCTTCGATCGA', circular=True)
        res = self.app.post('/genomes/1/fragments', data=json.dumps(data),
                            content_type='application/json')
        self.assertEquals(res.status, '201 CREATED')
        self.assertEquals(json.loads(res.data), {
            "circular": True,
            "id": 1,
            "length": len(data['sequence']),
            "name": data['name'],
            "uri": "/fragments/1"
        })

    def test_get_genome_returns_fragments(self):
        data = dict(name='chrI', sequence='AGCTAGCTTCGATCGA', circular=True)
        res = self.app.post('/genomes/1/fragments', data=json.dumps(data),
                            content_type='application/json')
        res = self.app.get('/genomes/1')
        self.assertEquals(json.loads(res.data)['fragments'], [{
            "circular": True,
            "id": 1,
            "length": len(data['sequence']),
            "name": data['name'],
            "uri": "/fragments/1"
        }])

    def test_can_use_uri_from_add_fragment_to_fetch_fragment(self):
        data = dict(name='chrI', sequence='AGCTAGCTTCGATCGA')
        res = self.app.post('/genomes/1/fragments', data=json.dumps(data),
                            content_type='application/json')
        re2 = self.app.get(json.loads(res.data)['uri'])
        self.assertEquals(re2.status, '200 OK')
        self.assertEquals(json.loads(res.data), json.loads(re2.data))


class FragmentTest(unittest.TestCase):

    def setUp(self):
        self.db_fd, edge_srv.app.config['DATABASE'] = tempfile.mkstemp()
        edge_srv.app.config['TESTING'] = True
        self.app = edge_srv.app.test_client()
        edge_srv.init_db()
        self.sequence = 'AGCTAGCTTCGATCGA'
        self.name = 'foo'
        res = self.app.post('/fragments',
                            data=json.dumps(dict(name=self.name, sequence=self.sequence)),
                            content_type='application/json')

    def tearDown(self):
        os.close(self.db_fd)
        os.unlink(edge_srv.app.config['DATABASE'])

    def test_get_all_non_genomic_fragments(self):
        res = self.app.get('/fragments')
        self.assertEquals(res.status, '200 OK')
        self.assertEquals(json.loads(res.data), [{
            "id": 1,
            "name": self.name,
            "length": len(self.sequence),
            "uri": '/fragments/1',
            "circular": False,
        }])

    def test_get_fragment(self):
        res = self.app.get('/fragments/1')
        self.assertEquals(res.status, '200 OK')
        self.assertEquals(json.loads(res.data), {
            "id": 1,
            "name": self.name,
            "length": len(self.sequence),
            "uri": '/fragments/1',
            "circular": False,
        })

    def test_get_fragment_sequence(self):
        res = self.app.get('/fragments/1/sequence')
        self.assertEquals(res.status, '200 OK')
        self.assertEquals(json.loads(res.data), {
            'sequence': self.sequence,
            'first_bp': 1,
            'last_bp': len(self.sequence)
        })

    def test_get_fragment_sequence_by_position(self):
        res = self.app.get('/fragments/1/sequence?f=3&l=10')
        self.assertEquals(res.status, '200 OK')
        self.assertEquals(json.loads(res.data), {
            'sequence': self.sequence[2:10],
            'first_bp': 3,
            'last_bp': 10
        })

    def test_returns_404_if_fragment_does_not_exist(self):
        res = self.app.get('/fragments/100')
        self.assertEquals(res.status, '404 NOT FOUND')

    def test_add_annotation_on_forward_strand(self):
        data = dict(first_bp=2, last_bp=9, name='proC', type='promoter', strand=1)
        res = self.app.post('/fragments/1/annotations', data=json.dumps(data),
                            content_type='application/json')
        self.assertEquals(res.status, '201 CREATED')
        self.assertEquals(json.loads(res.data), {})
        res = self.app.get('/fragments/1/annotations')
        self.assertEquals(res.status, '200 OK')
        self.assertEquals(json.loads(res.data), [{
            "first_bp": 2,
            "last_bp": 9,
            "name": 'proC',
            "type": 'promoter',
            "strand": 1,
            "annotation_full_length": 8,
            "annotation_first_bp": 1,
            "annotation_last_bp": 8,
        }])

    def test_add_annotation_on_reverse_strand(self):
        data = dict(first_bp=3, last_bp=10, name='proC', type='promoter', strand=-1)
        res = self.app.post('/fragments/1/annotations', data=json.dumps(data),
                            content_type='application/json')
        self.assertEquals(res.status, '201 CREATED')
        self.assertEquals(json.loads(res.data), {})
        res = self.app.get('/fragments/1/annotations')
        self.assertEquals(res.status, '200 OK')
        self.assertEquals(json.loads(res.data), [{
            "first_bp": 3,
            "last_bp": 10,
            "name": 'proC',
            "type": 'promoter',
            "strand": -1,
            "annotation_full_length": 8,
            "annotation_first_bp": 1,
            "annotation_last_bp": 8,
        }])

    def test_get_multiple_annotations(self):
        data = dict(first_bp=2, last_bp=9, name='proC', type='promoter', strand=1)
        res = self.app.post('/fragments/1/annotations', data=json.dumps(data),
                            content_type='application/json')
        data = dict(first_bp=3, last_bp=10, name='proD', type='promoter', strand=-1)
        res = self.app.post('/fragments/1/annotations', data=json.dumps(data),
                            content_type='application/json')
        res = self.app.get('/fragments/1/annotations')
        self.assertEquals(res.status, '200 OK')
        self.assertEquals(json.loads(res.data), [{
            "first_bp": 2,
            "last_bp": 9,
            "name": 'proC',
            "type": 'promoter',
            "strand": 1,
            "annotation_full_length": 8,
            "annotation_first_bp": 1,
            "annotation_last_bp": 8,
        }, {
            "first_bp": 3,
            "last_bp": 10,
            "name": 'proD',
            "type": 'promoter',
            "strand": -1,
            "annotation_full_length": 8,
            "annotation_first_bp": 1,
            "annotation_last_bp": 8,
        }])

    def test_get_annotations_by_region(self):
        data = dict(first_bp=2, last_bp=8, name='proC', type='promoter', strand=1)
        res = self.app.post('/fragments/1/annotations', data=json.dumps(data),
                            content_type='application/json')
        data = dict(first_bp=10, last_bp=13, name='proD', type='promoter', strand=-1)
        res = self.app.post('/fragments/1/annotations', data=json.dumps(data),
                            content_type='application/json')
        res = self.app.get('/fragments/1/annotations')
        self.assertEquals(res.status, '200 OK')
        self.assertEquals(json.loads(res.data), [{
            "first_bp": 2,
            "last_bp": 8,
            "name": 'proC',
            "type": 'promoter',
            "strand": 1,
            "annotation_full_length": 7,
            "annotation_first_bp": 1,
            "annotation_last_bp": 7,
        }, {
            "first_bp": 10,
            "last_bp": 13,
            "name": 'proD',
            "type": 'promoter',
            "strand": -1,
            "annotation_full_length": 4,
            "annotation_first_bp": 1,
            "annotation_last_bp": 4,
        }])
        res = self.app.get('/fragments/1/annotations?f=8&l=14')
        self.assertEquals(res.status, '200 OK')
        self.assertEquals(json.loads(res.data), [{
            "first_bp": 2,
            "last_bp": 8,
            "name": 'proC',
            "type": 'promoter',
            "strand": 1,
            "annotation_full_length": 7,
            "annotation_first_bp": 1,
            "annotation_last_bp": 7,
        }, {
            "first_bp": 10,
            "last_bp": 13,
            "name": 'proD',
            "type": 'promoter',
            "strand": -1,
            "annotation_full_length": 4,
            "annotation_first_bp": 1,
            "annotation_last_bp": 4,
        }])
        res = self.app.get('/fragments/1/annotations?f=9&l=14')
        self.assertEquals(res.status, '200 OK')
        self.assertEquals(json.loads(res.data), [{
            "first_bp": 10,
            "last_bp": 13,
            "name": 'proD',
            "type": 'promoter',
            "strand": -1,
            "annotation_full_length": 4,
            "annotation_first_bp": 1,
            "annotation_last_bp": 4,
        }])


class GenomeAnnotationsTest(unittest.TestCase):

    def setUp(self):
        self.db_fd, edge_srv.app.config['DATABASE'] = tempfile.mkstemp()
        edge_srv.app.config['TESTING'] = True
        self.app = edge_srv.app.test_client()
        edge_srv.init_db()
        res = self.app.post('/genomes', data=json.dumps(dict(name='foo', notes='bar')),
                            content_type='application/json')
        data = dict(name='chrI', sequence='AGCTAGCTTCGATCGA')
        res = self.app.post('/genomes/1/fragments', data=json.dumps(data),
                            content_type='application/json')
        self.fragment_data = json.loads(res.data)
        data = dict(first_bp=2, last_bp=9, name='proC', type='promoter', strand=1)
        res = self.app.post('/fragments/1/annotations', data=json.dumps(data),
                            content_type='application/json')

    def tearDown(self):
        os.close(self.db_fd)
        os.unlink(edge_srv.app.config['DATABASE'])

    def test_find_annotation(self):
        res = self.app.get('/genomes/1/annotations?q=proC')
        self.assertEquals(res.status, '200 OK')
        self.assertEquals(json.loads(res.data), [[
            self.fragment_data, [
                {"first_bp": 2,
                 "last_bp": 9,
                 "name": 'proC',
                 "type": 'promoter',
                 "strand": 1,
                 "annotation_full_length": 8,
                 "annotation_first_bp": 1,
                 "annotation_last_bp": 8}]
        ]])


class GenomeFragmentTest(unittest.TestCase):

    def setUp(self):
        self.db_fd, edge_srv.app.config['DATABASE'] = tempfile.mkstemp()
        edge_srv.app.config['TESTING'] = True
        self.app = edge_srv.app.test_client()
        edge_srv.init_db()

        uri = '/genomes/1'
        res = self.app.post('/genomes', data=json.dumps(dict(name='foo', notes='bar')),
                            content_type='application/json')
        self.sequence = 'AGCTAGCTTCGATCGA'
        self.fragment_name = 'chrI'
        data = dict(name=self.fragment_name, sequence=self.sequence)
        res = self.app.post(uri+'/fragments', data=json.dumps(data),
                            content_type='application/json')
        self.uri = uri+json.loads(res.data)['uri']

    def tearDown(self):
        os.close(self.db_fd)
        os.unlink(edge_srv.app.config['DATABASE'])

    def test_insert_bases(self):
        data = dict(name='chrI_m', op='insert_bases', before_bp=3, sequence='GATACA')
        res = self.app.put(self.uri, data=json.dumps(data),
                           content_type='application/json')
        self.assertEquals(res.status, '201 CREATED')
        print json.loads(res.data)
        self.assertEquals(json.loads(res.data), {
            "fragments": [{
                'circular': False,
                'id': 2,
                'length': len(self.sequence)+6,
                'name': self.fragment_name,
                'uri': '/fragments/2',
                'changes': [[1, len(self.sequence)+6]],
            }],
            "id": 2,
            "name": "chrI_m",
            "notes": None,
            "parent_id": 1,
            "uri": "/genomes/2"
        })
        res = self.app.get('/fragments/2/sequence')
        self.assertEquals(json.loads(res.data)['sequence'],
                          self.sequence[0:2]+'GATACA'+self.sequence[2:])

    def test_remove_bases(self):
        data = dict(name='chrI_m', op='remove_bases', before_bp=3, length=4)
        res = self.app.put(self.uri, data=json.dumps(data),
                           content_type='application/json')
        self.assertEquals(res.status, '201 CREATED')
        self.assertEquals(json.loads(res.data), {
            "fragments": [{
                'circular': False,
                'id': 2,
                'length': len(self.sequence)-4,
                'name': self.fragment_name,
                'uri': '/fragments/2',
                'changes': [[1, len(self.sequence)-4]],
            }],
            "id": 2,
            "name": "chrI_m",
            "notes": None,
            "parent_id": 1,
            "uri": "/genomes/2"
        })
        res = self.app.get('/fragments/2/sequence')
        self.assertEquals(json.loads(res.data)['sequence'],
                          self.sequence[0:2]+self.sequence[6:])

    def test_replace_bases(self):
        data = dict(name='chrI_m', op='replace_bases', before_bp=3, length=4, sequence='GATACA')
        res = self.app.put(self.uri, data=json.dumps(data),
                           content_type='application/json')
        self.assertEquals(res.status, '201 CREATED')
        self.assertEquals(json.loads(res.data), {
            "fragments": [{
                'circular': False,
                'id': 2,
                'length': len(self.sequence)+2,
                'name': self.fragment_name,
                'uri': '/fragments/2',
                'changes': [[1, len(self.sequence)+2]],
            }],
            "id": 2,
            "name": "chrI_m",
            "notes": None,
            "parent_id": 1,
            "uri": "/genomes/2"
        })
        res = self.app.get('/fragments/2/sequence')
        self.assertEquals(json.loads(res.data)['sequence'],
                          self.sequence[0:2]+'GATACA'+self.sequence[6:])

    def test_insert_fragment(self):
        sequence = 'GATACA'
        name = 'bar'
        res = self.app.post('/fragments', data=json.dumps(dict(name=name, sequence=sequence)),
                            content_type='application/json')
        fragment_id = json.loads(res.data)['id']

        data = dict(name='chrI_m', op='insert_fragment', before_bp=3, fragment_id=fragment_id)
        res = self.app.put(self.uri, data=json.dumps(data),
                           content_type='application/json')
        self.assertEquals(res.status, '201 CREATED')
        self.assertEquals(json.loads(res.data), {
            "fragments": [{
                'circular': False,
                'id': 3,
                'length': len(self.sequence)+6,
                'name': self.fragment_name,
                'uri': '/fragments/3',
                'changes': [[1, len(self.sequence)+6]],
            }],
            "id": 2,
            "name": "chrI_m",
            "notes": None,
            "parent_id": 1,
            "uri": "/genomes/2"
        })
        res = self.app.get('/fragments/3/sequence')
        self.assertEquals(json.loads(res.data)['sequence'],
                          self.sequence[0:2]+'GATACA'+self.sequence[2:])

    def test_replace_with_fragment(self):
        sequence = 'GATACA'
        name = 'bar'
        res = self.app.post('/fragments', data=json.dumps(dict(name=name, sequence=sequence)),
                            content_type='application/json')
        fragment_id = json.loads(res.data)['id']

        data = dict(name='chrI_m', op='replace_with_fragment', before_bp=3, length=4,
                    fragment_id=fragment_id)
        res = self.app.put(self.uri, data=json.dumps(data),
                           content_type='application/json')
        self.assertEquals(res.status, '201 CREATED')
        self.assertEquals(json.loads(res.data), {
            "fragments": [{
                'circular': False,
                'id': 3,
                'length': len(self.sequence)+2,
                'name': self.fragment_name,
                'uri': '/fragments/3',
                'changes': [[1, len(self.sequence)+2]],
            }],
            "id": 2,
            "name": "chrI_m",
            "notes": None,
            "parent_id": 1,
            "uri": "/genomes/2"
        })
        res = self.app.get('/fragments/3/sequence')
        self.assertEquals(json.loads(res.data)['sequence'],
                          self.sequence[0:2]+'GATACA'+self.sequence[6:])


if __name__ == '__main__':
    unittest.main()
