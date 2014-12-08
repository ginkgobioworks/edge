import os
import json
import re
import tempfile
from django.test import TestCase
from edge.models import Operation


class GenomeModifyTest(TestCase):

    def setUp(self):
        res = self.client.post('/edge/genomes/', data=json.dumps(dict(name='foo', notes='bar')),
                               content_type='application/json')
        self.genome_uri = json.loads(res.content)['uri']
        m = re.match(r'^/edge/genomes/(\d+)/$', self.genome_uri)
        self.genome_id = int(m.group(1))
        self.genome_name = 'foo'

        self.sequence = 'AGCTAGCTTCGATCGA'
        self.fragment_name = 'chrI'
        data = dict(name=self.fragment_name, sequence=self.sequence)
        res = self.client.post(self.genome_uri+'fragments/', data=json.dumps(data),
                               content_type='application/json')

    def test_modify_creates_operation(self):
        self.assertEquals(Operation.objects.count(), 0)
        data = dict(genome_name='foo-bar', notes='blah', operation='CRISPR',
                    name='m1', circular=False, sequence='AGCT')
        res = self.client.post(self.genome_uri+'modify/', data=json.dumps(data),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)
        self.assertEquals(Operation.objects.count(), 1)

        child_uri = json.loads(res.content)['uri']
        m = re.match(r'^/edge/genomes/(\d+)/$', child_uri)
        child_id = int(m.group(1))
        self.assertNotEquals(child_id, self.genome_id)

        fragment_uri = json.loads(res.content)['fragments'][0]['uri']
        m = re.match(r'^/edge/fragments/(\d+)/$', fragment_uri)
        fragment_id = int(m.group(1))

        op_fragment_id = Operation.objects.all()[0].fragment_id
        op_fragment_uri = json.loads(res.content)['operations'][0]['fragment']['uri']

        expected = {
            "fragments": [{
                'circular': False,
                'id': fragment_id,
                'length': len(self.sequence),
                'name': self.fragment_name,
                'uri': fragment_uri,
                "parent_id": None
            }],
            "id": child_id,
            "name": "foo-bar",
            "notes": "blah",
            "parent_id": self.genome_id,
            "parent_name": self.genome_name,
            "uri": child_uri,
            "operations": [{
                "type": 'CRISPR',
                "notes": None,
                "fragment": {
                    'id': op_fragment_id,
                    'name': 'm1',
                    'circular': False,
                    'length': 4,
                    'uri': op_fragment_uri,
                    "parent_id": None
                }
            }]
        }
        self.assertEquals(json.loads(res.content), expected)
        res = self.client.get(fragment_uri+'sequence/')
        self.assertEquals(json.loads(res.content)['sequence'], self.sequence)

    def test_modify_requires_operation(self):
        data = dict(genome_name='foo-bar', notes='blah', operation='CRISPR',
                    name='m1', circular=False, sequence='AGCT')
        res = self.client.post(self.genome_uri+'modify/', data=json.dumps(data),
                               content_type='application/json')
        self.assertEquals(res.status_code, 201)

        data = dict(genome_name='foo-bar', notes='blah',
                    name='m1', circular=False, sequence='AGCT')
        res = self.client.post(self.genome_uri+'modify/', data=json.dumps(data),
                               content_type='application/json')
        self.assertEquals(res.status_code, 400)
