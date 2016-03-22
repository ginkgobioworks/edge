# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Fragment'
        db.create_table(u'edge_fragment', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('circular', self.gf('django.db.models.fields.BooleanField')()),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=256)),
            ('parent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['edge.Fragment'], null=True)),
            ('start_chunk', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['edge.Chunk'], null=True)),
        ))
        db.send_create_signal('edge', ['Fragment'])

        # Adding model 'Chunk'
        db.create_table(u'edge_chunk', (
            ('id', self.gf('django.db.models.fields.BigIntegerField')(primary_key=True)),
            ('initial_fragment', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['edge.Fragment'])),
            ('sequence', self.gf('django.db.models.fields.TextField')(null=True)),
        ))
        db.send_create_signal('edge', ['Chunk'])

        # Adding model 'Edge'
        db.create_table(u'edge_edge', (
            ('id', self.gf('django.db.models.fields.BigIntegerField')(primary_key=True)),
            ('from_chunk', self.gf('django.db.models.fields.related.ForeignKey')(related_name='out_edges', to=orm['edge.Chunk'])),
            ('fragment', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['edge.Fragment'])),
            ('to_chunk', self.gf('django.db.models.fields.related.ForeignKey')(related_name='in_edges', null=True, to=orm['edge.Chunk'])),
        ))
        db.send_create_signal('edge', ['Edge'])

        # Adding model 'Feature'
        db.create_table(u'edge_feature', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('type', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('strand', self.gf('django.db.models.fields.IntegerField')(null=True)),
            ('length', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal('edge', ['Feature'])

        # Adding model 'Chunk_Feature'
        db.create_table(u'edge_chunk_feature', (
            ('id', self.gf('django.db.models.fields.BigIntegerField')(primary_key=True)),
            ('chunk', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['edge.Chunk'])),
            ('feature', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['edge.Feature'])),
            ('feature_base_first', self.gf('django.db.models.fields.IntegerField')()),
            ('feature_base_last', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal('edge', ['Chunk_Feature'])

        # Adding model 'Fragment_Chunk_Location'
        db.create_table(u'edge_fragment_chunk_location', (
            ('id', self.gf('django.db.models.fields.BigIntegerField')(primary_key=True)),
            ('fragment', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['edge.Fragment'])),
            ('chunk', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['edge.Chunk'])),
            ('base_first', self.gf('django.db.models.fields.IntegerField')()),
            ('base_last', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal('edge', ['Fragment_Chunk_Location'])

        # Adding unique constraint on 'Fragment_Chunk_Location', fields ['fragment', 'chunk']
        db.create_unique(u'edge_fragment_chunk_location', ['fragment_id', 'chunk_id'])

        # Adding model 'Genome'
        db.create_table(u'edge_genome', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('name', self.gf('django.db.models.fields.CharField')(max_length=256)),
            ('parent', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['edge.Genome'], null=True)),
            ('notes', self.gf('django.db.models.fields.TextField')(null=True)),
        ))
        db.send_create_signal('edge', ['Genome'])

        # Adding model 'Genome_Fragment'
        db.create_table(u'edge_genome_fragment', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('genome', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['edge.Genome'])),
            ('fragment', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['edge.Fragment'])),
            ('inherited', self.gf('django.db.models.fields.BooleanField')()),
        ))
        db.send_create_signal('edge', ['Genome_Fragment'])


    def backwards(self, orm):
        # Removing unique constraint on 'Fragment_Chunk_Location', fields ['fragment', 'chunk']
        db.delete_unique(u'edge_fragment_chunk_location', ['fragment_id', 'chunk_id'])

        # Deleting model 'Fragment'
        db.delete_table(u'edge_fragment')

        # Deleting model 'Chunk'
        db.delete_table(u'edge_chunk')

        # Deleting model 'Edge'
        db.delete_table(u'edge_edge')

        # Deleting model 'Feature'
        db.delete_table(u'edge_feature')

        # Deleting model 'Chunk_Feature'
        db.delete_table(u'edge_chunk_feature')

        # Deleting model 'Fragment_Chunk_Location'
        db.delete_table(u'edge_fragment_chunk_location')

        # Deleting model 'Genome'
        db.delete_table(u'edge_genome')

        # Deleting model 'Genome_Fragment'
        db.delete_table(u'edge_genome_fragment')


    models = {
        'edge.chunk': {
            'Meta': {'object_name': 'Chunk'},
            'id': ('django.db.models.fields.BigIntegerField', [], {'primary_key': 'True'}),
            'initial_fragment': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Fragment']"}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True'})
        },
        'edge.chunk_feature': {
            'Meta': {'object_name': 'Chunk_Feature'},
            'chunk': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Chunk']"}),
            'feature': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Feature']"}),
            'feature_base_first': ('django.db.models.fields.IntegerField', [], {}),
            'feature_base_last': ('django.db.models.fields.IntegerField', [], {}),
            'id': ('django.db.models.fields.BigIntegerField', [], {'primary_key': 'True'})
        },
        'edge.edge': {
            'Meta': {'object_name': 'Edge'},
            'fragment': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Fragment']"}),
            'from_chunk': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'out_edges'", 'to': "orm['edge.Chunk']"}),
            'id': ('django.db.models.fields.BigIntegerField', [], {'primary_key': 'True'}),
            'to_chunk': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'in_edges'", 'null': 'True', 'to': "orm['edge.Chunk']"})
        },
        'edge.feature': {
            'Meta': {'object_name': 'Feature'},
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'length': ('django.db.models.fields.IntegerField', [], {}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'strand': ('django.db.models.fields.IntegerField', [], {'null': 'True'}),
            'type': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        'edge.fragment': {
            'Meta': {'object_name': 'Fragment'},
            'circular': ('django.db.models.fields.BooleanField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Fragment']", 'null': 'True'}),
            'start_chunk': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Chunk']", 'null': 'True'})
        },
        'edge.fragment_chunk_location': {
            'Meta': {'unique_together': "(('fragment', 'chunk'),)", 'object_name': 'Fragment_Chunk_Location'},
            'base_first': ('django.db.models.fields.IntegerField', [], {}),
            'base_last': ('django.db.models.fields.IntegerField', [], {}),
            'chunk': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Chunk']"}),
            'fragment': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Fragment']"}),
            'id': ('django.db.models.fields.BigIntegerField', [], {'primary_key': 'True'})
        },
        'edge.genome': {
            'Meta': {'object_name': 'Genome'},
            'fragments': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['edge.Fragment']", 'through': "orm['edge.Genome_Fragment']", 'symmetrical': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'notes': ('django.db.models.fields.TextField', [], {'null': 'True'}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Genome']", 'null': 'True'})
        },
        'edge.genome_fragment': {
            'Meta': {'object_name': 'Genome_Fragment'},
            'fragment': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Fragment']"}),
            'genome': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Genome']"}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'inherited': ('django.db.models.fields.BooleanField', [], {})
        }
    }

    complete_apps = ['edge']