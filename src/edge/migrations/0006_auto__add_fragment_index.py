# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Fragment_Index'
        db.create_table(u'edge_fragment_index', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('fragment', self.gf('django.db.models.fields.related.OneToOneField')(to=orm['edge.Fragment'], unique=True)),
            ('fresh', self.gf('django.db.models.fields.BooleanField')()),
            ('updated_on', self.gf('django.db.models.fields.DateTimeField')(null=True)),
        ))
        db.send_create_signal('edge', ['Fragment_Index'])


    def backwards(self, orm):
        # Deleting model 'Fragment_Index'
        db.delete_table(u'edge_fragment_index')


    models = {
        'edge.chunk': {
            'Meta': {'object_name': 'Chunk'},
            'id': ('django.db.models.fields.BigIntegerField', [], {'primary_key': 'True'}),
            'initial_fragment': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Fragment']", 'on_delete': 'models.PROTECT'}),
            'sequence': ('django.db.models.fields.TextField', [], {'null': 'True'})
        },
        'edge.chunk_feature': {
            'Meta': {'object_name': 'Chunk_Feature'},
            'chunk': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Chunk']", 'on_delete': 'models.PROTECT'}),
            'feature': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Feature']", 'on_delete': 'models.PROTECT'}),
            'feature_base_first': ('django.db.models.fields.IntegerField', [], {}),
            'feature_base_last': ('django.db.models.fields.IntegerField', [], {}),
            'id': ('django.db.models.fields.BigIntegerField', [], {'primary_key': 'True'})
        },
        'edge.edge': {
            'Meta': {'object_name': 'Edge'},
            'fragment': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Fragment']", 'on_delete': 'models.PROTECT'}),
            'from_chunk': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'out_edges'", 'on_delete': 'models.PROTECT', 'to': "orm['edge.Chunk']"}),
            'id': ('django.db.models.fields.BigIntegerField', [], {'primary_key': 'True'}),
            'to_chunk': ('django.db.models.fields.related.ForeignKey', [], {'related_name': "'in_edges'", 'null': 'True', 'on_delete': 'models.PROTECT', 'to': "orm['edge.Chunk']"})
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
            'created_on': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'null': 'True', 'blank': 'True'}),
            'est_length': ('django.db.models.fields.IntegerField', [], {'null': 'True', 'blank': 'True'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Fragment']", 'null': 'True', 'on_delete': 'models.PROTECT'}),
            'start_chunk': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Chunk']", 'null': 'True', 'on_delete': 'models.PROTECT'})
        },
        'edge.fragment_chunk_location': {
            'Meta': {'unique_together': "(('fragment', 'chunk'),)", 'object_name': 'Fragment_Chunk_Location', 'index_together': "(('fragment', 'base_last'), ('fragment', 'base_first'))"},
            'base_first': ('django.db.models.fields.IntegerField', [], {}),
            'base_last': ('django.db.models.fields.IntegerField', [], {}),
            'chunk': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Chunk']", 'on_delete': 'models.PROTECT'}),
            'fragment': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Fragment']", 'on_delete': 'models.PROTECT'}),
            'id': ('django.db.models.fields.BigIntegerField', [], {'primary_key': 'True'})
        },
        'edge.fragment_index': {
            'Meta': {'object_name': 'Fragment_Index'},
            'fragment': ('django.db.models.fields.related.OneToOneField', [], {'to': "orm['edge.Fragment']", 'unique': 'True'}),
            'fresh': ('django.db.models.fields.BooleanField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'updated_on': ('django.db.models.fields.DateTimeField', [], {'null': 'True'})
        },
        'edge.genome': {
            'Meta': {'object_name': 'Genome'},
            'created_on': ('django.db.models.fields.DateTimeField', [], {'auto_now_add': 'True', 'null': 'True', 'blank': 'True'}),
            'fragments': ('django.db.models.fields.related.ManyToManyField', [], {'to': "orm['edge.Fragment']", 'through': "orm['edge.Genome_Fragment']", 'symmetrical': 'False'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '256'}),
            'notes': ('django.db.models.fields.TextField', [], {'null': 'True', 'blank': 'True'}),
            'parent': ('django.db.models.fields.related.ForeignKey', [], {'to': "orm['edge.Genome']", 'null': 'True', 'on_delete': 'models.PROTECT'})
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