# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import migrations, models
import edge.models.genome_updater
import edge.models.fragment_updater
import django.db.models.deletion
import edge.models.fragment_writer
import edge.models.fragment_annotator


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Chunk',
            fields=[
                ('id', models.BigIntegerField(serialize=False, primary_key=True)),
                ('sequence', models.TextField(null=True)),
            ],
        ),
        migrations.CreateModel(
            name='Chunk_Feature',
            fields=[
                ('id', models.BigIntegerField(serialize=False, primary_key=True)),
                ('feature_base_first', models.IntegerField()),
                ('feature_base_last', models.IntegerField()),
                ('chunk', models.ForeignKey(to='edge.Chunk', on_delete=django.db.models.deletion.PROTECT)),
            ],
        ),
        migrations.CreateModel(
            name='Edge',
            fields=[
                ('id', models.BigIntegerField(serialize=False, primary_key=True)),
            ],
        ),
        migrations.CreateModel(
            name='Feature',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=100)),
                ('type', models.CharField(max_length=100)),
                ('strand', models.IntegerField(null=True)),
                ('length', models.IntegerField()),
                ('_qualifiers', models.TextField(null=True, db_column=b'qualifiers')),
            ],
        ),
        migrations.CreateModel(
            name='Fragment',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('circular', models.BooleanField()),
                ('name', models.CharField(max_length=256)),
                ('est_length', models.IntegerField(null=True, verbose_name=b'Estimated length', blank=True)),
                ('created_on', models.DateTimeField(auto_now_add=True, verbose_name=b'Created', null=True)),
                ('active', models.BooleanField(default=True)),
                ('parent', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='edge.Fragment', null=True)),
                ('start_chunk', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='edge.Chunk', null=True)),
            ],
        ),
        migrations.CreateModel(
            name='Fragment_Chunk_Location',
            fields=[
                ('id', models.BigIntegerField(serialize=False, primary_key=True)),
                ('base_first', models.IntegerField()),
                ('base_last', models.IntegerField()),
                ('chunk', models.ForeignKey(to='edge.Chunk', on_delete=django.db.models.deletion.PROTECT)),
                ('fragment', models.ForeignKey(to='edge.Fragment', on_delete=django.db.models.deletion.PROTECT)),
            ],
        ),
        migrations.CreateModel(
            name='Fragment_Index',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('fresh', models.BooleanField()),
                ('updated_on', models.DateTimeField(null=True, verbose_name=b'Updated')),
                ('fragment', models.OneToOneField(to='edge.Fragment')),
            ],
        ),
        migrations.CreateModel(
            name='Genome',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('name', models.CharField(max_length=256)),
                ('notes', models.TextField(null=True, blank=True)),
                ('created_on', models.DateTimeField(auto_now_add=True, verbose_name=b'Created', null=True)),
                ('active', models.BooleanField(default=True)),
                ('blastdb', models.TextField(null=True, blank=True)),
            ],
            bases=(edge.models.genome_updater.Genome_Updater, models.Model),
        ),
        migrations.CreateModel(
            name='Genome_Fragment',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('inherited', models.BooleanField()),
                ('fragment', models.ForeignKey(to='edge.Fragment')),
                ('genome', models.ForeignKey(to='edge.Genome')),
            ],
        ),
        migrations.CreateModel(
            name='Operation',
            fields=[
                ('id', models.AutoField(verbose_name='ID', serialize=False, auto_created=True, primary_key=True)),
                ('type', models.IntegerField(choices=[(1, b'Homologous Recombination'), (2, b'CRISPR-Cas9 WT (Double Stranded Break)'), (3, b'PCR Product Sequence Verification')])),
                ('params', models.TextField(null=True, blank=True)),
                ('genome', models.ForeignKey(to='edge.Genome')),
            ],
        ),
        migrations.AddField(
            model_name='genome',
            name='fragments',
            field=models.ManyToManyField(to='edge.Fragment', through='edge.Genome_Fragment'),
        ),
        migrations.AddField(
            model_name='genome',
            name='parent',
            field=models.ForeignKey(related_name='children', on_delete=django.db.models.deletion.PROTECT, to='edge.Genome', null=True),
        ),
        migrations.AddField(
            model_name='feature',
            name='operation',
            field=models.ForeignKey(to='edge.Operation', null=True),
        ),
        migrations.AddField(
            model_name='edge',
            name='fragment',
            field=models.ForeignKey(to='edge.Fragment', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='edge',
            name='from_chunk',
            field=models.ForeignKey(related_name='out_edges', on_delete=django.db.models.deletion.PROTECT, to='edge.Chunk'),
        ),
        migrations.AddField(
            model_name='edge',
            name='to_chunk',
            field=models.ForeignKey(related_name='in_edges', on_delete=django.db.models.deletion.PROTECT, to='edge.Chunk', null=True),
        ),
        migrations.AddField(
            model_name='chunk_feature',
            name='feature',
            field=models.ForeignKey(to='edge.Feature', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.AddField(
            model_name='chunk',
            name='initial_fragment',
            field=models.ForeignKey(to='edge.Fragment', on_delete=django.db.models.deletion.PROTECT),
        ),
        migrations.CreateModel(
            name='Indexed_Fragment',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=(edge.models.fragment_writer.Fragment_Writer, edge.models.fragment_annotator.Fragment_Annotator, edge.models.fragment_updater.Fragment_Updater, 'edge.fragment'),
        ),
        migrations.CreateModel(
            name='Indexed_Genome',
            fields=[
            ],
            options={
                'proxy': True,
            },
            bases=('edge.genome',),
        ),
        migrations.AlterUniqueTogether(
            name='fragment_chunk_location',
            unique_together=set([('fragment', 'chunk')]),
        ),
        migrations.AlterIndexTogether(
            name='fragment_chunk_location',
            index_together=set([('fragment', 'base_first'), ('fragment', 'base_last')]),
        ),
    ]
