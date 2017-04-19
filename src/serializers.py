from edge.models import *
from rest_framework import serializers

class FragmentSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Fragment
        fields = ('url', 'name', 'parent', 'circular', 'est_length', 'created_on')

class GenomeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Genome
        fields = ('url', 'name', 'fragments', 'created_on') 