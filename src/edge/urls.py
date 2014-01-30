from django.conf.urls import *
from edge.views import *

urlpatterns = patterns('',
  url('^fragments/$', FragmentListView.as_view(), name='fragment_list'),
  url('^fragments/(?P<fragment_id>\d+)/$',
      FragmentView.as_view(), name='fragment'),
  url('^fragments/(?P<fragment_id>\d+)/sequence/$',
      FragmentSequenceView.as_view(), name='fragment_sequence'),
  url('^fragments/(?P<fragment_id>\d+)/annotations/$',
      FragmentAnnotationsView.as_view(), name='fragment_annotations'),
  url('^genomes/$',
      GenomeListView.as_view(), name='genome_list'),
  url('^genomes/(?P<genome_id>\d+)/annotations/$',
      GenomeAnnotationsView.as_view(), name='genome_annotations'),
  url('^genomes/(?P<genome_id>\d+)/fragments/$',
      GenomeFragmentListView.as_view(), name='genome_fragment_list'),
  url('^genomes/(?P<genome_id>\d+)/fragments/(?P<fragment_id>\d+)/$',
      GenomeFragmentView.as_view(), name='genome_fragment'),
  url('^genomes/(?P<genome_id>\d+)/$',
      GenomeView.as_view(), name='genome'),
)
