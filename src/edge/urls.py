from django.conf.urls import url
from django.views.generic.base import TemplateView

from edge.views import (
    FragmentView,
    FragmentListView,
    FragmentSequenceView,
    FragmentAnnotationsView,
    GenomeView,
    GenomeListView,
    GenomeAnnotationsView,
    GenomeFragmentListView,
    GenomeBlastView,
    GenomePcrView,
    GenomeRecombinationView,
    GenomeCrisprDSBView,
    genome_export,
    genome_import,
)

urlpatterns = [

    # UI: previously we set / to redirect to a static page, permanently. That
    # was a bad idea. First, we really want to use a template to take advantage
    # of the staticfile facility from Django for managing asset versions.
    # Second, the permanent redirect is cached in browser permanently, so now
    # we have to add a redirect from the static page to the /ui/ URL, for any
    # browser that have the redirect memorized.

    url(r'^$', TemplateView.as_view(template_name='edge/edge.html'), name='edge-ui'),
    url(r'^ui/?$', TemplateView.as_view(template_name='edge/edge.html')),

    # APIs

    url(r'^import_genome/$', genome_import, name='import'),

    url(r'^genomes/$', GenomeListView.as_view(), name='genome_list'),
    url(r'^genomes/(?P<genome_id>\d+)/$', GenomeView.as_view(), name='genome'),

    url('^fragments/$', FragmentListView.as_view(), name='fragment_list'),
    url(r'^fragments/(?P<fragment_id>\d+)/$', FragmentView.as_view(), name='fragment'),

    url(r'^fragments/(?P<fragment_id>\d+)/sequence/$', FragmentSequenceView.as_view()),
    url(r'^fragments/(?P<fragment_id>\d+)/annotations/$', FragmentAnnotationsView.as_view()),

    url(r'^genomes/(?P<genome_id>\d+)/annotations/$', GenomeAnnotationsView.as_view()),
    url(r'^genomes/(?P<genome_id>\d+)/fragments/$', GenomeFragmentListView.as_view()),
    url(r'^genomes/(?P<genome_id>\d+)/blast/$', GenomeBlastView.as_view()),
    url(r'^genomes/(?P<genome_id>\d+)/pcr/$', GenomePcrView.as_view()),
    url(r'^genomes/(?P<genome_id>\d+)/recombination/$', GenomeRecombinationView.as_view()),
    url(r'^genomes/(?P<genome_id>\d+)/crispr/dsb/$', GenomeCrisprDSBView.as_view()),

    url(r'^genomes/(?P<genome_id>\d+)/export/$', genome_export),
]
