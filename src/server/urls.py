from django.conf.urls import patterns, include, url
from rest_framework import routers
from django.contrib import admin
from edge import views

admin.autodiscover()

router = routers.DefaultRouter()
router.register(r'fragments', views.FragmentViewSet)
router.register(r'genomes', views.GenomeViewSet)

urlpatterns = patterns('',
    url(r'^', include(router.urls)),
    url(r'^genomes-api/(?P<genome_id>[0-9]+)/(?P<option>[\b(export|import)\b]*)[/]*$', views.GenomeAPI),
    url(r'^admin/', include(admin.site.urls)),
    url(r'^edge/', include('edge.urls')),
    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework'))
)