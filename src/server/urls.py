from django.conf.urls import patterns, include, url
from django.contrib import admin
from edge import views

admin.autodiscover()

urlpatterns = patterns('',
    url(r'^admin/', include(admin.site.urls)),
    url(r'^edge/', include('edge.urls')),
)
