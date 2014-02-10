from django_assets import Bundle, register

js = Bundle('edge/js/*.js', filters='jsmin', output='edge/static/edge/edge.min.js')
register('edge_js', js)

css = Bundle('edge/css/app.css',
             output='edge/static/edge/edge.bundle.css')
register('edge_css', css)

