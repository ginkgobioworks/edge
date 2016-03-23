from django_assets import Bundle, register

js = Bundle('edge/src/js/*.js', filters='jsmin', output='edge/edge.min.js')
register('edge_js', js)

js = Bundle('edge/src/js/*.js', output='edge/edge.js')
register('edge', js)

css = Bundle('edge/src/css/app.css', output='edge/edge.css')
register('edge_css', css)

jst = Bundle('edge/src/partials/*.html', filters='jst', output='edge/edge_jst.js')
register('edge_jst', jst)
