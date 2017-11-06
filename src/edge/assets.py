from django_assets import Bundle, register

js = Bundle('edge/static/edge/src/js/*.js',
            output='edge/static/edge/assets/edge.js')
register('edge_js', js)

jsmin = Bundle(js, filters='jsmin',
               output='edge/static/edge/assets/edge.min.js')
register('edge_jsmin', jsmin)

css = Bundle('edge/static/edge/src/css/app.css',
             output='edge/static/edge/assets/edge.css')
register('edge_css', css)

jst = Bundle('edge/static/edge/src/partials/*.html', filters='jst',
             output='edge/static/edge/assets/edge_jst.js')
register('edge_jst', jst)
