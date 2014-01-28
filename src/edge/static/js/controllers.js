'use strict';

function edgeAnnotationDisplayName(annotation) {
    var s = annotation['name']
    if (s === undefined || s == '') { s = annotation['type']; }
    if (annotation['annotation_first_bp'] != 1 ||
        annotation['annotation_last_bp'] != annotation['annotation_full_length']) {
        s += '[';
        if (annotation['annotation_first_bp'] != 1) {
            s += annotation['annotation_first_bp'];
        }
        s += ':';
        if (annotation['annotation_last_bp'] != annotation['annotation_full_length']) {
            s += annotation['annotation_last_bp']+'/'+annotation['annotation_full_length'];
        }
        s += ']';
    }
    return s;
}

function edgeAnnotationSummaryCSS(annotation) {
    var css = ['annotation-summary'];
    if (annotation['type'] == 'gene') { css.push('btn-danger'); }
    else if (annotation['type'] == 'pseudogene') { css.push('btn-warning'); }
    else if (annotation['type'] == 'promoter') { css.push('btn-primary'); }
    else if (annotation['type'] == 'RBS') { css.push('btn-success'); }
    else { css.push('btn-info'); }
    css.push('btn btn-xs');
    if (annotation['strand'] > 0) { css.push('annotation-strand-fwd'); }
    else { css.push('annotation-strand-rev'); }
    return css.join(' ');
}

function edgeAnnotationZoomCSS(annotation) {
    var css = ['annotation-zoom'];
    if (annotation['type'] == 'gene') { css.push('btn-danger'); }
    else if (annotation['type'] == 'pseudogene') { css.push('btn-warning'); }
    else if (annotation['type'] == 'promoter') { css.push('btn-primary'); }
    else if (annotation['type'] == 'RBS') { css.push('btn-success'); }
    else { css.push('btn-info'); }
    css.push('btn btn-xs');
    if (annotation['strand'] > 0) { css.push('annotation-strand-fwd'); }
    else { css.push('annotation-strand-rev'); }
    return css.join(' ');
}

function edgeFetchChanges($http, fragment_id, changed_loc, cb) {
    var url = '/fragments/'+fragment_id+'/annotations?f='+changed_loc[0]+'&l='+changed_loc[1];
    $http.get(url).success(function(data) {
        var annotations = [];
        data.forEach(function(annotation) {
            var s = edgeAnnotationDisplayName(annotation);
            annotations.push(s);
        });
        cb(annotations.join(', '));
    });
}

function GenomeListController($scope, $http) {
    function fetchAndSaveParent(genome) {
        var url = '/genomes/'+genome['parent_id'];
        $http.get(url).success(function(data) {
            genome['parent'] = data['name'];
        });
    }

    $http.get('/genomes').success(function(data) {
        $scope.genomes = data;
        data.forEach(function(genome) {
            if (genome['parent_id']) { fetchAndSaveParent(genome); }
        });
    });

    $scope.orderProp = 'id';
}

function FragmentListController($scope, $http) {
    $http.get('/fragments').success(function(data) {
        $scope.fragments = data;
    });

    $scope.orderProp = 'id';
    $scope.add_fragment_error = undefined;

    $scope.addFragment = function(fragment) {
        var data = JSON.stringify(fragment);
        $http.post('/fragments', data).
            success(function(data) {
                $scope.fragments.push(data);
                $scope.add_fragment_error = undefined;
            }).
            error(function(data, status, headers, config) {
                $scope.add_fragment_error = data;
            });
    }
}

function GenomeDetailController($scope, $routeParams, $http) {
    $scope.genomeId = $routeParams.genomeId;
    $scope.changes = [];
    $scope.parent = undefined;

    function fetchParent(parent_id) {
        var url = '/genomes/'+parent_id;
        $http.get(url).success(function(data) {
            $scope.parent = '<a href="#/genomes/'+data['id']+'">Genome '+data['id']+': '+data['name']+'</a>';
        });
    }

    /* driver */
    $http.get('/genomes/'+$scope.genomeId).success(function(data) {
        $scope.genome = data;
        if (data['parent_id']) { fetchParent(data['parent_id']); }
        data['fragments'].forEach(function(fragment) {
            if (fragment['changes']) {
                fragment['changes'].forEach(function(changed_loc) {
                    edgeFetchChanges($http, fragment['id'], changed_loc, function(desc) {
                        $scope.changes.push(fragment['name']+': '+desc);
                    });
                });
            }
        });
    });
}

function FragmentController($scope, $routeParams, $http) {
    var DEFAULT_ZOOM = 5000;

    $scope.featureTypes = edgeFeatureTypes;
    $scope.fragmentId = $routeParams.fragmentId;
    $scope.fragment = undefined;
    $scope.annotations = []
    $scope.display_summary = true;
    $scope.summary_annotations = [];
    $scope.zoom = {};

    function layoutAnnotations(annotations, first_bp, last_bp) {
        var full_length = last_bp-first_bp+1;
        var display = [];
        var layers = [];

        annotations.forEach(function(a) {
            var d = {};
            d['annotation'] = a;
            d['css'] = edgeAnnotationZoomCSS(a);

            var left = (a['first_bp']-first_bp)*100.0/full_length;
            d['left'] = ''+left+'%';
            var width = (a['last_bp']-a['first_bp']+1)*100.0/full_length;
            d['width'] = ''+width+'%';

            /* find a layer to put the feature */
            var found_layer_i = undefined;
            for (var i=0; i<layers.length; i++) {
                var layer = layers[i];
                var overlap = false;
                /* anything in this layer overlap? */
                layer.forEach(function(layer_a) {
                    if (layer_a['last_bp']+100 >= a['first_bp']) { overlap = true; } 
                });
                if (overlap == false) {
                    found_layer_i = i;
                    layer.push(a);
                    break;
                }
            }
            if (found_layer_i === undefined) { 
                found_layer_i = layers.length;
                layers.push([a]);
            }
            d['css'] += ' annotation-layer-'+found_layer_i;
            display.push(d);
        });
        return display;
    }

    $scope.fetchSequence = function() {
        if ($scope.display_summary == false) {
            var f = $scope.zoom['first_bp'];
            var l = $scope.zoom['last_bp'];
            $http.get('/fragments/'+$scope.fragmentId+'/sequence?f='+f+'&l='+l).success(function(data) {
                $scope.zoom['sequence'] = data;
                $scope.zoom['has_sequence'] = true;
            });
        }
    }

    $scope.zoomRefresh = function(first_bp, last_bp) {
        var zoom_annotations = [];
        var min_bp = undefined;
        var max_bp = undefined;
        $scope.annotations.forEach(function(a) {
            if (a['last_bp'] > first_bp && a['first_bp'] < last_bp) {
                zoom_annotations.push(a);
                if (min_bp === undefined || a['first_bp'] < min_bp) {
                    min_bp = a['first_bp'];
                }
                if (max_bp === undefined || a['last_bp'] > max_bp) {
                    max_bp = a['last_bp'];
                }
            }
        });
        if (zoom_annotations.length > 0) {
            zoom_annotations.sort(function(a,b) { return a['first_bp']-b['first_bp']; });
            $scope.zoom['first_bp'] = Math.min(min_bp, first_bp);
            $scope.zoom['last_bp'] = Math.max(max_bp, last_bp);
        }
        else {
            $scope.zoom['first_bp'] = first_bp;
            $scope.zoom['last_bp'] = last_bp;
        }
        $scope.zoom['annotations'] = zoom_annotations;
        var display = layoutAnnotations(zoom_annotations, $scope.zoom['first_bp'], $scope.zoom['last_bp']);
        $scope.zoom['display'] = display;
        $scope.zoom['has_sequence'] = false;
    }

    $scope.zoomAt = function(annotation) {
        // reset
        $scope.zoom = {};
        $scope.display_summary = false;
        // get DEFAULT_ZOOM bps upstream and downstream from annotation
        var first_bp = annotation['first_bp']-DEFAULT_ZOOM;
        if (first_bp < 1) { first_bp = 1; }
        var last_bp = annotation['last_bp']+DEFAULT_ZOOM;
        if (last_bp > $scope.fragment['length']) { last_bp = $scope.fragment['length']; }
        $scope.zoomRefresh(first_bp, last_bp);
    }

    $scope.zoomMoveRight = function() {
        var cur_zoom = $scope.zoom['last_bp']-$scope.zoom['first_bp']+1;
        var first_bp = $scope.zoom['last_bp'];
        var last_bp = first_bp+cur_zoom;
        if (last_bp > $scope.fragment['length']) { last_bp = $scope.fragment['length']; }
        if (first_bp < last_bp) {
            $scope.zoomRefresh(first_bp, last_bp);
        }
    }

    $scope.zoomMoveLeft = function() {
        var cur_zoom = $scope.zoom['last_bp']-$scope.zoom['first_bp']+1;
        if ($scope.zoom['first_bp'] > 1) {
            var last_bp = $scope.zoom['first_bp'];
            var first_bp = last_bp-cur_zoom;
            if (first_bp < 1) { first_bp = 1; }
            if (first_bp < last_bp) {
                $scope.zoomRefresh(first_bp, last_bp);
            }
        }
    }

    $scope.zoomOut = function() {
        var first_bp = $scope.zoom['first_bp']-DEFAULT_ZOOM;
        var last_bp = $scope.zoom['last_bp']+DEFAULT_ZOOM;
        if (first_bp < 1) { first_bp = 1; }
        if (last_bp > $scope.fragment['length']) { last_bp = $scope.fragment['length']; }
        $scope.zoomRefresh(first_bp, last_bp);
    }

    $scope.zoomIn = function() {
        if ($scope.zoom['last_bp']-$scope.zoom['first_bp'] > 2*DEFAULT_ZOOM) {
            var first_bp = $scope.zoom['first_bp']+DEFAULT_ZOOM;
            var last_bp = $scope.zoom['last_bp']-DEFAULT_ZOOM;
            $scope.zoomRefresh(first_bp, last_bp);
        }
    }

    $scope.showSummary = function() {
        $scope.zoom = {
            'first_bp': 1,
            'last_bp': $scope.fragment['length'],
            'annotations': $scope.annotations
        };
        $scope.display_summary = true;
    }

    $scope.annotate_error = undefined;
    $scope.addAnnotation = function(annotation) {
        var data = JSON.stringify(annotation);
        $http.post('/fragments/'+$scope.fragmentId+'/annotations', data).
            success(function(data) {
                var len = annotation['last_bp']-annotation['first_bp']+1;
                var new_a = {
                    first_bp: annotation['first_bp'],
                    last_bp: annotation['last_bp'],
                    name: annotation['name'],
                    type: annotation['type'],
                    strand: annotation['strand'],
                    annotation_full_length: len,
                    annotation_first_bp: 1,
                    annotation_last_bp: len
                };
                new_a['display_name'] = edgeAnnotationDisplayName(new_a);
                new_a['display_css'] = edgeAnnotationSummaryCSS(new_a);
                $scope.annotations.push(new_a);
                $scope.annotate_error = undefined;
                $scope.zoomRefresh($scope.zoom['first_bp'], $scope.zoom['last_bp']);
            }).
            error(function(data, status, headers, config) {
                $scope.annotate_error = data;
            });
    }

    $http.get('/fragments/'+$scope.fragmentId).success(function(fragment) {
        $scope.fragment = fragment;
    });

    $http.get('/fragments/'+$scope.fragmentId+'/annotations').success(function(annotations) {
        $scope.annotations = annotations;
        annotations.forEach(function(annotation) {
            annotation['display_name'] = edgeAnnotationDisplayName(annotation);
            annotation['display_css'] = edgeAnnotationSummaryCSS(annotation);
            var gap = 0;
            // scale up gap based on number of annotations, so we don't have
            // too many annotations to display in summary mode
            gap += parseInt(annotations.length/200)*500;

            var last_annotation = undefined;
            if ($scope.summary_annotations.length > 0) {
                last_annotation = $scope.summary_annotations[$scope.summary_annotations.length-1];
            }
            if (last_annotation === undefined || last_annotation['last_bp']+gap < annotation['first_bp']) {
                $scope.summary_annotations.push(annotation);
            }
        });
        $scope.showSummary();
    });

    $scope.query = undefined;
    $scope.annotationOrderProp = 'first_bp';
}

function GenomeFragmentController($scope, $routeParams, $injector, $http) {
    $injector.invoke(FragmentController, this, { $scope: $scope });
    $scope.genomeId = $routeParams.genomeId;
    $scope.genome = undefined;
    $scope.changes = [];

    $http.get('/genomes/'+$scope.genomeId).success(function(genome) {
        $scope.genome = genome;
        genome['fragments'].forEach(function(fragment) {
            if (fragment['id'] == $scope.fragmentId) {
                if (fragment['changes']) {
                    fragment['changes'].forEach(function(changed_loc) {
                        edgeFetchChanges($http, fragment['id'], changed_loc, function(desc) {
                            $scope.changes.push(desc);
                        });
                    });
                }
            }
        });
    });
}

function GenomeOpController($scope, $http, $location) {
    $scope.op_error = undefined;

    $scope.doOp = function(action, op) {
        op['op'] = action;
        var data = JSON.stringify(op);
        $http.put('/genomes/'+$scope.genomeId+'/fragments/'+$scope.fragmentId, data).
            success(function(genome) {
                // got a new genome back
                $scope.op_error = undefined;
                var url = '/genomes/'+genome.id;
                $location.path(url);
            }).
            error(function(data, status, headers, config) {
                $scope.op_error = data;
            });
    }
}

function GenomeOpWithFragmentController($scope, $http) {
    $http.get('/fragments').success(function(data) {
        $scope.fragments = data;
    });

    $scope.orderProp = 'name';
}

