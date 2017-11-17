'use strict';

function edgeAnnotationDisplayName(annotation) {
    var s = annotation['name']
    if (s === undefined || s == '') { s = annotation['type']; }
    if (annotation['feature_base_first'] != 1 ||
        annotation['feature_base_last'] != annotation['feature_full_length']) {
        s += '[';
        if (annotation['feature_base_first'] != 1) {
            s += annotation['feature_base_first'];
        }
        s += ':';
        if (annotation['feature_base_last'] != annotation['feature_full_length']) {
            s += annotation['feature_base_last']+'/'+annotation['feature_full_length'];
        }
        s += ']';
    }
    return s;
}

function edgeAnnotationColor(annotation) {
    if (annotation['type'] == 'gene') { return 'danger'; }
    else if (annotation['type'] == 'pseudogene') { return 'warning'; }
    else if (annotation['type'] == 'promoter') { return 'primary'; }
    else if (annotation['type'] == 'RBS') { return 'success'; }
    return 'info';
}

function edgeAnnotationSummaryCSS(annotation) {
    var css = ['annotation-summary'];
    css.push('btn-'+edgeAnnotationColor(annotation));
    css.push('btn btn-xs');
    if (annotation['strand'] > 0) { css.push('annotation-strand-fwd'); }
    else { css.push('annotation-strand-rev'); }
    return css.join(' ');
}

function edgeAnnotationZoomCSS(annotation) {
    var css = ['annotation-zoom'];
    css.push('btn-'+edgeAnnotationColor(annotation));
    css.push('btn btn-xs');
    if (annotation['strand'] > 0) { css.push('annotation-strand-fwd'); }
    else { css.push('annotation-strand-rev'); }
    return css.join(' ');
}

function PaginateController($scope, $http, $timeout) {
    $scope.query = '';
    $scope.curPos = 0;
    $scope.pageSize = 20;
    $scope.hasPrev = false;
    $scope.hasMore = false;
    $scope.delayPromise = undefined;

    $scope.fetch = function () {
        var base = $scope.getBaseURL();
        var ps = $scope.pageSize+1;
        var url = base+'?q='+encodeURIComponent($scope.query)+'&s='+$scope.curPos+'&p='+ps;
        $http.get(url).success(function(data) {
            $scope.processData(data.slice(0, ps-1));
            if (data.length < ps) { $scope.hasMore = false; }
            else { $scope.hasMore = true; }
        });
    }

    $scope.nextPage = function () {
        $scope.curPos += $scope.pageSize;
        $scope.hasPrev = true;
        $scope.fetch();
    };

    $scope.prevPage = function () {
        $scope.curPos -= $scope.pageSize;
        if ($scope.curPos <= 0) { $scope.curPos = 0; $scope.hasPrev = false; }
        $scope.fetch();
    };

    $scope.delayFetch = function() {
        if ($scope.delayPromise !== undefined) {
            $timeout.cancel($scope.delayPromise);
            $scope.delayPromise = undefined;
        }
        $scope.delayPromise = $timeout(function() {
            $scope.fetch();
            $scope.delayPromise = undefined;
        }, 300);
    };

}

function ImportController($scope, $http, $injector) {
    $scope.getBaseURL = function() { return '/edge/'; }
    $scope.addGenome = function(){
        var file = $scope.gffFile;
        $scope.addGenomeStatus = {
            "error": null,
            "pending": null,
            "results": null,
        }
        if (!file) {
            $scope.addGenomeStatus.error = "Need to choose GFF File";
            return;
        }
        $scope.addGenomeStatus.pending = true;
        var fd = new FormData();
        var gffName = $scope.genome.name;
        fd.append(gffName.toString(), file);
        $http.post('/edge/import_genome/', fd, {
            transformRequest: angular.identity,
            headers: { 'Content-Type': undefined },
        }).
            success(function(data) {
                $scope.addGenomeStatus.pending = null;
                $scope.addGenomeStatus.results = data;
                console.log('data is', data);
            }).
            error(function(data, status, headers, config) {
                $scope.addGenomeStatus.error = data;
                $scope.addGenomeStatus.pending = null;
                $scope.addGenomeStatus.results = null;
            });
    };
}

function GenomeListController($scope, $injector) {
    $injector.invoke(PaginateController, this, { $scope: $scope });
    $scope.getBaseURL = function() { return '/edge/genomes/'; }
    $scope.processData = function(data) { $scope.genomes = data; }
    $scope.fetch();
}

function FragmentListController($scope, $http, $injector) {
    $injector.invoke(PaginateController, this, { $scope: $scope });
    $scope.getBaseURL = function() { return '/edge/fragments/'; }
    $scope.processData = function(data) { $scope.fragments = data; }
    $scope.fetch();

    $scope.add_fragment_error = undefined;
    $scope.addFragment = function(fragment) {
        var data = JSON.stringify(fragment);
        $http.post('/edge/fragments/', data).
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
        var url = '/edge/genomes/'+parent_id+'/';
        $http.get(url).success(function(data) {
            $scope.parent = '<a href="#/genomes/'+data['id']+'">Genome '+data['id']+': '+data['name']+'</a>';
        });
    }

    /* driver */
    $http.get('/edge/genomes/'+$scope.genomeId+'/').success(function(data) {
        $scope.genome = data;
        if (data['parent_id']) { fetchParent(data['parent_id']); }
    });
}

function FragmentControllerBase($scope, $routeParams, $http, $location) {
    var DEFAULT_ZOOM = 2000;

    $scope.featureTypes = edgeFeatureTypes;
    $scope.fragmentId = $routeParams.fragmentId;
    $scope.fragment = undefined;
    $scope.fetchedAnnotations = false;
    $scope.annotations = []
    $scope.display_summary = true;
    $scope.summary_annotations = [];
    $scope.zoom = {};

    function layoutAnnotations(annotations, base_first, base_last) {
        var full_length = base_last-base_first+1;
        var display = [];
        var layers = [];

        annotations.forEach(function(a) {
            var d = {};
            d['annotation'] = a;
            d['title'] = a.display_name+' ('+a.type+')';
            d['formatted_qualifiers'] = a.formatted_qualifiers;
            d['css'] = edgeAnnotationZoomCSS(a);
            var sliced = false;
            var slice = [1, a['base_last']-a['base_first']+1];

            var a_first = a['base_first'];
            if (a_first < base_first) {
                a_first = base_first;
                sliced = true;
                slice[0] = base_first-a['base_first']+1;
            }
            var a_last = a['base_last'];
            if (a_last > base_last) {
                a_last = base_last; 
                sliced = true;
                slice[1] = base_last-a['base_first']+1;
            }
            if (sliced) {
                d['title'] += ' ('+slice[0]+'-'+slice[1]+')';
            }

            var left = (a_first-base_first)*100.0/full_length;
            d['left'] = ''+left+'%';
            var width = (a_last-a_first+1)*100.0/full_length;
            d['width'] = ''+width+'%';

            /* find a layer to put the feature */
            var found_layer_i = undefined;
            for (var i=0; i<layers.length; i++) {
                var layer = layers[i];
                var overlap = false;
                /* anything in this layer overlap? */
                layer.forEach(function(layer_a) {
                    if (layer_a['base_last']+100 >= a['base_first']) { overlap = true; } 
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
            var f = $scope.zoom['base_first'];
            var l = $scope.zoom['base_last'];
            $http.get('/edge/fragments/'+$scope.fragmentId+'/sequence/?f='+f+'&l='+l).success(function(data) {
                $scope.zoom['sequence'] = data;
                $scope.zoom['has_sequence'] = true;
                if ($scope.zoom['sequence_viewer'] === undefined) {
                    jQuery('#sequence-viewer').empty();
                    $scope.zoom['sequence_viewer'] = SvJs(jQuery, '#sequence-viewer');
                }
                // $scope.zoom['sequence_viewer'].setSequence(data.sequence, f);
                var annotations = $scope.zoom['annotations'];
                $scope.zoom['sequence_viewer'].setSequenceWithAnnotations(data.sequence, annotations, f);
            });
        }
    }

    $scope.fetchAllSequence = function() {
        $http.get('/edge/fragments/'+$scope.fragmentId+'/sequence/').success(function(data) {
            $scope.sequence = data;
        });
    };

    $scope.zoomRefresh = function(base_first, base_last) {
        var zoom_annotations = [];
        var min_bp = undefined;
        var max_bp = undefined;
        $scope.annotations.forEach(function(a) {
            if (a['base_last'] > base_first && a['base_first'] < base_last) {
                zoom_annotations.push(a);
                if (min_bp === undefined || a['base_first'] < min_bp) {
                    min_bp = a['base_first'];
                }
                if (max_bp === undefined || a['base_last'] > max_bp) {
                    max_bp = a['base_last'];
                }
            }
        });
        if (zoom_annotations.length > 0 && false) {
            zoom_annotations.sort(function(a,b) { return a['base_first']-b['base_first']; });
            $scope.zoom['base_first'] = Math.min(min_bp, base_first);
            $scope.zoom['base_last'] = Math.max(max_bp, base_last);
        }
        else {
            $scope.zoom['base_first'] = base_first;
            $scope.zoom['base_last'] = base_last;
        }
        $scope.zoom['annotations'] = zoom_annotations;
        var display = layoutAnnotations(zoom_annotations, $scope.zoom['base_first'], $scope.zoom['base_last']);
        $scope.zoom['display'] = display;
        $scope.zoom['has_sequence'] = false;
        if ($scope.zoom['sequence_viewer'] !== undefined) {
            $scope.zoom['sequence_viewer'].clear(); 
        }
    }

    $scope.zoomAt = function(annotation) {
        // get DEFAULT_ZOOM bps upstream and downstream from annotation
        var base_first = annotation['base_first']-DEFAULT_ZOOM;
        if (base_first < 1) { base_first = 1; }
        var base_last = annotation['base_last']+DEFAULT_ZOOM;
        if (base_last > $scope.fragment['length']) { base_last = $scope.fragment['length']; }

        $scope.goto(base_first, base_last, true);
    }

    $scope.zoomAtExact = function(annotation) {
        var base_first = annotation['base_first'];
        var base_last = annotation['base_last'];
        if (base_last > $scope.fragment['length']) { base_last = $scope.fragment['length']; }

        $scope.goto(base_first, base_last, true);
    }

    $scope.zoomMoveRight = function() {
        var cur_zoom = $scope.zoom['base_last']-$scope.zoom['base_first']+1;
        var base_first = $scope.zoom['base_last'];
        var base_last = base_first+cur_zoom;
        if (base_last > $scope.fragment['length']) { base_last = $scope.fragment['length']; }

        if (base_first < base_last) {
            $scope.goto(base_first, base_last, false);
        }
    }

    $scope.zoomMoveLeft = function() {
        var cur_zoom = $scope.zoom['base_last']-$scope.zoom['base_first']+1;
        if ($scope.zoom['base_first'] > 1) {
            var base_last = $scope.zoom['base_first'];
            var base_first = base_last-cur_zoom;
            if (base_first < 1) { base_first = 1; }
            if (base_first < base_last)
                $scope.goto(base_first, base_last, false);
        }
    }

    $scope.zoomOut = function() {
        var base_first = $scope.zoom['base_first']-DEFAULT_ZOOM;
        var base_last = $scope.zoom['base_last']+DEFAULT_ZOOM;
        if (base_first < 1) { base_first = 1; }
        if (base_last > $scope.fragment['length']) { base_last = $scope.fragment['length']; }

        $scope.goto(base_first, base_last, false);
    }

    $scope.zoomIn = function() {
        if ($scope.zoom['base_last']-$scope.zoom['base_first'] > 2*DEFAULT_ZOOM) {
            var base_first = $scope.zoom['base_first']+DEFAULT_ZOOM;
            var base_last = $scope.zoom['base_last']-DEFAULT_ZOOM;

            $scope.goto(base_first, base_last, false);
        }
    }

    $scope.goto = function(base_first, base_last, reset) {
        if (reset) {
            $scope.zoom = {};
            $scope.display_summary = false;
        }
        $scope.zoomRefresh(base_first, base_last);
    }

    $scope.updateBaseInput = function() {
        $scope.zoom['base_first'] = parseInt($scope.zoom['base_first']);
        $scope.zoom['base_last'] = parseInt($scope.zoom['base_last']);
    }

    $scope.userGoto = function() {
        $scope.goto($scope.zoom['base_first'], $scope.zoom['base_last'], true);
    }

    $scope.showSummary = function() {
        $scope.zoom = {
            'base_first': 1,
            'base_last': $scope.fragment['length'],
            'annotations': $scope.annotations
        };
        $scope.display_summary = true;
    }

    $scope.annotate_error = undefined;
    $scope.addAnnotation = function(annotation) {
        var data = JSON.stringify(annotation);
        $http.post('/edge/fragments/'+$scope.fragmentId+'/annotations/', data).
            success(function(data) {
                var len = annotation['base_last']-annotation['base_first']+1;
                var new_a = {
                    base_first: annotation['base_first'],
                    base_last: annotation['base_last'],
                    name: annotation['name'],
                    type: annotation['type'],
                    strand: annotation['strand'],
                    feature_full_length: len,
                    feature_base_first: 1,
                    feature_base_last: len
                };
                new_a['display_name'] = edgeAnnotationDisplayName(new_a);
                new_a['display_css'] = edgeAnnotationSummaryCSS(new_a);
                new_a['qualifiers'] = {};
                new_a['formatted_qualifiers'] = "";
                $scope.annotations.push(new_a);
                $scope.annotate_error = undefined;

                $scope.goto($scope.zoom['base_first'], $scope.zoom['base_last'], false);
            }).
            error(function(data, status, headers, config) {
                $scope.annotate_error = data;
            });
    }

    function tabify() {
        jQuery('#fragment-tab a').click(function (e) {
            e.preventDefault();
            jQuery(this).tab('show');
        });
        jQuery('#fragment-tab a:first').tab('show');
    }

    $http.get('/edge/fragments/'+$scope.fragmentId+'/').success(function(fragment) {
        $scope.fragment = fragment;
        // do this here, after templates have been rendered
        tabify();

        $http.get('/edge/fragments/'+$scope.fragmentId+'/annotations/').success(function(annotations) {
            $scope.annotations = annotations;
            $scope.fetchedAnnotations = true;
            annotations.forEach(function(annotation) {
                annotation['display_name'] = edgeAnnotationDisplayName(annotation);
                annotation['display_css'] = edgeAnnotationSummaryCSS(annotation);
                annotation['formatted_qualifiers'] = [];
                for (var field in annotation.qualifiers) {
                    annotation['formatted_qualifiers'].push(''+field+': '+annotation.qualifiers[field]);
                }
                annotation['formatted_qualifiers'].sort();
                annotation['formatted_qualifiers'] = annotation['formatted_qualifiers'].join('\n');

                var gap = 0;
                // scale up gap based on number of annotations, so we don't have
                // too many annotations to display in summary mode
                gap += parseInt(annotations.length/200)*500;

                var last_annotation = undefined;
                if ($scope.summary_annotations.length > 0) {
                    last_annotation = $scope.summary_annotations[$scope.summary_annotations.length-1];
                }
                if (last_annotation === undefined || last_annotation['base_last']+gap < annotation['base_first']) {
                    $scope.summary_annotations.push(annotation);
                }
            });

            $scope.postIndexCallbacks.forEach(function(cb) { cb(); })


            $scope.showSummary();

            var bp = $location.search().bp;
            if (bp) {
                bp = bp.split(',').slice(0,2);
                bp[0] = parseInt(bp[0]);
                bp[1] = parseInt(bp[1]);
                if (bp[0]) {
                    if (!bp[1]) {
                        bp[1] = bp[0]+DEFAULT_ZOOM;
                        bp[0] = bp[0]-DEFAULT_ZOOM;
                     }
                    $scope.goto(bp[0], bp[1], true);
                }
            }

        });
    });

    $scope.postIndexCallbacks = [];

    // fetch fragment again after indexing to get up to date length
    function refreshFragment() {
        $http.get('/edge/fragments/'+$scope.fragmentId+'/').success(function(fragment) {
            $scope.fragment = fragment;
        });
    }
    $scope.postIndexCallbacks.push(refreshFragment);

    $scope.query = undefined;
    $scope.annotationOrderProp = 'base_first';
}

function FragmentController($scope, $routeParams, $injector, $http) {
    $injector.invoke(FragmentControllerBase, this, { $scope: $scope });
    $scope.genomes = undefined;

    function findGenome() {
        $http.get('/edge/genomes/?f='+$scope.fragmentId).success(function(genomes) {
            if (genomes.length) { $scope.genomes = genomes; }
        });
    }

    findGenome();
    $scope.postIndexCallbacks.push(findGenome);
}

function GenomeFragmentController($scope, $routeParams, $injector, $http) {
    $injector.invoke(FragmentControllerBase, this, { $scope: $scope });
    $scope.genomeId = $routeParams.genomeId;
    $scope.genome = undefined;
    $scope.changes_and_locs = [];

    function fetchGenome() {
        $http.get('/edge/genomes/'+$scope.genomeId+'/').success(function(genome) {
            $scope.genome = genome;
        });
    }

    fetchGenome();
    $scope.postIndexCallbacks.push(fetchGenome);
}
