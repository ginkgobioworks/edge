var app = angular.module('edge', ['ngRoute', 'ngSanitize']).
    config(['$routeProvider', function($routeProvider) {
        $routeProvider.
            when('/import',
                 {template: JST['import'],
                  controller: ImportController}).
            when('/genomes',
                 {template: JST['genome-list'],
                  controller: GenomeListController}).
            when('/genomes/:genomeId',
                 {template: JST['genome-detail'],
                  controller: GenomeDetailController}).
            when('/genomes/:genomeId/recombination',
                 {template: JST['genome-recombination'],
                  controller: GenomeRecombinationController}).
            when('/genomes/:genomeId/blast',
                 {template: JST['genome-blast'],
                  controller: GenomeBlastController}).
            when('/genomes/:genomeId/pcr',
                 {template: JST['genome-pcr'],
                  controller: GenomePcrController}).
            when('/genomes/:genomeId/fragments/:fragmentId',
                 {template: JST['genome-fragment'],
                  controller: GenomeFragmentController}).
            when('/fragments',
                 {template: JST['fragment-list'],
                  controller: FragmentListController}).
            when('/fragments/:fragmentId',
                 {template: JST['fragment-detail'],
                  controller: FragmentController}).
            otherwise({redirectTo: '/genomes'});
    }]);

app.directive('partial', function($compile) {
    var linker = function(scope, element, attrs) {
        element.html(JST[attrs.template]());
        $compile(element.contents())(scope);
    };
    return {
        link: linker,
        restrict: 'E'
    }
});

app.directive('fileModel', ['$parse', function ($parse) {
    return {
        restrict: 'A',
        link: function(scope, element, attrs) {
            var model = $parse(attrs.fileModel);
            var modelSetter = model.assign;
            element.bind('change', function() {
                scope.$apply(function() {
                    modelSetter(scope, element[0].files[0]);
                });
            });
        }
    };
}]);

