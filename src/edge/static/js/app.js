angular.module('edge', []).
    config(['$routeProvider', function($routeProvider) {
        $routeProvider.
            when('/genomes',
                 {templateUrl: 'partials/genome-list.html',
                  controller: GenomeListController}).
            when('/genomes/:genomeId',
                 {templateUrl: 'partials/genome-detail.html',
                  controller: GenomeDetailController}).
            when('/genomes/:genomeId/fragments/:fragmentId',
                 {templateUrl: 'partials/genome-fragment.html',
                  controller: GenomeFragmentController}).
            when('/fragments',
                 {templateUrl: 'partials/fragment-list.html',
                  controller: FragmentListController}).
            when('/fragments/:fragmentId',
                 {templateUrl: 'partials/fragment-detail.html',
                  controller: FragmentController}).
            otherwise({redirectTo: '/genomes'});
    }])
