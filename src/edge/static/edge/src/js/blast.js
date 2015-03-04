function GenomeBlastController($scope, $routeParams, $http) {
    $scope.genomeId = $routeParams.genomeId;
    $scope.genome = undefined;
    $scope.results = undefined;
    $scope.error = undefined;

    $scope.Blast = function(query, program) {
        $scope.waiting = true;
        var data = JSON.stringify({query: query, program: program});
        $http.post('/edge/genomes/'+$scope.genomeId+'/blast/', data)
             .success(function(results) {
                 $scope.waiting = false;
                 $scope.results = results;
             })
             .error(function(data, status, headers, config) {
                 $scope.waiting = false;
                 $scope.errors = data;
             });
    };

    $http.get('/edge/genomes/'+$scope.genomeId+'/').success(function(data) { $scope.genome = data; });
}
