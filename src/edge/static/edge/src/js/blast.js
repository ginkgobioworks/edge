function GenomeBlastController($scope, $routeParams, $http) {
    $scope.genomeId = $routeParams.genomeId;
    $scope.genome = undefined;
    $scope.results = undefined;
    $scope.error = undefined;

    $scope.Blast = function(query, program) {
        var data = JSON.stringify({query: query, program: program});
        $http.post('/edge/genomes/'+$scope.genomeId+'/blast/', data)
             .success(function(results) {
                 $scope.results = results;
             })
             .error(function(data, status, headers, config) { $scope.errors = data; });
    };

    $http.get('/edge/genomes/'+$scope.genomeId+'/').success(function(data) { $scope.genome = data; });
}
