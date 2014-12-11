function GenomePcrController($scope, $routeParams, $http) {
    $scope.genomeId = $routeParams.genomeId;
    $scope.genome = undefined;
    $scope.results = undefined;
    $scope.error = undefined;

    $scope.Pcr = function(primer_a, primer_b) {
        var data = JSON.stringify({primers: [primer_a, primer_b]});
        $http.post('/edge/genomes/'+$scope.genomeId+'/pcr/', data)
             .success(function(results) {
                 $scope.results = results;
             })
             .error(function(data, status, headers, config) { $scope.errors = data; });
    };

    $http.get('/edge/genomes/'+$scope.genomeId+'/').success(function(data) { $scope.genome = data; });
}
