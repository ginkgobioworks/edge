function GenomeRecombinationController($scope, $routeParams, $http, $location) {
    $scope.genomeId = $routeParams.genomeId;
    $scope.cassette = "";
    $scope.homology_arm_length = 30;
    $scope.regions = undefined;
    $scope.errors = undefined;
    $scope.new_genome_name = undefined;

    $http.get('/edge/genomes/'+$scope.genomeId+'/').success(function(data) { $scope.genome = data; });

    $scope.FindRegions = function() {
        $scope.cassette = $scope.cassette.replace(/\s+/g,'');
        var data = JSON.stringify({cassette: $scope.cassette,
                                   homology_arm_length: $scope.homology_arm_length,
                                   create: false});
        $http.post('/edge/genomes/'+$scope.genomeId+'/recombination/', data)
             .success(function(data) {
                          console.log(data);
                          $scope.regions = data;
                      })
             .error(function(data, status, headers, config) { $scope.errors = data; });
    };

    $scope.Recombine = function() {
        $scope.cassette = $scope.cassette.replace(/\s+/g,'');
        var data = JSON.stringify({cassette: $scope.cassette,
                                   homology_arm_length: $scope.homology_arm_length,
                                   genome_name: $scope.new_genome_name,
                                   create: true});
        $http.post('/edge/genomes/'+$scope.genomeId+'/recombination/', data)
             .success(function(genome) {
                 console.log(genome);
                 var url = '/genomes/'+genome.id+'/';
                 $location.path(url);
             })
             .error(function(data, status, headers, config) { $scope.errors = data; });
    };
}
