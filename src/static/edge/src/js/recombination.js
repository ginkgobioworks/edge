function GenomeRecombinationController($scope, $routeParams, $http, $location) {
    $scope.genomeId = $routeParams.genomeId;
    $scope.cassette = "";
    $scope.homology_arm_length = 30;
    $scope.primer3_opts = {
      PRIMER_OPT_SIZE: 20,
      PRIMER_MIN_SIZE: 18,
      PRIMER_MAX_SIZE: 25,
      PRIMER_OPT_TM: 60.0,
      PRIMER_MIN_TM: 57.0,
      PRIMER_MAX_TM: 63.0,
      PRIMER_MIN_GC: 20.0,
      PRIMER_MAX_GC: 80.0,
      PRIMER_SALT_MONOVALENT: 50.0,
      PRIMER_DNA_CONC: 5000.0,
    };
    $scope.regions = undefined;
    $scope.errors = undefined;
    $scope.new_genome_name = undefined;

    $http.get('/edge/genomes/'+$scope.genomeId+'/').success(function(data) { $scope.genome = data; });

    $scope.FindRegions = function() {
        $scope.regions = undefined;
        $scope.waiting = true;
        $scope.cassette = $scope.cassette.replace(/\s+/g,'');
        var data = JSON.stringify({cassette: $scope.cassette,
                                   homology_arm_length: $scope.homology_arm_length,
                                   create: false,
                                   design_primers: true,
                                   primer3_opts: $scope.primer3_opts});
        $http.post('/edge/genomes/'+$scope.genomeId+'/recombination/', data)
             .success(function(data) {
                 $scope.waiting = false;
                 $scope.regions = data;
             })
             .error(function(data, status, headers, config) { $scope.errors = data; });
    };

    $scope.Recombine = function() {
        $scope.waiting = true;
        $scope.cassette = $scope.cassette.replace(/\s+/g,'');
        var data = JSON.stringify({cassette: $scope.cassette,
                                   homology_arm_length: $scope.homology_arm_length,
                                   genome_name: $scope.new_genome_name,
                                   create: true});
        $http.post('/edge/genomes/'+$scope.genomeId+'/recombination/', data)
             .success(function(genome) {
                 $scope.waiting = false;
                 var url = '/genomes/'+genome.id+'/';
                 $location.path(url);
             })
             .error(function(data, status, headers, config) { $scope.errors = data; });
    };
}
