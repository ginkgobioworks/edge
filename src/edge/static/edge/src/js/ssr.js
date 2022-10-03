function GenomeSSRController($scope, $routeParams, $http, $location) {
    $scope.genomeId = $routeParams.genomeId;
    $scope.reaction = undefined
    $scope.new_genome_name = undefined;
    $scope.donor = "";
    $scope.is_donor_circular = false;
    $scope.events = undefined;
    $scope.waiting = false;
    $scope.previewed_json = undefined;
    $scope.errors = undefined;
  
    $http.get("/edge/genomes/" + $scope.genomeId + "/").success(function(data) {
        $scope.genome = data;
    });

    $scope.PreviewSSR = function() {
        console.log($scope.new_genome_name);
        if (($scope.reaction === undefined) || ($scope.new_genome_name === undefined)) {
            $scope.previewed_json = undefined;
            $scope.errors = "Please specify (1) reaction and (2) new genome name";
        } else {
            $scope.waiting = true;
            $scope.donor = $scope.donor.replace(/\s+/g, "");
            var previewed_json = {
                donor: $scope.donor,
                is_donor_circular: $scope.is_donor_circular,
                reaction: $scope.reaction,
                genome_name: $scope.new_genome_name,
                create: false
            }
            var data = JSON.stringify(previewed_json);
            $http
                .post("/edge/genomes/" + $scope.genomeId + "/ssr/", data)
                .success(function(data) {
                    $scope.previewed_json = previewed_json;
                    $scope.events = data;
                    $scope.waiting = false;
                    $scope.errors = undefined;
                })
                .error(function(data, status, headers, config) {
                    $scope.previewed_json = undefined;
                    $scope.events = undefined;
                    $scope.waiting = false;
                    $scope.errors = data;
                });
        };
    };

    $scope.SSR = function() {
        $scope.waiting = true;
        $scope.previewed_json.create = true;
        var data = JSON.stringify($scope.previewed_json);
        $http
            .post("/edge/genomes/" + $scope.genomeId + "/ssr/", data)
            .success(function(genome) {
                $scope.waiting = false;
                var url = "/genomes/" + genome.id + "/";
                $location.path(url);
            })
            .error(function(data, status, headers, config) {
                $scope.waiting = false;
                $scope.errors = data;
            });
    };
}
