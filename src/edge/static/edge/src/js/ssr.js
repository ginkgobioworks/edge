function GenomeSSRController($scope, $routeParams, $http, $location) {
    $scope.genomeId = $routeParams.genomeId;
    $scope.reaction = undefined
    $scope.new_genome_name = "";
    $scope.donor = "";
    $scope.is_donor_circular = false;
    $scope.checked = undefined;
    $scope.errors = undefined;
  
    $http.get("/edge/genomes/" + $scope.genomeId + "/").success(function(data) {
        $scope.genome = data;
    });

    $scope.CheckSSR = function() {
        if (($scope.donor === "") || ($scope.reaction === undefined) || ($scope.new_genome_name === "")) {
            $scope.checked = undefined;
            $scope.errors = "Please specify (1) SSR operation type, (2) new genome name, and (3) donor sequence";
        } else {
            $scope.waiting = true;
            $scope.donor = $scope.donor.replace(/\s+/g, "");
            var data = JSON.stringify({
                donor: $scope.donor,
                is_donor_circular: $scope.is_donor_circular,
                reaction: $scope.reaction,
                genome_name: $scope.genome_name,
                create: false
            });
            $http
                .post("/edge/genomes/" + $scope.genomeId + "/ssr/", data)
                .success(function(data) {
                    $scope.waiting = false;
                    $scope.checked = true;
                    $scope.errors = undefined;
                })
                .error(function(data, status, headers, config) {
                    $scope.waiting = false;
                    $scope.checked = false;
                    $scope.errors = data;
                });
        };
    };

    $scope.SSR = function() {
        $scope.waiting = true;
        $scope.donor = $scope.donor.replace(/\s+/g, "");
        var data = JSON.stringify({
            donor: $scope.donor,
            is_donor_circular: $scope.is_donor_circular,
            reaction: $scope.reaction,
            genome_name: $scope.genome_name,
            create: true
        });
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
