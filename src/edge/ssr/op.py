import json
from edge.models import Operation
from edge.ssr.crelox import CreLoxReaction
from edge.ssr.flp import FlpReaction


class SSROp(object):
    JSON_ARGUMENT_REACTION_CRELOX = "crelox"
    JSON_ARGUMENT_REACTION_FLP = "flp"

    SUPPORTED = [
        JSON_ARGUMENT_REACTION_CRELOX,
        JSON_ARGUMENT_REACTION_FLP
    ]

    @staticmethod
    def get_reaction(
        genome,
        donor,
        is_donor_circular,
        reaction_name
    ):
        if reaction_name.lower() not in SSROp.SUPPORTED:
            raise Exception("Unknown SSR reaction, expecting one of %s" % (SSROp.SUPPORTED,))

        if reaction_name.lower() == SSROp.JSON_ARGUMENT_REACTION_CRELOX:
            reaction = CreLoxReaction(genome, donor, is_donor_circular)
        elif reaction_name.lower() == SSROp.JSON_ARGUMENT_REACTION_FLP:
            reaction = FlpReaction(genome, donor, is_donor_circular)
        else:
            assert False

        return reaction

    @staticmethod
    def check(
        genome,
        donor,
        is_donor_circular,
        reaction_name,
        **kwargs
    ):
        reaction = SSROp.get_reaction(genome, donor, is_donor_circular, reaction_name)
        reaction.group_into_events()
        return reaction.events

    @staticmethod
    def get_operation(
        donor,
        is_donor_circular,
        reaction_name,
        **kwargs
    ):
        params = dict(donor=donor, is_donor_circular=is_donor_circular, reaction_name=reaction_name)
        op = Operation(type=Operation.SSR[0], params=json.dumps(params))
        return op

    @staticmethod
    def perform(
        genome,
        donor,
        is_donor_circular,
        reaction_name,
        genome_name,
        notes,
        annotations=None,
    ):
        reaction = SSROp.get_reaction(genome, donor, is_donor_circular, reaction_name)
        new_genome = reaction.run_reaction(genome_name, notes=notes)
        return new_genome
