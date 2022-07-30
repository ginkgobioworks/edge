from edge.models import Operation
from edge.ssr.crelox import CreLoxReaction


class SSROp(object):
    @staticmethod
    def check(
        genome,
        donor,
        is_donor_circular,
        reaction,
        **kwargs
    ):
        # XXX

    @staticmethod
    def get_operation(
        donor,
        is_donor_circular,
        reaction
    ):
        params = dict(donor=donor, is_donor_circular=is_donor_circular, reaction=reaction)
        op = Operation(type=Operation.SSR[0], params=json.dumps(params))
        return op

    @staticmethod
    def perform(
        genome,
        donor,
        is_donor_circular,
        reaction,
        genome_name,
        notes,
        annotations=None,
    ):
        # XXX
        return new_genome
