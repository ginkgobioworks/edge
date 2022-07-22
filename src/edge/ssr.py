
def _c(s):
  return s.replace(" ", "").lower()


class Sites(object):
    """
    Cre/Lox sites
    """

    loxP     = _c("ATAACTTCGTATA GCATACAT TATACGAAGTTAT")
    lox66    = _c("ATAACTTCGTATA GCATACAT TATACGAACGGTA")
    lox71    = _c("TACCGTTCGTATA GCATACAT TATACGAAGTTAT")
    lox72    = _c("TACCGTTCGTATA GCATACAT TATACGAACGGTA")

    lox5171  = _c("ATAACTTCGTATA GtAcACAT TATACGAAGTTAT")
    lox2272  = _c("ATAACTTCGTATA GgATACtT TATACGAAGTTAT")

    loxm2    = _c("ATAACTTCGTATA TGGTTTCT TATACGAAGTTAT")
    loxm2_66 = _c("ATAACTTCGTATA TGGTTTCT TATACGAACGGTA")
    loxm2_71 = _c("TACCGTTCGTATA TGGTTTCT TATACGAAGTTAT")
    loxm2_72 = _c("TACCGTTCGTATA TGGTTTCT TATACGAACGGTA")

    """
    Flp/FRT sites
    """

    Frt  = _c("GAAGTTCCTATTC tctagaaa GTATAGGAACTTC")
    Frt1 = _c("GAAGTTCCTATTC tctagata GTATAGGAACTTC")
    Frt2 = _c("GAAGTTCCTATTC tctactta GTATAGGAACTTC")
    Frt3 = _c("GAAGTTCCTATTC ttcaaata GTATAGGAACTTC")
    Frt4 = _c("GAAGTTCCTATTC tctagaag GTATAGGAACTTC")
    Frt5 = _c("GAAGTTCCTATTC ttcaaaag GTATAGGAACTTC")

    """
    PhiC31 sites
    """

    PhiC31_attB_tt = _c("TGCGGGTGCCAGGGCGTGCCC tt GGGCTCCCCGGGCGCGTACTCC")
    PhiC31_attP_tt = _c("GTGCCCCAACTGGGGTAACCT tt GAGTTCTCTCAGTTGGGGG")

    """
    Bxb1 sites
    """

    Bxb1_attB =     _c("TCGGCCGGCTTGTCGACGACG gcggtctc CGTCGTCAGGATCATCCGGGC")
    Bxb1_attP = _c("GTCGTGGTTTGTCTGGTCAACCACC gcggtctc AGTGGTGTACGGTACAAACCCCGAC")


class SiteLocation(objects):

    def __init__(self, site, fragment, start, end, direction):
        self.site = site
        self.fragment = fragment
        self.start = start
        self.end = end
        self.direction = direction

        # after conversion from start/end to chunks
        self.prev_chunk = None
        self.first_chunk = None
        self.last_chunk = None
        self.next_chunk = None

        XXX define for both fragment and payload


class Recombination(object):
    """
    Generic class to model a recombination event.
    """

    def possible_locations(self):
        XXX
        XXX also disambiguity

    def generate_events(self, locations, event_cls):
        location_sets = self.possible_locations(locations)
        events = []
        for locs in location_sets:
            events.append(event_cls(self, locs))
        return events


class Integration(Recombination):

    def __init__(self, site_insert, site_genome, recombined_site_genome_left, recombined_site_genome_right):
        self.site_insert = site_insert
        self.site_genome = site_genome
        self.recombined_site_genome_left = recombined_site_genome_left
        self.recombined_site_genome_right = recombined_site_genome_right

    def required_genome_sites(self):
        return [self.site_genome]

    def required_insert_sites(self):
        return [self.site_insert]

    def events(self, locations):
        return self.generate_events(locations, IntegrationEvent)


class Excision(Recombination):
    def __init__(self, site_left, site_right, recombined_site):
        self.site_left = site_left
        self.site_right = site_right
        self.recombined_site = recombined_site

    def required_genome_sites(self):
        return [self.site_left, self.site_right]

    def required_insert_sites(self):
        return []

    def events(self, locations):
        return self.generate_events(locations, ExcisionEvent)


class Inversion(Recombination):
    def __init__(self, site_left, site_right, recombined_site_left, recombined_site_right):
        self.site_left = site_left
        self.site_right = site_right
        self.recombined_site_left = recombined_site_left
        self.recombined_site_right = recombined_site_right

    def required_genome_sites(self):
        return [self.site_left, self.site_right]

    def required_insert_sites(self):
        return []

    def events(self, locations):
        return self.generate_events(locations, InversionEvent)


class RMCE(Recombination):
    def __init__(self, site_left_insert, site_right_insert, site_left_genome, site_right_genome,
                 recombined_site_left_genome, recombined_site_right_genome):
        self.site_left_insert = site_left_insert
        self.site_right_insert = site_right_insert
        self.site_left_genome = site_left_genome
        self.site_right_genome = site_right_genome
        self.recombined_site_left_genome = recombined_site_left_genome
        self.recombined_site_right_genome = recombined_site_right_genome

    def required_genome_sites(self):
        return [self.site_left_genome, self.site_right_genome]

    def required_insert_sites(self):
        return [self.site_left_insert, self.site_right_insert]

    def events(self, locations):
        return self.generate_events(locations, RMCEEvent)


class Event(object):

    def __init__(self, recombination, locations):
        self.recombination = recombination
        self.locations = locations

    @property
    def fragment_id(self):
        XXX


class IntegrationEvent(Event):

    def run(self, new_fragment, payload, is_payload_circular):
        payload = self.rotate_payload()
        assert payload.index(recombination.site_insert) == 0

        XXX
        remove old site?
        create new chunk with payload
        flank with two sites


class ExcisionEvent(Event):

    def run(self, new_fragment, payload, is_payload_circular):
        # TODO
        pass


class InversionEvent(Event):

    def run(self, new_fragment, payload, is_payload_circular):
        # TODO
        pass


class RMCEEvent(Event):

    def run(self, new_fragment, payload, is_payload_circular):
        # TODO
        pass


class Reaction(object):

    @staticmethod
    def allowed():
        return []

    def __init__(self, parent_genome, payload, is_payload_circular):
        self.parent_genome = parent_genome
        self.payload = payload
        self.is_payload_circular = is_payload_circular
        self.locations = {}
        self.events = None
        self.errors = []

    def __ensure_locations(self):
        XXX

    def group_into_events(self):
	# TODO if we attempt an integration with loxP on loxP, how do we return
	# that as an error not just silently ignore?

        self.__ensure_locations(self)
        self.events = []

        for recombination in self.allowed():
            self.events.extend(recombination.events(self.locations))

    def run_reaction(self):
        self.__ensure_locations()
        if self.errors:
            return

        self.group_into_events()
        if self.errors:
            return
        if len(self.events) == 0:
            return

        new_genome = self.parent_genome.update()
        new_genome.name = genome_name
        new_genome.notes = notes
        new_genome.save()

        old_to_new_fragment_dict = {}
        for event in self.events:
            if event.fragment_id not in old_to_new_fragment_dict:
                with new_genome.update_fragment_by_fragment_id(
                    event.fragment_id, new_fragment=True
                ) as f:
                    old_to_new_fragment_dict[event.fragment_id] = f
            new_fragment = old_to_new_fragment_dict[event.fragment_id]
            event.run(new_fragment)
        
        return new_genome


class CreLoxReaction(Reaction):

    @staticmethod
    def allowed():
        return [
            Integration(Sites.lox66, Sites.lox71, Sites.lox72, Sites.loxP),
            Integration(Sites.lox71, Sites.lox66, Sites.loxP, Sites.lox72),

            Excision(Sites.loxP, Sites.loxP, Sites.loxP),
            Excision(Sites.lox66, Sites.loxP, Sites.loxP),
            Excision(Sites.loxP, Sites.lox66, Sites.lox66),
            Excision(Sites.lox71, Sites.loxP, Sites.lox71),
            Excision(Sites.loxP, Sites.lox71, Sites.loxP),

            Inversion(Sites.loxP, rc(Sites.loxP), Sites.loxP, rc(Sites.loxP)),
            Inversion(Sites.lox66, rc(Sites.loxP), Sites.loxP, rc(Sites.lox66)),
            Inversion(Sites.loxP, rc(Sites.lox66), Sites.lox66, rc(Sites.loxP)),
            Inversion(Sites.lox71, rc(Sites.loxP), Sites.lox71, rc(Sites.loxP)),
            Inversion(Sites.loxP, rc(Sites.lox71), Sites.loxP, rc(Sites.lox71)),

            RMCE(Sites.loxP, Sites.lox2272, Sites.loxP, Sites.lox2272, Sites.loxP, Sites.lox2272),
            RMCE(Sites.lox66, Sites.lox2272, Sites.lox71, Sites.lox2272, Sites.lox72, Sites.lox2272),
        ]
