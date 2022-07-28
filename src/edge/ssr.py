from Bio.Seq import Seq


def rc(s):
    return str(Seq(s).reverse_complement())


def _c(s):
    return s.replace(" ", "").lower()


def find_indices(big, small):
    """
    Returns 0 based indices of all occurances of small in big.
    """

    res = []
    offset = 0
    while len(big) >= len(small):
        if small in big:
            res.append(offset + big.index(small))
            offset += big.index(small) + len(small)
            big = big[big.index(small) + len(small):]
        else:
            return res
    return res


class SiteLocation(object):

    def __init__(self, site, fragment, insert, start_0based):
        self.site = site
        self.fragment = fragment
        self.insert = insert
        self.start_0based = start_0based  # 0-indexed
        self.adjusted_start_0based = self.start_0based

    @property
    def fragment_id(self):
        return self.fragment.id if self.fragment is not None else None

    @property
    def on_insert(self):
        return self.insert is not None


class Recombination(object):
    """
    Generic class to model a recombination event.
    """

    def required_genome_sites(self):
        return []

    def required_insert_sites(self):
        return []

    def __init__(self):
        self.errors = []

    def find_matching_locations(self, all_locations, required_sites, on_insert):
        # XXX disambiguity
        # XXX does not allow multiple combination of sites on insert for integration
        # XXX handle circularity for sites on insert

        candidate_single_site_locations = all_locations
        candidate_single_site_locations = \
          [loc for loc in all_locations if on_insert == loc.on_insert and loc.site in required_sites]
        candidate_single_site_locations = sorted(
            candidate_single_site_locations,
            key=lambda loc: (loc.fragment_id, loc.start_0based)
        )

        matching_locations = []
        for i, single_site_loc in enumerate(candidate_single_site_locations):
            if single_site_loc.site == required_sites[0]:
                if [loc.site for loc in candidate_single_site_locations[i:i+len(required_sites)]] == list(required_sites):
                    matching_locations.append(candidate_single_site_locations[i:i+len(required_sites)])
        return matching_locations 

    def possible_locations(self, all_locations):
        required_insert_sites = self.required_insert_sites()
        insert_locations = []

        if len(required_insert_sites):
            locations_on_insert = []
            locations_on_insert.extend(self.find_matching_locations(all_locations, required_insert_sites, True))
            locations_on_insert.extend(self.find_matching_locations(all_locations, [rc(s) for s in required_insert_sites][::-1], True))

            if len(locations_on_insert) == 0:
                self.errors.append("Requires site(s) %s on insert, but did not find any" % (required_insert_sites,))
                return []
            elif len(locations_on_insert) > 1:
                self.errors.append("Requires one site or one set of sites %s on insert, found multiple" % (required_insert_sites,))
                return []
          
            insert_locations = locations_on_insert[0] 

        possible_locations = []

        required_genome_sites = self.required_genome_sites()
        locations_on_genome = []
        locations_on_genome.extend(self.find_matching_locations(all_locations, required_genome_sites, False))
        locations_on_genome.extend(self.find_matching_locations(all_locations, [rc(s) for s in required_genome_sites][::-1], False))
        for locs in locations_on_genome:
            possible_locations.append(locs+insert_locations)

        return possible_locations

    def generate_events(self, locations, event_cls):
        location_sets = self.possible_locations(locations)
        events = []
        for locs in location_sets:
            events.append(event_cls(self, locs))
        return events


class Integration(Recombination):

    def __init__(self, site_insert, site_genome, recombined_site_left_genome, recombined_site_right_genome):
        super(Integration, self).__init__()
        self.site_insert = site_insert
        self.site_genome = site_genome
        self.recombined_site_left_genome = recombined_site_left_genome
        self.recombined_site_right_genome = recombined_site_right_genome

    def required_genome_sites(self):
        return [self.site_genome]

    def required_insert_sites(self):
        return [self.site_insert]

    def events(self, locations):
        return self.generate_events(locations, IntegrationEvent)


class Excision(Recombination):
    def __init__(self, site_left, site_right, recombined_site):
        super(Excision, self).__init__()
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
        super(Inversion, self).__init__()
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
        super(RMCE, self).__init__()
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
        for loc in self.locations:
            if loc.fragment_id:
                return loc.fragment_id
        assert False

    def is_reversed(self):
        first_loc = self.genomic_locations[0]

        if first_loc.site == self.recombination.required_genome_sites()[0]:
            return False
        if rc(first_loc.site) == self.recombination.required_genome_sites()[-1]:
            return True
        assert False

    @property
    def genomic_locations(self):
        return sorted([l for l in self.locations if l.fragment_id is not None], key=lambda l: l.start_0based)

    @property
    def genomic_start_0based(self):
        return self.genomic_locations[0].start_0based

    @property
    def genomic_adjusted_start_0based(self):
        return self.genomic_locations[0].adjusted_start_0based

    def adjust_genomic_start_0based(self, shift):
        for loc in self.genomic_locations:
            loc.adjusted_start_0based += shift


class IntegrationEvent(Event):

    def get_integrated_aligned_with_site_direction(self, insert):
        site_insert = self.recombination.site_insert
        if site_insert not in insert:
           insert = rc(insert)
        assert insert.index(site_insert) >= 0
        insert = insert*2
        bps = find_indices(insert, site_insert)
        assert len(bps) == 2
        return insert[bps[0]+len(site_insert):bps[1]]

    def run(self, new_fragment, insert, is_insert_circular):
        integrated = self.get_integrated_aligned_with_site_direction(insert)
        bps_to_replace = len(self.site_genome)

        if self.is_reversed():
            integrated = rc(integrated)
            new_site_left = rc(self.recombination.recombined_site_right_genome)
            new_site_right = rc(self.recombination.recombined_site_left_genome)
        else:
            new_site_left = self.recombination.recombined_site_left_genome
            new_site_right = self.recombination.recombined_site_right_genome

        integrated = new_site_left+integrated+new_site_right
        new_fragment.replace_bases(
            genomic_locations[0].adjusted_start_0based+1,
            bps_to_replace,
            integrated
        )

        return len(integrated)-bps_to_replace


class ExcisionEvent(Event):

    def run(self, new_fragment, insert, is_insert_circular):
        genomic_locations = self.genomic_locations
        assert len(genomic_locations) == 2

        if not self.is_reversed():
            assert genomic_locations[0].site == self.recombination.site_left
            assert genomic_locations[1].site == self.recombination.site_right
            new_site = self.recombination.recombined_site
            bps_to_remove = genomic_locations[1].adjusted_start_0based -\
                            genomic_locations[0].adjusted_start_0based +\
                            len(self.recombination.site_right)
        else:
            assert genomic_locations[0].site == rc(self.recombination.site_right)
            assert genomic_locations[1].site == rc(self.recombination.site_left)
            new_site = rc(self.recombination.recombined_site)
            bps_to_remove = genomic_locations[1].adjusted_start_0based -\
                            genomic_locations[0].adjusted_start_0based +\
                            len(self.recombination.site_left)

        new_fragment.replace_bases(
            genomic_locations[0].adjusted_start_0based+1,
            bps_to_remove,
            new_site
        )

        return len(new_site)-bps_to_remove


class InversionEvent(Event):

    def run(self, new_fragment, insert, is_insert_circular):
        genomic_locations = self.genomic_locations
        assert len(genomic_locations) == 2

        if not self.is_reversed():
            assert genomic_locations[0].site == self.recombination.site_left
            assert genomic_locations[1].site == self.recombination.site_right
            old_site_left = self.recombination.site_left
            old_site_right = self.recombination.site_right
            new_site_left = self.recombination.recombined_site_left
            new_site_right = self.recombination.recombined_site_right

            bps_to_replace = genomic_locations[1].adjusted_start_0based -\
                             genomic_locations[0].adjusted_start_0based +\
                             len(self.recombination.site_right)
        else:
            assert genomic_locations[0].site == rc(self.recombination.site_right)
            assert genomic_locations[1].site == rc(self.recombination.site_left)
            old_site_left = rc(self.recombination.site_right)
            old_site_right = rc(self.recombination.site_left)
            new_site_left = rc(self.recombination.recombined_site_right)
            new_site_right = rc(self.recombination.recombined_site_left)

            bps_to_replace = genomic_locations[1].adjusted_start_0based -\
                             genomic_locations[0].adjusted_start_0based +\
                             len(self.recombination.site_left)

        sequence_to_replace = new_fragment.get_sequence(bp_lo=genomic_locations[0].adjusted_start_0based+len(old_site_left)+1,
                                                        bp_hi=genomic_locations[1].adjusted_start_0based+1-1)
        new_sequence = new_site_left+rc(sequence_to_replace)+new_site_right
        assert len(old_site_left)+len(sequence_to_replace)+len(old_site_right) == bps_to_replace

        new_fragment.replace_bases(
            genomic_locations[0].adjusted_start_0based+1,
            bps_to_replace,
            new_sequence
        )

        return len(new_sequence)-(len(old_site_left)+len(sequence_to_replace)+len(old_site_right))


class RMCEEvent(Event):

    def get_integrated_aligned_with_site_direction(self, insert):
        site_left_insert = self.recombination.site_left_insert
        site_right_insert = self.recombination.site_right_insert
        if site_left_insert not in insert:
           insert = rc(insert)

        assert insert.index(site_left_insert) == 1
        assert insert.index(site_right_insert) == 1
        return insert[insert.index(site_left_insert)+len(site_left_insert):insert.index(site_right_insert)]

    def run(self, new_fragment, insert, is_insert_circular):
        integrated = self.get_integrated_aligned_with_site_direction(insert)

        if not self.is_reversed():
            assert genomic_locations[0].site == self.recombination.site_left_genome
            assert genomic_locations[1].site == self.recombination.site_right_genome
            old_site_left = self.recombination.site_left_genome
            old_site_right = self.recombination.site_right_genome
            new_site_left = self.recombination.recombined_site_left_genome
            new_site_right = self.recombination.recombined_site_right_genome

            bps_to_replace = genomic_locations[1].adjusted_start_0based -\
                             genomic_locations[0].adjusted_start_0based +\
                             len(self.recombination.site_right_genome)
        else:
            integrated = rc(integrated)

            assert genomic_locations[0].site == rc(self.recombination.site_right_genome)
            assert genomic_locations[1].site == rc(self.recombination.site_left_genome)
            old_site_left = rc(self.recombination.site_right_genome)
            old_site_right = rc(self.recombination.site_left_genome)
            new_site_left = rc(self.recombination.recombined_site_right_genome)
            new_site_right = rc(self.recombination.recombined_site_left_genome)

            bps_to_replace = genomic_locations[1].adjusted_start_0based -\
                             genomic_locations[0].adjusted_start_0based +\
                             len(self.recombination.site_left_genome)

        integrated = new_site_left+integrated+new_site_right
        new_fragment.replace_bases(
            genomic_locations[0].adjusted_start_0based+1,
            bps_to_replace,
            integrated
        )

        return len(integrated)-bps_to_replace


def add_reverse_sites(sites):
    new_sites = []
    for s in sites:
        new_sites.append(s)
        new_sites.append(rc(s))
    return set(new_sites)


def find_site_locations_on_sequence(sequence, is_circular, sites, fragment_obj=None):
    sites = add_reverse_sites(sites)
    template = sequence
    if is_circular:
        template = sequence*2
    template = template.lower()

    locations = []
    for site in sites:
        locations.extend(
            find_query_locations_on_duplicated_template(
                template, len(sequence), site, fragment_obj=fragment_obj
            )
        )
    return locations


def find_query_locations_on_duplicated_template(template, sequence_len, site, fragment_obj=None):
    indices = find_indices(template, site)
    indices = [i for i in indices if i < sequence_len]
    return [
      SiteLocation(
          site,
          fragment_obj,
          template[0:sequence_len] if fragment_obj is None else None,
          i
      )
      for i in indices
    ]


class Reaction(object):

    @staticmethod
    def allowed():
        return []

    def __init__(self, parent_genome, insert, is_insert_circular, parent_fragments=None):
        self.parent_genome = parent_genome
        if parent_fragments is not None:
            self.parent_fragments = parent_fragments
        else:
            self.parent_fragments = [f.indexed_fragment() for f in self.parent_genome.fragments.all()]

        self.insert = insert
        self.is_insert_circular = is_insert_circular
        self.locations = None
        self.events = None
        self.errors = []

    @classmethod
    def __sites_on_genome(cls):
        sites = []
        for recombination in cls.allowed():
            sites.extend(recombination.required_genome_sites())
        return set(sites)

    @classmethod
    def __sites_on_insert(cls):
        sites = []
        for recombination in cls.allowed():
            sites.extend(recombination.required_insert_sites())
        return set(sites)

    def determine_site_locations(self):
	# FIXME this is a very slow implementation, we can make this faster by
	# using cached results

        locations = []

        if self.insert:
            locations.extend(
                find_site_locations_on_sequence(
                    self.insert,
                    self.is_insert_circular,
                    self.__sites_on_insert()
                )
            )

        sites_on_genome = self.__sites_on_genome()
        for fragment in self.parent_fragments:
            locations.extend(
                find_site_locations_on_sequence(
                    fragment.sequence,
                    fragment.circular,
                    sites_on_genome,
                    fragment_obj=fragment
                )
            )

        self.locations = locations

    def group_into_events(self):
	# TODO if we attempt an integration with loxP on loxP, how do we return
	# that as an error not just silently ignore?

        self.determine_site_locations()
        self.events = []

        for recombination in self.allowed():
            self.events.extend(recombination.events(self.locations))
            self.errors.extend(recombination.errors)

    def run_reaction(self, new_genome_name):
        self.group_into_events()
        if self.errors:
            return
        if len(self.events) == 0:
            return

        new_genome = self.parent_genome.update()
        new_genome.name = new_genome_name
        new_genome.save()

        old_to_new_fragment_dict = {}

	# needs to run through events by fragment and start bp, and after each
	# event shift the coordinates of the remaining events

        events = sorted(self.events, key=lambda e: (e.fragment_id, e.genomic_adjusted_start_0based))
        while len(events) > 0:
            event = events[0]

            if event.fragment_id not in old_to_new_fragment_dict:
                with new_genome.update_fragment_by_fragment_id(
                    event.fragment_id, new_fragment=True
                ) as f:
                    old_to_new_fragment_dict[event.fragment_id] = f

            new_fragment = old_to_new_fragment_dict[event.fragment_id]
            bp_shift = event.run(new_fragment, self.insert, self.is_insert_circular)

            for future_event in events[1:]:
                if future_event.fragment_id == event.fragment_id:
                    future_event.adjust_genomic_start_0based(bp_shift)
            events = events[1:]


        return new_genome


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
