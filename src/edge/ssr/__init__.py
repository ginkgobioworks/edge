from Bio.Seq import Seq


# sites separated by more than these number of bps will not be considered as
# part of the same combination
SITE_SEPARATION_MAX = 50000


def rc(s):
    return str(Seq(s).reverse_complement())


def lower_no_whitespace(s):
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

    def __init__(self, site, fragment, insert, is_fragment_or_insert_circular, start_0based):
        self.site = site
        self.fragment = fragment
        self.insert = insert
        self.start_0based = start_0based  # 0-indexed
        self.used = False
        self.is_fragment_or_insert_circular = is_fragment_or_insert_circular

    def use(self):
        self.used = True

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

    def find_matching_locations(self, all_locations, required_sites, on_insert, errors):
        all_candidate_single_site_locations = \
          [loc for loc in all_locations
           if on_insert == loc.on_insert and loc.site in required_sites and loc.used is False]

        candidate_single_site_locations_by_fragment = {}
        for loc in all_candidate_single_site_locations:
            if loc.fragment_id not in candidate_single_site_locations_by_fragment:
                candidate_single_site_locations_by_fragment[loc.fragment_id] = []
            candidate_single_site_locations_by_fragment[loc.fragment_id].append(loc)

        matching_locations = []

        for fragment_id, candidate_single_site_locations in \
                candidate_single_site_locations_by_fragment.items():

            candidate_single_site_locations = sorted(
                candidate_single_site_locations,
                key=lambda loc: loc.start_0based
            )

            original_candidate_single_site_locations_length = len(candidate_single_site_locations)
            if len(candidate_single_site_locations) > 0 and \
               candidate_single_site_locations[0].is_fragment_or_insert_circular is True:
                candidate_single_site_locations = \
                    candidate_single_site_locations + candidate_single_site_locations

            for i, single_site_loc in enumerate(candidate_single_site_locations):
                if i >= original_candidate_single_site_locations_length:
                    break
                if single_site_loc.site == required_sites[0]:
                    next_locs = candidate_single_site_locations[i:i + len(required_sites)]
                    next_sites = [loc.site for loc in next_locs]
                    site_separation = next_locs[-1].start_0based - next_locs[0].start_0based
                    if next_sites == list(required_sites) and site_separation < SITE_SEPARATION_MAX:

                        has_conflict = False
                        for loc in next_locs:
                            if loc.used:
                                errors.append("Found %s, can trigger multiple events" % loc.site)
                                has_conflict = True
                            loc.use()

                        if not has_conflict:
                            matching_locations.append(next_locs)

        return matching_locations

    def possible_locations(self, all_locations, errors):
        required_insert_sites = self.required_insert_sites()
        insert_locations = []

        if len(required_insert_sites):
            locations_on_insert = []
            reverse_of_required_insert_sites = [rc(s) for s in required_insert_sites][::-1]
            locations_on_insert.extend(
                self.find_matching_locations(all_locations, required_insert_sites, True, errors)
            )
            if required_insert_sites != reverse_of_required_insert_sites:
                locations_on_insert.extend(
                    self.find_matching_locations(
                        all_locations, reverse_of_required_insert_sites, True, errors
                    )
                )

            if len(locations_on_insert) == 0:
                errors.append(
                    "Required site(s) %s on insert missing" % (required_insert_sites,)
                )
                return []
            elif len(locations_on_insert) > 1:
                errors.append(
                    "Requires one site or one set of sites %s on insert, found multiple"
                    % (required_insert_sites,)
                )
                return []

            insert_locations = locations_on_insert[0]

        possible_locations = []

        required_genome_sites = self.required_genome_sites()
        reverse_of_required_genome_sites = [rc(s) for s in required_genome_sites][::-1]
        locations_on_genome = []
        locations_on_genome.extend(
            self.find_matching_locations(all_locations, required_genome_sites, False, errors)
        )
        if required_genome_sites != reverse_of_required_genome_sites:
            locations_on_genome.extend(
                self.find_matching_locations(
                    all_locations, reverse_of_required_genome_sites, False, errors
                )
            )
        for locs in locations_on_genome:
            possible_locations.append(locs + insert_locations)

        return possible_locations

    def generate_events(self, locations, event_cls, errors):
        location_sets = self.possible_locations(locations, errors)
        events = []
        for locs in location_sets:
            events.append(event_cls(self, locs))
        return events


class Integration(Recombination):

    def __init__(self, site_insert, site_genome,
                 recombined_site_left_genome, recombined_site_right_genome):
        super(Integration, self).__init__()
        self.site_insert = site_insert
        self.site_genome = site_genome
        self.recombined_site_left_genome = recombined_site_left_genome
        self.recombined_site_right_genome = recombined_site_right_genome

    def required_genome_sites(self):
        return [self.site_genome]

    def required_insert_sites(self):
        return [self.site_insert]

    def events(self, locations, errors):
        return self.generate_events(locations, IntegrationEvent, errors)


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

    def events(self, locations, errors):
        return self.generate_events(locations, ExcisionEvent, errors)


class Inversion(Recombination):
    def __init__(self, site_left, site_right, recombined_site_left, recombined_site_right):
        super(Inversion, self).__init__()
        assert site_left not in site_right
        assert site_right not in site_left
        self.site_left = site_left
        self.site_right = site_right
        self.recombined_site_left = recombined_site_left
        self.recombined_site_right = recombined_site_right

    def required_genome_sites(self):
        return [self.site_left, self.site_right]

    def required_insert_sites(self):
        return []

    def events(self, locations, errors):
        return self.generate_events(locations, InversionEvent, errors)


class RMCE(Recombination):
    def __init__(self, site_left_insert, site_right_insert, site_left_genome, site_right_genome,
                 recombined_site_left_genome, recombined_site_right_genome):
        super(RMCE, self).__init__()
        assert site_left_insert not in site_right_insert
        assert site_right_insert not in site_left_insert
        assert site_left_genome not in site_right_genome
        assert site_right_genome not in site_left_genome

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

    def events(self, locations, errors):
        return self.generate_events(locations, RMCEEvent, errors)


class Event(object):

    def __init__(self, recombination, locations):
        self.recombination = recombination
        self.locations = locations

    @property
    def fragment_id(self):
        for loc in self.genomic_locations:
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
        return sorted([l for l in self.locations
                       if l.fragment_id is not None], key=lambda l: l.start_0based)

    @property
    def genomic_start_0based(self):
        return self.genomic_locations[0].start_0based

    def to_dict(self):
        return dict(
            recombination=dict(
                type=self.recombination.__class__.__name__,
                sites=self.recombination.__dict__
            ),
            genomic_locations=[
                dict(
                    site=loc.site,
                    fragment_id=loc.fragment_id,
                    fragment_name=loc.fragment.name,
                    start=loc.start_0based + 1
                ) for loc in self.genomic_locations
            ]
        )


class IntegrationEvent(Event):

    def get_integrated_aligned_with_site_direction(self, insert, is_insert_circular):
        insert = insert * 3
        site_insert = self.recombination.site_insert
        if site_insert not in insert:
            insert = rc(insert)
        assert insert.index(site_insert) >= 0
        bps = find_indices(insert, site_insert)
        assert len(bps) >= 2
        return insert[bps[0] + len(site_insert):bps[1]]

    def run(self, new_fragment, insert, is_insert_circular):
        genomic_locations = self.genomic_locations
        integrated = self.get_integrated_aligned_with_site_direction(insert, is_insert_circular)
        bps_to_replace = len(self.recombination.site_genome)

        if self.is_reversed():
            integrated = rc(integrated)
            new_site_left = rc(self.recombination.recombined_site_right_genome)
            new_site_right = rc(self.recombination.recombined_site_left_genome)
        else:
            new_site_left = self.recombination.recombined_site_left_genome
            new_site_right = self.recombination.recombined_site_right_genome

        integrated = new_site_left + integrated + new_site_right

        new_fragment.replace_bases(
            genomic_locations[0].start_0based + 1,
            bps_to_replace,
            integrated
        )

        return len(integrated) - bps_to_replace


class ExcisionEvent(Event):

    def run(self, new_fragment, insert, is_insert_circular):
        genomic_locations = self.genomic_locations
        assert len(genomic_locations) == 2

        if not self.is_reversed():
            assert genomic_locations[0].site == self.recombination.site_left
            assert genomic_locations[1].site == self.recombination.site_right
            new_site = self.recombination.recombined_site
            bps_to_remove = \
                genomic_locations[1].start_0based - \
                genomic_locations[0].start_0based + \
                len(self.recombination.site_right)
        else:
            assert genomic_locations[0].site == rc(self.recombination.site_right)
            assert genomic_locations[1].site == rc(self.recombination.site_left)
            new_site = rc(self.recombination.recombined_site)
            bps_to_remove = \
                genomic_locations[1].start_0based - \
                genomic_locations[0].start_0based + \
                len(self.recombination.site_left)

        new_fragment.replace_bases(
            genomic_locations[0].start_0based + 1,
            bps_to_remove,
            new_site
        )

        return len(new_site) - bps_to_remove


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

            bps_to_replace = \
                genomic_locations[1].start_0based - \
                genomic_locations[0].start_0based + \
                len(self.recombination.site_right)
        else:
            assert genomic_locations[0].site == rc(self.recombination.site_right)
            assert genomic_locations[1].site == rc(self.recombination.site_left)
            old_site_left = rc(self.recombination.site_right)
            old_site_right = rc(self.recombination.site_left)
            new_site_left = rc(self.recombination.recombined_site_right)
            new_site_right = rc(self.recombination.recombined_site_left)

            bps_to_replace = \
                genomic_locations[1].start_0based - \
                genomic_locations[0].start_0based + \
                len(self.recombination.site_left)

        sequence_to_replace = new_fragment.get_sequence(
            bp_lo=genomic_locations[0].start_0based + len(old_site_left) + 1,
            bp_hi=genomic_locations[1].start_0based + 1 - 1
        )
        new_sequence = new_site_left + rc(sequence_to_replace) + new_site_right
        assert \
            len(old_site_left) + len(sequence_to_replace) + len(old_site_right) == bps_to_replace

        new_fragment.replace_bases(
            genomic_locations[0].start_0based + 1,
            bps_to_replace,
            new_sequence
        )

        return len(new_sequence) - \
            (len(old_site_left) + len(sequence_to_replace) + len(old_site_right))


class RMCEEvent(Event):

    def get_integrated_aligned_with_site_direction(self, insert):
        insert = insert * 2
        site_left_insert = self.recombination.site_left_insert
        site_right_insert = self.recombination.site_right_insert
        if site_left_insert not in insert:
            insert = rc(insert)
        # below, to handle when sites are across circular boundary
        left_trimmed_insert = insert[insert.index(site_left_insert) + len(site_left_insert):]
        return left_trimmed_insert[:left_trimmed_insert.index(site_right_insert)]

    def run(self, new_fragment, insert, is_insert_circular):
        genomic_locations = self.genomic_locations
        integrated = self.get_integrated_aligned_with_site_direction(insert)

        if not self.is_reversed():
            assert genomic_locations[0].site == self.recombination.site_left_genome
            assert genomic_locations[1].site == self.recombination.site_right_genome
            new_site_left = self.recombination.recombined_site_left_genome
            new_site_right = self.recombination.recombined_site_right_genome

            bps_to_replace = genomic_locations[1].start_0based - \
                genomic_locations[0].start_0based + \
                len(self.recombination.site_right_genome)
        else:
            integrated = rc(integrated)

            assert genomic_locations[0].site == rc(self.recombination.site_right_genome)
            assert genomic_locations[1].site == rc(self.recombination.site_left_genome)
            new_site_left = rc(self.recombination.recombined_site_right_genome)
            new_site_right = rc(self.recombination.recombined_site_left_genome)

            bps_to_replace = genomic_locations[1].start_0based - \
                genomic_locations[0].start_0based + \
                len(self.recombination.site_left_genome)

        integrated = new_site_left + integrated + new_site_right
        new_fragment.replace_bases(
            genomic_locations[0].start_0based + 1,
            bps_to_replace,
            integrated
        )

        return len(integrated) - bps_to_replace


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
        template = sequence * 2
    template = template.lower()

    locations = []
    for site in sites:
        locations.extend(
            find_query_locations_on_duplicated_template(
                template, len(sequence), is_circular, site, fragment_obj=fragment_obj
            )
        )

    uniqdict = {(l.fragment_id, l.start_0based):l for l in locations}
    locations = list(uniqdict.values())
    return locations


def find_query_locations_on_duplicated_template(
        template, sequence_len, is_fragment_or_insert_circular, site,
        fragment_obj=None):

    indices = find_indices(template, site)
    indices = [i for i in indices if i < sequence_len]
    return [
      SiteLocation(
          site,
          fragment_obj,
          template[0:sequence_len] if fragment_obj is None else None,
          is_fragment_or_insert_circular,
          i,
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
            self.parent_fragments = [
                f.indexed_fragment() for f in self.parent_genome.fragments.all()
            ]

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

        if self.events is None:
            self.determine_site_locations()
            self.events = []
            self.errors = []

            # note that we iterate through allowed() in order - earlier
            # recombination definitions take precedence
            for recombination in self.allowed():
                self.events.extend(recombination.events(self.locations, self.errors))

    def run_reaction(self, new_genome_name, notes=None):
        self.group_into_events()

        if len(self.events) == 0:
            print("errors", self.errors)
            return

        event_insert_sites = [
            site
            for event in self.events
            for site in event.recombination.required_insert_sites()
        ]
        if len(event_insert_sites) == 0 and self.insert:
            print("errors", self.errors)
            print("has insert, but no events use insert")
            return

        new_genome = self.parent_genome.update()
        if new_genome_name:
            new_genome.name = new_genome_name
        else:
            new_genome.name = self.parent_genome.name + " SSR modified"
        new_genome.notes = notes
        new_genome.save()

        old_to_new_fragment_dict = {}

    # needs to run through events by fragment and start bp in reverse order
    # on the fragment, so we don't have to adjust coordinates of events
    # that have not been run yet.

        events = sorted(self.events,
                        key=lambda e: (e.fragment_id, -e.genomic_start_0based))
        while len(events) > 0:
            event = events[0]

            if event.fragment_id not in old_to_new_fragment_dict:
                with new_genome.update_fragment_by_fragment_id(
                    event.fragment_id, new_fragment=True
                ) as f:
                    old_to_new_fragment_dict[event.fragment_id] = f

            new_fragment = old_to_new_fragment_dict[event.fragment_id]
            event.run(new_fragment, self.insert, self.is_insert_circular)

            events = events[1:]

        print("new genome in run_reaction", new_genome)
        return new_genome
