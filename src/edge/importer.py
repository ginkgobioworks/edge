import time

from BCBio import GFF
from django.db import connection

from edge.models import Fragment, Fragment_Chunk_Location


class GFFImporter(object):
    def __init__(self, genome, gff_fasta_fn):
        self.__genome = genome
        self.__gff_fasta_fn = gff_fasta_fn

    def do_import(self):
        in_file = self.__gff_fasta_fn
        in_handle = open(in_file)

        # In DEBUG=True mode, Django keeps list of queries and blows up memory
        # usage when doing a big import. The following line disables this
        # logging.
        connection.use_debug_cursor = False

        for rec in GFF.parse(in_handle):
            if self.__genome.fragments.filter(name=rec.id).count() > 0:
                print("skipping %s, already imported" % rec.id)
            else:
                f = GFFFragmentImporter(rec).do_import()
                self.__genome.genome_fragment_set.create(fragment=f, inherited=False)

        # Be nice and turn debug cursor back on
        connection.use_debug_cursor = True
        in_handle.close()


class GFFFragmentImporter(object):
    def __init__(self, gff_rec):
        self.__rec = gff_rec
        self.__sequence = None
        self.__features = None
        self.__fclocs = None
        self.__subfeatures_dict = {}

    def do_import(self):
        self.parse_gff()
        t0 = time.time()
        f = self.build_fragment()
        print("build fragment: %.4f" % (time.time() - t0,))
        t0 = time.time()
        self.annotate(f)
        print("annotate: %.4f" % (time.time() - t0,))
        return f

    def parse_gff(self):
        name_fields = (
            "name",
            "Name",
            "gene",
            "locus",
            "locus_tag",
            "product",
            "protein_id",
        )

        self.__sequence = str(self.__rec.seq)
        seqlen = len(self.__sequence)
        print("%s: %s" % (self.__rec.id, seqlen))

        features = []
        for feature in self.__rec.features:
            # skip features that cover the entire sequence
            if feature.location.start == 0 and feature.location.end == seqlen:
                continue

            # get name
            name = feature.id
            if name == "":
                name = feature.type
            for field in name_fields:
                if field in feature.qualifiers:
                    v = feature.qualifiers[field]
                    if len(v) > 0:
                        name = v[0]
                        break
            name = name[0:100]

            # get qualifiers
            qualifiers = {}
            for field in feature.qualifiers:
                v = feature.qualifiers[field]
                if len(v) > 0:
                    qualifiers[field] = v

            # start in Genbank format is start after, so +1 here
            features.append(
                (
                    int(feature.location.start) + 1,
                    int(feature.location.end),
                    name,
                    feature.type,
                    feature.strand,
                    qualifiers,
                )
            )

            feature_name = name
            # add sub features for chunking for CDS only
            self.__subfeatures_dict[feature_name] = []
            for sub in feature.sub_features:
                # change name for sub feature
                if hasattr(sub, 'qualifiers') and 'Name' in sub.qualifiers:
                    subfeature_name = sub.qualifiers['Name'][0]
                elif sub.id != '':
                    subfeature_name = sub.id
                else:
                    subfeature_name = feature_name

                # check that the type is right
                if sub.type.upper() in ['MRNA', 'CDS', 'EXON', 'INTRON']:
                    sub_tup = (
                            int(sub.location.start) + 1,
                            int(sub.location.end),
                            subfeature_name,
                            sub.type,
                            sub.strand,
                            qualifiers,
                        )

                    # if it has no id, it belongs to the feature
                    # otherwise, mark it as its own feature
                    if subfeature_name == feature_name:
                        self.__subfeatures_dict[feature_name].append(sub_tup)
                    else:
                        if subfeature_name in self.__subfeatures_dict:
                            self.__subfeatures_dict[subfeature_name].append(sub_tup)
                        else:
                            features.append(sub_tup)
                            self.__subfeatures_dict[subfeature_name] = [sub_tup]

                    for sub_sub in sub.sub_features:
                        # change name for sub sub feature
                        if hasattr(sub, 'qualifiers') and 'Name' in sub.qualifiers:
                            subsubfeature_name = sub_sub.qualifiers['Name'][0]
                        elif sub_sub.id != '':
                            subsubfeature_name = sub_sub.id
                        else:
                            subsubfeature_name = subfeature_name

                        sub_sub_tup = (
                                int(sub_sub.location.start) + 1,
                                int(sub_sub.location.end),
                                subsubfeature_name,
                                sub_sub.type,
                                sub_sub.strand,
                                qualifiers,
                            )

                        # if it has no id and the sub feature has no id, it belongs to the feature
                        # if it has no id and the sub feature has id, it belongs to the sub feature
                        # otherwise, mark it as its own feature
                        if subsubfeature_name == feature_name:
                            self.__subfeatures_dict[feature_name].append(sub_sub_tup)
                        elif subsubfeature_name == subfeature_name:
                            self.__subfeatures_dict[subfeature_name].append(sub_sub_tup)
                        else:
                            if subsubfeature_name in self.__subfeatures_dict:
                                self.__subfeatures_dict[subsubfeature_name].append(sub_sub_tup)
                            else:
                                features.append(sub_sub_tup)
                                self.__subfeatures_dict[subsubfeature_name] = [sub_sub_tup]

        self.__features = features

        # update features made from only subfeatures
        for feature in features:
            if self.__subfeatures_dict[feature[2]] != []:
                features.remove(feature)
                new_start = self.__subfeatures_dict[feature[2]][0][0]
                new_end = self.__subfeatures_dict[feature[2]][-1][1]
                new_feature = (new_start, new_end, feature[2], feature[3], feature[4], feature[5])
                features.append(new_feature)

    def build_fragment(self):
        # pre-chunk the fragment sequence at feature start and end locations.
        # there should be no need to further divide any chunk during import.
        starts_and_ends = []
        for feature in self.__features:
            name = feature[2]
            starts_and_ends.append(feature[0])
            starts_and_ends.append(feature[1] + 1)
            for subfeature in self.__subfeatures_dict[name]:
                starts_and_ends.append(subfeature[0])
                starts_and_ends.append(subfeature[1] + 1)
        break_points = sorted(list(set(starts_and_ends)))

        cur_len = 0
        chunk_sizes = []
        seq_len = len(self.__sequence)
        for i, bp in enumerate(break_points):
            if i == 0:
                if bp > 1:
                    chunk_sizes.append(break_points[i] - 1)
                    cur_len += chunk_sizes[-1]
            else:
                chunk_sizes.append(break_points[i] - break_points[i - 1])
                cur_len += chunk_sizes[-1]

        if cur_len < seq_len:
            chunk_sizes.append(seq_len - cur_len)

        new_fragment = Fragment(
            name=self.__rec.id, circular=False, parent=None, start_chunk=None
        )
        new_fragment.save()
        new_fragment = new_fragment.indexed_fragment()

        # divide chunks bigger than a certain threshold to smaller chunks, to
        # allow insertion of sequence into database. e.g. MySQL has a packet
        # size that prevents chunks that are too large from being inserted.
        chunk_size_limit = 1000000
        new_chunk_sizes = []
        for original_chunk_size in chunk_sizes:
            if original_chunk_size < chunk_size_limit:
                new_chunk_sizes.append(original_chunk_size)
            else:
                divided_chunks = []
                while original_chunk_size > 0:
                    divided_chunks.append(min(original_chunk_size, chunk_size_limit))
                    original_chunk_size -= chunk_size_limit
                new_chunk_sizes.extend(divided_chunks)
        chunk_sizes = new_chunk_sizes
        print("%d chunks" % (len(chunk_sizes),))

        prev = None
        fragment_len = 0
        for chunk_size in chunk_sizes:
            t0 = time.time()
            prev = new_fragment._append_to_fragment(
                prev,
                fragment_len,
                self.__sequence[fragment_len : fragment_len + chunk_size],
            )
            fragment_len += chunk_size
            print("add chunk to fragment: %.4f\r" % (time.time() - t0,), end="")

        return new_fragment

    def annotate(self, fragment):
        self.__fclocs = {
            c.base_first: c
            for c in Fragment_Chunk_Location.objects.select_related("chunk").filter(
                fragment=fragment
            )
        }

        for feature in self.__features:
            t0 = time.time()
            f_start, f_end, f_name, f_type, f_strand, f_qualifiers = feature
            # print('  %s %s: %s-%s %s' % (f_type, f_name, f_start, f_end, f_strand))
            new_feature = self._annotate_feature(
                fragment, f_start, f_end, f_name, f_type, f_strand, f_qualifiers
            )
            feature_base_first = 1
            sorted_subfeatures = sorted(
                self.__subfeatures_dict[f_name],
                key=lambda x: f_strand * x[0]
            )
            for subfeature in sorted_subfeatures:
                sf_start, sf_end, sf_name, sf_type, sf_strand, sf_qualifiers = subfeature
                self._annotate_feature(
                    fragment, sf_start, sf_end, sf_name, sf_type, sf_strand, sf_qualifiers,
                    feature=new_feature, feature_base_first=feature_base_first
                )
                feature_base_first += sf_end - sf_start + 1
            print("annotate feature: %.4f\r" % (time.time() - t0,), end="")
        print("\nfinished annotating feature")

    def _annotate_feature(
        self, fragment, first_base1, last_base1, name, type,
        strand, qualifiers, feature=None, feature_base_first=1
    ):
        if fragment.circular and last_base1 < first_base1:
            # has to figure out the total length from last chunk
            length = len(self.__sequence) - first_base1 + 1 + last_base1
        else:
            length = last_base1 - first_base1 + 1
            if length <= 0:
                raise Exception("Annotation must have length one or more")

        if first_base1 not in self.__fclocs or (
            (last_base1 < len(self.__sequence) and last_base1 + 1 not in self.__fclocs)) or (
            last_base1 > len(self.__sequence)
        ):
            raise Exception(
                "Cannot find appropriate sequence for feature: %s, start %s, end %s"
                % (name, first_base1, last_base1)
            )

        bases = []
        for key in self.__fclocs:
            fcloc = self.__fclocs[key]
            if fcloc.base_first >= first_base1 and fcloc.base_last <= last_base1:
                bases.append((fcloc.base_first, fcloc.base_last))
        bases.sort(key=lambda x: strand * x[0])

        length = 0
        for first_base1, last_base1 in bases:
            length += fragment.bp_covered_length(first_base1, last_base1)

        new_feature = fragment._add_feature(
            name, type, length, strand, qualifiers
        ) if feature is None else feature

        if feature is not None or self.__subfeatures_dict[name] == []:
            for first_base1, last_base1 in bases:
                region_length = fragment.bp_covered_length(first_base1, last_base1)
                fragment.annotate_chunk(
                    new_feature, feature_base_first, first_base1, last_base1
                )
                feature_base_first += region_length

        return new_feature
