import json

from edge.models.chunk import Chunk_Feature, Feature


class Fragment_Annotator(object):
    """
    Mixin for annotating a fragment.
    """

    def _add_feature(
        self, name, type, length, strand, qualifiers=None,
        operation=None, force_create=False
    ):
        if strand not in (1, -1, None):
            raise Exception("Strand must be 1, -1, or None")
        qualifiers = {} if qualifiers is None else qualifiers

        if force_create:
            f = Feature(
                name=name, type=type, length=length, strand=strand,
                operation=operation
            )
            f.set_qualifiers(qualifiers)
            f.save()
        else:
            f, _ = Feature.objects.get_or_create(
                name=name, type=type, length=length, strand=strand,
                operation=operation, _qualifiers=json.dumps(qualifiers)
            )

        return f

    def annotate_chunk(self, feature, feature_base_first, first_base1, last_base1):
        _, annotation_start = self._find_and_split_before(first_base1)
        annotation_end, _ = self._find_and_split_before(last_base1 + 1)

        # did two splits, so must reload annotation_start in case that got split
        annotation_start = annotation_start.reload()

        # now, starting with chunk annotation_start, walk through chunks until
        # we hit annotation_end, and add annotation for each chunk
        chunk = annotation_start
        a_i = 0
        cfs = []
        while True:
            fc = self.fragment_chunk(chunk)
            if feature.strand is None or feature.strand > 0:
                cf = self._create_chunk_annotation(
                    chunk,
                    feature,
                    feature_base_first + a_i,
                    feature_base_first + a_i + chunk.length - 1
                )
                if cf is not None:
                    cfs.append(cf)
            else:
                subfeature_length = self.bp_covered_length(first_base1, last_base1)
                feature_base_last = feature_base_first + subfeature_length - 1
                cf = self._create_chunk_annotation(
                    chunk,
                    feature,
                    feature_base_last - a_i - chunk.length + 1,
                    feature_base_last - a_i
                )
                if cf is not None:
                    cfs.append(cf)
            a_i += chunk.length
            if chunk.id == annotation_end.id:
                break
            chunk = fc.next_chunk
            if chunk is None:
                chunk = self.start_chunk

        # generate objects ahead of time, then bulk create at end of call
        Chunk_Feature.bulk_create(cfs)

    def bp_covered_length(self, first_base1, last_base1):
        if self.circular and last_base1 < first_base1:
            # has to figure out the total length from last chunk
            length = self.length - first_base1 + 1 + last_base1
        else:
            length = last_base1 - first_base1 + 1
            if length <= 0:
                raise Exception("Annotation must have length one or more")
        return length

    def annotate_chunks(
        self, bases, name, type, strand, qualifiers=None, operation=None,
        force_new_annotation=False
    ):
        length = 0
        for first_base1, last_base1 in bases:
            length += self.bp_covered_length(first_base1, last_base1)

        new_feature = self._add_feature(
            name, type, length, strand, qualifiers=qualifiers,
            operation=operation, force_create=force_new_annotation
        )
        feature_base_first = 1
        for first_base1, last_base1 in bases:
            region_length = self.bp_covered_length(first_base1, last_base1)
            self.annotate_chunk(
                new_feature, feature_base_first, first_base1, last_base1
            )
            feature_base_first += region_length

        return new_feature

    def annotate(
        self,
        first_base1,
        last_base1,
        name,
        type,
        strand,
        qualifiers=None,
        operation=None,
        force_new_annotation=False
    ):
        return self.annotate_chunks(
            [(first_base1, last_base1)],
            name,
            type,
            strand,
            qualifiers=qualifiers,
            operation=operation,
            force_new_annotation=force_new_annotation
        )
