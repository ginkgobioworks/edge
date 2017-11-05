from edge.models.chunk import Feature


class Fragment_Annotator(object):
    """
    Mixin for annotating a fragment.
    """

    def _add_feature(self, name, type, length, strand, qualifiers=None, operation=None):
        if strand not in (1, -1, None):
            raise Exception('Strand must be 1, -1, or None')
        qualifiers = {} if qualifiers is None else qualifiers
        f = Feature(name=name, type=type, length=length, strand=strand, operation=operation)
        f.set_qualifiers(qualifiers)
        f.save()
        return f

    def annotate(self, first_base1, last_base1, name, type, strand,
                 qualifiers=None, operation=None):
        if self.circular and last_base1 < first_base1:
            # has to figure out the total length from last chunk
            length = self.length - first_base1 + 1 + last_base1
        else:
            length = last_base1 - first_base1 + 1
            if length <= 0:
                raise Exception('Annotation must have length one or more')

        prev_chunk, annotation_start = self._find_and_split_before(first_base1)
        annotation_end, next_chunk = self._find_and_split_before(last_base1 + 1)

        # did two splits, so must reload annotation_start in case that got split
        annotation_start = annotation_start.reload()

        new_feature = self._add_feature(name, type, length, strand, qualifiers, operation)

        # now, starting with chunk annotation_start, walk through chunks until
        # we hit annotation_end, and add annotation for each chunk
        chunk = annotation_start
        a_i = 1
        while True:
            fc = self.fragment_chunk(chunk)
            self._annotate_chunk(chunk, new_feature, a_i, a_i + len(chunk.sequence) - 1)
            a_i += len(chunk.sequence)
            if chunk.id == annotation_end.id:
                break
            chunk = fc.next_chunk
            if chunk is None:
                chunk = self.start_chunk

        return new_feature
