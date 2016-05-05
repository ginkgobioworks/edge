from contextlib import contextmanager
from edge.models.fragment import Fragment


class Genome_Updater(object):
    """
    Mixin with helpers for updating genome.
    """

    @contextmanager
    def annotate_fragment_by_name(self, name):
        f = [x for x in self.fragments.all() if x.name == name]
        if len(f) != 1:
            raise Exception('Zero or more than one fragments have name %s' % (name,))
        u = f[0].indexed_fragment()
        yield u

    @contextmanager
    def annotate_fragment_by_fragment_id(self, fragment_id):
        f = [x for x in self.fragments.all() if x.id == fragment_id]
        if len(f) != 1:
            raise Exception('Zero or more than one fragments have ID %s' % (fragment_id,))
        u = f[0].indexed_fragment()
        yield u

    @contextmanager
    def update_fragment_by_name(self, name, new_name=None):
        if self.parent is None:
            raise Exception('Cannot update fragment without a parent genome. Try editing instead.')
        f = [x for x in self.fragments.filter(name=name)]
        if len(f) != 1:
            raise Exception('Zero or more than one fragments have name %s' % (name,))
        new_name = name if new_name is None else new_name
        u = f[0].indexed_fragment().update(new_name)
        yield u
        self._add_updated_fragment(u)

    @contextmanager
    def update_fragment_by_fragment_id(self, fragment_id, new_name=None, new_fragment=True):
        if self.parent is None:
            raise Exception('Cannot update fragment without a parent genome. Try editing instead.')
        f = [x for x in self.fragments.filter(id=fragment_id)]
        if len(f) != 1:
            raise Exception('Zero or more than one fragments have ID %s' % (fragment_id,))
        new_name = f[0].name if new_name is None else new_name
        u = f[0].indexed_fragment()
        if new_fragment is True:
            u = u.update(new_name)
        yield u
        if new_fragment is True:
            self._add_updated_fragment(u)

    def add_fragment(self, name, sequence, circular=False):
        if len(sequence) == 0:
            raise Exception('Cannot create a fragment of length zero')
        new_fragment = Fragment.create_with_sequence(name=name,
                                                     sequence=sequence,
                                                     circular=circular)
        self.genome_fragment_set.create(fragment=new_fragment, inherited=False)
        return new_fragment

    def _add_updated_fragment(self, fragment):
        existing_fragment_ids = [f.id for f in self.fragments.all()]
        if fragment.parent_id in existing_fragment_ids:
            gf = self.genome_fragment_set.get(fragment=fragment.parent)
            gf.fragment = fragment
            gf.inherited = False
            gf.save()
        else:
            raise Exception('Fragment parent not part of the genome')
