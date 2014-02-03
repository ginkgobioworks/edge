from contextlib import contextmanager
from edge.models import *


class Genome_Updater(Genome):
    class Meta:
        app_label = "edge"
        proxy = True

    def save(self, *args, **kwargs):
        super(Genome_Updater, self).save(*args, **kwargs)
        return Genome.objects.get(pk=self.pk)

    @contextmanager
    def annotate_fragment_by_name(self, name):
        f = [x for x in self.fragments.all() if x.name == name]
        if len(f) != 1:
            raise Exception('Zero or more than one fragments have name %s' % (name,))
        u = f[0].annotate()
        yield u
        u.save()

    @contextmanager
    def annotate_fragment_by_fragment_id(self, fragment_id):
        f = [x for x in self.fragments.all() if x.id == fragment_id]
        if len(f) != 1:
            raise Exception('Zero or more than one fragments have ID %s' % (fragment_id,))
        u = f[0].annotate()
        yield u
        u.save()

    @contextmanager
    def update_fragment_by_name(self, name):
        if self.parent is None:
            raise Exception('Cannot update fragment without a parent genome. Try editing instead.')
        f = [x for x in self.fragments.all() if x.name == name]
        if len(f) != 1:
            raise Exception('Zero or more than one fragments have name %s' % (name,))
        u = f[0].update(name)
        yield u
        u.save()
        self._add_updated_fragment(u)

    @contextmanager
    def update_fragment_by_fragment_id(self, fragment_id):
        if self.parent is None:
            raise Exception('Cannot update fragment without a parent genome. Try editing instead.')
        f = [x for x in self.fragments.all() if x.id == fragment_id]
        if len(f) != 1:
            raise Exception('Zero or more than one fragments have ID %s' % (fragment_id,))
        u = f[0].update(f[0].name)
        yield u
        u.save()
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

    def import_gff(self, gff_fasta_fn):
        from edge.importer import GFFImporter
        GFFImporter(self, gff_fasta_fn).do_import()
