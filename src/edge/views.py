import json
import random
from django.core.urlresolvers import reverse
from django.http import HttpResponse
from django.views.generic.base import View
from django.shortcuts import get_object_or_404
from django.db.models import Q
from edge.models import *
from edge.tasks import build_genome_blastdb


def get_genome_or_404(pk):
    return get_object_or_404(Genome, pk=pk)


def get_fragment_or_404(pk):
    return get_object_or_404(Fragment, pk=pk)


class ViewBase(View):
    def get(self, request, *args, **kwargs):
        res = self.on_get(request, *args, **kwargs)
        return HttpResponse(json.dumps(res), content_type='application/json')

    def put(self, request, *args, **kwargs):
        res, status = self.on_put(request, *args, **kwargs)
        return HttpResponse(json.dumps(res), status=status, content_type='application/json')

    def post(self, request, *args, **kwargs):
        res, status = self.on_post(request, *args, **kwargs)
        return HttpResponse(json.dumps(res), status=status, content_type='application/json')


class RequestParser(object):

    def __init__(self):
        self.__args = []

    def add_argument(self, name, field_type, required=False, default=None, location='get'):
        if type(field_type) not in (list, tuple):
            if field_type == str:
                field_type = [str, unicode]
            else:
                field_type = [field_type]
        self.__args.append((name, field_type, required, default, location))

    def parse_args(self, request):
        json_payload = None
        args = {}
        for name, field_type, required, default, location in self.__args:
            if location == 'get':
                d = request.GET
            elif location == 'post':
                d = request.POST
            else:
                if json_payload is None:
                    json_payload = json.loads(request.body)
                    d = json_payload
            if name not in d and required:
                raise Exception('Missing required field "%s"' % (name,))
            if name not in d:
                args[name] = default
            else:
                v = d[name]
                if type(v) not in field_type:
                    if int in field_type:
                        v = int(v)
                    elif float in field_type:
                        v = float(v)
                    else:
                        raise Exception('Field should be of type "%s", got "%s"' %
                                        (field_type, type(v)))
                args[name] = v
        return args


fragment_parser = RequestParser()
fragment_parser.add_argument('name', field_type=str, required=True, location='json')
fragment_parser.add_argument('sequence', field_type=str, required=True, location='json')
fragment_parser.add_argument('circular', field_type=bool, default=False, location='json')


class FragmentView(ViewBase):

    @staticmethod
    def to_dict(fragment, compute_length=True):
        length = fragment.est_length
        if compute_length is True and fragment.has_location_index:
            length = fragment.indexed_fragment().length
        return dict(id=fragment.id,
                    uri=reverse('fragment', kwargs=dict(fragment_id=fragment.id)),
                    name=fragment.name,
                    circular=fragment.circular,
                    parent_id=fragment.parent.id if fragment.parent else None,
                    length=length)

    def on_get(self, request, fragment_id):
        fragment = get_fragment_or_404(fragment_id)
        return FragmentView.to_dict(fragment)


class FragmentSequenceView(ViewBase):

    def on_get(self, request, fragment_id):
        q_parser = RequestParser()
        q_parser.add_argument('f', field_type=int, location='get')
        q_parser.add_argument('l', field_type=int, location='get')
        args = q_parser.parse_args(request)
        f = args['f']
        l = args['l']

        fragment = get_fragment_or_404(fragment_id)
        s = fragment.indexed_fragment().get_sequence(bp_lo=f, bp_hi=l)
        if f is None:
            f = 1
        if l is None:
            l = f+len(s)-1
        return {'sequence': s, 'base_first': f, 'base_last': l}


class FragmentAnnotationsView(ViewBase):

    @staticmethod
    def to_dict(annotation):
        return dict(base_first=annotation.base_first, base_last=annotation.base_last,
                    name=annotation.feature.name,
                    type=annotation.feature.type,
                    strand=annotation.feature.strand,
                    qualifiers=annotation.feature.qualifiers,
                    feature_full_length=annotation.feature.length,
                    feature_base_first=annotation.feature_base_first,
                    feature_base_last=annotation.feature_base_last)

    def on_get(self, request, fragment_id):
        q_parser = RequestParser()
        q_parser.add_argument('f', field_type=int, location='get')
        q_parser.add_argument('l', field_type=int, location='get')
        q_parser.add_argument('m', field_type=int, location='get')
        args = q_parser.parse_args(request)
        f = args['f']
        l = args['l']
        m = args['m']

        fragment = get_fragment_or_404(fragment_id)
        annotations = fragment.indexed_fragment().annotations(bp_lo=f, bp_hi=l)
        if m is not None and len(annotations) > m:
            to_return = []
            while len(to_return) < m:
                i = random.randint(0, len(annotations)-1)
                to_return.append(annotations[i])
                new_a = annotations[0:i]+annotations[i+1:]
                annotations = new_a
            annotations = to_return
        return [FragmentAnnotationsView.to_dict(annotation) for annotation in annotations]

    def on_post(self, request, fragment_id):
        annotation_parser = RequestParser()
        annotation_parser.add_argument('base_first', field_type=int, required=True, location='json')
        annotation_parser.add_argument('base_last', field_type=int, required=True, location='json')
        annotation_parser.add_argument('name', field_type=str, required=True, location='json')
        annotation_parser.add_argument('type', field_type=str, required=True, location='json')
        annotation_parser.add_argument('strand', field_type=int, required=True, location='json')

        args = annotation_parser.parse_args(request)
        fragment = get_fragment_or_404(fragment_id)
        fragment = fragment.indexed_fragment()
        fragment.annotate(first_base1=args['base_first'],
                          last_base1=args['base_last'],
                          name=args['name'],
                          type=args['type'],
                          strand=args['strand'])
        return {}, 201


class FragmentListView(ViewBase):

    def on_get(self, request):
        q_parser = RequestParser()
        q_parser.add_argument('q', field_type=str, location='get')
        q_parser.add_argument('s', field_type=int, location='get', default=0)
        q_parser.add_argument('p', field_type=int, location='get', default=100)
        args = q_parser.parse_args(request)
        s = args['s']
        p = args['p']
        q = args['q']
        p = 200 if p > 200 else p
        if q is not None and q.strip() != '':
            fragments = Fragment.non_genomic_fragments(Q(name__icontains=q), s, s+p)
        else:
            fragments = Fragment.non_genomic_fragments(None, s, s+p)
        return [FragmentView.to_dict(fragment) for fragment in fragments]

    def on_post(self, request):
        args = fragment_parser.parse_args(request)
        fragment = Fragment.create_with_sequence(name=args['name'],
                                                 sequence=args['sequence'],
                                                 circular=args['circular'])
        return FragmentView.to_dict(fragment), 201


class GenomeView(ViewBase):

    @staticmethod
    def to_dict(genome, include_changes=True, compute_length=True, include_fragments=True):
        changes = None
        if include_changes is True:
            if genome.has_location_index:
                genome = genome.indexed_genome()
                changes = genome.changed_locations_by_fragment()
                changes = {f.id: v for f, v in changes.iteritems()}
        fragments = None
        if include_fragments:
            fragments = []
            for f in genome.fragments.all():
                d = FragmentView.to_dict(f, compute_length=compute_length)
                if changes is not None and f.id in changes:
                    d['changes'] = changes[f.id]
                fragments.append(d)

        return dict(id=genome.id,
                    uri=reverse('genome', kwargs=dict(genome_id=genome.id)),
                    name=genome.name,
                    notes=genome.notes,
                    parent_id=genome.parent_id,
                    parent_name=genome.parent.name if genome.parent is not None else '',
                    fragments=fragments)

    def on_get(self, request, genome_id):
        genome = get_genome_or_404(genome_id)
        return GenomeView.to_dict(genome)


class GenomeAnnotationsView(ViewBase):

    def on_get(self, request, genome_id):
        genome = get_genome_or_404(genome_id)
        q_parser = RequestParser()
        q_parser.add_argument('q', field_type=str, required=True)
        args = q_parser.parse_args(request)

        res = []
        fragment_annotations = genome.indexed_genome().find_annotation_by_name(args['q'])
        for fragment_id in fragment_annotations:
            fragment = get_fragment_or_404(fragment_id)
            annotations = fragment_annotations[fragment_id]
            d = FragmentView.to_dict(fragment)
            a = [FragmentAnnotationsView.to_dict(x) for x in annotations]
            res.append((d, a))

        return res


class GenomeFragmentListView(ViewBase):

    def on_post(self, request, genome_id):  # adding new fragment
        args = fragment_parser.parse_args(request)
        genome = get_genome_or_404(genome_id)
        fragment = None
        fragment = genome.add_fragment(name=args['name'], sequence=args['sequence'],
                                       circular=args['circular'])
        return FragmentView.to_dict(fragment), 201


class GenomeFragmentView(ViewBase):

    def insert_bases(self, request, genome_id, fragment_id):
        op_parser = RequestParser()
        op_parser.add_argument('name', field_type=str, required=True, location='json')
        op_parser.add_argument('before_bp', field_type=int, required=True, location='json')
        op_parser.add_argument('sequence', field_type=str, required=True, location='json')
        args = op_parser.parse_args(request)

        genome = get_genome_or_404(genome_id)
        fragment = get_fragment_or_404(fragment_id)

        u = genome.update(name=args['name'])
        with u.update_fragment_by_fragment_id(fragment.id) as f:
            f.insert_bases(args['before_bp'], args['sequence'])

        build_genome_blastdb.delay(u.id)
        return GenomeView.to_dict(u), 201

    def remove_bases(self, request, genome_id, fragment_id):
        op_parser = RequestParser()
        op_parser.add_argument('name', field_type=str, required=True, location='json')
        op_parser.add_argument('before_bp', field_type=int, required=True, location='json')
        op_parser.add_argument('length', field_type=int, required=True, location='json')
        args = op_parser.parse_args(request)

        genome = get_genome_or_404(genome_id)
        fragment = get_fragment_or_404(fragment_id)

        u = genome.update(name=args['name'])
        with u.update_fragment_by_fragment_id(fragment.id) as f:
            f.remove_bases(args['before_bp'], args['length'])

        build_genome_blastdb.delay(u.id)
        return GenomeView.to_dict(u), 201

    def replace_bases(self, request, genome_id, fragment_id):
        op_parser = RequestParser()
        op_parser.add_argument('name', field_type=str, required=True, location='json')
        op_parser.add_argument('before_bp', field_type=int, required=True, location='json')
        op_parser.add_argument('length', field_type=int, required=True, location='json')
        op_parser.add_argument('sequence', field_type=str, required=True, location='json')
        args = op_parser.parse_args(request)

        genome = get_genome_or_404(genome_id)
        fragment = get_fragment_or_404(fragment_id)

        u = genome.update(name=args['name'])
        with u.update_fragment_by_fragment_id(fragment.id) as f:
            f.replace_bases(args['before_bp'], args['length'], args['sequence'])

        build_genome_blastdb.delay(u.id)
        return GenomeView.to_dict(u), 201

    def insert_fragment(self, request, genome_id, fragment_id):
        op_parser = RequestParser()
        op_parser.add_argument('name', field_type=str, required=True, location='json')
        op_parser.add_argument('before_bp', field_type=int, required=True, location='json')
        op_parser.add_argument('fragment_id', field_type=int, required=True, location='json')
        args = op_parser.parse_args(request)

        genome = get_genome_or_404(genome_id)
        fragment = get_fragment_or_404(fragment_id)
        new_fragment = get_fragment_or_404(args['fragment_id'])

        u = genome.update(name=args['name'])
        with u.update_fragment_by_fragment_id(fragment.id) as f:
            f.insert_fragment(args['before_bp'], new_fragment)

        build_genome_blastdb.delay(u.id)
        return GenomeView.to_dict(u), 201

    def replace_with_fragment(self, request, genome_id, fragment_id):
        op_parser = RequestParser()
        op_parser.add_argument('name', field_type=str, required=True, location='json')
        op_parser.add_argument('before_bp', field_type=int, required=True, location='json')
        op_parser.add_argument('length', field_type=int, required=True, location='json')
        op_parser.add_argument('fragment_id', field_type=int, required=True, location='json')
        args = op_parser.parse_args(request)

        genome = get_genome_or_404(genome_id)
        fragment = get_fragment_or_404(fragment_id)
        new_fragment = get_fragment_or_404(args['fragment_id'])

        u = genome.update(name=args['name'])
        with u.update_fragment_by_fragment_id(fragment.id) as f:
            f.replace_with_fragment(args['before_bp'], args['length'], new_fragment)

        build_genome_blastdb.delay(u.id)
        return GenomeView.to_dict(u), 201

    def on_put(self, request, genome_id, fragment_id):  # updating genome and fragment
        op_parser = RequestParser()
        op_parser.add_argument('op', field_type=str, required=True, location='json')
        args = op_parser.parse_args(request)

        if args['op'] == 'insert_bases':
            return self.insert_bases(request, genome_id, fragment_id)
        elif args['op'] == 'remove_bases':
            return self.remove_bases(request, genome_id, fragment_id)
        elif args['op'] == 'replace_bases':
            return self.replace_bases(request, genome_id, fragment_id)
        elif args['op'] == 'insert_fragment':
            return self.insert_fragment(request, genome_id, fragment_id)
        elif args['op'] == 'replace_with_fragment':
            return self.replace_with_fragment(request, genome_id, fragment_id)


class GenomeListView(ViewBase):

    def on_get(self, request):
        if 'f' in request.GET:
            fragment_ids = []
            for f in request.GET.getlist('f'):
                try:
                    fragment_ids.append(int(f))
                except ValueError:
                    return []
            if len(fragment_ids) == 0:
                return []

            sql_joins = 5
            q = Genome.objects.filter(genome_fragment__fragment_id=fragment_ids[0])
            for i in range(1, sql_joins):
                if i < len(fragment_ids):
                    q = q.filter(genome_fragment__fragment_id=fragment_ids[i])
            candidates = list(q)
            genomes = []
            for g in candidates:
                x = Genome.objects.raw("""
                  SELECT edge_genome.id,
                         GROUP_CONCAT(edge_genome_fragment.fragment_id) as fragment_id_list
                    FROM edge_genome
                    JOIN edge_genome_fragment ON edge_genome_fragment.genome_id = edge_genome.id
                   WHERE edge_genome.id = %s
                """, [g.id])
                x = list(x)[0]
                id_list = [int(n) for n in x.fragment_id_list.split(',')]
                if set(id_list) == set(fragment_ids):
                    genomes.append(g)

        else:
            q_parser = RequestParser()
            q_parser.add_argument('q', field_type=str, location='get')
            q_parser.add_argument('s', field_type=int, location='get', default=0)
            q_parser.add_argument('p', field_type=int, location='get', default=100)
            args = q_parser.parse_args(request)
            s = args['s']
            p = args['p']
            q = args['q']
            p = 200 if p > 200 else p
            if q is not None and q.strip() != '':
                where = Q(name__icontains=q)
                try:
                    id = int(q)
                except:
                    pass
                else:
                    where = where | Q(id=q)
                genomes = Genome.objects.filter(active=True)\
                                        .filter(where)\
                                        .order_by('-id')[s:s+p]
            else:
                genomes = Genome.objects.filter(active=True).order_by('-id')[s:s+p]
        return [GenomeView.to_dict(genome,
                                   include_changes=False,
                                   compute_length=False,
                                   include_fragments=False)
                for genome in genomes]

    def on_post(self, request):
        genome_parser = RequestParser()
        genome_parser.add_argument('name', field_type=str, required=True, location='json')
        genome_parser.add_argument('notes', field_type=str, location='json')

        args = genome_parser.parse_args(request)
        genome = Genome.create(name=args['name'], notes=args['notes'])
        return GenomeView.to_dict(genome), 201


class GenomeBlastView(ViewBase):

    def on_post(self, request, genome_id):
        from edge.blast import blast_genome
        genome = get_genome_or_404(genome_id)

        parser = RequestParser()
        parser.add_argument('query', field_type=str, required=True, location='json')
        parser.add_argument('program', field_type=str, required=True, location='json')

        args = parser.parse_args(request)
        results = blast_genome(genome, args['program'], args['query'])
        results = [r.to_dict() for r in results]
        return results, 200


class GenomePcrView(ViewBase):

    def on_post(self, request, genome_id):
        from edge.pcr import pcr_from_genome
        genome = get_genome_or_404(genome_id)

        parser = RequestParser()
        parser.add_argument('primers', field_type=list, required=True, location='json')

        args = parser.parse_args(request)
        primers = args['primers']
        if len(primers) != 2:
            raise Exception('Expecting two primers, got %s' % (primers,))

        r = pcr_from_genome(genome, primers[0], primers[1])
        r = (r[0], [b.to_dict() for b in r[1]], [b.to_dict() for b in r[2]], r[3])
        return r, 200


class GenomeRecombinationView(ViewBase):

    def on_post(self, request, genome_id):
        from edge.recombine import find_swap_region, recombine
        genome = get_genome_or_404(genome_id)

        parser = RequestParser()
        parser.add_argument('cassette', field_type=str, required=True, location='json')
        parser.add_argument('homology_arm_length', field_type=int, required=True, location='json')
        parser.add_argument('create', field_type=bool, required=True, location='json')
        parser.add_argument('name', field_type=str, required=False, default=None, location='json')

        args = parser.parse_args(request)
        cassette = args['cassette']
        arm_length = args['homology_arm_length']
        create = args['create']
        name = args['name']

        if create is False:
            r = find_swap_region(genome, cassette, arm_length)
            return [x.to_dict() for x in r], 200
        else:
            c = recombine(genome, cassette, arm_length, name)
            if c is None:
                return None, 400
            else:
                build_genome_blastdb.delay(c.id)
                return GenomeView.to_dict(c), 201
