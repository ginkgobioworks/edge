import os
import re
import json
import random
import tempfile

from django.urls import reverse
from django.http import HttpResponse, JsonResponse
from django.views.generic.base import View
from django.views.decorators.http import require_http_methods
from django.shortcuts import get_object_or_404
from django.db.models import Q
from django.db import transaction

from edge.models import Fragment, Genome, Operation
from edge.io import IO
from edge import import_gff
from edge.tasks import build_genome_blastdb


def genome_export(request, genome_id):
    get_genome_or_404(genome_id)
    io = IO(Genome.objects.get(pk=genome_id))
    response = HttpResponse(content_type='text/plain')
    response['Content-Disposition'] = "attachment; filename=\"g%s.gff\"" % genome_id
    io.to_gff_file(response)
    return response


@require_http_methods(['POST'])
def genome_import(request):
    res = {
      'imported_genomes': [],
    }

    for name in request.FILES:
        with tempfile.NamedTemporaryFile(mode='w+b', delete=False) as gff:
            for chuck in request.FILES.get(name).chunks():
                gff.write(chuck)
        g = import_gff(name, gff.name)
        os.unlink(gff.name)

        blastdb_task = build_genome_blastdb.delay(g.id)

        res['imported_genomes'].append({
            'id': g.id,
            'name': g.name,
            'blastdb_task_id': blastdb_task.id,
        })

    return JsonResponse(res)


def schedule_building_blast_db(genome_id, countdown=None):
    # scheduling building genome DB in the future, so a) our transaction has a
    # chance to commit and b) we avoid an immediate followup operation building
    # db on demand at the same time as this delayed building of database
    countdown = 30 if countdown is None else countdown
    build_genome_blastdb.apply_async((genome_id, ), countdown=countdown)


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
                field_type = [bytes, str]
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
        ll = args['l']

        fragment = get_fragment_or_404(fragment_id)
        s = fragment.indexed_fragment().get_sequence(bp_lo=f, bp_hi=ll)
        if f is None:
            f = 1
        if ll is None:
            ll = f + len(s) - 1
        return {'sequence': s, 'base_first': f, 'base_last': ll}


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
        ll = args['l']
        m = args['m']

        fragment = get_fragment_or_404(fragment_id)
        annotations = fragment.indexed_fragment().annotations(bp_lo=f, bp_hi=ll)
        if m is not None and len(annotations) > m:
            to_return = []
            while len(to_return) < m:
                i = random.randint(0, len(annotations) - 1)
                to_return.append(annotations[i])
                new_a = annotations[0:i] + annotations[i + 1:]
                annotations = new_a
            annotations = to_return
        return [FragmentAnnotationsView.to_dict(annotation) for annotation in annotations]

    @transaction.atomic()
    def on_post(self, request, fragment_id):
        annotation_parser = RequestParser()
        annotation_parser.add_argument('base_first', field_type=int, required=True, location='json')
        annotation_parser.add_argument('base_last', field_type=int, required=True, location='json')
        annotation_parser.add_argument('name', field_type=str, required=True, location='json')
        annotation_parser.add_argument('type', field_type=str, required=True, location='json')
        annotation_parser.add_argument('strand', field_type=int, required=True, location='json')
        annotation_parser.add_argument('qualifiers', field_type=dict, required=False,
                                       default=None, location='json')

        args = annotation_parser.parse_args(request)
        fragment = get_fragment_or_404(fragment_id)
        fragment = fragment.indexed_fragment()
        fragment.annotate(first_base1=args['base_first'],
                          last_base1=args['base_last'],
                          name=args['name'],
                          type=args['type'],
                          strand=args['strand'],
                          qualifiers=args['qualifiers'])
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
            fragments = Fragment.user_defined_fragments(Q(name__icontains=q), s, s + p)
        else:
            fragments = Fragment.user_defined_fragments(None, s, s + p)
        return [FragmentView.to_dict(fragment) for fragment in fragments]

    @transaction.atomic()
    def on_post(self, request):
        args = fragment_parser.parse_args(request)
        fragment = Fragment.create_with_sequence(name=args['name'],
                                                 sequence=args['sequence'],
                                                 circular=args['circular'])
        return FragmentView.to_dict(fragment), 201


class GenomeView(ViewBase):

    @staticmethod
    def op_to_dict(genome, op):
        choices = Operation._meta.get_field('type').choices
        type_str = [t[1] for t in choices if t[0] == op.type]
        if len(type_str) > 0:
            type_str = type_str[0]
        else:
            type_str = ''

        annotations = {}
        if genome.has_location_index:
            genome = genome.indexed_genome()
            for feature in op.feature_set.all():
                fragment_annotations = genome.find_annotation_by_feature(feature)
                for fragment_id, fragment_annotations in fragment_annotations.items():
                    if fragment_id not in annotations:
                        annotations[fragment_id] = []
                    annotations[fragment_id].extend(fragment_annotations)

        annotation_list = []
        for fragment_id in annotations:
            f = genome.fragments.filter(id=fragment_id)[0]
            v = annotations[fragment_id]
            v = [FragmentAnnotationsView.to_dict(x) for x in v]
            for a in v:
                a['fragment_id'] = fragment_id
                a['fragment_name'] = f.name
                annotation_list.append(a)

        d = dict(type=type_str, params=json.loads(op.params), annotations=annotation_list)
        return d

    @staticmethod
    def to_dict(genome, compute_length=True, include_fragments=True, include_operations=True):
        operations = []
        if include_operations:
            for op in genome.operation_set.order_by('id'):
                d = GenomeView.op_to_dict(genome, op)
                operations.append(d)

        fragments = None
        if include_fragments:
            fragments = []
            for f in genome.fragments.all():
                d = FragmentView.to_dict(f, compute_length=compute_length)
                fragments.append(d)

        d = dict(id=genome.id,
                 uri=reverse('genome', kwargs=dict(genome_id=genome.id)),
                 name=genome.name,
                 notes=genome.notes,
                 parent_id=genome.parent_id,
                 parent_name=genome.parent.name if genome.parent is not None else '',
                 fragments=fragments)

        if len(operations):
            d['operations'] = operations

        return d

    def on_get(self, request, genome_id):
        genome = get_genome_or_404(genome_id)
        return GenomeView.to_dict(genome)


class GenomeAnnotationsView(ViewBase):

    def on_get(self, request, genome_id):
        genome = get_genome_or_404(genome_id)
        q_parser = RequestParser()
        q_parser.add_argument('q', field_type=str, required=True)
        q_parser.add_argument('field', field_type=str, required=False, default=None)
        args = q_parser.parse_args(request)
        field = args['field']

        res = []
        if field is None:
            fragment_annotations = genome.indexed_genome().find_annotation_by_name(args['q'])
        else:
            fragment_annotations = genome.indexed_genome().find_annotation_by_qualifier(
                args['q'], fields=[field])
        for fragment_id in fragment_annotations:
            fragment = get_fragment_or_404(fragment_id)
            annotations = fragment_annotations[fragment_id]
            d = FragmentView.to_dict(fragment)
            a = [FragmentAnnotationsView.to_dict(x) for x in annotations]
            res.append((d, a))

        return res


class GenomeFragmentListView(ViewBase):

    @transaction.atomic()
    def on_post(self, request, genome_id):  # adding new fragment
        args = fragment_parser.parse_args(request)
        genome = get_genome_or_404(genome_id)
        fragment = None
        fragment = genome.add_fragment(name=args['name'], sequence=args['sequence'],
                                       circular=args['circular'])
        return FragmentView.to_dict(fragment), 201


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
                    int(q)  # See if q can be converted to an int
                except BaseException:
                    pass
                else:
                    where = where | Q(id=q)
                genomes = Genome.objects.filter(active=True)\
                                        .filter(where)\
                                        .order_by('-id')[s:s + p]
            else:
                genomes = Genome.objects.filter(active=True).order_by('-id')[s:s + p]
        return [GenomeView.to_dict(genome,
                                   compute_length=False,
                                   include_fragments=False,
                                   include_operations=False)
                for genome in genomes]

    @transaction.atomic()
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
        from edge.blastdb import check_and_build_genome_db

        genome = get_genome_or_404(genome_id)
        check_and_build_genome_db(genome)

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
        from edge.blastdb import check_and_build_genome_db

        genome = get_genome_or_404(genome_id)
        check_and_build_genome_db(genome)

        parser = RequestParser()
        parser.add_argument('primers', field_type=list, required=True, location='json')

        args = parser.parse_args(request)
        primers = args['primers']
        if len(primers) != 2:
            raise Exception('Expecting two primers, got %s' % (primers,))

        r = pcr_from_genome(genome, primers[0], primers[1])
        r = (r[0], [b.to_dict() for b in r[1]], [b.to_dict() for b in r[2]], r[3])
        return r, 200


class GenomeOperationViewBase(ViewBase):

    @transaction.atomic()
    def on_post(self, request, genome_id):
        from edge.blastdb import check_and_build_genome_db

        genome = get_genome_or_404(genome_id)
        check_and_build_genome_db(genome)

        # always require a 'create' argument
        parser = RequestParser()
        parser.add_argument('create', field_type=bool, required=True, location='json')
        args = parser.parse_args(request)
        create = args['create']

        errors = []
        parsed = self.parse_arguments(request, errors)
        if parsed is None:
            return dict(errors=" ".join(errors)), 400

        args, op_class = parsed
        if create is False:
            r = op_class.check(genome, **args)
            if r is None:
                return None, 400
            return [x.to_dict() for x in r], 200

        else:
            op = op_class.get_operation(**args)

            # find another a child genome with same operation
            child = None
            for existing_child in genome.children.all():
                if existing_child.operation_set.count() == 1 and\
                   existing_child.operation_set.all()[0].type == op.type and\
                   existing_child.operation_set.all()[0].params == op.params:
                    child = existing_child

            if child is None:
                child = op_class.perform(genome, **args)
                if child:
                    print(f'Generated child genome {child.id} from parent genome {genome_id}')
                    schedule_building_blast_db(child.id)
                    return GenomeView.to_dict(child), 201
            else:  # found existing child, update genome name and set to active
                if "genome_name" in args:
                    child.name = args["genome_name"]
                child.active = True
                child.save()

            if child is None:
                return None, 400
            else:
                return GenomeView.to_dict(child), 200


class GenomeCrisprDSBView(GenomeOperationViewBase):

    def parse_arguments(self, request, errors):
        from edge.crispr import CrisprOp

        parser = RequestParser()
        parser.add_argument('guide', field_type=str, required=True,
                            default=None, location='json')
        parser.add_argument('pam', field_type=str, required=True,
                            default=None, location='json')
        parser.add_argument('genome_name', field_type=str, required=False,
                            default=None, location='json')
        parser.add_argument('notes', field_type=str, required=False,
                            default=None, location='json')

        args = parser.parse_args(request)
        guide = args['guide']
        pam = args['pam']
        genome_name = args['genome_name']
        notes = args['notes']

        return (dict(guide=guide, pam=pam, genome_name=genome_name, notes=notes), CrisprOp)


class GenomeRecombinationView(GenomeOperationViewBase):
    DEFAULT_HA_LENGTH = 30

    @staticmethod
    def validate_annotations(cassette, annotations):
        """
        If user supplied annotations for the donor sequence, to make sure we
        can precisely use the annotation and not leave any ambiguity, we
        require the donor dna to not contain overhangs and backbone
        modifications. This method returns True if that's the case, and
        required fields for annotations exist.
        """

        if re.match(r'^[A-Za-z]+$', cassette) and \
           len([a for a in annotations
                if "base_first" not in a
                or "base_last" not in a
                or "name" not in a
                or "type" not in a
                or "strand" not in a]) == 0:
            return True

        # if there are no supplied annotations, we don't care what format the
        # donor sequence is
        return len(annotations) == 0

    def parse_arguments(self, request, errors):
        from edge.recombine import RecombineOp

        parser = RequestParser()
        parser.add_argument('cassette', field_type=str, required=True, location='json')
        parser.add_argument('homology_arm_length', field_type=int, required=False, location='json')
        parser.add_argument('genome_name', field_type=str, required=False,
                            default=None, location='json')
        parser.add_argument('cassette_name', field_type=str, required=False,
                            default=None, location='json')
        parser.add_argument('notes', field_type=str, required=False,
                            default=None, location='json')
        parser.add_argument('design_primers', field_type=bool, required=False,
                            default=False, location='json')
        parser.add_argument('primer3_opts', field_type=dict, required=False,
                            default=None, location='json')
        parser.add_argument('annotations', field_type=list, required=False,
                            default=None, location='json')

        args = parser.parse_args(request)
        cassette = args['cassette'].strip()
        homology_arm_length = args['homology_arm_length']
        genome_name = args['genome_name']
        cassette_name = args['cassette_name']
        notes = args['notes']
        design_primers = args['design_primers']
        primer3_opts = args['primer3_opts']
        annotations = args['annotations']

        if primer3_opts is None:
            primer3_opts = {}
        if annotations is None:
            annotations = []
        elif GenomeRecombinationView.validate_annotations(cassette, annotations) is False:
            errors.append(
                "Annotations failed validation: \
please make sure donor sequence does not have overhangs \
and annotation array elements have all the required fields."
            )
            return None

        if homology_arm_length is None:
            homology_arm_length = GenomeRecombinationView.DEFAULT_HA_LENGTH

        return (dict(cassette=cassette, homology_arm_length=homology_arm_length,
                     genome_name=genome_name, cassette_name=cassette_name,
                     notes=notes, design_primers=design_primers,
                     primer3_opts=primer3_opts, annotations=annotations),
                RecombineOp)
