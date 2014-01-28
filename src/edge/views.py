from flask import Flask, g
from flask.ext.restful import reqparse, abort, Api, Resource

#
# Setup App and API
#

app = Flask(__name__)
api = Api(app)

app.config.update(dict(
    DATABASE='/tmp/world.db',
    DEBUG=True,
    USERNAME='',
    PASSWORD='',
))


#
# Handle DB connection
#

def connect_db():
    from edge.connector import Connector
    connector = getattr(g, '_connector', None)
    if connector is None:
        connector = g._connector = Connector.open_db(app.config['DATABASE'])
    return connector


@app.teardown_appcontext
def close_connection(exception):
    connector = getattr(g, '_connector', None)
    if connector is not None:
        connector.close()


def init_db():
    from edge.connector import Connector
    Connector.create_db(app.config['DATABASE'])


#
# Actual App/API code
#

def get_genome_or_404(genome_id):
    from edge.genome import Genome, GenomeNotFound
    connector = connect_db()
    try:
        return Genome(connector, genome_id)
    except GenomeNotFound as e:
        abort(404, message=str(e))


def get_fragment_or_404(fragment_id):
    from edge.fragment import Fragment_Operator, FragmentNotFound
    connector = connect_db()
    try:
        return Fragment_Operator(connector, fragment_id)
    except FragmentNotFound as e:
        abort(404, message=str(e))


fragment_parser = reqparse.RequestParser()
fragment_parser.add_argument('name', type=str, required=True, location='json')
fragment_parser.add_argument('sequence', type=str, required=True, location='json')
fragment_parser.add_argument('circular', type=bool, default=False, location='json')


class FragmentResource(Resource):

    @staticmethod
    def to_dict(fragment):
        return dict(id=fragment.id,
                    uri=api.url_for(FragmentResource, fragment_id=fragment.id),
                    name=fragment.name,
                    circular=fragment.circular,
                    length=fragment.length)

    def get(self, fragment_id):
        fragment = get_fragment_or_404(fragment_id)
        return FragmentResource.to_dict(fragment)


class FragmentSequenceResource(Resource):

    def get(self, fragment_id):
        q_parser = reqparse.RequestParser()
        q_parser.add_argument('f', type=int, location='args')
        q_parser.add_argument('l', type=int, location='args')
        args = q_parser.parse_args()
        f = args['f'] if 'f' in args and args['f'] is not None else None
        l = args['l'] if 'l' in args and args['l'] is not None else None

        fragment = get_fragment_or_404(fragment_id)
        s = fragment.get_sequence(bp_lo=f, bp_hi=l)
        if f is None:
            f = 1
        if l is None:
            l = f+len(s)-1
        return {'sequence': s, 'first_bp': f, 'last_bp': l}


class FragmentAnnotationsResource(Resource):

    @staticmethod
    def to_dict(annotation):
        return dict(first_bp=annotation.first_bp,
                    last_bp=annotation.last_bp,
                    name=annotation.name, type=annotation.type, strand=annotation.strand,
                    annotation_full_length=annotation.annotation_full_length,
                    annotation_first_bp=annotation.annotation_first_bp,
                    annotation_last_bp=annotation.annotation_last_bp)

    def get(self, fragment_id):
        q_parser = reqparse.RequestParser()
        q_parser.add_argument('f', type=int, location='args')
        q_parser.add_argument('l', type=int, location='args')
        args = q_parser.parse_args()
        f = args['f'] if 'f' in args and args['f'] is not None else None
        l = args['l'] if 'l' in args and args['l'] is not None else None

        fragment = get_fragment_or_404(fragment_id)
        return [FragmentAnnotationsResource.to_dict(annotation)
                for annotation in fragment.annotations(bp_lo=f, bp_hi=l)]

    def post(self, fragment_id):
        annotation_parser = reqparse.RequestParser()
        annotation_parser.add_argument('first_bp', type=int, required=True, location='json')
        annotation_parser.add_argument('last_bp', type=int, required=True, location='json')
        annotation_parser.add_argument('name', type=str, required=True, location='json')
        annotation_parser.add_argument('type', type=str, required=True, location='json')
        annotation_parser.add_argument('strand', type=int, required=True, location='json')

        args = annotation_parser.parse_args()
        fragment = get_fragment_or_404(fragment_id)
        with fragment.annotate() as u:
            u.annotate(first_base1=args['first_bp'],
                       last_base1=args['last_bp'],
                       annotation_name=args['name'],
                       annotation_type=args['type'],
                       strand=args['strand'])
        return {}, 201


class FragmentListResource(Resource):

    def get(self):
        connector = connect_db()
        fragments = connector.non_genomic_fragments()
        return [FragmentResource.to_dict(fragment) for fragment in fragments]

    def post(self):
        args = fragment_parser.parse_args()
        connector = connect_db()
        fragment = connector.create_fragment_with_sequence(name=args['name'],
                                                           sequence=args['sequence'],
                                                           circular=args['circular'])
        return FragmentResource.to_dict(fragment), 201


class GenomeResource(Resource):

    @staticmethod
    def to_dict(genome):
        changes = genome.changed_locations_by_fragment()
        changes = {f.id: v for f, v in changes.iteritems()}
        fragments = []
        for f in genome.fragments():
            d = FragmentResource.to_dict(f)
            if f.id in changes:
                d['changes'] = changes[f.id]
            fragments.append(d)

        return dict(id=genome.id,
                    uri=api.url_for(GenomeResource, genome_id=genome.id),
                    name=genome.name,
                    notes=genome.notes,
                    parent_id=genome.parent_id,
                    fragments=fragments)

    def get(self, genome_id):
        genome = get_genome_or_404(genome_id)
        return GenomeResource.to_dict(genome)


class GenomeAnnotationsResource(Resource):

    def get(self, genome_id):
        genome = get_genome_or_404(genome_id)
        q_parser = reqparse.RequestParser()
        q_parser.add_argument('q', type=str, required=True)
        args = q_parser.parse_args()

        res = []
        fragment_annotations = genome.find_annotation(args['q'])
        for fragment_id in fragment_annotations:
            fragment = get_fragment_or_404(fragment_id)
            annotations = fragment_annotations[fragment_id]
            d = FragmentResource.to_dict(fragment)
            a = [FragmentAnnotationsResource.to_dict(x) for x in annotations]
            res.append((d, a))

        return res


class GenomeFragmentListResource(Resource):

    def post(self, genome_id):  # adding new fragment
        args = fragment_parser.parse_args()
        genome = get_genome_or_404(genome_id)
        fragment = None
        with genome._edit() as u:
            fragment = \
                u.add_fragment(name=args['name'], sequence=args['sequence'],
                               circular=args['circular'], annotate=False)
        return FragmentResource.to_dict(fragment), 201


class GenomeFragmentResource(Resource):

    def insert_bases(self, genome_id, fragment_id):
        op_parser = reqparse.RequestParser()
        op_parser.add_argument('name', type=str, required=True, location='json')
        op_parser.add_argument('before_bp', type=int, required=True, location='json')
        op_parser.add_argument('sequence', type=str, required=True, location='json')
        args = op_parser.parse_args()

        genome = get_genome_or_404(genome_id)
        fragment = get_fragment_or_404(fragment_id)

        with genome.update(name=args['name']) as u:
            with u.update_fragment_by_fragment_id(fragment.id) as f:
                f.insert_bases(args['before_bp'], args['sequence'])
        g2 = genome.last_updated()
        return GenomeResource.to_dict(g2), 201

    def remove_bases(self, genome_id, fragment_id):
        op_parser = reqparse.RequestParser()
        op_parser.add_argument('name', type=str, required=True, location='json')
        op_parser.add_argument('before_bp', type=int, required=True, location='json')
        op_parser.add_argument('length', type=int, required=True, location='json')
        args = op_parser.parse_args()

        genome = get_genome_or_404(genome_id)
        fragment = get_fragment_or_404(fragment_id)

        with genome.update(name=args['name']) as u:
            with u.update_fragment_by_fragment_id(fragment.id) as f:
                f.remove_bases(args['before_bp'], args['length'])
        g2 = genome.last_updated()
        return GenomeResource.to_dict(g2), 201

    def replace_bases(self, genome_id, fragment_id):
        op_parser = reqparse.RequestParser()
        op_parser.add_argument('name', type=str, required=True, location='json')
        op_parser.add_argument('before_bp', type=int, required=True, location='json')
        op_parser.add_argument('length', type=int, required=True, location='json')
        op_parser.add_argument('sequence', type=str, required=True, location='json')
        args = op_parser.parse_args()

        genome = get_genome_or_404(genome_id)
        fragment = get_fragment_or_404(fragment_id)

        with genome.update(name=args['name']) as u:
            with u.update_fragment_by_fragment_id(fragment.id) as f:
                f.replace_bases(args['before_bp'], args['length'], args['sequence'])
        g2 = genome.last_updated()
        return GenomeResource.to_dict(g2), 201

    def insert_fragment(self, genome_id, fragment_id):
        op_parser = reqparse.RequestParser()
        op_parser.add_argument('name', type=str, required=True, location='json')
        op_parser.add_argument('before_bp', type=int, required=True, location='json')
        op_parser.add_argument('fragment_id', type=int, required=True, location='json')
        args = op_parser.parse_args()

        genome = get_genome_or_404(genome_id)
        fragment = get_fragment_or_404(fragment_id)
        new_fragment = get_fragment_or_404(args['fragment_id'])

        with genome.update(name=args['name']) as u:
            with u.update_fragment_by_fragment_id(fragment.id) as f:
                f.insert_fragment(args['before_bp'], new_fragment)
        g2 = genome.last_updated()
        return GenomeResource.to_dict(g2), 201

    def replace_with_fragment(self, genome_id, fragment_id):
        op_parser = reqparse.RequestParser()
        op_parser.add_argument('name', type=str, required=True, location='json')
        op_parser.add_argument('before_bp', type=int, required=True, location='json')
        op_parser.add_argument('length', type=int, required=True, location='json')
        op_parser.add_argument('fragment_id', type=int, required=True, location='json')
        args = op_parser.parse_args()

        genome = get_genome_or_404(genome_id)
        fragment = get_fragment_or_404(fragment_id)
        new_fragment = get_fragment_or_404(args['fragment_id'])

        with genome.update(name=args['name']) as u:
            with u.update_fragment_by_fragment_id(fragment.id) as f:
                f.replace_with_fragment(args['before_bp'], args['length'], new_fragment)
        g2 = genome.last_updated()
        return GenomeResource.to_dict(g2), 201

    def put(self, genome_id, fragment_id):  # updating genome and fragment
        op_parser = reqparse.RequestParser()
        op_parser.add_argument('op', type=str, required=True, location='json')
        args = op_parser.parse_args()

        if args['op'] == 'insert_bases':
            return self.insert_bases(genome_id, fragment_id)
        elif args['op'] == 'remove_bases':
            return self.remove_bases(genome_id, fragment_id)
        elif args['op'] == 'replace_bases':
            return self.replace_bases(genome_id, fragment_id)
        elif args['op'] == 'insert_fragment':
            return self.insert_fragment(genome_id, fragment_id)
        elif args['op'] == 'replace_with_fragment':
            return self.replace_with_fragment(genome_id, fragment_id)


class GenomeListResource(Resource):

    def get(self):
        connector = connect_db()
        genomes = connector.genomes()
        return [GenomeResource.to_dict(genome) for genome in genomes]

    def post(self):
        genome_parser = reqparse.RequestParser()
        genome_parser.add_argument('name', type=str, required=True, location='json')
        genome_parser.add_argument('notes', type=str, location='json')

        args = genome_parser.parse_args()
        connector = connect_db()
        genome = connector.create_genome(name=args['name'], notes=args['notes'])
        return GenomeResource.to_dict(genome), 201


#
# Setup routing
#

api.add_resource(FragmentListResource, '/fragments')
api.add_resource(FragmentResource, '/fragments/<int:fragment_id>')
api.add_resource(FragmentSequenceResource, '/fragments/<int:fragment_id>/sequence')
api.add_resource(FragmentAnnotationsResource, '/fragments/<int:fragment_id>/annotations')
api.add_resource(GenomeListResource, '/genomes')
api.add_resource(GenomeResource, '/genomes/<int:genome_id>')
api.add_resource(GenomeAnnotationsResource, '/genomes/<int:genome_id>/annotations')
api.add_resource(GenomeFragmentListResource, '/genomes/<int:genome_id>/fragments')
api.add_resource(GenomeFragmentResource, '/genomes/<int:genome_id>/fragments/<int:fragment_id>')

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
