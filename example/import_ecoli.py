from edge.connector import Connector

c = Connector.open_db('/tmp/world.db')
g = c.create_genome('E. coli MG1655')
with g._edit() as u:
  u.import_gff('./ecoli-mg1656.gff')

