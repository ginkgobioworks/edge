from edge.connector import Connector

c = Connector.open_db('/tmp/world.db')
g = c.create_genome('Saccharomyces cerevisiae')
with g._edit() as u:
  u.import_gff('./yeast.gff')

