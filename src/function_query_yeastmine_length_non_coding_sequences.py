#!/usr/bin/env python

# This is an automatically generated script to run your query
# to use it you will require the intermine python client.
# To install the client, run the following command from a terminal:
#
#     sudo easy_install intermine
#
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/

# The line below will be needed if you are running this script with python 2.
# Python 3 will ignore it.
from __future__ import print_function

# The following two lines will be needed in every python script:

from intermine.webservice import Service
service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("IntergenicRegion")

# Type constraints should come early - before all mentions of the paths they constrain
query.add_constraint("chromosomeLocation.locatedOn", "Chromosome")

# The view specifies the output columns
query.add_view(
    "primaryIdentifier", "length",
    "chromosomeLocation.locatedOn.primaryIdentifier"
)

# Uncomment and edit the line below (the default) to select a custom sort order:
# query.add_sort_order("IntergenicRegion.primaryIdentifier", "ASC")

# You can edit the constraint values below
query.add_constraint("organism.shortName", "=", "S. cerevisiae", code="A")
query.add_constraint("adjacentGenes.featureType", "!=", "not in systematic sequence of S288C", code="D")
query.add_constraint("chromosomeLocation.locatedOn.primaryIdentifier", "=", "chrI", code="B")
query.add_constraint("length", ">=", "100", code="C")

# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A and B and C and D")

# Outer Joins
# (display properties of these relations if they exist,
# but also show objects without these relationships)
query.outerjoin("adjacentGenes")

for row in query.rows():
    print(row["primaryIdentifier"], row["length"], \
        row["chromosomeLocation.locatedOn.primaryIdentifier"])

