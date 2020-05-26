'''
This script calculates the number of pages needed to download
TROPOMI data given a XML query file (argument 1) and a number
of results per page (argument 2).
'''

import xml.etree.ElementTree as ET
import sys

# Get the file that contains the query and the number of
# results downloaded per page as arguments
query_file=sys.argv[1]
results_per_page=int(sys.argv[2])

# Open the xml file
root = ET.parse(query_file).getroot()

# The information we want is contained in the subtitle tag
for item in root.findall('{http://www.w3.org/2005/Atom}subtitle'):
    query_results = item.text

# Get the total results from that text
total_results = int(query_results.split(' ')[5])

# And then calculate the total number of pages
# needed to download the number of results
total_pages = round(total_results/results_per_page)

# Return the total number of pages
print(total_pages)
