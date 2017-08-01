import itertools
import string
import os
from xml.etree import cElementTree as ET


# Generate atom element names
filename = os.path.join(os.getcwd(), 'test.txt')

file_obj  = open("test.txt", 'w')

for i in range(0,3070):
    file_obj.write("C\n")

file_obj.close()

# Original file parser CODE
parser = ET.XMLParser(encoding="utf-8")
root = ET.fromstring(search_mapping_filename, parser=parser)
root = ET.parse(search_mapping_filename).getroot()
