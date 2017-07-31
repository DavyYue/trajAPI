import itertools
import string
import os
from xml.etree import cElementTree as ET

filename = os.path.join(os.getcwd(), 'test.txt')

file_obj  = open("test.txt", 'w')

for i in range(0,3070):
    file_obj.write("C\n")

file_obj.close()
