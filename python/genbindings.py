#! /usr/bin/env python2

import sys

from pybindgen import FileCodeSink
from pybindgen.gccxmlparser import ModuleParser

def my_module_gen():
    module_parser = ModuleParser('mpc', '::')
    gccxml_options = dict(include_paths=[sys.argv[2]])
    includes=[]
    for f in sys.argv[3:]:
        includes.append('"'+f.replace(sys.argv[2]+"/", '')+'"')
    module_parser.parse(sys.argv[3:], includes=includes, pygen_sink=FileCodeSink(open(sys.argv[1], "w")), gccxml_options=gccxml_options)

if __name__ == '__main__':
    my_module_gen()
