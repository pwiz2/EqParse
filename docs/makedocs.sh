#!/bin/bash
#eval "ls"
eval "make html"
eval "make pdf"

EPIDOCCMD="epydoc --html ../src/eqparse -o eqparse_epydoc --graph all --name EqParse --inheritance grouped"
eval $(EPIDOCCMD)
