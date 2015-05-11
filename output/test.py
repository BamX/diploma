#!/usr/local/bin/python
import sys

test = open('test_view.csv', 'r')
view = open('view.csv', 'r')

isOk = True
for t, v in zip([l.split(',') for l in test], [l.split(',') for l in view]):
        if t[0] != v[0]:
            print "Wrong time grid %s vs %s" % (t[0], v[0])
        deltas = [abs(float(tView) - float(vView)) for (tView, vView) in zip(t[1:], v[1:])]
        if len(filter(lambda delta: delta > 0.011, deltas)) > 0:
            print t[0], deltas
            isOk = False

if isOk:
    print "OK"
