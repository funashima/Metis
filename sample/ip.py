#!/usr/bin/env python3
for i in range(3):
    for j in range(3):
        if i == j:
            continue
        k = 3 - (i + j)
        print('i= {0}, j={1} -> k={2}'.format(i, j, k))
