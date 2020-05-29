#!/usr/bin/env python

# Authors:Xiaohong Zhuang, Yalun Yu, Dr.Jeffery B. Klauda
# Date:9/27/2016
# Calculate the number of clusters and number of lipids in each cluster


import numpy as np


def number_corr(data, portion=0.9):
    length = data.shape[0]
    # the step here is unit-less
    cnt = [1.0]
    for delta_t in range(1, int(portion * length)):
        print(delta_t)
        num_pairs = length - delta_t
        production = data[:-delta_t] * data[delta_t:]
        sum_prod = production.sum()
        average = sum_prod / num_pairs
        cnt.append(average)
    return cnt


a = number_corr(np.array([1, 0, 1, 0, 1, 0, 1, 0, 1, 0]), portion=0.9)
print(a)