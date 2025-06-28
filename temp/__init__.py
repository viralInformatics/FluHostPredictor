# !/usr/bin/python 
# -*- coding: utf-8 -*-
# Author:lihuiru
# Created on 2023/11/15 11:16
import time
import functools


def stat_time(func):
    @functools.wraps(func)
    def wrapper():
        start_time = time.time()
        func()
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time: {elapsed_time} seconds")

    return wrapper