# -*- coding: utf-8 -*-
"""
Created on Mon May 20 23:14:02 2019

@author: Ryan
"""

class wave:
    ''' A class for EM waves '''
    def __init__(self,kvec=None,efield=None,medium=None):
        self.kvec = None if kvec is None else kvec
        self.efield = None if efield is None else efield
        self.medium = None if medium is None else medium
        
    def poynting(self):
        return print('Need to add a Poynting attribute')


class surface:
    ''' A class for planar interfaces '''
    def __init__(self,normal=None,into=None,out=None):
        self.normal = None if normal is None else normal
        self.into = None if into is None else into
        self.out = None if out is None else out
        