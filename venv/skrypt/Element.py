import numpy as np

class Element:
    def __init__(self, source_1, source_2, source_3, source_4):
        self.ID = []
        self.ID.append(source_1)
        self.ID.append(source_2)
        self.ID.append(source_3)
        self.ID.append(source_4)

    def __repr__(self):
        return "<Element %s>" % (self.ID)

