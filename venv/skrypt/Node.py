class Node:
    def __init__(self, source_x, source_y, source_t, border):
        self.x = source_x
        self.y = source_y
        self.temp = source_t
        self.border = border # border == 1 => jest brzegowy
    def __repr__(self):
        return "<Node x: %s y: %s t: %s border: %s>" % (round(self.x,5), round(self.y,5), self.temp, self.border)