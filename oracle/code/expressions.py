"""Finite-field expression DAG used by coordinate and conversion circuits."""


class Expr:
    """Immutable finite-field straight-line program; integer references."""

    def __init__(self, p, ninputs):
        self.p = p
        self.ninputs = ninputs
        self.nodes = []
        self.constants = {}

    def const(self, v):
        v %= self.p
        if v not in self.constants:
            self.constants[v] = self.op("const", v)
        return self.constants[v]

    def op(self, kind, *args):
        i = self.ninputs + len(self.nodes)
        self.nodes.append((kind, args))
        return i

    def add(self, a, b):
        return self.op("add", a, b)

    def sub(self, a, b):
        return self.op("sub", a, b)

    def mul(self, a, b):
        return self.sq(a) if a == b else self.op("mul", a, b)

    def sq(self, a):
        return self.op("square", a)

    def scale(self, a, k):
        return self.op("scale", a, k % self.p)

    def inv(self, a):
        return self.op("inv", a)

    def power(self, a, e):
        if e == 0:
            return self.const(1)
        acc = a
        for bit in bin(e)[3:]:
            acc = self.sq(acc)
            if bit == "1":
                acc = self.mul(acc, a)
        return acc

    def evaluate(self, inputs):
        vals = list(inputs)
        p = self.p
        for op, a in self.nodes:
            if op == "const":
                v = a[0]
            elif op == "add":
                v = vals[a[0]] + vals[a[1]]
            elif op == "sub":
                v = vals[a[0]] - vals[a[1]]
            elif op == "mul":
                v = vals[a[0]] * vals[a[1]]
            elif op == "square":
                v = vals[a[0]] ** 2
            elif op == "scale":
                v = vals[a[0]] * a[1]
            elif op == "inv":
                v = pow(vals[a[0]], -1, p) if vals[a[0]] else 0
            vals.append(v % p)
        return vals
