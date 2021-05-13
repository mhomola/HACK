from common_constants import Constants

class Example(Constants):
    def define_altitude(self, h):
        self.ISA_calculator(h_input=h)
    def compute_something(self):
        return self.p_0 / self.p



x = Example()
x.define_altitude(h=20000)
r = x.compute_something()
print(x.p_0, '/', x.p, ' = ', r)