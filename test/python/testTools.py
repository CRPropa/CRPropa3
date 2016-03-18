from crpropa import *

class CountingModule(Module):

	def __init__(self):
		Module.__init__(self)
		self.count = 0

	def process(self, c):
		self.count += 1

f = ParticleFilter([-1,1])
a = CountingModule()
r = CountingModule()
f.onAccept(a)
f.onReject(r)

c = Candidate()

for id in [-1 , 1, 6, 9, -19, 23, 100010001]:
  c.current.setId(id)
  f.process(c)
print (a.count, r.count)
