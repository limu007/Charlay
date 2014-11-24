from scanner import labin
import numpy as nu
self=labin.avantes()
self.setup([1.1,2.5],50,250)
self.adrive(8,4)
dd=nu.array(self.measure())

self.dark=dd
data=[self.measure(prep_task="self.rate(1,1212)") for i in range(5)]
import cPickle
#cPickle.dump(nu.array([self.pixtable]+data),open("oikus.dat","wb"))