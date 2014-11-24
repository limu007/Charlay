--- extra.py	(original)
+++ extra.py	(refactored)
@@ -69,20 +69,20 @@
         vstd=sata.std()*nsig
     if len(sata)/10<ndiv: 
         ndiv=max([int(len(sata)/10),10])
-        print 'reducing no. of bins for histog. to %i'%ndiv
+        print('reducing no. of bins for histog. to %i'%ndiv)
     vstp=2*vstd/ndiv
     vmin=vmax-vstd
     vmax=vmax+vstd
     c=arange(vmin,vmax,vstp)
-    if loud: print 'histog. from %4.2f to %4.2f'%(vmin,vmax)
+    if loud: print('histog. from %4.2f to %4.2f'%(vmin,vmax))
     a,b=histogram(data,c)
     if loud>0:
-        if loud>2: print b,a
-        else: print 'histog. sum %i [%i above]'%(sum(a),sum(a>min_count))
+        if loud>2: print(b,a)
+        else: print('histog. sum %i [%i above]'%(sum(a),sum(a>min_count)))
     if sum(a>min_count)<10: bord=2
     else: bord=3
     if sum(a>min_count)<min_bins:
-        print 'too few bins to fit a gausian'
+        print('too few bins to fit a gausian')
         if rep==2: return b[argmax(a[bord:-bord])+bord],(vmax-vmin)/2
         else: return vstd
     #if rep>3: return get_width(b[a>min_count],a[a>min_count],bord=bord,rep=rep),b[a>min_count],a[a>min_count]
@@ -97,7 +97,7 @@
     for i in range(nclip):
         avg=val[sel].mean()
         std=val[sel].std()
-        if loud>0: print "selecting %f +- %f"%(avg,std*nsig)
+        if loud>0: print("selecting %f +- %f"%(avg,std*nsig))
         sel=abs(val-avg)<std*nsig
         if sum(sel)==len(sel): break # no sense of repeating
     if rep==1: return val[sel].min(),val[sel].max()
@@ -109,7 +109,7 @@
     dlen=len(data)
     idi=arange(dlen)[sele]
     if len(idi)==0: 
-        if loud>0:print 'no good bins'
+        if loud>0:print('no good bins')
         return []
     neidiff=idi[1:]-idi[:-1]
     #gstep=pidi[1:]-pidi[:-1]
@@ -120,8 +120,8 @@
     gbegs.insert(0,idi[0])
     gends=list(idi[gidi]+1)
     if gidi[-1]!=len(idi)-1: gends.append(idi[-1]+1)
-    gtis=zip(gbegs,gends)
-    if loud>0: print '%i good bins, %i intervals'%(sum(sele),len(gtis))
+    gtis=list(zip(gbegs,gends))
+    if loud>0: print('%i good bins, %i intervals'%(sum(sele),len(gtis)))
     return [a for a in gtis if a[1]-a[0]>min_wid] #good intervals longer than ...
 
 def mask_data(data,bias,sigma=None):
@@ -191,22 +191,22 @@
             step=bin*median(xbin[1:]-xbin[:-1])
         labs=((xbin-xbin[0])/step).astype(int32)
         if sum(labs<0)>0:
-            print 'x-axis not rising'
+            print('x-axis not rising')
             return [],[],[]
         bin=int(len(xbin)/labs[-1])
-        if loud: print 'using step %.3f, grouping in average by %i bins'%(step,bin) 
+        if loud: print('using step %.3f, grouping in average by %i bins'%(step,bin)) 
     else: labs=(arange(len(data))/bin).astype(int32)
     if sigma!=None: labs[sele==False]=-1
     cnts=zeros((max(labs)+1,))
-    if loud>1: print "labels [%i-%i], length %i"%(min(labs),max(labs),len(cnts))
+    if loud>1: print("labels [%i-%i], length %i"%(min(labs),max(labs),len(cnts)))
     for l in labs[sele]:
         cnts[l]+=1
     idx=arange(len(cnts))[cnts>=bin*frac_min]
     if len(idx)<1:
-        print "all bins empty"
+        print("all bins empty")
         return None
     else:
-        if loud>0: print "%i bin(s) empty"%(len(cnts)-len(idx))
+        if loud>0: print("%i bin(s) empty"%(len(cnts)-len(idx)))
     #print 'max. index %i'%(labs[-1])
     avgs=ndimage.mean(data,labs,idx)
     if errs!=None: 
@@ -240,7 +240,7 @@
     iend=0
     if xbin!=None: xmark=xbin[0]
     else: 
-        xbin=range(len(data))
+        xbin=list(range(len(data)))
         xmark=0
     for i in range(ndiv):
         if step!=None:
@@ -300,7 +300,7 @@
     global chi2
     assert type(r)==int
     if len(y)<r+1:
-        print 'too few data to fit a polynom'
+        print('too few data to fit a polynom')
         if full==0: return None    
         else: return None,None,None,None
     arrmat=ones((size(y),r+1),float64)
@@ -317,7 +317,7 @@
         resid/=resid.mean()
         if lim>0:
             sel=resid<lim
-            print "using %i points [%.4f]"%(sum(sel),sum(resid[sel]**2))
+            print("using %i points [%.4f]"%(sum(sel),sum(resid[sel]**2)))
             return polyfit(x[sel],y[sel],r,w,clip=clip-1,lim=lim,full=full)
         if w==None: w=ones(x.shape)
         return polyfit(x,y,r,w/resid,clip=clip-1,full=full)
@@ -378,9 +378,9 @@
     old = n.array( a.shape )
     ndims = len( a.shape )
     if len( newdims ) != ndims:
-        print "[congrid] dimensions error. " \
+        print("[congrid] dimensions error. " \
             "This routine currently only support " \
-            "rebinning to the same number of dimensions."
+            "rebinning to the same number of dimensions.")
         return None
     newdims = n.asarray( newdims, dtype=float )
     dimlist = []
@@ -407,7 +407,7 @@
         mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
         newa = mint( dimlist[-1] )
 
-        trorder = [ndims - 1] + range( ndims - 1 )
+        trorder = [ndims - 1] + list(range( ndims - 1))
         for i in range( ndims - 2, -1, -1 ):
             newa = newa.transpose( trorder )
 
@@ -425,7 +425,7 @@
         nslices = [ slice(0,j) for j in list(newdims) ]
         newcoords = n.mgrid[nslices]
 
-        newcoords_dims = range(n.rank(newcoords))
+        newcoords_dims = list(range(n.rank(newcoords)))
         #make first index last
         newcoords_dims.append(newcoords_dims.pop(0))
         newcoords_tr = newcoords.transpose(newcoords_dims)
@@ -441,9 +441,9 @@
         newa = scipy.ndimage.map_coordinates(a, newcoords)
         return newa
     else:
-        print "Congrid error: Unrecognized interpolation type.\n", \
+        print("Congrid error: Unrecognized interpolation type.\n", \
             "Currently only \'neighbour\', \'nearest\',\'linear\',", \
-            "and \'spline\' are supported."
+            "and \'spline\' are supported.")
         return None
     
 
@@ -475,7 +475,7 @@
         if type(errs)==list: errs=errs[sel]
         #xbase=[data[xids[0]]-errs*fact,data[xids[1]]+errs*fact]
     step=(xbase[1]-xbase[0])/ndiv
-    print 'range %5.2f-%5.2f step:%5.2f'%(xbase[0],xbase[1],step)
+    print('range %5.2f-%5.2f step:%5.2f'%(xbase[0],xbase[1],step))
     if iterable(errs)>0:
         #errm=(max(errs)+min(errs))/2.
         errm=errs[errs>0].mean()
@@ -523,7 +523,7 @@
             try:
                 vals[pos:pos+alen]+=kern[nos:nos+alen]
             except:
-                print "casting %i:+%i to %i:+%i (%i / %i)"%(nos,alen,pos,alen,len(vals),len(kern))
+                print("casting %i:+%i to %i:+%i (%i / %i)"%(nos,alen,pos,alen,len(vals),len(kern)))
     return vals,axis
 
 def fill_k2de(xdata,xerrs,ydata,yerrs,fact=2.,form='gauss',xdiv=100,xbase=None,ydiv=100,ybase=None,xlog=False,ylog=False):
@@ -546,7 +546,7 @@
         ydiv=(ybase[1]-ybase[0])//ystep
         ybase[1]=ybase[0]+ystep*ydiv
     else: ystep=(ybase[1]-ybase[0])/ydiv
-    print 'step:%5.2f x %5.2f'%(xstep,ystep)
+    print('step:%5.2f x %5.2f'%(xstep,ystep))
     if iterable(xerrs)>0:
         #xerrm=(max(xerrs)+min(xerrs))/2.
         xerrm=xerrs[xerrs>0].mean()
@@ -593,7 +593,7 @@
     elif xerrm==0: xvals=exp(-xern**2/2/xerrs**2)
     elif yerrm==0: yvals=exp(-yern**2/2/yerrs**2)
     else:
-        print 'mean errors:%5.2f x %5.2f'%(xerrm,yerrm)
+        print('mean errors:%5.2f x %5.2f'%(xerrm,yerrm))
     xaxis=arange(xbase[0],xbase[1],xstep)
     yaxis=arange(ybase[0],ybase[1],ystep)
     vals=zeros((xdiv,ydiv),dtype=float64)
@@ -648,10 +648,10 @@
     smin=arep.argmin()
     spos=[smin//len(grid[0]),smin%len(grid[0])]
     amin=arep[spos[0],spos[1]]
-    print 'min. value %.2f at %i,%i : %.3f,%.3f'%(amin,spos[0],spos[1],grid[0][spos[1]],grid[1][spos[0]])
+    print('min. value %.2f at %i,%i : %.3f,%.3f'%(amin,spos[0],spos[1],grid[0][spos[1]],grid[1][spos[0]]))
     sbord=((arep[0,0]+arep[-1,-1])/2.-amin)*shallow
     if sum(arep-amin>sbord)<arep.size*bord_frac:
-        print 'too flat map'
+        print('too flat map')
         return arep
     if blofit=='origauss':
         from numpy import dot
@@ -720,13 +720,13 @@
         if ybnd[0]+grid[1][0]<yaxis[0] or ybnd[1]+grid[1][-1]>yaxis[-1]: #need to limit the range
             sele=(yval+grid[1][0]>yaxis[0])*(yval+grid[1][-1]<=yaxis[-1])
             if sum(sele)<2:
-                print 'for slope %.3f: out of range'%(pars[0]+dx)
+                print('for slope %.3f: out of range'%(pars[0]+dx))
                 cnt.append(0)
                 rep.append([0.]*len(grid[1]))
                 continue
             yval=yval[sele]
             if loud>1:
-                print 'for slope %.3f: going from %.3f to %.3f: %i'%(pars[0]+dx,yval[0],yval[-1],sum(sele))
+                print('for slope %.3f: going from %.3f to %.3f: %i'%(pars[0]+dx,yval[0],yval[-1],sum(sele)))
             rep.append([kmap[xpos[sele],((yval-yaxis[0]+dy)/ystep).astype(int)].sum() for dy in grid[1]])
             cnt.append(len(yval))
         else:
@@ -802,7 +802,7 @@
     if ylog: selv*=yarr>0
     n=sum(selv)
     if n<=2:
-        print 'not enough data'
+        print('not enough data')
         return
     xcor=xarr[selv]
     ycor=yarr[selv]
@@ -853,7 +853,7 @@
             wei/=wei.sum()
             arrmat[:,0]*=wei
             ycor2=ycor*wei
-            if loud>2: print "weight:%s"%wei
+            if loud>2: print("weight:%s"%wei)
         else: ycor2=ycor
         for i in range(r):
             arrmat[:,i+1]=arrmat[:,i]*xcor
@@ -866,7 +866,7 @@
     if dsum>0:
         eslop=dsum/sqrt(xxcor)
         ebias=dsum*sqrt(1./n+xavg**2/xxcor)
-    print 'correlation (%5.2f) : %5.2f(+-%5.2f) *x + %5.2f(+-%5.2f)'%(dsum,slop,eslop,bias,ebias)
+    print('correlation (%5.2f) : %5.2f(+-%5.2f) *x + %5.2f(+-%5.2f)'%(dsum,slop,eslop,bias,ebias))
     
     if grid==None:
         if loud==2: return corr_dist(xcor,ycor,slop,bias,perp,xerr,yerr)
