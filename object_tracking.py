#!/usr/local/sci/bin/python2.7

import numpy as np
from scipy import interpolate
import scipy.ndimage as ndimage
import datetime
from os.path import isfile, isdir
import matplotlib.pyplot as plt
import warnings

###################################################################
# Initiate StormS class with object properties. Can be adjusted to store additional object properties.
# Future release could relate StormS and write_storms functions for easier management and user adaptation.
###################################################################

class StormS():
    def __init__(self, jj,StormLabels,var, xmat,ymat,newwas, newumat, newvmat, num_dt, misval, doradar, under_threshold, extra_thresh=[], storm_history=False, string=None, rarray=[],azarray=[]):
        if string == None: # initialiase to default values
            C=np.where(StormLabels==jj)
            self.storm=int(jj)
            self.area=int(np.size(C,1))
            if under_threshold:
                self.extreme=np.min(var[C])
            else:
                self.extreme=np.max(var[C])
            self.meanvar=np.mean(var[C])
            if len(extra_thresh) > 0:
                self.extra_area=[(var[C] < a).sum() for a in extra_thresh]
            self.centroidx=np.mean(xmat[C])
            self.centroidy=np.mean(ymat[C])
            self.centroidx0=misval
            self.centroidy0=misval
            self.boxleft = np.min(xmat[C])
            self.boxup = np.max(ymat[C])
            self.boxwidth = np.max(xmat[C])-np.min(xmat[C])
            self.boxheight = np.max(ymat[C])-np.min(ymat[C])
            self.life=1
            if storm_history:
                self.was=int(jj)
                self.dx=np.mean(newumat[C])/num_dt
                self.dy=np.mean(newvmat[C])/num_dt
            else: # First image to be considered, so no dx or dy or previous label
                self.was=newwas
                self.dx=0
                self.dy=0
            self.parent=[misval]
            self.child=misval
            self.wasdist=misval
            self.accreted=[misval]
            if doradar:
                self.rangel=np.min(rarray[C])
                self.rangeu=np.max(rarray[C])
                if np.min(rarray[C])==0:
                    self.azimuthl=0.
                    self.azimuthu=360.
                elif (np.min(azarray[C])==0) & (np.max(azarray[C])>180):
                    azxy = np.where((xmat==np.round(self.centroidx)) & (ymat==np.round(self.centroidy)))
                    azoffset = np.fmod(np.round(azarray[azxy])+180.,360.)
                    aznotind = np.where((StormLabels==jj) & (azarray<azoffset))
                    if np.size(aznotind)==0:
                        self.azimuthl=0.
                        self.azimuthu=360.
                    elif np.max(azarray[aznotind])>azoffset-1.:
                        self.azimuthl=0.
                        self.azimuthu=360.
                    else:
                        azleftind=np.where((StormLabels==jj) & (azarray>=azoffset))
                        if np.size(azleftind)==0:
                            self.azimuthl = np.min(azarray[aznotind])
                            self.azimuthu = np.max(azarray[aznotind])
                        else:
                            self.azimuthl = np.min(azarray[azleftind])
                            self.azimuthu = np.max(azarray[aznotind])
                else:
                    self.azimuthl = np.min(azarray[C])
                    self.azimuthu = np.max(azarray[C])
        else: # Use string which is line in save file to initialise
            self.storm=int(jj)
            self.area=int([d for d in string.split() if d.startswith('area=')][0].replace('area=', ''))
            self.extreme=float([d for d in string.split() if d.startswith('extreme=')][0].replace('extreme=', ''))
            self.meanvar=float([d for d in string.split() if d.startswith('meanv=')][0].replace('meanv=', ''))
            if len(extra_thresh) > 0:
                self.extra_area=[int([d for d in string.split() if d.startswith('area<'+str(e))][0].split('=')[-1]) for e in extra_thresh]
            self.centroidx=float([d for d in string.split() if d.startswith('centroid=')][0].replace('centroid=', '').split(',')[0])
            self.centroidy=float([d for d in string.split() if d.startswith('centroid=')][0].replace('centroid=', '').split(',')[1])
            self.centroidx0=float([d for d in string.split() if d.startswith('centroid0=')][0].replace('centroid0=', '').split(',')[0])
            self.centroidy0=float([d for d in string.split() if d.startswith('centroid0=')][0].replace('centroid0=', '').split(',')[1])
            self.life=int([d for d in string.split() if d.startswith('life=')][0].replace('life=', ''))
            self.was=int(string.split()[1])
            self.dx=float([d for d in string.split() if d.startswith('dx=')][0].replace('dx=', ''))
            self.dy=float([d for d in string.split() if d.startswith('dy=')][0].replace('dy=', ''))
            self.parent=[int(p) for p in [d for d in string.split() if d.startswith('parent=')][0].replace('parent=', '').split(',')]
            self.child=[int(p) for p in [d for d in string.split() if d.startswith('child=')][0].replace('child=', '').split(',')]
            self.accreted=[int(p) for p in [d for d in string.split() if d.startswith('accreted=')][0].replace('accreted=', '').split(',')]
            box = [d for d in string.split() if d.startswith('box=')][0].replace('box=', '').split(',')
            self.boxleft=float(box[0])
            self.boxup=float(box[1])
            self.boxheight=float(box[2])
            self.boxwidth=float(box[3])
            if doradar:
                self.range=float([d for d in string.split() if d.startswith('range=')][0].replace('range=', ''))
                self.azimuthl=float([d for d in string.split() if d.startswith('azimuth=')][0].replace('azimuth=', '')[0])
                self.azimuthu=([d for d in string.split() if d.startswith('azimuth=')][0].replace('azimuth=', '')[1])


    def inherit_properties(self, jj, OldStormData, kindex, QuvL, StormLabels, qhist, lapthresh, misval, single_overlap=False):
        self.was = OldStormData[kindex].was
        self.life=OldStormData[kindex].life+1
        self.wasdist=np.size(np.where((QuvL==kindex+1) & (StormLabels==jj)),1)
        self.centroidx0 = OldStormData[kindex].centroidx0
        self.centroidy0 = OldStormData[kindex].centroidy0
        qind = kindex+1
        # Code below is only required when multiple clouds overlap
        if not(single_overlap):
            alllaps=np.where(qhist[1:] >= lapthresh)
            for kkind in range(np.size(alllaps,1)):
                allindex = np.squeeze(alllaps[0][kkind])
                if allindex==kindex:
                   continue
                if self.accreted[-1]==misval:
                   self.accreted[-1]=OldStormData[allindex].was
                else:
                   self.accreted.append(OldStormData[allindex].was)

###################################################################
# TRACKING ALGORITHM
# 1. Correlate previous and current time step to find (dx,dy) displacements.
# 2. Propagate features from previous time step to current time step using (dx,dy) displacements. 
# 3. Iterate through objects to check for overlap and inherit object properties.
# 4. Iterate through objects to check for splitting and merging events.
###################################################################

def track_storms(OldStormData, var, newwas, StormLabels, OldStormLabels, xmat, ymat, fftpixels, dd_tolerance, halosq, squarehalf, oldbt, newbt, num_dt, lapthresh, misval, doradar, under_threshold, IMAGES_DIR, write_file_ID, flagplot, rarray=[],azarray=[]):
    
    ###############################################################
    # PARAMETERS FOR FUTURE FUNCTIONALITY
    ###################################################################
    tukey_window=1
    extra_thresh=[]
    
    ###############################################################
    # START TRACKING!!
    ###################################################################
    
    newumat=0
    newvmat=0
    wasarray=0*StormLabels # set up array of zeros.
    lifearray=0*StormLabels
    numstorms = StormLabels.max()
    print('numstorms = ', numstorms)
    StormData = []

    if len(OldStormData) == 0:
       waslabels=[]
       for ns in range(numstorms):
           jj = ns+1 # First storm is labelled 1, but python indeces start at 0.
           C=np.where(StormLabels==jj)
           StormData += [StormS(jj,StormLabels,var, xmat,ymat,newwas, 0,0, num_dt, misval, doradar, under_threshold, extra_thresh=extra_thresh, storm_history=False, string=None, rarray=rarray, azarray=azarray)]
           wasarray[C]=newwas
           newwas=newwas+1
           lifearray[C]=1
           waslabels.append(StormData[ns].was)
           StormData[ns].centroidx0 = StormData[ns].centroidx+0.0
           StormData[ns].centroidy0 = StormData[ns].centroidy+0.0
    ###################################################
    # OldStormData & StormData ARE NOT EMPTY, SO USE FFT TO GET VELOCITIES
    # AND UPDATE UVLABEL IN OldStormData ACCORDINGLY
    # Estimate velocities using squares within domain
    ###################################################
    elif np.max(OldStormLabels) > 0 and np.max(StormLabels) > 0:
       xint,yint = np.meshgrid(range(xmat[0,0]+squarehalf,xmat[0,-1],squarehalf),range(ymat[0,0]+squarehalf,ymat[-1,0],squarehalf))
       buu = np.full(xint.shape, np.NaN)
       bvv = np.full(xint.shape, np.NaN)
       bww = np.full(xint.shape, np.NaN)
       for corx in range(0,int(np.size(xint,0))): 
           if flagplot==True:
               nij=-3
               #fig, axs = plt.subplots(np.size(xint,1),3, figsize=(6,2*np.size(xint,1)), facecolor='w', edgecolor='k')
               fig, axs = plt.subplots(int(0.5*np.size(xint,1))+1,6, figsize=(6,np.size(xint,1)), facecolor='w', edgecolor='k')
               axs = axs.ravel()
           for cory in range(0,int(np.size(xint,1))):
               if flagplot==True:
                   nij=nij+3
               oldsquare=oldbt[(squarehalf)*corx:(squarehalf)*corx+2*squarehalf,(squarehalf)*cory:(squarehalf)*cory+2*squarehalf]
               newsquare=newbt[(squarehalf)*corx:(squarehalf)*corx+2*squarehalf,(squarehalf)*cory:(squarehalf)*cory+2*squarehalf]
               if np.sum(oldsquare)<fftpixels or np.sum(newsquare)<fftpixels: # if there are too few storms, don't try to derive motion vectors. 
                  buu[corx,cory]=np.NaN
                  bvv[corx,cory]=np.NaN
                  bww[corx,cory]=np.NaN
               else:
                  dx,dy,amplitude,corrval = ffttrack(oldsquare,newsquare,tukey_window)
                  buu[corx,cory]=dx
                  bvv[corx,cory]=dy ## indices are upside down so need minus to get real-world dy-velocity
                  bww[corx,cory]=amplitude
                  if flagplot==True:
                      #print(corx)
                      #print(cory)
                      #print(nij)
                      axs[nij].pcolormesh(oldsquare)
                      axs[nij].set_title(str(int(np.sum(oldsquare))))
                      axs[nij+1].pcolormesh(newsquare)
                      axs[nij+1].set_title(str(int(np.sum(newsquare))))
                      axs[nij+2].pcolormesh(corrval)
                      axs[nij+2].set_title('(' + str(dx) + ',' + str(dy) + ')')
           if flagplot==True:
               plt.savefig(IMAGES_DIR + 'Correlations_' + write_file_ID + '_' + str(corx) +'.png')
               plt.close()
       
       # CHECK NEIGHBOURING VALUES FOR SMOOTHNESS
       # Ignore warnings about mean over empty array in this section
       with warnings.catch_warnings():
           warnings.simplefilter("ignore", category=RuntimeWarning)
           for corx in range(0,int(np.size(xint,0))):
               for cory in range(0,int(np.size(xint,1))):
                   bu_nb=np.nan
                   bv_nb=np.nan
                   if np.isnan(buu[corx,cory]) and np.isnan(bvv[corx,cory]):
                       continue
                   if corx==0:
                       if cory==0:
                           bu_nb = np.nanmean([buu[0,1],buu[1,0],buu[1,1]])
                           bv_nb = np.nanmean([bvv[0,1],bvv[1,0],bvv[1,1]])
                       elif cory==int(np.size(xint,1))-1:
                           bu_nb = np.nanmean([buu[0,cory-1],buu[1,cory],buu[1,cory-1]])
                           bv_nb = np.nanmean([bvv[0,cory-1],bvv[1,cory],bvv[1,cory-1]])
                       else:
                           bu_nb = np.nanmean([buu[0,cory+1],buu[0,cory-1],buu[1,cory-1],buu[1,cory],buu[1,cory+1]])
                           bv_nb = np.nanmean([bvv[0,cory+1],bvv[0,cory-1],bvv[1,cory-1],bvv[1,cory],bvv[1,cory+1]])
                   elif corx==int(np.size(xint,0))-1:
                       if cory==0:
                           bu_nb = np.nanmean([buu[corx,1],buu[corx-1,0],buu[corx-1,1]])
                           bv_nb = np.nanmean([bvv[corx,1],bvv[corx-1,0],bvv[corx-1,1]])
                       elif cory==int(np.size(xint,1))-1:
                           bu_nb = np.nanmean([buu[corx,cory-1],buu[corx-1,cory],buu[corx-1,cory-1]])
                           bv_nb = np.nanmean([bvv[corx,cory-1],bvv[corx-1,cory],bvv[corx-1,cory-1]])
                       else:
                           bu_nb = np.nanmean([buu[corx,cory+1],buu[corx,cory-1],buu[corx-1,cory-1],buu[corx-1,cory],buu[corx-1,cory+1]])
                           bv_nb = np.nanmean([bvv[corx,cory+1],bvv[corx,cory-1],bvv[corx-1,cory-1],bvv[corx-1,cory],bvv[corx-1,cory+1]])
                   elif cory==int(np.size(xint,1))-1:
                       bu_nb = np.nanmean([buu[corx,cory-1],buu[corx-1,cory],buu[corx-1,cory-1],buu[corx+1,cory-1],buu[corx+1,cory]])
                       bv_nb = np.nanmean([bvv[corx,cory-1],bvv[corx-1,cory],bvv[corx-1,cory-1],bvv[corx+1,cory-1],bvv[corx+1,cory]])
                   else:
                       bu_nb = np.nanmean([buu[corx,cory+1],buu[corx,cory-1],buu[corx-1,cory-1],buu[corx-1,cory],buu[corx-1,cory+1],buu[corx+1,cory-1],buu[corx+1,cory],buu[corx+1,cory+1]])
                       bv_nb = np.nanmean([bvv[corx,cory+1],bvv[corx,cory-1],bvv[corx-1,cory-1],bvv[corx-1,cory],bvv[corx-1,cory+1],bvv[corx+1,cory-1],bvv[corx+1,cory],bvv[corx+1,cory+1]])
                   if np.abs(buu[corx,cory]-bu_nb)>dd_tolerance*num_dt:
                       buu[corx,cory]=np.nan
                   if np.abs(bvv[corx,cory]-bv_nb)>dd_tolerance*num_dt:
                       bvv[corx,cory]=np.nan		   
		   
       ## ACTUAL DISPLACEMENT
       # Interpolate these displacements onto the full grid
       newumat = interpolate_speeds(xint, yint, xmat, ymat, buu)
       newvmat = interpolate_speeds(xint, yint, xmat, ymat, bvv)
       if len(newumat)==0:
           newumat=np.zeros(np.shape(OldStormLabels))
       if len(newvmat)==0:
           newvmat=np.zeros(np.shape(OldStormLabels))
       
       # Assign displacement to each of the old storms. 
       newlabel = np.zeros(OldStormLabels.shape)
       for ns in range(len(OldStormData)):
           jj = OldStormData[ns].storm
           labelind = np.where(OldStormLabels == jj)
           dx = np.mean(newumat[labelind])
           dy = np.mean(newvmat[labelind])
           if dx == 0.0 and dy == 0.0:
              newlabel[labelind] = jj
           else:
              for ii in range(np.size(labelind,1)):
                  newyind=labelind[1][ii]+int(np.around(dx))
                  newxind=labelind[0][ii]+int(np.around(dy))
                  if newxind>np.size(newlabel,0)-1 or newyind>np.size(newlabel,1)-1 or newxind<0 or newyind<0:
                     continue
                  elif newlabel[newxind,newyind]>0:
                     nq=int(newlabel[newxind,newyind]-1)
                     olddist=(xmat[newxind,newyind]-OldStormData[nq].centroidx)**2+(ymat[newxind,newyind]-OldStormData[nq].centroidy)**2
                     newdist=(xmat[newxind,newyind]-OldStormData[ns].centroidx)**2+(ymat[newxind,newyind]-OldStormData[ns].centroidy)**2
                     if newdist < olddist:
                        newlabel[newxind,newyind]=jj
                  else:
                     newlabel[newxind,newyind]=jj
       QuvL = newlabel
       AdvectedStorms=np.zeros([len(OldStormData),3])
       for ns in range(len(OldStormData)):
           jj = OldStormData[ns].storm
           centrind=np.where(QuvL==jj)
           if np.size(centrind,1) == 0:
              continue
           else:
              AdvectedStorms[ns][0] = np.mean(xmat[centrind])
              AdvectedStorms[ns][1] = np.mean(ymat[centrind])
              AdvectedStorms[ns][2] = int(np.size(centrind,1))
    
       ###################################################
       # NOW LOOP THROUGH StormData AND CHECK FOR OVERLAP WITH
       # ADVECTED OldStormData STORMS
       ###################################################
       wasnum = np.zeros(len(StormData))
       qbins = range(int(np.max(OldStormLabels))+2)
       qarea = np.ones([int(np.max(OldStormLabels))+1])
       qlife = np.ones([int(np.max(OldStormLabels))+1])
       for qq in range(len(OldStormData)):
           if AdvectedStorms[qq,2] > 0:
              qarea[qq+1]=AdvectedStorms[qq,2]
           qlife[qq+1]=OldStormData[qq].life
       for ns in range(numstorms):
           jj = ns + 1 # first storm is labelled 1, but python indeces start at 0.
           C=np.where(StormLabels==jj)
           StormData += [StormS(jj,StormLabels,var, xmat,ymat,newwas, newumat, newvmat, num_dt, misval, doradar, under_threshold, extra_thresh=extra_thresh, storm_history=True, string=None, rarray=rarray, azarray=azarray)]
           wasarray[C]=int(jj)
           lifearray[C]=1
           
           ###################################################
           # CHECK OVERLAP WITH QHIST
           # IF NO OVERLAP, THEN
           # GENERATE (halo) km RADIUS AROUND CENTROID
           # CHECK FOR OVERLAP WITHIN (halo) km OF CENTROID
           ###################################################
           qhist = (np.histogram(QuvL[np.where(StormLabels==jj)], qbins))[0][:]/float(StormData[ns].area) + (np.histogram(QuvL[np.where(StormLabels==jj)],qbins))[0][:]/qarea[:]
           #if nt==35 and jj==131:
           #   raise ValueError("Check storm 131 (or 120)")
               
           if np.max(qhist[1:]) < lapthresh:
              newblob = 0*xmat
              blobind = np.where((xmat-StormData[ns].centroidx)**2+(ymat-StormData[ns].centroidy)**2 < halosq)
              newblob[blobind] = newblob[blobind]+1
              qhist = (np.histogram(QuvL[np.where(newblob==1)],qbins))[0][:]/float(StormData[ns].area)+(np.histogram(QuvL[np.where(newblob==1)],qbins))[0][:]/qarea[:]
           ###################################################
           # IF OVERLAP, THEN
           # - INHERIT "WAS"
           # - UPDATE "LIFE" AND "TRACK" AND "WASDIST"
           # - INHERIT "dx" AND "dy" (ONLY UPDATE IF SINGLE OVERLAP)
           ###################################################
           if np.max(qhist[1:]) >= lapthresh:
              numlaps = np.where(qhist[1:] >= lapthresh)
              ###################################################
              # IF MORE THAN ONE GOOD OVERLAP
              # KEEP PROPERTIES OF STORM WITH LARGEST OVERLAP
              # IF MORE THAN ONE LARGEST, KEEP NEAREST IN CENTROID
              ###################################################
              if np.size(numlaps,1) > 1:
                 lapdist = np.zeros([np.size(numlaps,1)])
                 sectlap = np.zeros([np.size(numlaps,1)])
                 for kkind in range(np.size(numlaps,1)):
                     qindex = np.squeeze(numlaps[0][kkind])
                     lapdist[kkind] = np.sqrt((StormData[ns].centroidx-AdvectedStorms[qindex,0])**2+(StormData[ns].centroidy-AdvectedStorms[qindex,1])**2)
                     sectlap[kkind] = np.size(np.where((QuvL==qindex+1) & (StormLabels==jj)),1)
                 kmax = np.where(sectlap==np.max(sectlap))
                 if np.size(kmax,1) > 1:
                    kkmax = kmax[0][np.where(lapdist[kmax[0][:]]==np.min(lapdist[kmax[0][:]]))]
                    if np.size(kkmax) > 1:
                        kkmax = kmax[0][kkmax[0]]
                 else:
                    kkmax = kmax[0][0]  
                 kindex = np.squeeze(numlaps[0][kkmax])
                 StormData[ns].inherit_properties(jj, OldStormData, kindex, QuvL, StormLabels, qhist, lapthresh, misval, single_overlap=False)
                 wasarray[C] = OldStormData[kindex].was
                 lifearray[C]=StormData[ns].life
              ###################################################
              # SINGLE OVERLAP
              ###################################################
              else:
                 zindex = np.squeeze(numlaps[0][0])
                 lapdist = np.sqrt((StormData[ns].centroidx-OldStormData[zindex].centroidx)**2+(StormData[ns].centroidy-OldStormData[zindex].centroidy)**2)
                 StormData[ns].inherit_properties(jj, OldStormData, zindex, QuvL, StormLabels, qhist, lapthresh, misval, single_overlap=True)
                 wasarray[C] = OldStormData[zindex].was
                 lifearray[C]=StormData[ns].life

           ###################################################
           # IF NO OVERLAP, THEN (NEW STORM)
           # - "WAS" SET TO CURRENT MAX LABEL +1
           # - UPDATE "LIFE" AND "TRACK" AND "WASDIST" FOR A NEW STORM
           ###################################################
           else:
              StormData[ns].was=newwas
              wasarray[C] = newwas
              StormData[ns].life=1
              lifearray[C] = 1
              newwas=newwas+1
              StormData[ns].centroidx0 = StormData[ns].centroidx+0.0
              StormData[ns].centroidy0 = StormData[ns].centroidy+0.0
       wasnum=np.array([StormData[ns].was for ns in range(len(StormData))])
       ###################################################
       # QUICK SANITY CHECK
       # ACCRETED SHOULD NEVER BE A VALUE
       # SIMILAR TO EXISTING STORM ID
       ###################################################
       for ns in range(len(StormData)):
           jj=StormData[ns].storm
           if StormData[ns].accreted[-1]==misval:
              continue
           else:
              for acnum in range(np.size(StormData[ns].accreted)):
                  acind=np.where((wasnum-StormData[ns].accreted[acnum])==0)
                  if np.size(acind,1) > 0:
                     StormData[ns].accreted[acnum]=misval
              #acnew=np.where(StormData[ns].accreted > misval) 
              acnew = [aci for aci in StormData[ns].accreted if aci > misval]
              if np.size(acnew) > 0:
                 for acindex in range(np.size(acnew)):
                    StormData[ns].accreted[acnum]=acnew[acindex]
              else:
                 StormData[ns].accreted=[misval]
       ###################################################
       # TRACKING MERGING BREAKING
       # MULTIPLE STORMS AT T (StormData) MAY HAVE SAME LABEL "WAS"
       # FIND STORM WITH LARGEST OVERLAP AT T+1 WITH ADVECTED q(T)
       # THIS IS THE "PARENT" STORM,
       # "PARENT" VECTOR WITH INDICES OF NEW LABELS FOR "CHILD" STORMS
       # STORMS WITH SAME WAS BUT FUTHER FROM CENTROID ARE "CHILD", VALUE "PARENT"
       ###################################################
       for ns in range(len(StormData)):
           jj=StormData[ns].storm
           if StormData[ns].wasdist==[misval]:
              continue
           wasind=np.where(wasnum==wasnum[ns])
           wasseplength=0
           for kkind in range(np.size(wasind)):
              if StormData[wasind[0][kkind]].wasdist==[misval]:
                 continue
              else:
                 wasseplength=wasseplength+1
           wassep=np.zeros(wasseplength)
           if np.size(wasind) > 1:
              kkval=0
              for kkind in range(np.size(wasind)):
                 if StormData[wasind[0][kkind]].wasdist==[misval]:
                    continue
                 else:
                    wassep[kkval]=StormData[wasind[0][kkind]].wasdist
                    kkval=kkval+1
           else:
              wassep=StormData[wasind[0][0]].wasdist
           ########################################
           # WASSEP NOW CONTAINS ALL NON-ZERO OVERLAP VALUES
    	   # FIND THE MAXIMUM (THIS WILL BE THE PARENT)
    	   # ALL OTHER STORMS WILL BE THE CHILDREN
    	   #########################################
           kmax=np.where(wassep==np.max(wassep))
           kkmax=np.min(kmax)
           children=[]
           for kkind in range(np.size(wassep)):
              if not kkind == kkmax:
                 StormData[wasind[0][kkind]].child=StormData[wasind[0][kkmax]].was
                 StormData[wasind[0][kkind]].was=newwas
                 wasarray[np.where(StormLabels==wasind[0][kkind]+1)] = newwas
                 StormData[wasind[0][kkind]].life=StormData[wasind[0][kkmax]].life
                 lifearray[np.where(StormLabels==wasind[0][kkind]+1)] = StormData[wasind[0][kkmax]].life
                 newwas=newwas+1
                 wasnum[wasind[0][kkind]]=StormData[wasind[0][kkind]].was
                 children.append(StormData[wasind[0][kkind]].was)
                 StormData[wasind[0][kkind]].wasdist=misval
                 StormData[wasind[0][kkind]].centroidx0 = StormData[wasind[0][kkind]].centroidx + 0.0
                 StormData[wasind[0][kkind]].centroidy0 = StormData[wasind[0][kkind]].centroidy + 0.0
           ###################################################
           # UPDATE PARENT STORM WITH CHILDREN
           ###################################################
           if np.size(children) > 0:
               StormData[wasind[0][kkmax]].parent = children
    
    
    return StormData, newwas, StormLabels, newumat, newvmat, wasarray, lifearray

###################################################
# interpolate_speeds used for (dx,dy) calculation where no objects are identified.
###################################################

def interpolate_speeds(xint, yint, xmat, ymat, buu):
    valid_mask = ~np.isnan(buu)
    coords = np.array(np.nonzero(valid_mask)).T
    values = buu[valid_mask]
    if np.size(values)>=4:
       it = interpolate.LinearNDInterpolator(coords, values, fill_value=0)
       filled = it(list(np.ndindex(buu.shape))).reshape(buu.shape)
       fu = interpolate.interp2d(xint[0,:],yint[:,0],filled,kind='cubic')
       newumat = fu(xmat[0,:],ymat[:,0])
    else:
       newumat = []
    return newumat

###################################################
# label_storms IS A FLOOD FILL ALGORITHM RETURNING
# AN ARRAY OF UNIQUELY LABELLED STORMS
# THIS FUNCTION IS ESSENTIAL AND SHOULD ONLY BE ALTERED 
# BY EXPERIENCED USERS
###################################################

def label_storms(bt,minarea,threshold,struct,under_threshold):
    binbt=np.zeros_like(bt)
    if under_threshold:
        binbt[np.where(bt<threshold)]=1
    else:
        binbt[np.where(bt>threshold)]=1
    id_regions, num_ids = ndimage.label(binbt, structure=struct)
    id_sizes = np.array(ndimage.sum(binbt, id_regions, range(num_ids+1)))
    area_mask = (id_sizes < minarea)
    binbt[area_mask[id_regions]] = 0
    id_regions, num_ids = ndimage.label(binbt, structure=struct)
    print('num_ids = ',num_ids)
    
    return id_regions

##############################################################
# ffttrack  
##############################################################
# [dx, dy, amp] = ffttrack(s1, s2, method)
# Input: 
# s1 = oldsquare
# s2 = newsquare
# method = 1 for TUKEY WINDOW (TAPERED COSINE)
# Output:
# dx = distance in x-direction from previous cell
# dy = distance in y-direction from previous cell
# amp = amplitude
# ffv = full output, only needed for testing (if plotting with flagplot)
##############################################################

def ffttrack(s1,s2,method):

    leno = max(np.size(s1,0),np.size(s1,1))

    if method == 1:
       alpha = max(0.1, 10.0/leno)
       xhan = np.array(np.arange(0.5,leno+0.5))
       hann1 = np.ones([np.size(xhan)])
       hann1[np.where(xhan<alpha*leno/2.)] = 0.5*(1+np.cos(np.pi*(2*xhan[np.where(xhan<alpha*leno/2.)]/(alpha*leno)-1)))
       hann1[np.where(xhan>leno*(1-alpha/2.))] = 0.5*(1+np.cos(np.pi*(2*xhan[np.where(xhan>leno*(1-alpha/2.))]/(alpha*leno)-2./alpha+1)))
       hann2 = hann1.conj().transpose()*hann1
    else:
       xhan = np.array(np.arange(0.5,leno+0.5))
       hann1 = np.ones([np.size(xhan)])
       hann2 = hann1.conj().transpose()*hann1

## FIND CONVOLUTION S1, S2 USING FFT

    b1 = s1*hann2
    b2 = s2*hann2

    m1 = b1-np.mean(b1)
    m2 = b2-np.mean(b2)

    normval = np.sqrt(np.sum(m1**2)*np.sum(m2**2))
	
    #ffv = signal.fftconvolve(s1,s2,mode='same')
    ffv = np.real(np.fft.ifft2(np.fft.fft2(m2)*(np.fft.fft2(m1)).conj()))

    val = np.max(ffv)
    ind = np.where(ffv==val)

    #print 'max ffv and ind -> ',val, ind
    dx = ind[1][0]
    dy = ind[0][0]

## 1hour -> 25km(leno/2) ; 5mins -> 2km(leno/10) : 10mins -> 4km(leno/5)
    cv = leno/2 # Org. from Thorld = 25km
    #cv = leno/2 # For 200m grids = 20km
    if dx > cv:
       dx = dx - leno # Org. from Thorld
       #dx = dx - (cv*(dx/cv))
    if dy > cv:
       dy = dy - leno # Org. from Thorld
       #dy = dy - (cv*(dy/cv))
    amp = val/normval

    return dx, dy, amp, ffv

###################################################
# write_storms produces TXT file for analysis of tracked object properties
###################################################
    
def write_storms(file_ID, init_time, now_time, label_method, squarelength, rafraction, newwas, StormData, doradar, misval, IMAGES_DIR):
   if not(isdir(IMAGES_DIR)): makedirs(IMAGES_DIR)
   #print("IMAGES_DIR + file_ID +'.txt'=", IMAGES_DIR + file_ID +'.txt')
   fw = open(IMAGES_DIR + 'history_' + file_ID +'.txt','w')
   fw.write('missing_value=' + str(misval) + '\r\n')
   fw.write('Start date and time=' + init_time.strftime('%d/%m/%y-%H%M') + '\r\n')
   fw.write('Current date and time=' + now_time.strftime('%d/%m/%y-%H%M') + '\r\n')
   fw.write('Label method=' + label_method + '\r\n')
   fw.write('Squarelength=' + str(squarelength) + '\r\n')
   fw.write('Rafraction=' + str(rafraction) + '\r\n')
   fw.write('total number of tracked storms=' + str(newwas-1) + '\r\n')
   for ns in range(len(StormData)):
       fw.write('storm ' + str(StormData[ns].was))
#       fw.write(' label=' + str(StormData[ns].storm)) # Matches storm to label in mask. Actually no need for this as it is the same as it matches the order of the storms.
       fw.write(' area=' + str(StormData[ns].area))
       fw.write(' centroid=' + str(round(StormData[ns].centroidx,2)) + ',' + str(round(StormData[ns].centroidy,2)))
       fw.write(' centroid0=' + str(round(StormData[ns].centroidx0,2)) + ',' + str(round(StormData[ns].centroidy0,2)))
       fw.write(' box=' + str(StormData[ns].boxleft) + ',' + str(StormData[ns].boxup) + ',' + str(StormData[ns].boxwidth) + ',' + str(StormData[ns].boxheight))
       fw.write(' life=' + str(StormData[ns].life))
       fw.write(' dx=' + str(round(StormData[ns].dx,2)) + ' dy=' + str(round(StormData[ns].dy,2)))

       if doradar:
           fw.write(' range=' + str(round(StormData[ns].rangel,2)) + ',' + str(round(StormData[ns].rangeu,2)))
           fw.write(' azimuth=' + str(round(StormData[ns].azimuthl,2)) + ',' + str(round(StormData[ns].azimuthu,2)))
       fw.write(' meanv=' + str(round(StormData[ns].meanvar,2)))
       fw.write(' extreme=' + str(round(StormData[ns].extreme,2)) + ' accreted=')
       if np.size(StormData[ns].accreted) > 1:
           for acind in range(np.size(StormData[ns].accreted)-1):
               fw.write(str(StormData[ns].accreted[acind]) + ',')
           fw.write(str(StormData[ns].accreted[np.size(StormData[ns].accreted)-1]) + ' parent=')
       else:
           fw.write(str(StormData[ns].accreted[-1]) + ' parent=')
       fw.write(str(StormData[ns].child) + ' child=')
       if np.size(StormData[ns].parent) > 1:   
           for acind in range(np.size(StormData[ns].parent)-1):
               fw.write(str(StormData[ns].parent[acind]) + ',')
           fw.write(str(StormData[ns].parent[np.size(StormData[ns].parent)-1]) + '\r\n')
       else:
           fw.write(str(StormData[ns].parent[-1]) + '\r\n')
   fw.close()
