"""
   A module for SBAS application in python.
   The first version was developed by Wanpeng, Feng, @NRCan, 2016-03-21
   
"""
import numpy as np, pSAR
import scipy.stats
import datetime, opt
###############################################################################
def jd2datetime(jd):
    #
    return pSAR.jdutil.jd_to_datetime(jd).strftime("%Y%m%dT%H%M%s")
    
def dates2jd(in_date_str):
    #
    '''
    in_date_str is a format like 20150101T023830
    
    '''
    inyear = float(in_date_str[0:4])
    inmonth= float(in_date_str[4:6])
    inday  = float(in_date_str[6:8])
    in_fraction_day = float(in_date_str[9:11])/24.+\
                      float(in_date_str[11:13])/(60.*24.) + \
                      float(in_date_str[13:])/(60.*24.*60.)
    return pSAR.jdutil.date_to_jd(inyear,inmonth,inday+in_fraction_day)
    
###############################################################################    
def points2file(vec,outfile):
    #
    # dims = vec.shape
    # ncol = dims[1]
    # 
    np.savetxt(outfile,vec)
#################################################################
def checkevent(infile,eventdate):
    """
      #
      Check if input file covers the specific date
      the date should be in YYYYMMDD format...
      #
    """
    master,slave = name2dateinfo(infile)
    if master < eventdate and slave > eventdate:
       return True
    else:
       return False
#################################################################
def name2dateinfo(infile):
    """
       File name should be in geo_XXXXXXXX_XXXXXXXX_TRACK.phs
       #
    """
    file_names = infile.split('.')[0]
    master = file_names[4:12]
    slave  = file_names[13:21]
    #
    return int(master),int(slave)
    
###############################################################################
def nw_infs2udates(in_infs):
    #
    # return unique dates from a set of interferograms
    #
    in_infs = np.array(in_infs)
    in_infs = np.reshape(in_infs,in_infs.shape[0]*in_infs.shape[1])
    return np.unique(in_infs)
###############################################################################
def nw_read_datalist(inlist):
    #
    """
    inlist is a list file recording all data in an ascii file
    format of each line is simple like:
    #
    <fullpath_of_file> <baseline> <master date> <slave date>
    """
    datapath=[]
    fid = open(inlist,'r')
    #
    dinfo = []
    for cline in fid:
        #
        cline = cline.split('\n')[0]
        tmp = cline.split()
        datapath.append(tmp[0])
        dinfo.append([float(tmp[1]),float(tmp[2]),float(tmp[3])])
    #
    datainfo = np.array(dinfo)
     #
    return datapath, datainfo
###############################################################################
def nw_min_lag_check(udates,master,slave):
    #
    flag1 = udates > master
    flag2 = udates < slave
    flag  = np.logical_and(flag1,flag2)
    if sum(flag) > 0:
        return False
    else:
        return True
###############################################################################
def nw_list2minlagpairs(in_list):
    #
    dpath,dinfo = nw_read_datalist(in_list)
    dpath       = np.array(dpath)
    #
    indMatrix,udate   = nw_dateinfo2uniquedate(dinfo[:,1:3])
    selinds = np.zeros(len(dpath),dtype=bool)
    #
    for ni in range(len(dpath)):
        #
        flag = nw_min_lag_check(udate,dinfo[ni,1],dinfo[ni,2])
        if flag:
            selinds[ni] = True
    #
    return dpath[selinds],dinfo[selinds,:]
###############################################################################
def nw_synsign(in_pair_1,in_pair_2,expected_pair):
    #
    # in_pair_1 = np.array(in_pair_1)
    # in_pair_2 = np.array(in_pair_2)
    # expected_pair = np.array(expected_pair)
    #
    master_sign = None
    slave_sign  = None
    ###########################################################################
    in_index_1 = np.array([1,1])
    in_index_2 = np.array([1,1])
    if in_pair_1[0] == expected_pair[0]:
        master_sign = in_index_1[0]
    if in_pair_1[1] == expected_pair[0]:
        master_sign = in_index_1[1] * -1
    #
    if in_pair_1[0] == expected_pair[1]:
        master_sign = in_index_1[0] * -1
    if in_pair_1[1] == expected_pair[1]:
        master_sign = in_index_1[1]
    # 
    # slave
    if in_pair_2[0] == expected_pair[0]:
        slave_sign = in_index_2[0]
    if in_pair_2[1] == expected_pair[0]:
        slave_sign = in_index_2[1] * -1
    if in_pair_2[0] == expected_pair[1]:
        slave_sign = in_index_2[0] * -1
    if in_pair_2[1] == expected_pair[1]:
        slave_sign = in_index_2[1]
    return [master_sign,slave_sign]
#
#############################################################################
def nw_locatepairsbydate(ref_pairs,in_date):
    # flag for locating the place 
    flag = np.copy(ref_pairs) * 0.
    flag[ref_pairs[:,0] == in_date,0] = 1
    flag[ref_pairs[:,1] == in_date,1] = 1
    # 
    # flag is N x 2 array 
    # 
    return flag[:,0] + flag[:,1], flag 
#
###############################################################################
def nw_pairindex(in_pairs,master,slave):
    #
    flag1,mflag = nw_locatepairsbydate(in_pairs,master)
    flag2,sflag = nw_locatepairsbydate(in_pairs,slave)
    flag = flag1 + flag2
    #
    if sum(flag==2) < 1:
        return None
    else:
        return np.where(flag==2)

#
###############################################################################
#
def nw_infs2newpair(ref_pairs,master,slave):
    #
    # return index of existing pairs for generating a new pair with master and slave
    #
    udates = nw_infs2udates(ref_pairs)
    #
    if len(np.logical_and(udates==master,udates==slave)) == 0:
        return False,np.array([None,None,None,None])
    else:
        #
        outindex = []
        outflag  = False
        flag1,flag_master = nw_locatepairsbydate(ref_pairs,master)
        flag2,flag_slave  = nw_locatepairsbydate(ref_pairs,slave)
        mpairs = ref_pairs[flag1>0,:]
        spairs = ref_pairs[flag2>0,:]
        #
        for inmaster in range(mpairs.shape[0]):
            #
            for inslave in range(spairs.shape[0]):
                
                cdates = np.array([mpairs[inmaster,0],mpairs[inmaster,1],\
                          spairs[inslave, 0],spairs[inslave,1]])
                #
                cudates = np.unique(cdates)
                flag_m  = sum(cdates == master)
                flag_s  = sum(cdates == slave)
                #
                if ( len(cudates) == 3 and flag_m == 1 and flag_s == 1 ):
                    #
                    outindex.append([mpairs[inmaster,0],mpairs[inmaster,1],\
                                     spairs[inslave, 0],spairs[inslave, 1]])
                    outflag = True
                #
        #
        return outflag, np.array(outindex)
        
#    
###############################################################################
#
def nw_syn2third(date1,date2):
    #
    # here, date1 and data2 are time information for pair1 and pair2
    # helpinformation was provided by Wanpeng Feng, @NRCan, 2017-02-17
    #
    if date1[0] == date2[0]:
       if date1[1] < date2[1]:
          sign = [-1,1]
       else:
          sign = [1,-1]
    else:
       if date1[0] > date2[0] and date1[1] == date2[1]:
          sign = [-1,1]
       elif date1[0] < date2[0] and date1[1] == date2[1]:
          sign = [1,-1]
       else:
          print(" ERROR: they have to have a common date")
          sign = [0,0]
    #
    return sign      
###############################################################################
def nw_phaseloop2weight(data):
    #
    data = np.sqrt(data ** 2)
    mean_data = np.mean(data[np.where(data!=0)])
    rdata = data / mean_data
    ind   = np.argsort(rdata[np.where(rdata!=0)])
    update_mean = np.mean(data[np.where(data!=0)][ind[0:int(ind.shape[0]/4)]])
    #print(mean_data,mean_data2)
    return data / update_mean
#################################################################
def nw_phaseloop_outname(pair1,pair2,pair3):
    #print(type(pair1))
    pair = []
    pair.append(pair1)
    pair.append(pair2)
    pair.append(pair3)
    pair = np.array(pair)
    pair = np.unique(pair)
    outname = str(int(pair[0])) + '-' + str(int(pair[1])) + '-' + str(int(pair[2])) 
    return outname
    
###############################################################################
def nw_udates2fullinfs(udates):
    #
    # return a full list of potential interferograms from SAR acquisitions
    # by Wanpeng Feng, @Ottawa, 2017-02-19
    #
    ndate    = len(udates)
    fullinfs = np.zeros([ ndate * (ndate - 1) / 2 ,2])
    counter  = 0
    for nk in range(len(udates)):
        for nj in range(len(udates)):
            master = pSAR.sbas.dates2jd(str(int(udates[nk]))+"T010101")
            slave  = pSAR.sbas.dates2jd(str(int(udates[nj]))+"T010101")
            #
            if master < slave:
                fullinfs[counter,:] = [udates[nk],udates[nj]]
                counter += 1
    #
    return fullinfs
#
###############################################################################
#
def nw_phaseloop_sign(pair1,pair2,pair3):
    #
    p1_m = pair1[0]
    p1_s = pair1[1]
    p2_m = pair2[0]
    p2_s = pair2[1]
    p3_m = pair3[0]
    p3_s = pair3[1]
    #
    signs = [1,1,1]
    if p1_m == p2_m:
       if p1_s > p2_s:
          signs[0] = -1
       else:
          signs[1] = -1
    #
    if p1_m == p3_m:
       if p1_s > p3_s:
          signs[0] = -1
       else:
          signs[2] = -1
    #
    if p2_m == p3_m:
       if p2_s > p3_s:
          signs[1] = -1
       else:
          signs[2] = -1
    #
    return signs
#################################################################
def sbas_invmatrix(indMatrix,baseinfo,dem=True,ndate=10,wavelength=5.54,r=693000,inc=32):
    #
    #
    ninf = indMatrix.shape[0] 
    A    = np.zeros([ninf+1,ndate])
    for ni in np.arange(ninf):
        A[ni,indMatrix[ni,0]] = -1
        A[ni,indMatrix[ni,1]] =  1
    #
    if dem:
       ba = 4 * np.pi * baseinfo / (wavelength * r * np.sin(inc/180.*np.pi))  
       MA = np.zeros([ninf+1,ndate+1])
       MA[:,0:-1] = A
       MA[0:-1,-1] = ba
       A = MA
       # A = np.concatenate((A,ba),axis=1)
    #
    A[ninf,0] = 1 
    return A

#################################################################
def nw_read_uniquedate(inlist):
    #
    datapath,datainfo = nw_read_datalist(inlist)
    #
    indMatrix,udate   = nw_dateinfo2uniquedate(datainfo[:,1:3])
    #
    return indMatrix, udate
#################################################################
def nw_dateinfo2uniquedate(datas):
    #
    # master and slave time mastrix to unique dates
    #
    # print(datas)
    dims = datas.shape
    dates= datas.reshape([dims[0]*dims[1],1])
    udate= np.unique(dates).astype(int)
    indMatrix = datas.astype(int) * 0
    #
    for cind in range(dims[0]):
        ind1 = np.where(udate==datas[cind,0])[0]
        ind2 = np.where(udate==datas[cind,1])[0] 
        indMatrix[cind,:] = [ind1,ind2] 
    #
    return indMatrix, udate
################################################################
def nw_checkpair(inpair,indMatrix):
    #
    outind = -1
    dims = indMatrix.shape
    for cind in range(dims[0]):
        mid = indMatrix[cind,0]
        sid = indMatrix[cind,1]
        if mid == inpair[0] and sid == inpair[1]:
           outind = cind
           break
    #
    return outind
###############################################################      
def nw_dateinfo2loopind(indMatrix):
    #
    dims = indMatrix.shape
    allind = indMatrix.reshape([dims[0]*dims[1],1]);
    uinds  = np.unique(allind).astype(int)
    ndate  = uinds.shape
    counter= 0
    #
    for mind in uinds:
        for sind in range(mind+1,ndate[0]):
            for tind in range(sind+1,ndate[0]):
                pair1 = [mind,sind]
                pair2 = [mind,tind]
                pair3 = [sind,tind]
                #
                outind1 = nw_checkpair(pair1,indMatrix)
                outind2 = nw_checkpair(pair2,indMatrix)
                outind3 = nw_checkpair(pair3,indMatrix)
                if outind1 > -1 and outind2 > -1 and outind3 > -1:
                   counter += 1
                   if counter == 1: 
                      outloopid = np.unique(np.array([outind1,outind2,outind3])).reshape((1,3))
                   else:
                      outloopid = np.concatenate((outloopid,np.unique(np.array([outind1,outind2,outind3])).reshape(1,3)))
    #
    # ndims = outloopid.shape
    # 
    return outloopid
#################################################################
def demcorrelatedsign_removal(phs,mlon,mlat,dem,maskphs=None,x0=[1.,1.,1.],\
                              maxfev=1000000,iterations=1,scale=1):
    #
    if maskphs is None:
       maskphs = np.copy(phs)
    #
    subphs  = maskphs[::scale,::scale]
    subdem  = dem[::scale,::scale]
    sublon  = mlon[::scale,::scale]
    sublat  = mlat[::scale,::scale]
    #
    #
    for ni in range(iterations):
      #
      demvec = subdem[np.where(subphs!=0)]
      phsvec = subphs[np.where(subphs!=0)]
      phsvec = phsvec[np.where(demvec>0)]
      demvec = demvec[np.where(demvec>0)]
      indem,inphs = opt.powlaw_sortdata(demvec,phsvec)
      coef,res    = opt.est_powlawcoef(indem,inphs,opt.powlaw_errfun_1,x0,maxfev=maxfev)
      simphs      = opt.powlaw_fitfun_1(coef,subdem,iskeepcons=True)
      tmpphs      = subphs - simphs
      tmpphs[subphs==0] = 0
      #
      print(" NO %d: %f %f %f" % (ni,coef[0],coef[1],coef[2]))
      corphs,simph= orbitalCor(tmpphs,sublon,sublat)
      #
      tmphs      = subphs - simph
      tmphs[subphs==0] = 0
      subphs     = np.copy(tmphs)
      #
      x0          = coef
    # 
    #
    demphs      = opt.powlaw_fitfun_1(coef,dem,iskeepcons=False)
    masksub     = maskphs - demphs
    masksub[maskphs==0] = 0
    corphs,orbphs  = orbitalCor(masksub,mlon,mlat)
    #orbphs = np.zeros(demphs.shape)
    outphs         = phs - demphs - orbphs
    outphs[phs==0] = 0
    
    #
    return outphs,demphs,orbphs,coef
###############################################################################
def orbitalCor_theisen(phs,mlon,mlat,scale=1,npara=3):
    #
    corphs,sim = orbitalCor(phs,mlon,mlat,dem='NULL',thirddep='NULL',scale=1,npara=3)
    res        = np.abs(phs - sim)
    medres     = np.median(res)
    res_shift  = medres - np.mean(res)
    flag       = res < res_shift+medres
    corphs,sim = orbitalCor(phs[flag],mlon[flag],mlat[flag],dem='NULL',thirddep='NULL',scale=1,npara=3)
#    
###############################################################################
#    
def orbitalCor(phs,mlon,mlat,dem='NULL',thirddep='NULL',scale=1,npara=3):
    #
    # Make sure dem is an float array if it is given...
    # print(phs.shape,mlon.shape,mlat.shape,dem.shape)
    # 
    # ---
    # npara can be 3,4,6
    # added by Wanpeng Feng, @NRCan, 2017-03-09
    #
    phs[np.isnan(phs)] = 0
    #
    cphs  = phs[::scale,::scale]
    cmlon = mlon[::scale,::scale]
    cmlat = mlat[::scale,::scale]
    #
    if not isinstance(dem,str):
       #
       cdem = dem[::scale,::scale]
       cdem[np.isnan(cdem)] = 0
       #
    else: 
       cdem = np.array([])
    #
    if not isinstance(thirddep,str):
       cthirddep = thirddep[::scale,::scale]
    else:
       cthirddep = np.array([])
    #
    cmlon = cmlon[np.where(cphs != 0)]
    cmlat = cmlat[np.where(cphs != 0)]
    #
    # extract dem value if not zero
    #
    if len(cdem)>0:
       cdem = cdem[np.where(cphs != 0)]
    if len(cthirddep)>0:
       cthirddep = cthirddep[np.where(cphs != 0)]
    #
    cphs = cphs[np.where(cphs != 0)]
    cons = cphs * 0 + 1
    #
    # create A matrix for orbital error correction
    #
    if npara == 3:
       if len(cdem) == 0:
          A = np.array([cons,cmlon,cmlat])
       else:
          if len(cthirddep) == 0:
             print( " Processes: Generating A...")
             A = np.array([cons,cmlon,cmlat,cdem])
          else:
             A = np.array([cons,cmlon,cmlat,cdem,cthirddep])
    elif npara == 4:
       lonlatm = cmlon * cmlat
       if len(cdem) == 0:
          A = np.array([cons,cmlon,cmlat,lonlatm])
       else:
          if len(cthirddep) == 0:
             A = np.array([cons,cmlon,cmlat,lonlatm,cdem])
          else:
             A = np.array([cons,cmlon,cmlat,lonlatm,cdem,cthirddep])
    elif npara == 6:
        lonlatm = cmlon * cmlat
        #
        if len(cdem) == 0:
           A = np.array([cons,cmlon,cmlat,lonlatm,cmlon**2,cmlat**2])
        else:
           A = np.array([cons,cmlon,cmlat,lonlatm,cmlon**2,cmlat**2,cdem])
    #
    # solve a large linear problem with numpy least square method
    A = A.T
    
    cof,res,nrnk, s = np.linalg.lstsq(A,cphs)
    #
    if npara == 3: 
       sim = cof[0] + mlon * cof[1] + mlat * cof[2] 
       nodem = 3
    elif npara == 4:
       sim = cof[0] + mlon * cof[1] + mlat * cof[2] + mlon * mlat * cof[3] 
       nodem = 4
    elif npara == 6:
       sim = cof[0] + mlon * cof[1] + mlat * cof[2] + mlon * mlat * cof[3] +\
                      mlon**2 * cof[4] + mlat ** 2 * cof[5]
       nodem = 6
    #
    # print(cof[nodem])
    demsim = 0
    thirdsim = 0
    if len(cdem) > 0:
       demsim = dem * cof[nodem]
       nodem = nodem + 1
    #
    if len(cthirddep)>0:
       thirdsim = thirddep * cof[nodem]
    #
    sim = sim + demsim + thirdsim
    #
     
    corphs = phs - sim
    corphs[phs==0] = 0
    # sim[np.where(phs==0)] = 0
    return corphs,sim
################################################################### 
def dates2dates(dates):
    #
    dims = dates.shape
    #
    outdates = np.array([datetime.datetime(int(str(dates[i])[0:4]),\
                                  int(str(dates[i])[4:6]),\
                                  int(str(dates[i])[6:8])) for i in range(dims[0])])
    return outdates
###################################################################
def newdate(indate,shiftdays):
    #
    # return a new data by given a shift of few days...
    # for example "20150101" 
    #     with a shift 0f  1 can reach "20150102"
    #     with a shift of -1 can reach "20141231"
    #
    indate = str(int(indate))
    my = int(indate[0:4])
    mm = int(indate[4:6])
    md = int(indate[6:8]) 
    dform = datetime.datetime(my,mm,md)
    newd  = dform + datetime.timedelta(days=shiftdays)
    return int(newd.strftime('%Y%m%d'))
    #print(newd)
###################################################################
def diff2dates(mdate,sdate):
    #
    # mdate and sdate are both in string format...
    # like '20150101'
    #
    my = int(mdate[0:4])
    mm = int(mdate[4:6])
    md = int(mdate[6:8])
    sy = int(sdate[0:4])
    sm = int(sdate[4:6])
    sd = int(sdate[6:8])
    # print(my,mm,md)
    # print(sy,sm,sd)
    d1 = datetime.datetime(my,mm,md)
    d2 = datetime.datetime(sy,sm,sd)
    dd = d2 - d1
    return dd.days
#################################################################
def tinfo2uniquedate(infpairs):
    #
    tinfo = infpairs.ravel()
    dates = np.unique(tinfo)
    return dates
#################################################################
def nw_coefA(loopmatrix,pairInfo,weights):
    #
    maxweig = np.max(weights)
    #
    nloops  = loopmatrix.shape[0]
    A1      = np.zeros([nloops, pairInfo.shape[0]])
    A2      = np.zeros([pairInfo.shape[0],pairInfo.shape[0]])
    #
    for ni in range(nloops):
       sign = nw_phaseloop_sign(pairInfo[loopmatrix[ni,0],:],pairInfo[loopmatrix[ni,1],:],pairInfo[loopmatrix[ni,2],:])
       A1[ni,loopmatrix[ni,0]] = sign[0]
       A1[ni,loopmatrix[ni,1]] = sign[1]
       A1[ni,loopmatrix[ni,2]] = sign[2] 
    #
    A1 = A1 * maxweig
    #
    for ni in range(pairInfo.shape[0]):
       A2[ni,ni] = 1 * weights[ni]
    #
    A = np.concatenate((A1,A2),axis=0)
    #
    return A,A1,A2
#################################################################
#
def tinfo2datematrix(infpairs):
    #
    # infpairs is an array saving time information of master and slave for each interferogram
    #
    dates  = tinfo2uniquedate(infpairs)
    npair  = infpairs.shape[0]
    ndate  = dates.shape[0]
    tinfomatrix = np.zeros((npair,ndate))
    for pair_ind in range(0,npair):
        for date_ind in range(0,ndate):
            #
            m_ind = np.where(dates == infpairs[pair_ind,0])
            s_ind = np.where(dates == infpairs[pair_ind,1])
            tinfomatrix[pair_ind,m_ind] = 1
            tinfomatrix[pair_ind,s_ind] = -1
    #
    return tinfomatrix
   
##################################################################
def nw_phscor_lsq(A1,A2,D1,D2,weigs,iterations=5):
    #
    ci  = 0
    WA2 = np.copy(A2)
    WD2 = np.copy(D2)
    for ni in range(weigs.shape[0]):
        WA2[ni,:] = WA2[ni,:] * weigs[ni]
    #
    WD2 = WD2 * weigs
    D = np.concatenate((D1,WD2))
    #
    A = np.concatenate((A1*np.max(weigs),WA2))
    #
    while ci < iterations:
       #
       corphs,res,rnk,s = np.linalg.lstsq(A,D)
       simD = np.dot(A2,corphs)
       resD = np.abs(D2 - simD)
       index = np.where(resD == np.max(resD))
       D2[index] = simD[index]
       WD2 = np.copy(D2)
       WD2 = WD2 * weigs
       D = np.concatenate((D1,WD2))
       # print(" NO: %d with RES: %f " %(ci,res))
       ci += 1
    #
    return corphs
##################################################################
def sbas_lsq(A,D,iterations=3,percent=0.99):
    #
    # iterative linear least-square method  
    # 
    ci = 0
    while ci < iterations:
       #
       bperp,res,rnk,s = np.linalg.lstsq(A,D)
       simD = np.dot(A,bperp)
       resD = np.abs(D - simD)
       #
       D[np.where(resD>=np.max(resD)*percent)] = simD[np.where(resD>=np.max(resD)*percent)]
       #
       ci += 1
    #
    # print(resD)
    #
    return bperp,res
#
##################################################################
#
def tinfo2baseline(infpairs,baseline,refindex=1):
    #
    # Calculate approximate baseline for each aquisition....
    tinfomatrix = tinfo2datematrix(infpairs)
    #
    presetbase  = tinfomatrix[0,:] * 0
    presetbase[refindex-1] = 1
    zerobase   = np.zeros(1) 
    #
    # Create coefficient matrix 
    A = np.concatenate((tinfomatrix,[presetbase]),axis=0)
    D = np.concatenate((baseline,zerobase),axis=0)
    # bperp,res,rnk,s = np.linalg.lstsq(A,D) 
    bperp,res = sbas_lsq(A,D,iterations=5)
    #
    return bperp
##################################################################
def optimalmaster(pair_dates,baseline,Tc=10000.,Bc=20000.):
    #
    # Tc a critical Timperal baseline,   days
    # Bc a critical Spatial baseline,    meters
    #
    # return unique dates from a series of InSAR pairs
    #
    udates = tinfo2uniquedate(pair_dates)
    #
    refindex = 1
    #
    pberp    = tinfo2baseline(pair_dates,baseline,refindex=refindex)
    #
    tinfo    = udates * 0
    binfo    = udates * 0
    dims     = udates.shape
    optimalVal = tinfo * 0
    Tc         = float(Tc)
    Bc         = float(Bc)
    #
    for index in range(0,dims[0]):
        #
        tmptime = tinfo * 0
        tmpbase = tinfo * 0
        tmpopts = tinfo * 0
        #
        for cindex in range(0,dims[0]):
            # print(np.absolute(diff2dates(str(udates[cindex]), str(udates[index]))) / Tc)
            tmptime[cindex] = 1. - np.absolute(diff2dates(str(udates[cindex]), str(udates[index]))) / Tc
            tmpbase[cindex] = 1. - np.absolute(pberp[cindex] - pberp[index]) / Bc
            tmpopts[cindex] = tmptime[cindex] * tmpbase[cindex]
        #
        # print(tmptime)
        #
        tinfo[index]      = scipy.stats.gmean(tmptime[np.where(tmptime>0)]) 
        binfo[index]      = scipy.stats.gmean(tmpbase[np.where(tmpbase>0)])
        optimalVal[index] = scipy.stats.gmean(tmpopts[np.where(tmpopts>0)])
    #
    master = udates[np.where(optimalVal==np.max(optimalVal))]
    #
    return master,tinfo,binfo,optimalVal,udates
