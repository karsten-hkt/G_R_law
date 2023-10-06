#! /usr/bin/env python
#
# A python module, which was re-organized from pDATA.py and pGAMMA.py 
# this module will be focused on RADARSAT-1/2 and RADARSAT Constellation 
# Developed by Wanpeng Feng, @NRCan, 2017-06-01
# wanpeng.feng@hotmail.com
#
###############################################################################
import numpy as np
import os
import sys
import glob
from   lxml import etree
import xml.dom.minidom
import zipfile
import shutil
from datetime import datetime
import pSAR
import pGAMMA
#
#
###############################################################################
def dir2newname(cdir):
    #
    outname = None
    cdate   = None
    cxml = glob.glob(cdir+'/product.xml')
    if len(cxml)>0:
        cdate = rs2date(cxml[0])
        beam  = rs2xml(cxml[0],inkeyword='beam')
        track = rs2relativeOrb(cxml[0])
        pdir  = rs2xml(cxml[0],inkeyword='pointing')
        outname = beam + '_' + str(track) + '_' + pdir
    return outname,cdate
#
###############################################################################
#
def rsatunzip(czip):
    #
    cxml = rsatzip2xml(czip)
    odir = os.path.dirname(czip)
    unzip_dir = os.path.basename(czip).split('.')[0]
    #
    if not os.path.exists(cxml):
         zip2xml(czip)
    #
    ctime = rs2timeinfo(cxml)
    cdate = ctime[0:8]
    slcs = glob.glob(odir+'/'+cdate+'*.slc')
    #
    if len(slcs) < 1:
        #
        goflag = True
    else:
        goflag = False
    #
    if goflag:
      #
      newunzipdir = os.path.join(odir,unzip_dir)
      if not os.path.exists(newunzipdir):
         os.makedirs(newunzipdir)
         #
         unzip_command = ' unzip %s -d %s > %s/log' % (czip,odir,newunzipdir)
         print(unzip_command)
         os.system(unzip_command)
    #
    if os.path.exists(newunzipdir):
        return True
    else:
        return False
###############################################################################
def rs2updateingorb(orb,orbdir=""):
    #
    # update from MDA automatically...
    updating_orb = 'pSAR_rs2orbit.py %s' % orb.lower()
    print(' pGAMMA: %s to %s' % (updating_orb,orbdir))
    os.system(updating_orb)
    #
    if os.path.exists(orbdir+'/'+orb.lower()):
        return orbdir+'/'+orb.lower()
    else:
        return None
#
###############################################################################
#
def rs2dir2orb(in_dir,orbdir=None):
    #
    orb = None
    xml = glob.glob(in_dir+'/*/product.xml')
    if len(xml)>0:
        orb = rs2xml(xml[0],inkeyword='orbitfile')
    #
    if 'RS2_ORB' in os.environ and orbdir is None:
        orbdir = os.environ['RS2_ORB']
    if (orb is not None and orbdir is not None):
        # 
        if 'PRED' in orb.upper():
            orb = orb.split('_')[0]+'_def.orb'
        #
        orbfile = glob.glob(orbdir+'/'+orb.lower())
        goflag = False
        if len(orbfile)>0:
            orbtmp = orbfile[0]
            if os.stat(orbtmp).st_size > 5:
                orb = orbtmp
                goflag = True
        if not goflag:
           orb = rs2updateingorb(orb,orbdir=orbdir)
    #
    return orb
#
###############################################################################
#
def slcs2slcpars(slcs):
    return [cslc+'.par' for cslc in slcs]
#
###############################################################################
#
def filterslcs(slcs):
    #
    slcs = np.array(slcs,dtype='str')
    nslcs = slcs.shape[0]
    #
    xmls,flag = slcs2productxmls(slcs)
    mtimes = np.zeros([xmls.shape[0],3])
    #
    #
    for ni in range(nslcs):
        #
        start_time = rs2xml(xmls[ni],inkeyword='start_time')
        end_time   = rs2xml(xmls[ni],inkeyword='end_time')
        Ts_time    = pSAR.ts.timestr2jd(start_time,fmt='%Y-%m-%dT%H:%M:%S.%fZ')
        Te_time    = pSAR.ts.timestr2jd(end_time,fmt='%Y-%m-%dT%H:%M:%S.%fZ')
        #
        diff_t     = (Ts_time - Te_time) * 24. * 60. * 60.
        #
        if diff_t > 0:
           mtimes[ni,:] = [Te_time,Ts_time,diff_t]
        else:
           mtimes[ni,:] = [Ts_time,Te_time,diff_t * -1.]
    #
    starttimes = mtimes[:,0]
    endtimes   = mtimes[:,1] 
    #
    index      = np.argsort(starttimes)
    endtimes   = endtimes[index]
    starttimes = starttimes[index]
    newslcs    = slcs[index]
    #
    #
    flag = np.zeros(nslcs)
    #
    for ni in range(nslcs):
        #
        if ni < 1:
           flag[ni] = 1
           latesttime = endtimes[ni]
        else:
           difftime = (endtimes[ni] - latesttime) * 3600. * 24.
           if difftime > 0.:
             flag[ni] = 1
             latesttime = np.vstack([latesttime,endtimes[ni]])
             latesttime = np.max(latesttime)
             #  
    #
    outslcs = newslcs[flag>0]
    return outslcs
#
###############################################################################
#    
def slcs2productxmls(slcs):
    #
    slcs = np.array(slcs,dtype='str')
    flag = np.zeros(slcs.shape[0])
    xmls = np.copy(slcs)
    for ni in range(slcs.shape[0]):
        cdir = os.path.dirname(slcs[ni])
        cxml = glob.glob(cdir+'/product.xml')
        if len(cxml) < 1:
            flag[ni] = -1
            xmls[ni] = 'NULL'
        else:
            xmls[ni] = cxml[0]
    return xmls,flag
#
###############################################################################
#
def dir2SLCcat(cdir,orbdir=None):
    #
    slcs = glob.glob(cdir+'/*/*.slc')
    print(' %s : %d ' % (cdir,len(slcs)))
    if len(slcs) == 1:
       #
       slclinks(slcs[0],cdir)
    if len(slcs) > 1:
       #
       # sorting SLCs based on the time information
       #
       slcs = filterslcs(slcs)
       #
       if slcs.shape[0] < 2:
          fullslc = os.path.abspath(slcs[0])
          newslc = cdir+'/'+os.path.basename(fullslc)
          fullslc_par = fullslc+'.par'
          newslc_par  = newslc +'.par'
          if not os.path.exists(newslc):
            os.symlink(fullslc,newslc)
          if not os.path.exists(newslc_par):
            os.symlink(fullslc_par,newslc_par)
       else:
         str_slcs = ' '.join(slcs)
         outslc = cdir+'/'+os.path.basename(slcs[0])
         #
         if not os.path.exists(outslc):
           orbitfile = rs2dir2orb(cdir,orbdir=orbdir)
           if not os.path.exists(orbitfile):
              orb_opt = ''
           else:
              orb_opt = ' -orb %s' % orbitfile   
           cmd = 'rs2_cat.sh %s -outslc %s %s' % (str_slcs,outslc,orb_opt)
           print(" pGAMMA: %s" % cmd)
           os.system(cmd)
           
def slclinks(cslc,outdir):
    #
    inslc = os.path.abspath(cslc)
    slcname = os.path.basename(inslc)
    newslc  = outdir+'/'+slcname
    newslc_par = newslc+'.par'
    inslc_par = inslc+'.par'
    if not os.path.exists(newslc):
       os.symlink(inslc,newslc)
    if not os.path.exists(newslc_par):
       os.symlink(inslc_par,newslc_par)
    return True
#
def invpol(inpol):
    if inpol.upper() == 'HH':
        return 'VV'
    else:
        return 'HH'
#
def rs2_dir2SLC(cdir,lutxml='Beta',pol='HH',orbdir=None):
    '''
    Create SLC from a local folder. 
    '''
    pol = pol.upper()
    #
    productxml= glob.glob(cdir+'/product.xml')
    geotiff   = glob.glob(cdir+'/imagery_%s.tif' % pol)
    clutxml   = glob.glob(cdir+'/lut%s.xml' % lutxml)
    #
    if len(productxml) < 1:
       print(' No product.xml found in %s' % cdir)
       return False
    if len(geotiff) < 1:
       print(' No imagery_%s.tif found in %s' % (pol,cdir))
       pol = invpol(pol)
       geotiff   = glob.glob(cdir+'/imagery_%s.tif' % pol)
       if len(geotiff) < 1:
           print(' No imagery_%s.tif found in %s either' % (pol,cdir))
           return False
    #
    if len(clutxml) < 1:
       print(' No lut%s.xml found in %s' % (lutxml,cdir))
       return False
    #
    cdate = rs2date(productxml[0])
    orbfile = rs2xml(productxml[0],inkeyword='orbitfile')
    if orbdir is None or not os.path.exists(orbdir):
        orbdir = os.environ['RS2_ORB']
    fullorb = orbdir+'/'+orbfile.lower()
    #
    outslc = cdir+'/'+cdate+'_'+pol+'.slc'
    outslc_par = outslc+'.par'
    #
    cmd = 'par_RSAT2_SLC %s %s %s %s %s %s ' % \
           (productxml[0],clutxml[0],geotiff[0],pol,outslc_par,outslc)
    #      
    print(' ++++++++++++++')       
    print(' + pRSAT: %s' % cmd)
    print(' ++++++++++++++')
    if not os.path.exists(outslc):
       os.system(cmd)
    else:
       print(' %s found already. Delete it before re-run.' % outslc)
    #
    if os.path.exists(fullorb):
       pGAMMA.gamma_updateRS2orb(outslc_par,fullorb,nstat=64)
    #
    if (os.path.exists(outslc) and os.path.exists(outslc_par)):
        #
        return True
    else:
        return False
#
###############################################################################
#         
def rsatzip2xml(czip):
    #
    outxml = czip.split('.zip')[0]+'.xml'
    if not os.path.exists(outxml):
       zip2xml(czip)
    return outxml
#
###############################################################################
def rsatzip2local_dir(indir,disp=False):
    #
    zips = glob.glob(indir+'/*.zip')
    #
    newzips = []
    for czip in zips:
        cnewzip = rsatzip2local(czip)
        if disp:
           print("   %s : %s" % (czip,cnewzip))
        newzips.append(cnewzip)
    return cnewzip
#
#
def rsatzip2local(czip):
    #
    cxml = rsatzip2xml(czip)
    if not os.path.exists(cxml):
       zip2xml(czip)
    #
    mission     = rs2satellite(cxml)
    beam        = rs2xml(cxml,inkeyword='beam')
    tracknum    = rs2relativeOrb(cxml)
    cdate       = rs2date(cxml)
    pointingdir = rs2xml(cxml,inkeyword='pointing')
    #
    # local_top   = mission+'_RAW_DIR/'
    #
    second_dir  = mission+'_'+beam+'_'+str(tracknum)+'_'+pointingdir 
    zipnewhome  = second_dir + '/SLC/' + cdate
    #
    if not os.path.exists(zipnewhome):
       os.makedirs(zipnewhome)
    #
    czip = os.path.abspath(czip)
    newzip = zipnewhome+'/'+os.path.basename(czip)
    if not os.path.exists(newzip):
        os.symlink(czip,newzip)
    #
    return newzip
#
###############################################################################
def zipdir2local(zip_dir):
    '''
    To link zip(s) to local directory. A specific folder will be assigned based
    on the data information.
    by Wanpeng Feng, @CCRS/NRCan, 2017-07-14
    
    NOTE: rsatzip2local_dir is identical to this function....
    '''
    #
    outnames = []
    #
    zips = glob.glob(zip_dir+'/*.zip')
    for czip in zips:
        #
        newzip = rsatzip2local(czip)
        #print(" %s : %s" % (czip,newzip))
        outnames.append(newzip)
    #
    return outnames
#        
###############################################################################        
#        
def zip2xml(czip,disp=True):
    #
    '''
    Return the product xml into the local folder with a new name...
    by Wanpeng Feng, @CCRS/NRCan, 2017-07-14
    '''
    
    try:
      #
      zf            = zipfile.ZipFile(czip)
      #print(" OK reading: %s" % czip)
      zip_files     = zf.namelist()
      zip_dir       = os.path.dirname(czip)
      #
      for cname in zip_files:
        #
        root_name    = os.path.basename(cname)
        dir_name     = os.path.dirname(cname)
        cur_dir_name = os.getcwd()
        #
        if root_name == "product.xml":
           c_dir_name   = os.path.basename(dir_name);
           c_root_name,ext = os.path.splitext(c_dir_name)
           output_xml = c_root_name + ".xml"
           #  print output_manifest
           zf.extract(cname)
           out_xml = os.path.join(zip_dir,output_xml)
           #
           shutil.move(cname,out_xml)
           if cur_dir_name != dir_name:
              shutil.rmtree(dir_name,ignore_errors=False)
    except:
      print(" ERROR: %s" % czip)
      out_xml = None
    if disp:
       print(" %s : %s" % (czip,out_xml))
    #
    return out_xml
###############################################################################
#     
def bytestoutf8(t):
    '''
    To convert a bytes-like variable to utf8
    '''
    if isinstance(t,bytes):
        return t.decode("utf-8")
    return t
#
###############################################################################
def rs2inROI(xml,lons,lats):
    #
    # xml is product.xml in the released RADARSAT-2 data folder
    #
    lons = np.array(lons)
    lats = np.array(lats)
    #
    minlon,maxlon = np.min(lons),np.max(lons)
    minlat,maxlat = np.min(lats),np.max(lats)
    #
    corners         = rs2corners(xml)
    meanlon,meanlat = np.mean(corners[:,0]),np.mean(corners[:,1])
    #
    # print(meanlon,meanlat)
    #
    if (meanlon >= minlon and meanlon <= maxlon and \
        meanlat >= minlat and meanlat <= maxlat):
        return True
    else:
        return False
###############################################################################
def rs2timeinfo(inxml):
    #
    timeinfo = rs2xml(inxml,inkeyword='time')
    
    return timeinfo.replace("-","").replace(":","")
    #
#
###############################################################################
def rs2relativeOrb(cxml):
    #
    toltracks = 343.
    orbitname = rs2xml(cxml,inkeyword='orbitfile')
    cnum = float(orbitname.split('_')[0])
    return int((cnum / toltracks - \
               int(cnum / toltracks)) * 100000)
#
###############################################################################
#
def rs2date(inxml):
    times = rs2xml(inxml,inkeyword='time')
    return times.replace('-','')[0:8]
###############################################################################
def rs2xml(inxml,inkeyword='pointing'):
    #
    if inxml is None:
       helpstr = \
       '''
           output = rs2xml(inxlm,inkeyword='pointing')
           +++++++++++++++++++++++++++++++++++++++++++
           + Inkeyword:
             sensor,     satellite
             pointing,   atennapointing
             beam,       beam mode
             orbitfile,  orbitDataFile
             passdir,    passDirection
             time,       acquisition time information
           + Return:
             value corresponding to the given keyword from the xml file
       '''
       print(helpstr)
       return None
       sys.exit(-1)
    #
    output     = None
    conts      = []
    DOMTree    = xml.dom.minidom.parse(inxml)
    collection = DOMTree.documentElement
    #
    if inkeyword.upper() == "LINES":
       tags       = collection.getElementsByTagName('numberOfLines')
       # 
       if len(tags) > 0:
          output = '%s' % tags[0].childNodes[0].data
                              
    if inkeyword.upper() == "PRF":
       tags       = collection.getElementsByTagName('pulseRepetitionFrequency')
       # 
       if len(tags) > 0:
          output = '%s' % tags[0].childNodes[0].data
    #
    if (inkeyword.upper() == "TIME" or inkeyword.upper() == "START_TIME"):
       tags       = collection.getElementsByTagName('sarProcessingInformation')
       # output     = None
       if len(tags) > 0:
          conts = tags[0].getElementsByTagName('zeroDopplerTimeFirstLine')
       if len(conts) > 0:
          output = ("%s" % conts[0].childNodes[0].data)
    if (inkeyword.upper() == "END_TIME"):
       tags       = collection.getElementsByTagName('sarProcessingInformation')
       # output     = None
       if len(tags) > 0:
          conts = tags[0].getElementsByTagName('zeroDopplerTimeLastLine')
       if len(conts) > 0:
          output = ("%s" % conts[0].childNodes[0].data) 
    #
    if inkeyword.upper() == "SENSOR":
       tags       = collection.getElementsByTagName("sourceAttributes")
       # output     = None
       if len(tags) > 0:
          conts = tags[0].getElementsByTagName("satellite")
       if len(conts) > 0:
          output = ("%s" % conts[0].childNodes[0].data)
    #    
    if inkeyword.upper() == "POINTING":
       tags       = collection.getElementsByTagName("radarParameters")
       # output     = None
       if len(tags) > 0:
          conts = tags[0].getElementsByTagName("antennaPointing")
       if len(conts) > 0:
          output = ("%s" % conts[0].childNodes[0].data)
    #      
    if inkeyword.upper() == "BEAM":
       tags       = collection.getElementsByTagName("radarParameters")
       # output     = None
       if len(tags) > 0:
          conts = tags[0].getElementsByTagName("beams")
       if len(conts) > 0:
          output = ("%s" % conts[0].childNodes[0].data)      
    #
    if inkeyword.upper() == "PASSDIR":
       tags       = collection.getElementsByTagName("orbitInformation")
       # output     = None
       if len(tags) > 0:
          conts = tags[0].getElementsByTagName("passDirection")
       if len(conts) > 0:
          output = ("%s" % conts[0].childNodes[0].data)       
    if inkeyword.upper() == "ORBITFILE":
       tags       = collection.getElementsByTagName("orbitInformation")
       # output     = None
       if len(tags) > 0:
          conts = tags[0].getElementsByTagName("orbitDataFile")
       if len(conts) > 0:
          output = ("%s" % conts[0].childNodes[0].data) 
    return output
############################################################################### 
def rs2duration(inxml):
    #
    DOMTree    = xml.dom.minidom.parse(inxml)
    collection = DOMTree.documentElement
    tags       = collection.getElementsByTagName("zeroDopplerTimeFirstLine")
    t1 = tags[0].childNodes[0].data
    tags       = collection.getElementsByTagName("zeroDopplerTimeLastLine")
    t2 = tags[0].childNodes[0].data
    t1 = datetime.strptime(t1, '%Y-%m-%dT%H:%M:%S.%fZ')
    t2 = datetime.strptime(t2, '%Y-%m-%dT%H:%M:%S.%fZ')
    return (t1-t2).seconds
#
def rs2satellite(inxml):
    '''
    Return <mission> of the input data...
    by Wanpeng Feng, @CCRS/NRCan, 2017-07-14
    
    '''
    #
    DOMTree    = xml.dom.minidom.parse(inxml)
    collection = DOMTree.documentElement
    tags       = collection.getElementsByTagName("satellite")
    mission    = None
    if len(tags)>0:
        mission = tags[0].childNodes[0].data
    return mission
#
###############################################################################
#
def rs2tiepointsv2(inxml):
    '''
    A new version to return geolocation of rs2.xml after python 3.x is used.
    
    '''
    DOMTree    = xml.dom.minidom.parse(inxml)
    collection = DOMTree.documentElement
    tags       = collection.getElementsByTagName("imageTiePoint")
    outdata    = np.zeros([len(tags),4])
    for ni in range(len(tags)):
        line  = tags[ni].getElementsByTagName("line")[0].childNodes[0].data
        pixel = tags[ni].getElementsByTagName("pixel")[0].childNodes[0].data
        lat   = tags[ni].getElementsByTagName('latitude')[0].\
                    childNodes[0].data
        lon   = tags[ni].getElementsByTagName('longitude')[0].\
                    childNodes[0].data 
        outdata[ni,:] = [float(line),float(pixel),float(lat),float(lon)]   
    return outdata[:,0],outdata[:,1],outdata[:,3],outdata[:,2]
#
###############################################################################
#
def rs2tiepoints(inxml):
    '''
    A version with python2.x
    '''
    inTree = etree.parse(inxml)
    uri    = inTree.getroot().nsmap[None]
    NS     = '{'+uri+'}'
    querylist = {\
             'geolocationline'     : \
             './/%simageAttributes//%sgeographicInformation//'+\
             '%sgeolocationGrid//%simageTiePoint//%sline' % (NS,NS,NS,NS,NS),\
             'geolocationpixel'    : \
             './/%simageAttributes//%sgeographicInformation//'+\
             '%sgeolocationGrid//%simageTiePoint//%spixel' % (NS,NS,NS,NS,NS),\
             'geolocationlat'      : \
             './/%simageAttributes//%sgeographicInformation//'+\
             '%sgeolocationGrid//%simageTiePoint//%slatitude' %\
             (NS,NS,NS,NS,NS),\
             'geolocationlon'      : \
             './/%simageAttributes//%sgeographicInformation//'+\
             '%sgeolocationGrid//%simageTiePoint//%slongitude' \
             % (NS,NS,NS,NS,NS)}
    #
    for key in querylist.keys():
    #
      try:
         vars()[key];                     # check if non-existed
      except KeyError or NameError:
         vars()[key] = [];                # if non-exist, initialize

      for nodes in inTree.findall(querylist[key]):
        #
        if key == 'geolocationline':
            # vars()[key].append(nodes.text)
            vars()[key].append(int(float(nodes.text)))
            line  = vars()[key]
        elif key == 'geolocationpixel':
            #
            vars()[key].append(int(float(nodes.text)))
            pixel = vars()[key]
            #
        elif key == 'geolocationlat':
            vars()[key].append(float(nodes.text))
            lat   = vars()[key]
        elif key == 'geolocationlon':
            vars()[key].append(float(nodes.text))
            lon   = vars()[key]
    #
    return np.array(line),np.array(pixel),np.array(lon),np.array(lat)

###############################################################################    
def rs2corners(inxml):
    '''
    To return corners of a RADARSAT-2 SAR SLC data based on its product.xml
    
    '''
    #
    line,pixel,lon,lat = rs2tiepointsv2(inxml)
    #
    flag_line_min  = line  == 0
    flag_line_max  = line  == np.amax(line)
    flag_pixe_min  = pixel == 0
    flag_pixe_max  = pixel == np.amax(pixel)
    flag_up_left   = flag_line_min & flag_pixe_min
    flag_bot_left  = flag_line_min & flag_pixe_max
    flag_bot_right = flag_line_max & flag_pixe_max
    flag_up_right  = flag_line_max & flag_pixe_min
    lon1,lat1      = lon[flag_up_left],   lat[flag_up_left]
    lon2,lat2      = lon[flag_bot_left],  lat[flag_bot_left]
    lon3,lat3      = lon[flag_bot_right], lat[flag_bot_right]
    lon4,lat4      = lon[flag_up_right],  lat[flag_up_right]
    #
    lon1           = lon1[0]
    lat1           = lat1[0]
    lon2           = lon2[0]
    lat2           = lat2[0]
    lon3           = lon3[0]
    lat3           = lat3[0]
    lon4           = lon4[0]
    lat4           = lat4[0]
    #
    corners        = np.array([[lon1,lon2,lon3,lon4,lon1],\
                               [lat1,lat2,lat3,lat4,lat1]]).T
    return corners  

