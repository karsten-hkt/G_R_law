import os
def gmt_drawing(roi, faults_file_locate, meca_file, outgmt_name, tmp_time, out_gmt_locate):
    gmt_file = out_gmt_locate + outgmt_name + ('%s' % tmp_time) + '.gmt'  # gmtb值文件
    gmt_cpt = out_gmt_locate + outgmt_name + '.cpt'  # cpt色标
    gmt_quake_cpt = out_gmt_locate + 'quake' + '.cpt'
    goStr = open("%s.sh" % gmt_file, 'w')
    goStr.write('gmt begin %s png,pdf\n' % gmt_file)
    goStr.write('gmt makecpt -Cjet -T0.4/1.5/0.1 -Ic -H > %s\n' % gmt_cpt)
    goStr.write('gmt basemap -JM4i -R%s/%s/%s/%s -BENws -Bxaf+l"Long(deg)" -Byaf+l"Lat(deg)" \n' % (
    roi[0], roi[1], roi[2], roi[3]))
    goStr.write('gmt grdimage @earth_relief_30s -I+d -R%s/%s/%s/%s -C../gmt_draw_research_area/grid.cpt\n' % (
    roi[0], roi[1], roi[2], roi[3]))
    goStr.write('gmt plot %s -R%s/%s/%s/%s -JM4i -C%s -W0.25p -L \n' % (
    gmt_file, roi[0], roi[1], roi[2], roi[3], gmt_cpt))
    goStr.write('gmt makecpt -Cinferno -T0/30 -Ic -H > %s\n' %(gmt_quake_cpt))
    goStr.write('gmt colorbar -C%s -DjBR+w2c/0.3c+ml+o0.6c/0.7c+h -Bx+lDepth -By+lkm\n' %(gmt_quake_cpt))
    goStr.write('gmt meca %s -A -Sm0.5c -C%s \n' % (meca_file,gmt_quake_cpt))
    goStr.write('gmt plot %s -W1p,black -JM4i \n' % (faults_file_locate))
    goStr.write('gmt colorbar -C%s -DjBL+w3c/0.3c+ml+o0.3c/0.3c -Bx -By+lb \n' % (gmt_cpt))
    goStr.write('gmt end show')
    goStr.close()
    os.system("sh %s.sh" % (gmt_file))

def gmt_draw_cross_section(ext, m_up_3_name, epi_1, epi_locate, outgmt_name, tmp_time, out_gmt_locate):
    ## 绘图所需要的文件名
    gmt_file = out_gmt_locate + outgmt_name + ('%s' % tmp_time) + '.gmt'  # gmtb值文件
    gmt_cpt = out_gmt_locate + outgmt_name + '.cpt'  # cpt色标
    ## 底图大小
    ysize = 5 * ((ext[3] - ext[2]) / (ext[1] - ext[0]))
    ## 画图
    goStr = open("%s.sh" % gmt_file, 'w')
    goStr.write('gmt begin %s png,pdf\n' % gmt_file)
    goStr.write('gmt makecpt -Cjet -T0.4/1.7/0.1 -Ic -H> %s\n' % gmt_cpt)
    goStr.write(
        'gmt plot %s -R%f/%f/%f/%f -JX5i/-%fi -C%s -W0.01p,white -L -BENws \
        -Bxaf+l"Distance to Epicenter(km)" -Byaf+l"Dep(km)" --FONT_LABEL=8p,Helvetica,black\n' % (
        gmt_file, ext[0], ext[1], ext[2], ext[3], ysize, gmt_cpt))
    goStr.write('gmt plot %s -Ggray -W0.05p -Sc \n' % (m_up_3_name))
    goStr.write('echo 0 %f |gmt plot -Ggray -W0.05p -Sc0.1c -l"Mw>3"\n' % (epi_1[2]))
    goStr.write('echo 0 %f |gmt plot -Gyellow -W0.05p -Sa0.2c -l"Mw7.8"\n' % (epi_1[2]))
    goStr.write('gmt plot %s -Gblue -W0.05p -Sa0.2c -l"large_quake"\n' % (epi_locate))
    goStr.write(
        'gmt colorbar -DJBR+w0.8/0.1+o-0.6c/-0.95c -C%s -Bxaf -By+l"b" --FONT_ANNOT_PRIMARY=3p,Helvetica,black\n' % (
            gmt_cpt))
    #goStr.write('gmt contour %s -W0.1p -A1+d+f5p+o -C1 -l"slip(m)" --FONT_ANNOT_PRIMARY=4.5p,1,black\n ' % (slip_name))
    goStr.write('gmt legend -DjBL+o0.01c/0.01c --FONT_ANNOT_PRIMARY=4.5p,1,black\n')
    goStr.write('gmt end show')
    goStr.close()
    os.system('sh %s.sh' % gmt_file)

def gmt_draw_cross_section_3d(roi, epi_1, outgmt_name, tmp_time, out_gmt_locate, p):
    gmt_file = out_gmt_locate + outgmt_name + ('%s' % tmp_time) + '.gmt'  # gmtb值文件
    gmt_cpt = out_gmt_locate + outgmt_name + '.cpt'  # cpt色标
    goStr = open("%s.sh" % gmt_file, 'w')
    goStr.write('gmt begin %s png,pdf\n' % gmt_file)
    goStr.write('gmt set MAP_FRAME_TYPE plain\n')
    goStr.write('gmt set MAP_GRID_PEN 0p,gray,-\n')
    goStr.write('gmt basemap -R%s/%s/%s/%s/%s/%s -JX10c/10c \
                -JZ3c -Bxa1fg -Bya1fg -Bza10fg+l"Depth (km)" -BwSEnZ -p%s/%s\n'
                % (roi[0], roi[1], roi[2], roi[3], roi[4], roi[5], p[0], p[1]))
    goStr.write('gmt makecpt -Cjet -T0.4/1.7/0.1 -Ic -H> %s\n' % gmt_cpt)
    goStr.write('gmt plot3d %s -C%s -L -W0p,gray -p\n' %(gmt_file, gmt_cpt))
    goStr.write('echo %f %f %f |gmt plot3d -Gyellow -W0.05p -Sa0.2c\n' % (epi_1[0], epi_1[1], epi_1[2],))
    goStr.write('gmt set FONT_ANNOT_PRIMARY 7p FONT_LABEL 8p\n')
    goStr.write('gmt set FONT_ANNOT_PRIMARY 7p FONT_LABEL 8p\n')
    goStr.write('gmt colorbar -C%s -DjBL+w2.5c/0.3c+h+o1.3c/2c+ml -Ba0.5 -Bx+l"b value"\n'%(gmt_cpt))
    goStr.write('gmt end show\n')
    goStr.close()
    os.system('sh %s.sh' % gmt_file)

