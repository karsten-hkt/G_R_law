import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import scipy
import area_cut
import inpoly

## 需要地震序列(原始，不用做mag剔除的数据)，b值，震级所在行
def b_error_shi(aftershocks_file,Mc,b,mag_line):
    aftershocks_file_up_Mc = aftershocks_file[aftershocks_file[:,mag_line]>Mc]
    mag_mean = np.mean(aftershocks_file_up_Mc[:,mag_line])
    N = aftershocks_file_up_Mc.shape[0]
    mag_sigma_sum = 0
    for i in range(N):
        mag_sigma_sum = mag_sigma_sum + (aftershocks_file_up_Mc[i,mag_line]-mag_mean)**2
    if mag_sigma_sum == 0:
        return np.inf
    else:
        b_error = 2.30*b**2*math.sqrt(mag_sigma_sum/(N*(N-1)))
        return b_error

## 极大似然法计算b值
## 需要传入地震序列(原始，不用做mag剔除的数据)，完备震级，震级差，震级所在行
def b_calculate(aftershocks_file,Mc,delta_M,mag_line):
    ## 剔除完备震级以下的数据
    aftershocks_file_up_Mc = aftershocks_file[aftershocks_file[:,mag_line]>Mc]
    N = aftershocks_file_up_Mc.shape[0]
    if N > 0:
        M_mean = np.mean(aftershocks_file_up_Mc[:,mag_line])
        b = math.log10(math.exp(1))/(M_mean-Mc)
        b_error = b_error_shi(aftershocks_file,Mc,b,mag_line)
        return b,b_error
    else:
        return np.inf,np.inf

## 计算区域内部的b值——完备震级使用fmd方法
def get_b_Mc(INdata,minnum,mbin,mag_line):
    tmp_INdata = INdata
    if tmp_INdata.shape[0]<minnum:
        return np.NAN,np.NAN,np.NAN
    Mc = get_maxc(tmp_INdata[:,mag_line], mbin =mbin)
    b,b_error = b_calculate(tmp_INdata,Mc,mbin,mag_line)
    return b,b_error,Mc

## 计算完备震级的另外一种方式，
def get_maxc(mag, mbin):
    this_fmd = fmd(mag, mbin) # FMD
    maxc = this_fmd[0][np.argmax(this_fmd[1])] # Mag bin with highest no. of events
    return round(maxc, 2)+0.2

def fmd(mag, mbin):
    minmag = math.floor(min(mag/mbin)) * mbin # Lowest magnitude bin
    maxmag = math.ceil(max(mag/mbin)) * mbin # Highest magnitude bin
    mi = np.arange(minmag, maxmag + mbin, mbin) # Sequence of magnitude bins
    nbm = len(mi) # No. of magnitude bins
    cumnbmag = np.zeros(nbm) # Pre-allocate array for cumulative no. of events in mag bin and higher

    # Get cumulative no. of events in mag bin and higher
    for i in range(nbm):
        cumnbmag[i] = np.where(mag > mi[i] - mbin/2)[0].shape[0]

    # Get no. of events in each mag bin:
    nbmag = abs(np.diff(np.append(cumnbmag, 0)))

    return mi, nbmag, cumnbmag # Return magnitude bins, no. of events in bin, and cumulative no. of events
#

def a_calculate(Mc,b,GR_Mc_N_counts,delta_M):
    Mc_N_now = GR_Mc_N_counts[GR_Mc_N_counts[:,0]==np.round(Mc,1),1]
    a = np.log10(Mc_N_now) + b*Mc
    return a
#
def count_GR_Mc_N(aftershocks_file,Mc,delta_M,mag_line):
    mag_max = np.round(np.max(aftershocks_file[:,mag_line]),2)
    mag_min = np.round(np.floor(np.min(aftershocks_file[:,mag_line])),2)
    mag_total = int(np.ceil((mag_max-mag_min)/delta_M))
    GR_Mc_N_counts = np.zeros((mag_total,2))
    for i in range(mag_total):
        mag_now = mag_min+delta_M*i
        GR_Mc_N_counts[i,0] = np.round(mag_now,1)
        aftershocks_file_now = aftershocks_file[aftershocks_file[:,mag_line]>=mag_now]
        GR_Mc_N_counts[i,1] = aftershocks_file_now.shape[0]
    return GR_Mc_N_counts
#
def FMD_drawing(aftershocks_file,Mc,delta_M,mag_line,b,tmp_time,b_error):
    GR_Mc_N_counts = count_GR_Mc_N(aftershocks_file,Mc,delta_M,mag_line)
    a = a_calculate(Mc,b,GR_Mc_N_counts,delta_M)
    x_gr = GR_Mc_N_counts[:,0]
    x_gr = x_gr[x_gr[:]>=Mc]
    y_gr = a-b*x_gr
    plt.scatter(GR_Mc_N_counts[:,0],GR_Mc_N_counts[:,1],marker='o',s=6)
    plt.vlines([Mc],1,max(10**y_gr[:]),linestyles='dashed',color='black')
    plt.plot(x_gr,10**y_gr,label='%s)FMD'%(chr(97+tmp_time)))
    plt.text(Mc,max(10**y_gr[:]),r'%.2f$\pm$%.3f'%(b,b_error),fontsize=10)

def month_b_value_calculate_and_draw(data, minnum, delta_M, mag_line, color):
    data_b, data_b_error, data_Mc = get_b_Mc(data, minnum, delta_M, mag_line)
    GR_Mc_N_counts = count_GR_Mc_N(data, data_Mc, delta_M, mag_line)
    ## 计算a值
    data_a = a_calculate(data_Mc, data_b, GR_Mc_N_counts, delta_M)
    ## GR_曲线
    data_x_gr = GR_Mc_N_counts[:, 0]
    data_x_gr = data_x_gr[data_x_gr[:] >= data_Mc]
    data_y_gr = data_a - data_b * data_x_gr
    plt.scatter(GR_Mc_N_counts[:, 0], GR_Mc_N_counts[:, 1], marker='o', c=color, s=6, label='mag_sum')
    plt.vlines([data_Mc], 1, max(10 ** data_y_gr[:]), color=color, linestyles='dashed', label='Mc')
    plt.plot(data_x_gr, 10 ** data_y_gr, color=color, label='G-R line')
    plt.text(data_Mc, max(10 ** data_y_gr[:]), r'%.2f$\pm$%.3f' % (data_b, data_b_error), color=color, fontsize=15)
    return data_b, data_b_error, data_Mc

## 完备震级计算
def fmbass(aftershocks_file,delta,mag_line):
    min_mag = min(aftershocks_file[:,mag_line])
    max_mag = max(aftershocks_file[:,mag_line])
    n = math.floor((-(min_mag-delta/2)+max_mag+delta/2)/delta)
    xc = np.zeros(n-1)
    x = np.zeros(n)
    log_nc = np.zeros(n-1)
    mag_count = np.zeros(n-1)
    for i in range(n):
        x[i] = round(min_mag + delta*i,2)
    xc = x[0:n-1]

    sum_mag_count = np.zeros(n-1)
    sum_ = 0
    for i in range(n-1):
        low_now = round(xc[i]-delta/2,2)
        up_now = round(low_now+delta,2)
        mag_count_file = aftershocks_file[aftershocks_file[:,mag_line]>low_now,:]
        mag_count_file = mag_count_file[mag_count_file[:,mag_line]<=up_now,:]
        mag_count[i] = mag_count_file.shape[0]
        sum_ = sum_ + mag_count[i]
        sum_mag_count[i] = sum_
        log_nc[i] = math.log10(1/delta* (aftershocks_file.shape[0]-sum_mag_count[i]) * delta)

    log_n = np.zeros(n)
    mag_count_2 = np.zeros(n)
    for i in range(n):
        low_now = round(x[i]-delta/2,2)
        up_now = round(low_now+delta,2)
        mag_count_file = aftershocks_file[aftershocks_file[:,mag_line]>low_now,:]
        mag_count_file = mag_count_file[mag_count_file[:,mag_line]<=up_now,:]
        mag_count_2[i] = mag_count_file.shape[0]
        if mag_count_2[i] == 0:
            x[i] = np.inf
            log_n[i] = np.inf
        else:
            log_n[i] = math.log10((1/delta)*mag_count_2[i]*delta)
    x = x[x[:] != np.inf]
    log_n = log_n[log_n[:] != np.inf]
    sl = np.diff(log_n)/np.diff(x)#segment slopes
    xsl = x[2:x.shape[0]]
    niter = 3
    tau = np.zeros(niter)
    pva = np.zeros(niter)
    N = sl.shape[0]
    j = 0 #iterations
    k = 0 #discontinuities
    SA = np.zeros(N)
    sl_pd = pd.Series(sl.round(10),dtype='float64')#in order to use the pandas function rank
    while j<niter:
        sl_sort = sl_pd.rank(method='first')
        for i in range(0,N,1):
            SA[i] = abs(2*sum(sl_sort[0:i]) -i*(N+1))
        SA[0:N-1]=SA[1:N]
        SA[N-1]=0
        N_argSA_len = np.argsort(SA).shape[0]
        SA_judge = np.argsort(SA)[N_argSA_len-1]
        n1 = np.where(SA == SA[SA_judge])[0]
        xn1 = sl[0:n1[0]+1]
        xn2 = sl[n1[0]+1:]
        stat,p = scipy.stats.ranksums(xn1, xn2)
        if (n1[0]>2 and n1[0]<=N-2) and p<0.05:
            pva[j] = p
            tau[j] = n1[0]
            if (k>0):
                medsl1 = np.median(sl[1:n0+1])
                medsl2 = np.median(sl[n1[0]+1:])
                for i in range(0, n0, 1):
                    sl[i] = sl[i]+medsl1
                for i in range(n0, sl.shape[0], 1):
                    sl[i] = sl[i]+medsl2
            medsl1 = np.median(sl[1:n1[0]+1])
            medsl2 = np.median(sl[n1[0]+1:])
            for i in range(0, n1[0], 1):
                sl[i] = sl[i]-medsl1
            for i in range(n1[0], sl.shape[0], 1):
                sl[i] = sl[i]-medsl2
            n0 = n1[0]
            k = k+1
        else:
            pva[j] = p
            tau[j] = n1[0]
            break
        j = j+1
    pva_pd = pd.Series(pva.round(10),dtype='float64')
    ip = pva_pd.rank(method='first')
    m0= round(xsl[int(tau[int(ip[0]-1)]-1)],2) #use the p min
    p = pva[int(ip[0]-1)]
    return m0,p

## 计算区域内部的b值——完备震级使用fmbass方法
def eqs2bINpoly_fmbass(INdata,polys,minnum,mbin,lon_line,mag_line,dep_line,m_total):
    IN,ON = inpoly.inpoly2(INdata[:,[lon_line,dep_line]], polys)
    tmp_INdata = INdata[IN==1,:]
    if tmp_INdata.shape[0]<minnum:
        return np.NAN,np.NAN
    m0,p0 = fmbass(tmp_INdata, mbin, mag_line)
    if p0 > 0.05:
        Mc = m_total
    else:
        Mc = m0
    b,b_error = b_calculate(tmp_INdata,Mc,mbin,mag_line)
    return b,b_error

## 计算区域内部的b值——完备震级使用fmd方法
def eqs2bINpoly_maxc(INdata,polys,minnum,mbin,lon_line,mag_line,dep_line,refpolys,tmp_time):
    if len(INdata) == 0:
        return np.NAN,np.NAN
    IN,ON = inpoly.inpoly2(INdata[:,[lon_line,dep_line]], polys)
    tmp_INdata = INdata[IN==1,:]
    if tmp_INdata.shape[0]<minnum:
        return np.NAN,np.NAN
    Mc = get_maxc(tmp_INdata[:,mag_line], mbin =mbin)
    b,b_error = b_calculate(tmp_INdata,Mc,mbin,mag_line)
    #绘制震中处的FMD
    epi_1=[0,2.1]
    if area_cut.is_in_poly(epi_1, refpolys) and tmp_time<5:
        FMD_drawing(tmp_INdata,Mc,mbin,mag_line,b,tmp_time,b_error)
    return b,b_error

def FMD_drawing(aftershocks_file,Mc,delta_M,mag_line,b,tmp_time,b_error):
    GR_Mc_N_counts = count_GR_Mc_N(aftershocks_file,Mc,delta_M,mag_line)
    a = a_calculate(Mc,b,GR_Mc_N_counts,delta_M)
    x_gr = GR_Mc_N_counts[:,0]
    x_gr = x_gr[x_gr[:]>=Mc]
    y_gr = a-b*x_gr
    plt.scatter(GR_Mc_N_counts[:,0],GR_Mc_N_counts[:,1],marker='o',s=6)
    plt.vlines([Mc],1,max(10**y_gr[:]),linestyles='dashed',color='black')
    plt.plot(x_gr,10**y_gr,label='%s)FMD'%(chr(97+tmp_time)))
    if tmp_time==3:
        plt.text(Mc,max(10**y_gr[:])-17,r'%.2f$\pm$%.3f'%(b,b_error),fontsize=10)
    elif tmp_time==4:
        plt.text(Mc,max(10**y_gr[:])-20,r'%.2f$\pm$%.3f'%(b,b_error),fontsize=10)
    else:
        plt.text(Mc,max(10**y_gr[:]),r'%.2f$\pm$%.3f'%(b,b_error),fontsize=10)


if __name__ == '__main__':
    print('This is for b calculation')