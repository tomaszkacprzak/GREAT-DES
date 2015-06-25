import numpy as np
import pylab as pl
import os
import tktools
import tktools.fitting


def get_model_prediction(filename_table_bias,xp,yp,plots=False):

    w,w_cov = get_model(filename_table_bias,plots)
    xxp = np.concatenate([xp[:,None],yp[:,None]],axis=1)
    print 'getting prediction for %d data points' % xxp.shape[0]
    p,sigma = tktools.fitting.predict(xxp,w,w_cov,expand=basis)

    return p

def get_model(filename_table_bias,plots=False):

    bias_table = tktools.load(filename_table_bias)
    print bias_table['m'][:3]
    xt=bias_table['vsnr_mid']
    yt=bias_table['vpsf_mid']
    zt=bias_table['m']-1
    st=bias_table['std_m']
    ng=bias_table['n_gals']

    w,w_cov = fit_model(xt,yt,zt,st,ng,plots)

    return w, w_cov



def get_optimal_w(snr_mid, m, s, expand_basis, n_gals=1):

        list_chi2 = []
        list_res_sig = []
        list_w = []
        list_w_cov = []
        list_res_tot = []
        eps_grid=np.linspace(-15,15,1000)
        for eps in eps_grid:
            w, w_cov = tktools.fitting.fit(snr_mid, m, s,expand=expand_basis,eps=10**eps)
            p_mid, sigma = tktools.fitting.predict(snr_mid,w,w_cov,expand=expand_basis)

            # chi2 =  np.sum( n_gals* (((m-p_mid)**2)/(s**2)) ) / np.sum(n_gals)
            # chi2 =  np.mean( (((m-p_mid)**2)/(s**2)) ) 
            # chi2 =  np.sum( (((m-p_mid)**2)/(s**2)) ) / float( len(snr_mid) - basis(snr_mid).shape[1] - 1  )
            # chi2 =  np.sum( n_gals*(m-p_mid) ) / np.sum(n_gals)
            chi2 =  np.mean( (((m-p_mid)**2)/(s**2)) )  - np.sum( n_gals*(m-p_mid) )/ float(np.sum(n_gals))
            # chi2 =  np.mean( (((m-p_mid)**2)/(s**2)) ) 
            res_sig = np.sum( n_gals*(m-p_mid) )/ float(np.sum(n_gals))
            res_tot = np.sum( n_gals*(p_mid) )/ float(np.sum(n_gals))
            # print 'optimal w: eps=%2.2e chi2=%2.2e mean_res/sig=%2.4e' % (10**(eps),chi2,res_sig)
            list_chi2.append(chi2)
            list_w.append(w)
            list_w_cov.append(w_cov)
            list_res_sig.append(res_sig)
            list_res_tot.append(res_tot)

        arr_chi2 = np.array(list_chi2)

        select = np.argmin(np.abs(arr_chi2-1))

        w = list_w[select]
        w_cov = list_w_cov[select]
        eps = eps_grid[select]
        chi2 = list_chi2[select]
        res_sig = list_res_sig[select]
        res_tot = list_res_tot[select]

        print '(multiplicative) final optimal w: eps=%2.2e chi2=%2.2f res_sig=%2.5f res_tot=%2.4f' % (10**(eps),chi2,res_sig,res_tot)

        return w, w_cov


def basis(x):

    n_func = 22
    X = np.zeros([x.shape[0],n_func])

    x0=x[:,0]/1000.
    x1=(x[:,1]-1)/10

    X[:,0]  = 1/x0**2    * 1/x1**2
    X[:,1]  = 1/x0**3    * 1/x1**3
    X[:,2]  = 1/x0**3    * 1/x1**2
    X[:,3]  = 1/x0**2    * 1/x1**3
    X[:,4]  = 1/x0**1.5  * 1/x1**4 # not used
    X[:,5]  = 1/x0**2    * 1/x1**4 # not used
    X[:,6]  = 1/x0**2.5  * 1/x1**4 # not used
    X[:,7]  = 1/x0**2.5  * 1/x1**2.5
    X[:,8]  = 1/x0**2.5  * 1/x1**3
    X[:,9]  = 1/x0**3    * 1/x1**2.5
    X[:,10] = 1/x0**1.5  * 1/x1**1.5
    X[:,11] = 1/x0**2.5  * 1/x1**1.5
    X[:,12] = 1/x0**1.5  * 1/x1**2.5
    X[:,13] = 1/x0**1.5  * 1/x1**2
    X[:,14] = 1/x0**2    * 1/x1**1.5
    X[:,15] = 1/x0**1.25 * 1/x1**1.75
    X[:,16] = 1/x0**1.75 * 1/x1**1.25
    X[:,17] = 1/x0**4 * 1/x1**4
    X[:,19] = 1/x0**4    * 1/x1**1.5 # not used
    X[:,18] = 1/x0**5    * 1/x1**1.5 # not used
    X[:,20] = 1/x0**6    * 1/x1**1.5 # not used

    return X
    return X

def fit_model(xt,yt,zt,st,ng,plots=False):

    xx=np.concatenate([xt[:,None],yt[:,None]],axis=1)


    # w,C=tktools.fitting.fit(xx,zt,s=st,expand=basis)
    w,w_cov=get_optimal_w(xx,zt,st,n_gals=ng,expand_basis=basis)


    b,sigma = tktools.fitting.predict(xx,w,w_cov,expand=basis)

    r = b-zt

    print 'expected error=%2.4f' % (np.sum(r*ng)/float(np.sum(ng)))


    if plots:

        # pl.style.use('supermongo')

        n_upsample = 200
        
        X,Y=np.meshgrid(np.logspace(1,2,n_upsample),np.linspace(1.15,2,n_upsample))
        x=X.flatten()
        y=Y.flatten()

        xp = x
        yp = y

        xxp = np.concatenate([xp[:,None],yp[:,None]],axis=1)
        p,sigma = tktools.fitting.predict(xxp,w,w_cov,expand=basis)

        pl.figure()
        pl.scatter(xt,yt,c=zt,s=100,lw=0)
        pl.xscale('log')
        pl.title(r'multiplicative shear bias $m$')
        pl.xlabel(r'SNR')
        pl.xlabel(r'R_{rpp}/R_{rp}')
        pl.clim([-0.6,0.1])
        pl.xticks([10,12,15,20,30,50,100])
        pl.yticks([1.15,1.3,1.45,1.6,1.75,1.9])
        pl.colorbar()

        pl.figure()
        pl.scatter(xt,yt,c=st,s=100,lw=0)
        pl.xscale('log')
        pl.title('std_m')
        pl.colorbar()

        pl.figure()
        pl.scatter(xt,yt,c=b,s=100,lw=0)
        pl.xscale('log')
        pl.title('model')
        pl.clim([-0.6,0.1])
        pl.colorbar()

        pl.figure()
        pl.scatter(xt,yt,c=(zt-b)/st,s=100,lw=0)
        pl.xscale('log')
        pl.title('residual/sigma')
        pl.colorbar()

        pl.figure()
        pl.scatter(xt,yt,c=(zt-b),s=100,lw=0)
        pl.xscale('log')
        pl.title('residual')
        pl.clim([-0.1,0.1])
        pl.colorbar()

        pl.figure()
        pl.scatter(xp,yp,c=p,s=100,lw=0)
        pl.xscale('log')
        pl.xlim([0,200])
        pl.clim([-0.6,0.1])
        pl.colorbar()

        pl.figure()
        pl.pcolormesh(np.reshape(xp,[n_upsample,n_upsample]),np.reshape(yp,[n_upsample,n_upsample]),np.reshape(p,[n_upsample,n_upsample]),cmap='jet')
        pl.title(r'multiplicative bias $m$')
        pl.xscale('log')
        pl.xlabel(r'SNR')
        pl.ylabel(r'$R_{gpp}/R_{rp}$')
        pl.colorbar()
        pl.clim([-0.3,0.1])
        pl.xticks([10,12,15,20,30,50,100])
        pl.yticks([1.15,1.3,1.45,1.6,1.75,1.9])
        pl.xlim([15,100])
        pl.ylim([1.15,2])
        pl.gca().xaxis.set_major_formatter(pl.matplotlib.ticker.ScalarFormatter())

    return w, w_cov


if __name__=='__main__':

    filename_table_bias = os.path.join('case-19-final3','bias_table.fits')

    w,w_cov = get_model(filename_table_bias,plots=True)
    pl.show()
    import pdb; pdb.set_trace()
