"""
Arquivo com todas as funções auxiliares utilizadas no projeto
"""
import re
import regex
from scipy import arange, array, exp
import pandas as pd
import random
import copy


def atof(text):
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval


def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    float regex comes from https://stackoverflow.com/a/12643073/190597
    '''
    return [atof(c) for c in re.split(r'[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)', text)]


def flatten(x):
    new_x = ' '.join([data for ele in x for data in ele])
    rx = regex.compile(r':.+?\.txt(*SKIP)(*FAIL)|\s+')
    return rx.split(new_x)


def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(list(map(pointwise, array(xs))))

    return ufunclike


def drange(start, stop, step):
    r = start
    while r < stop:
         yield r
         r += step


def multiply(*args):
    product = 1
    for a in args:
        if a is None:
            return None
        else:
            product *= a
    return product


def efficiency(Pout, t, ind, topology='TwoLev'):
    Pout *= 1e3
    ind = ind.result_inductor[ind.result_inductor['Flag'] == "True"]
    ind = ind.dropna()
    t = t.result_transistor.dropna()
    fsw_vector = pd.merge(t,
                          ind, on=['Fsw'], how='inner')["Fsw"].unique()
    deltaI_vector = pd.merge(t,
                          ind, on=['deltaI'], how='inner')['deltaI'].unique()
    d = []
    if topology == 'TwoLev':
        for deltaI in deltaI_vector:
            aux_L_delta = ind[ind['deltaI'].isin([deltaI])]
            for fsw in fsw_vector:
                aux_L_fsw = aux_L_delta[aux_L_delta['Fsw'].isin([fsw])]
                mat_vector = aux_L_fsw['Material'].unique()
                for mat in mat_vector:
                    aux_L_mat = aux_L_fsw[aux_L_fsw['Material'].isin([mat])]
                    aux_L = 3 * aux_L_mat['Loss_total'].iloc[0]
                    aux_T = 6 * t[t[
                                t['deltaI'].isin([deltaI])]['Fsw'].isin([fsw])]['Ptot_tot'].iloc[0]
                    aux_total = aux_L + aux_T
                    d.append(
                        {
                            'deltaI': deltaI,
                            'Fsw': fsw,
                            'Material': mat,
                            'T_loss': aux_T,
                            'L_loss': aux_L,
                            'Total_loss': aux_total,
                            'Eff': Pout / (Pout + aux_total)
                        })
    return pd.DataFrame(d)


def vol_dens(Pout, hs, ind):
    Cr = 0.6
    ind = ind.result_inductor[ind.result_inductor['Flag'] == "True"]
    ind = ind.dropna()
    hs = hs.result_hs.dropna()
    fsw_vector = pd.merge(hs,
                          ind, on=['Fsw'], how='inner')["Fsw"].unique()
    deltaI_vector = pd.merge(hs,
                             ind, on=['deltaI'], how='inner')['deltaI'].unique()
    d = []
    for deltaI in deltaI_vector:
        aux_i_delta = ind[ind['deltaI'].isin([deltaI])]
        aux_hs_delta = hs[hs['deltaI'].isin([deltaI])]
        for fsw in fsw_vector:
            aux_i_fsw = aux_i_delta[aux_i_delta['Fsw'].isin([fsw])]
            aux_hs_fsw = aux_hs_delta[aux_hs_delta['Fsw'].isin([fsw])]
            pn_vector = aux_hs_fsw['Part_Number'].unique()
            mat_vector = aux_i_fsw['Material'].unique()
            for mat in mat_vector: # pode existir em uma freq, mas n em outra..
                aux_i_mat = aux_i_fsw[aux_i_fsw['Material'].isin([mat])]
                for pn in pn_vector: # pode existir em uma freq, mas n em outra..
                    aux_ind = 3 * aux_i_mat['Vol_total'].iloc[0]
                    aux_hs_pn = aux_hs_fsw[aux_hs_fsw['Part_Number'].isin([pn])]
                    aux_hs_vol = aux_hs_pn['Volume'].iloc[0]
                    aux_total = (aux_ind + aux_hs_vol) / Cr
                    d.append(
                        {
                            'deltaI': deltaI,
                            'Fsw': fsw,
                            'Material': mat,
                            'Part_Number': pn,
                            'HS_vol': aux_hs_vol,
                            'L_vol': aux_ind,
                            'Total_vol': aux_total,
                            'Dens_pot_vol': Pout / aux_total
                        })

    return pd.DataFrame(d)



def weight_dens():
    pass
"""
def eff_sample(Pout, L_sample, T_sample):
    deltaI = T_sample['deltaI'].iloc[0]
    fsw = T_sample['Fsw'].iloc[0]
    T_loss = 6 * T_sample['Ptot_tot'].iloc[0]
    L_loss = 3 * L_sample['Loss_total'].iloc[0]
    Total_loss = T_loss + L_loss
    d.append(
        {
            'deltaI': deltaI,
            'Fsw': fsw,
            'Material': L_sample['Material'].iloc[0],
            'T_loss': T_loss,
            'L_loss': L_loss,
            'Total_loss': Total_loss,
            'Eff': Pout / (Pout + Total_loss)
        })
    return pd.DataFrame(d)
"""


def choose(x, n=1):
    return x.loc[random.sample(list(x.index), n)]

def eff_sample(Pout, L_sample, T_sample):
    T_loss = 6 * T_sample['Ptot_tot'].iloc[0]
    L_loss = 3 * L_sample['Loss_total'].iloc[0]
    Total_loss = T_loss + L_loss
    return Pout / (Pout + Total_loss)


def vol_dens_sample(Pout, L_sample, HS_sample):
    Cr = 0.6
    L_vol = 3 * L_sample['Vol_total'].iloc[0]
    HS_vol = HS_sample['Volume'].iloc[0]
    Total_vol = (L_vol + HS_vol) / Cr
    return Pout / Total_vol


def fitness(Pout, L_sample, T_sample, HS_sample):
    return eff_sample(Pout, L_sample, T_sample), vol_dens_sample(Pout, L_sample, HS_sample)