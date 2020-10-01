"""
Função principal do programa de otimização multi objetiva de uma UPS de dupla conversão
Autor: Pedro Augusto de Castro e Castro
Ideias: salvar todos os arquivos em um bando de dados, data lake? Na nuvem?
Ler tais arquivos, salvar em variaveis (classes).
Definir uma forma de otmização com GA.
"""
import os
import re
import classes
import funcs
import copy
import collections
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd
from functools import reduce
from deap import base
from deap import creator
from deap import tools
import random
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.factory import get_problem
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
from pymoo.configuration import Configuration
Configuration.show_compile_hint = False
# Ler os arquivos dos indutores

# Ler os arquivos dos capacitores, ou gerar um novo arquivo de resultados

# Ler os arquivos dos dissipadores, ou gerar um novo arquivo de resultados


if __name__ == "__main__":
    # Escolha os parâmetros a serem analisados:
    Vout = 127  # V
    Vcc = 400  # V
    Pout = 10  # kW
    ripple = 0.3
    mod = 'SCT3030AL'

    transistor = classes.Transistor(mod=mod, Vout=Vout, Vcc=Vcc, Pout=Pout, ripple=ripple)
    hs = classes.HeatSink(transistor, mod=mod, Vout=Vout, Vcc=Vcc, Pout=Pout, ripple=ripple)
    inductor = classes.Inductors(Vout=Vout, Vcc=Vcc, Pout=Pout, ripple=ripple)
    #print([key for key, value in hs_inv.result_hs['info_hs'][0].items()])
    #volume = [[] * len(hs_inv.result_hs['info_hs'][0]['HS_len'])] * len(hs_inv.result_hs['Part_Number'])
    #transistor.result_transistor.info(verbose=True)
    #inductor.result_inductor.info(verbose=True)
    #hs.result_hs.info(verbose=True)
    #print(transistor.result_transistor)
    #print(inductor.result_inductor)
    #fsw = pd.merge(transistor.result_transistor, inductor.result_inductor, on=['Fsw'], how='inner')["Fsw"]

    #fsw = reduce(lambda left, right: pd.merge(left, right, on='Fsw'), dfs)['Fsw']

    #d = funcs.efficiency(Pout, transistor, inductor)
    #d.info(verbose=True)
    #df = df.drop(columns=df.corr().query("-0.05 < col < 0.05").index)
    #aux= aux[index]

    #eff = funcs.efficiency(Pout, transistor, inductor)
    #dens = funcs.vol_dens(Pout, hs, inductor)

    """
    for deltaI in deltaI_vector:
        for fsw in fsw_vector:
            for mat in mat_vector:
                for pn in pn_vector:
                    print(pn)
                    aux_hs = hs[hs['Fsw'].isin([fsw])]
                    aux_hs = aux_hs[aux_hs['deltaI'].isin([deltaI])]
                    aux_hs = aux_hs[aux_hs['Part_Number'].isin([pn])]
                    print(aux_hs['Volume'])
    """

    """

    hs = hs.result_hs.dropna()
    aux_hs = hs[hs['Fsw'].isin([16320])]
    print(aux_hs['Part_Number'].unique())
    aux_hs = aux_hs[aux_hs['deltaI'].isin([0.3])]
    aux_hs = aux_hs[aux_hs['Part_Number'].isin(['HS6835'])]
    #print(aux_hs['Volume'])
    """

    L_p = inductor.result_inductor.dropna()  # inductor population
    T_p = transistor.result_transistor.dropna()  # transistor population
    HS_p = hs.result_hs.dropna()  # hs population

    aux_T = T_p[T_p[T_p['deltaI'].isin([0.3])]['Fsw'].isin([9600])].iloc[0]

    aux_L_delta = L_p[L_p['deltaI'].isin([0.3])]
    aux_L_fsw = aux_L_delta[aux_L_delta['Fsw'].isin([9600])]
    aux_L = aux_L_fsw

    # para cada ponto (deltaI, fsw) vc terá um individuo
    #print(aux_L.sample()['Material'].iloc[0])
    """
    aux_hs_delta = HS_p[HS_p['deltaI'].isin([0.3])]
    aux_hs_fsw = aux_hs_delta[aux_hs_delta['Fsw'].isin([9600])]
    aux_hs = aux_hs_fsw
    """

    fsw_vector = hs.result_hs["Fsw"].unique()
    deltaI_vector = hs.result_hs['deltaI'].unique()
    d = []
    aux_hs = hs.result_hs.dropna()
    for deltaI in deltaI_vector:
        aux_hs_delta = aux_hs[aux_hs['deltaI'].isin([deltaI])]
        for fsw in fsw_vector:
            aux_hs_fsw = aux_hs_delta[aux_hs_delta['Fsw'].isin([fsw])]
            aux_hs_vol = aux_hs_fsw[aux_hs_fsw['Volume'].isin([min(aux_hs_fsw['Volume'])])]
            d.append(
                {
                    'deltaI': deltaI,
                    'Fsw': fsw,
                    'Part_Number': aux_hs_vol['Part_Number'],
                    'HS_vol': aux_hs_vol['Volume'],
                }
            )

    result = collections.defaultdict(list)
    for i in d:
        for k, v in i.items():
            result[k].append(v)

    hs_min = pd.DataFrame(result)
    hs_vol = []
    hs_fsw = []
    for k, fsw in zip(hs_min['HS_vol'], hs_min['Fsw']):
        hs_vol.append(k.iloc[0])
        hs_fsw.append(fsw)


    print(hs_vol)
    plt.plot(hs_fsw, hs_vol)
    plt.show()


    N_CYCLES = 4
    IND_SIZE = 10

    toolbox = base.Toolbox()

    toolbox.register('Fitness', funcs.fitness)

    creator.create("FitnessMulti", base.Fitness, weights=(1.0, 1.0))
    creator.create("Individual", pd.core.frame.DataFrame, fitness=creator.FitnessMulti)

    toolbox.register("attr_hs", funcs.choose, aux_hs, n=1)
    toolbox.register("attr_L", funcs.choose, aux_L, n=1)
    toolbox.register("attr_T", funcs.choose, aux_T, n=1)
    toolbox.register("individual", tools.initCycle, creator.Individual,
                     (toolbox.attr_hs, toolbox.attr_L, toolbox.attr_T), n=N_CYCLES)


    toolbox.register("population", tools.initRepeat, pd.core.frame.DataFrame, toolbox.individual)

    toolbox.register("evaluate", funcs.fitness)
    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("select", tools.selTournament, tournsize=3)


    #pop = toolbox.population(n=min(len(aux_hs), len(aux_L)))


    # um individuo é composto por: dissipador, transistor e indutor

    #toolbox.register("attr_transistor", transistor_population.sample())





    #aux_hs_vol = aux_hs['Volume']

    #print(hs['Part_Number'].loc['HS271019'])
    #print(hs[hs['Part_Number'].isin(['HS271019'])])

# Delta i: 0.3 // fsw: 16320.0 // mat: High Flux // pn: HS6835
"""
['HS6835' 'HS11450' 'HS12168' 'HS12764' 'HS14050' 'HS14569' 'HS15073'
 'HS15450' 'HS15559' 'HS15560' 'HS17232' 'HS19334' 'HS21575' 'HS26574'
 'HS271019' 'HS125135']"""
