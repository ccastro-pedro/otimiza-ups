"""
Arquivo com todas as classes utilizadas no projeto
"""
import os
import re
import sys
import funcs
from collections import defaultdict
from scipy.interpolate import interp1d
import copy
import pandas as pd


class Folders:
    def __init__(self, Vout, Vcc, Pout, ripple, mod=''):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        #self.__transistor_folder = dir_path + "\\Transistors\\SCT3030AL\\Inverter\\127V_400V_10kW\\Ripple=0.3"
        #self.__thermal_folder = dir_path + "\\Transistors\\SCT3030AL"
        #self.__hs_folder = dir_path + "\\Heatsinks"
        #self.__inductor_folder = dir_path + "\\Indutores\\127V_400V_10kW"
        #self.__capacitor_folder = dir_path + "\\Capacitors\\Database"

        self.__transistor_folder = dir_path + "\\Transistors\\" + mod + "\\Inverter\\" + \
                                   str(Vout) + "V_" + str(Vcc) + "V_" + str(Pout) + "kW\\Ripple=" + str(ripple)
        self.__thermal_folder = dir_path + "\\Transistors\\" + mod
        self.__hs_folder = dir_path + "\\Heatsinks"
        self.__inductor_folder = dir_path + "\\Indutores\\" \
                                 + str(Vout) + "V_" + str(Vcc) + "V_" + str(Pout) + "kW"
        self.__capacitor_folder = dir_path + "\\Capacitors\\Database"


class Transistor(Folders):
    def __init__(self, mod, Vout, Vcc, Pout, ripple):

        Folders.__init__(self, Vout, Vcc, Pout, ripple, mod)
        self.__transistor_folder = self._Folders__transistor_folder
        self.__keys = (
            "deltaI", "freq", "ma", "Ecc", "Vin", "L", "rL", "C", "Linear", "Rcarga", "fp", "Pcarga", "Npontos",
            "Nciclos", "Topologia", "Mod", "Path", "Fsw", "mf", "Tamb", "Ron", "Roff", "Pcond_D1", "Pon_D1", "Poff_D1",
            "Psw_D1", "Ptot_D1", "Pcond_Q1", "Pon_Q1", "Poff_Q1", "Psw_Q1", "Ptot_Q1", "Pcond_D2", "Pon_D2", "Poff_D2",
            "Psw_D2", "Ptot_D2", "Pcond_Q2", "Pon_Q2", "Poff_Q2", "Psw_Q2", "Ptot_Q2", "Pcond_D5", "Pon_D5", "Poff_D5",
            "Psw_D5", "Ptot_D5", "Pcond_tot", "Pon_tot", "Poff_tot", "Psw_tot", "Ptot_tot", "Tj_D1", "Tc_D1", "Tj_Q1",
            "Tc_Q1", "Tj_D2", "Tc_D2", "Tj_Q2", "Tc_Q2", "Tj_D5", "Tc_D5", "Ts"
        )
        self.__read_and_save_transistor()
        """
        try:
            for key in keys:
                setattr(self, "__" + key, kwargs.get(key))  # Isso cria atributos privados? Verificar
        except TypeError:
            print("Algum erro ocorreu na classe Transistores")
        """

    # Ler os arquivos dos transistores
    def __read_and_save_transistor(self):
        """
        Método para ler e salvar em variáveis as informações dos transistores
        Cria-se um objeto da classe Transistor que possui
        :return:
        """
        number_of_files = len([item for item in sorted(os.listdir(self.__transistor_folder), key=funcs.natural_keys)
                               if os.path.isfile(os.path.join(self.__transistor_folder, item))])
        aux_transistors = [] * number_of_files

        # Get the ripple values, based on the folder's name
        self.__ripple = os.path.basename(self.__transistor_folder).split('=')[1].lstrip()

        for file in [item for item in sorted(os.listdir(self.__transistor_folder), key=funcs.natural_keys) if
                     os.path.isfile(os.path.join(self.__transistor_folder, item))]:
            with open(self.__transistor_folder + "\\" + file, 'r') as transistor_file:
                # aux = [line.split(";")[1:-2] for line_num, line in enumerate(transistor_file) if re.match(".\)", line)]
                aux = [line.split(";")[1:-2] for line_num, line in enumerate(transistor_file) if re.match(".\)", line)]
            aux_transistors.append(aux)

        # Split the loaded variable into a dict of the name of the Rth and its value
        # self.all_transistors_dict_aux = [dict(zip(i[0].split(" "), i[1:])) for i in all_transistors]

        # self.info_transistors = dict(zip(self.__keys, [all_transistors[i][j][k]]))

        all_transistors = [[[self.__ripple] if j == 0 else aux_transistors[i][j - 1] for j in
                            range(len(aux_transistors[i]) + 1)] for i in range(len(aux_transistors))]

        # self.__info_transistors = dict([zip(self.__keys, info_tran)  for info_tran in all_transistors])
        # self.info_transistors = {key: funcs.flatten(info_tran[i]) for info_tran in all_transistors for key in self.__keys}
        # self.__info_transistors = dict.fromkeys(self.__keys)

        # [zip(self.__keys, flatten(info_tran)) for info_tran in all_transistors]

        # self.__info_transistors = [zip(self.__keys, funcs.flatten(info_tran)) for info_tran in all_transistors]

        info_transistors_aux = defaultdict(list)
        for info_tran in all_transistors:
            for i in range(len(self.__keys)):
                info_transistors_aux[self.__keys[i]].append(funcs.flatten(info_tran)[i])
        self.__info_transistors = {key: tuple(map(lambda x: funcs.atof(x), value)) for key, value in
                                   info_transistors_aux.items()}
        self.__result_transistor = pd.DataFrame(self.__info_transistors)

    @property
    def result_transistor(self):
        return self.__result_transistor

        # Agora é preciso passar um elemento de cada vez para a classe

        # Agora é preciso passar um elemento de cada vez para a classe


class HeatSink(Folders):
    #  Algumas constantes
    layer_thickness = 0.508  # [mm] # Espessura d
    factor = 5 * 10 ** -3  # [W/(mm.K)]
    # SW_case = 'TO-247'
    width = 16 / 100  # [dm]  Largura
    height = 21 / 100  # [dm]  Altura
    leg = 20 / 100  # [dm] % Tamanho
    area = width * height

    def __init__(self, transistor, mod, Vout, Vcc, Pout, ripple, temp_j_max=150, temp_amb=40, vent=0, hs_config=1):
        Folders.__init__(self, Vout, Vcc, Pout, ripple, mod)
        self.__temp_j_max = temp_j_max
        self.__temp_amb = temp_amb
        self.__vent = vent
        self.__hs_config = hs_config
        self.__Rth_cs = self.layer_thickness / (self.factor * self.area * 10000)
        self.__thermal_folder = self._Folders__thermal_folder
        self.__hs_folder = self._Folders__hs_folder
        try:
            self.__info_transistors = transistor.result_transistor
        except AttributeError as AE:
            print(f"ErrorMsg: {AE} \n You must pass a object of type Transistor() to the class HeatSink() ")
            sys.exit(1)

        # Chamada dos métodos privados
        self.__load_Rth()
        self.__load_heat_sink()
        self.__calc_Rth()
        self.__fix_len()

    def __load_Rth(self):
        with open(self.__thermal_folder + "\\dados_termicos.txt", 'r') as thermal_file:
            read_Rth_hs = [line.split(";")[:-2] for line_num, line in enumerate(thermal_file)]

        # Split the loaded variable into a dict with the name of the Rth and its value
        read_Rth_hs_dict_aux = [dict(zip(i[0].split(","), i[1:])) for i in read_Rth_hs]

        # Transform the list of dicts into a single dict
        info_Rth_hs_dict_aux = {}
        for d in read_Rth_hs_dict_aux:
            info_Rth_hs_dict_aux.update(d)

        # Transform all string numbers to float
        self.__info_Rth_hs_dict = dict(zip(info_Rth_hs_dict_aux.keys(),
                                           tuple(map(lambda x: funcs.atof(x), info_Rth_hs_dict_aux.values()))))

        for key, value in self.__info_Rth_hs_dict.items():
            if value == 0:
                self.__info_Rth_hs_dict[key] = self.__Rth_cs

        self.__info_Rth_hs_dict = pd.DataFrame(self.__info_Rth_hs_dict,
                                               index=[os.path.basename(os.path.normpath(self.__thermal_folder))])

    def __load_heat_sink(self):
        with open(self.__hs_folder + "\\" + 'dissipadores.txt', 'r') as hs_file:
            self.__dimension = hs_file.readline().split(":")[1].replace("\n", "")
            keys = tuple(map(lambda x: x.replace("\n", ""), hs_file.readline().split(" ")))
            hs_aux = tuple([line.split(" ") for line_num, line in enumerate(hs_file)])

        hs_dict_aux = {keys[j]: tuple([hs_aux[i][j] for i in range(len(hs_aux))]) for j in range(len(keys))}
        # Info of the heatsink's database
        self.__info_hs_db = dict(
            zip(hs_dict_aux.keys(), tuple(map(lambda x: tuple([funcs.atof(i) for i in x]), hs_dict_aux.values()))))
        self.__info_hs_db['DimensionalFactor'] = {
            'mm': 100,
            'cm': 10,
            'dm': 1,
            'inches': 3.937
        }[self.__dimension.lower().lstrip()]

        self.__info_hs_db['Width'] = tuple(
            map(lambda x: round(x / self.__info_hs_db['DimensionalFactor'], 4), self.__info_hs_db['Width']))
        self.__info_hs_db['Height'] = tuple(
            map(lambda x: round(x / self.__info_hs_db['DimensionalFactor'], 4), self.__info_hs_db['Height']))
        self.__info_hs_db['Peso'] = tuple(
            map(lambda x: round(x / 10, 4), self.__info_hs_db['Peso']))

        self.__info_hs_db = pd.DataFrame(self.__info_hs_db)

    def __calc_Rth(self):
        info_hs_aux = defaultdict(list)
        for i in range(len(self.__info_transistors['deltaI'])):
            info_hs_aux['deltaI'].append(self.__info_transistors['deltaI'][i])
            info_hs_aux['Fsw'].append(self.__info_transistors['Fsw'][i])

            temp_hs_Q1 = self.__temp_j_max - self.__info_transistors['Ptot_Q1'][i] * \
                         (self.info_Rth_hs_dict[[key for key in self.info_Rth_hs_dict
                                                 if 'jc' in key and 'Q1' in key][0]].iloc[0] +
                          self.info_Rth_hs_dict[[key for key in self.info_Rth_hs_dict
                                                 if 'cs' in key and 'Q1' in key][0]].iloc[0])

            temp_hs_D1 = self.__temp_j_max - self.__info_transistors['Ptot_D1'][i] * \
                         (self.info_Rth_hs_dict[[key for key in self.info_Rth_hs_dict
                                                 if 'jc' in key and 'D1' in key][0]].iloc[0] +
                          self.info_Rth_hs_dict[[key for key in self.info_Rth_hs_dict
                                                 if 'cs' in key and 'D1' in key][0]].iloc[0])

            info_hs_aux['T_hs'].append(min(temp_hs_Q1, temp_hs_D1))
            info_hs_aux['delta_T_hs'].append(info_hs_aux['T_hs'][i] - self.__temp_amb)
            # Limitar o delta de temperatura em 74ºC
            if info_hs_aux['delta_T_hs'][i] >= 74:
                info_hs_aux['delta_T_hs'][i] = 74
            if self.__hs_config == 1:
                # R_hs é a resistencia necessaria para cada caso (delta I, fsw)
                info_hs_aux['R_hs'].append((info_hs_aux['delta_T_hs'][i]
                                            / (6 * (self.__info_transistors['Ptot_Q1'][i] +
                                                    self.__info_transistors['Ptot_D1'][i]))))
        info_hs_aux['R_hs'] = self.__fix_temp(delta_T=info_hs_aux['delta_T_hs'], R_hs=info_hs_aux['R_hs'])
        info_hs_aux['R_hs'] = self.__fix_vent(R_hs=info_hs_aux['R_hs'])
        self.__info_hs = {key: list(value) for key, value in info_hs_aux.items()}
        self.__info_hs = pd.DataFrame(self.__info_hs)

    def __fix_temp(self, delta_T, R_hs):
        with open(self.__hs_folder + "\\" + 'DeltaTemp.txt', 'r') as delta_temp_file:
            temp_x_axis = delta_temp_file.readline()
            temp_y_axis = delta_temp_file.readline()
        # Filtrar as strings em branco
        temp_x_axis = list(filter(lambda item: type(item) is not str,
                                  list(map(lambda x: funcs.atof(x.rstrip()), temp_x_axis.split(" ")))))
        temp_y_axis = list(filter(lambda item: type(item) is not str,
                                  list(map(lambda x: funcs.atof(x.rstrip()), temp_y_axis.split(" ")))))
        return list(map(lambda x, y: x / y, R_hs,
                        list(map(lambda x: 1 + interp1d(temp_y_axis, temp_x_axis)(x), delta_T))))

    def __fix_vent(self, R_hs):
        with open(self.__hs_folder + "\\" + 'ventilacao.txt', 'r') as vent_file:
            vent_x_axis = vent_file.readline()
            vent_y_axis = vent_file.readline()
        # Filtrar as strings em branco e passar para float
        vent_x_axis = list(filter(lambda item: type(item) is not str,
                                  list(map(lambda x: funcs.atof(x.rstrip()), vent_x_axis.split(" ")))))
        vent_y_axis = list(filter(lambda item: type(item) is not str,
                                  list(map(lambda x: funcs.atof(x.rstrip()), vent_y_axis.split(" ")))))
        return list(map(lambda y: y * interp1d(vent_x_axis, vent_y_axis)(self.__vent)[()], R_hs))

    def __fix_len(self):
        with open(self.__hs_folder + "\\" + 'CorrigirComprimento.txt', 'r') as len_file:
            len_x_axis = len_file.readline()
            len_y_axis = len_file.readline()
        len_x_axis = list(filter(lambda item: type(item) is not str,
                                 list(map(lambda x: funcs.atof(x.rstrip()), len_x_axis.split(" ")))))
        len_y_axis = list(filter(lambda item: type(item) is not str,
                                 list(map(lambda x: funcs.atof(x.rstrip()), len_y_axis.split(" ")))))

        len_x_axis = list(map(lambda x: x / self.info_hs_db['DimensionalFactor'].iloc[0], len_x_axis))
        len_y_axis_Rth = {'Rth_len':
                              tuple([[Rth * factor for factor in len_y_axis] for Rth in self.info_hs_db['Rth']]),
                          'Part_number':
                              tuple([[Code] for Code in self.info_hs_db['Code']])}
        min_len = (self.height + self.leg)
        max_len = max(len_x_axis)
        # self.__info_hs['HS_len'] = [None] * len(self.__info_hs['R_hs'])
        self.__info_hs['HS_len'] = pd.Series([None for _ in range(len(self.__info_hs))], index=self.__info_hs.index)
        self.__result_hs = defaultdict(list)
        d = []
        aux = copy.deepcopy(self.__info_hs)
        for n in range(len(len_y_axis_Rth['Part_number'])):
            max_y = max(len_y_axis_Rth['Rth_len'][n])
            min_y = min(len_y_axis_Rth['Rth_len'][n])
            for i in range(len(self.__info_hs)):
                if aux['R_hs'][i] > max_y:
                    aux.at[i, 'R_hs'] = max_y
                    # aux['R_hs'][i] = max_y
                if aux['R_hs'][i] < min_y:
                    # self.__info_hs['R_hs'][i] = min_y
                    # aux['HS_len'][i] = aux['HS_len'][i].replace(aux['HS_len'][i], None)
                    aux.at[i, 'HS_len'] = None
                    # aux['HS_len'][i] = None  # esta assim no programa original
                else:
                    if self.__hs_config == 1:
                        aux.at[i, 'HS_len'] = \
                            interp1d(len_y_axis_Rth['Rth_len'][n], len_x_axis)(aux.at[i, 'R_hs'])[()]
                        flag = 0
                        if aux['HS_len'][i] is not None:
                            # for tam in [x * 0.0001 for x in range(self.__info_hs['HS_len'][i], max_len)]:
                            for tam in funcs.drange(aux.at[i, 'HS_len'], max_len, 0.0001):
                                if tam >= min_len:
                                    flag = 1
                                    aux.at[i, 'HS_len'] = tam
                                    break
                            if not flag:
                                aux.at[i, 'HS_len'] = None

                d.append(
                    {
                        'deltaI': copy.deepcopy(aux['deltaI'][i]),
                        'Fsw': copy.deepcopy(aux['Fsw'][i]),
                        'Part_Number': len_y_axis_Rth['Part_number'][n][0],
                        'Width': self.__info_hs_db['Width'][n],
                        'Height': self.__info_hs_db['Height'][n],
                        'Length': copy.deepcopy(aux['HS_len'][i]),
                        'Volume': funcs.multiply(self.__info_hs_db['Width'][n],
                                                 self.__info_hs_db['Height'][n], copy.deepcopy(aux['HS_len'][i])),
                        'Peso': funcs.multiply(copy.deepcopy(aux['HS_len'][i]), self.__info_hs_db['Peso'][n])
                    }
                )
        self.__result_hs = pd.DataFrame(d)

        # if self.__hs_config == 1 and self.__result_hs['Width'][n] < 6 * (self.width + 0.5 * self.width):
        # self.__result_hs['Length'][i] = None
        # for n in range(len(self.__result_hs['Part_Number'])):
        # for i in range(len(self.result_hs['info_hs'][n]['R_hs'])):

        # self.__result_hs = pd.DataFrame(self.__result_hs)

    @property
    def info_hs_db(self):
        return self.__info_hs_db

    @property
    def info_Rth_hs_dict(self):
        return self.__info_Rth_hs_dict

    @property
    def info_hs(self):
        return {key: tuple(value) for key, value in self.__info_hs.items()}

    @property
    def result_hs(self):
        return self.__result_hs


class Inductors(Folders):

    def __init__(self, Vout, Vcc, Pout, ripple):
        Folders.__init__(self, Vout, Vcc, Pout, ripple)
        self.__inductor_folder = self._Folders__inductor_folder
        self.__read_and_save_inductors()

    def __read_and_save_inductors(self):
        number_of_files = len([item for item in sorted(os.listdir(self.__inductor_folder), key=funcs.natural_keys)
                               if os.path.isfile(os.path.join(self.__inductor_folder, item)) and item.startswith(
                "Inverter")])
        d = []
        for file in [item for item in sorted(os.listdir(self.__inductor_folder), key=funcs.natural_keys) if
                     os.path.isfile(os.path.join(self.__inductor_folder, item)) and item.startswith("Inverter")]:
            with open(self.__inductor_folder + "\\" + file, 'r') as inductor_file:
                # aux = [line.split(";")[1:-2] for line_num, line in enumerate(transistor_file) if re.match(".\)", line)]
                aux = [line.split("=") for line_num, line in enumerate(inductor_file)]
                if len(aux) == 25:
                    d.append(
                        {
                            'Fsw': copy.deepcopy(funcs.atof(aux[0][1].rstrip().split(' ')[1])),
                            'L_ind': copy.deepcopy(funcs.atof(aux[1][1].rstrip().split(' ')[1])),
                            'deltaI': copy.deepcopy(funcs.atof(aux[2][1].rstrip().split(' ')[1])),
                            'Flag': copy.deepcopy(funcs.atof(aux[3][1].rstrip().split(' ')[1])),
                            'L_init': copy.deepcopy(funcs.atof(aux[4][1].rstrip().split(' ')[1])),
                            'L_sat': copy.deepcopy(funcs.atof(aux[5][1].rstrip().split(' ')[1])),
                            'Error': copy.deepcopy(funcs.atof(aux[6][1].rstrip().split(' ')[1])),
                            'Material': copy.deepcopy(aux[7][1].lstrip().rstrip()),
                            'Part_Number': copy.deepcopy(aux[8][1].lstrip().rstrip()),
                            'Number_coils': copy.deepcopy(funcs.atof(aux[9][1].rstrip().split(' ')[1])),
                            'I_dens': copy.deepcopy(funcs.atof(aux[10][1].rstrip().split(' ')[1])),
                            'Fill_Factor': copy.deepcopy(funcs.atof(aux[11][1].rstrip().split(' ')[1])),
                            'Number_layers': copy.deepcopy(funcs.atof(aux[12][1].rstrip().split(' ')[1])),
                            'Number_paralel': copy.deepcopy(funcs.atof(aux[13][1].rstrip().split(' ')[1])),
                            'Loss_total': copy.deepcopy(funcs.atof(aux[15][1].rstrip().split(' ')[1])),
                            'Loss_core': copy.deepcopy(funcs.atof(aux[16][1].rstrip().split(' ')[1])),
                            'Loss_copper': copy.deepcopy(funcs.atof(aux[17][1].rstrip().split(' ')[1])),
                            'Loss_copper_skin': copy.deepcopy(funcs.atof(aux[18][1].rstrip().split(' ')[1])),
                            'Delta_Temp': copy.deepcopy(funcs.atof(aux[19][1].rstrip().split(' ')[1])),
                            'R_60hz': copy.deepcopy(funcs.atof(aux[20][1].rstrip().split(' ')[1])),
                            'R_fsw': copy.deepcopy(funcs.atof(aux[21][1].rstrip().split(' ')[1])),
                            'Vol_total': copy.deepcopy(funcs.atof(aux[22][1].rstrip().split(' ')[1])),
                            'Weight_core': copy.deepcopy(funcs.atof(aux[23][1].rstrip().split(' ')[1])),
                            'Weight_copper': copy.deepcopy(funcs.atof(aux[24][1].rstrip().split(' ')[1]))
                        }
                    )
                # Preciso salvar todos os arquivos que nao existem?
                else:
                    d.append(
                        {
                            'Fsw': copy.deepcopy(funcs.atof(aux[0][1].rstrip().split(' ')[1])),
                            'L_ind': copy.deepcopy(funcs.atof(aux[1][1].rstrip().split(' ')[1])),
                            'deltaI': copy.deepcopy(funcs.atof(aux[2][1].rstrip().split(' ')[1])),
                            'Material': copy.deepcopy(aux[3][1].lstrip().rstrip()),
                            'Flag': copy.deepcopy(funcs.atof(aux[4][1].rstrip().split(' ')[1])),
                            'L_init': None,
                            'L_sat': None,
                            'Error': None,
                            'Part_Number': None,
                            'Number_coils': None,
                            'I_dens': None,
                            'Fill_Factor': None,
                            'Number_layers': None,
                            'Number_paralel': None,
                            'Diam_wires': None,
                            'Loss_total': None,
                            'Loss_core': None,
                            'Loss_copper': None,
                            'Loss_copper_skin': None,
                            'Delta_Temp': None,
                            'R_60hz': None,
                            'R_fsw': None,
                            'Vol_total': None,
                            'Weight_core': None,
                            'Weight_copper': None
                        }
                    )

        self.__result_inductor = pd.DataFrame(d)

    @property
    def result_inductor(self):
        return self.__result_inductor


class Capacitors(Folders):
    def __init__(self):
        Folders.__init__(self)
        self.__capacitor_folder = self._Folders__capacitor_folder

    def __read_and_save_capacitor(self):
        number_of_files = len([item for item in sorted(os.listdir(self.__capacitor_folder), key=funcs.natural_keys)
                               if os.path.isfile(os.path.join(self.__capacitor_folder, item))])
        aux_transistors = [] * number_of_files

        # Get the ripple values, based on the folder's name
        self.__ripple = os.path.basename(self.__transistor_folder).split('=')[1].lstrip()

        for file in [item for item in sorted(os.listdir(self.__transistor_folder), key=funcs.natural_keys) if
                     os.path.isfile(os.path.join(self.__transistor_folder, item))]:
            with open(self.__transistor_folder + "\\" + file, 'r') as transistor_file:
                # aux = [line.split(";")[1:-2] for line_num, line in enumerate(transistor_file) if re.match(".\)", line)]
                aux = [line.split(";")[1:-2] for line_num, line in enumerate(transistor_file) if re.match(".\)", line)]
            aux_transistors.append(aux)
