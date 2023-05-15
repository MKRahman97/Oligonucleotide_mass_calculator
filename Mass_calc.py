import pandas as pd
import numpy as np
import regex as re
from collections import Counter, OrderedDict
import itertools, sys, os
import time
import pprint as pp
from myconfig import *
pd.options.mode.chained_assignment = None  # default='warn'



start = time.time()


''''
Section below should not be changed
'''



mass_df = pd.read_csv('Element_masses.csv', header=None) #read element masses
mass_dict = dict(sorted(mass_df.values.tolist())) #convert to dictionary
nucleic_acid_df = pd.read_csv('Nucleic_acids.csv', header=None) #read nucleic acids
nucleic_acid = dict(sorted(nucleic_acid_df.values.tolist())) #convert to dictioanry

tmp_Array = [] # global array to temporarily append results for base loss terminal
tmp_Array_prec = [] # global array to temporarily append results for other ions loss of prec
tmp_Array_frag = [] # global array to temporarily append results for other ions loss of frag
tmp_Array_int = [] # global array to temporarily append results for base loss of internal fragments


def get_base(sequence):
    '''
    Converts sequence to a  list of bases
    '''
    base = re.findall('[a-z]*[A-Z]', sequence)  # Converts sequence into a list
    # print(base)
    return base


def sequence_builder(base, nucleic_acid, ring, sample, modification, phos, end_3, end_5):
    '''
    Build an oligonucleotide sequence
    '''
    seq_list = []
    for i in base:
        seq = nucleic_acid[i] + ring + sample
        if 'dU' in i or 'dA' in i or 'dG' in i or 'dC' in i:
            seq = seq + modification
        seq_list.append(seq)
    for i in range(len(base) - 1):
        seq_list.append(phos)
    seq_list.extend([end_3, end_5])
    return seq_list


def clean_up(sequence):
    '''
    Clean up the formulas so that it looks like CXHYNZ instead of the elements repeating themselves and convert it to a dict
    '''
    result = re.findall('[A-Z][a-z]?[-]?\d*', sequence)
    a = [re.match('([A-Z][a-z]?)([-]?\d*)', i).groups() for i in result]
    d = dict(a)
    d = {x: 0 for x in d}
    for i in a:
        preValue = d.get(i[0])
        if i[0] in d:
            try:
                d[i[0]] = preValue + int(i[1])
            except ValueError:
                d[i[0]] = preValue + 1
    return d

def clean_up_counter(sequence):
    '''
    Clean up the formulas so that it looks like CXHYNZ instead of the elements repeating themselves and convert it to a dict
    '''
    result = re.findall('[A-Z][a-z]?[-]?\d*', sequence)
    a = [re.match('([A-Z][a-z]?)([-]?\d*)', i).groups() for i in result]
    d = dict(a)
    d = {x: 0 for x in d}
    for i in a:
        preValue = d.get(i[0])
        if i[0] in d:
            try:
                d[i[0]] = preValue + int(i[1])
            except ValueError:
                d[i[0]] = preValue + 1
    d = Counter(d)
    return d

def seq_part(row):
    '''
    a3 - selects the first 3 bases and determines the loss from that
    '''
    base = get_base(sequence)
    part_seq = base[:row]
    part_seq = list(dict.fromkeys(part_seq))  # remove dup bases as u cant tell which one is lost
    part_seq = list(reversed(part_seq)) # reverse it so it starts from right to left
    return part_seq


def base_loss(row):
    i = 0
    d = {}
    for each_base in row['part_seq']:
        part_seq_d = Counter()
        part_seq_d[each_base] = 1 # select the base and set the num of bases to 1
        assing_i_lab = combine_elemental_form_r(part_seq_d)
        assing_i_lab_no_int = re.findall(r'[a-zA-Z]+', assing_i_lab)[0]
        elem_comp_i = get_base_formula(part_seq_d, nucleic_acid)
        prec_lab = row['Assignement'] + ' - B(%s)'%(assing_i_lab_no_int)
        new_comp = row['Elemental_composition_dict'] - elem_comp_i #- Counter({'H': 2}) # remove 2 H from base loss
        if 'a' in row['Assignement']:
            new_comp = row['Elemental_composition_dict'] - elem_comp_i - Counter({'H':1})
        # else:
        #     new_comp = elem_com_d_c -  elem_comp_i + Counter({'H': 1})
        new_comp_mass_d = get_mass(new_comp.copy(), mass_dict)
        new_comp_mass = sum(new_comp_mass_d.values())
        mass = round((new_comp_mass - row['Charge2'] * mass_dict['Hp']) / row['Charge2'],
                     6)  # abstraction of a charge not proton therefore -2 e- i.e. add one e- then -2e- or just -1e-
        elem_comp = combine_elemental_form(new_comp)

        d[i] = {'Assignement': prec_lab, 'Elemental composition': elem_comp,
                                 'Mass': mass, 'Charge': row['Charge']}
        i = i + 1
    df = pd.DataFrame.from_dict(d, "index")
    tmp_Array.append(df)


def other_ions_loss_prec(row):
    # Charge2 is the int of Charge
    i = 0
    d = {}
    for k in l_l:
        Assignment = '%s - %s' % (row['Assignement'], k)
        if k == 'H' or k == '2H' or k == '3H':
            k, k_i = swap(k)
            Assignment = '%s - %sH' % (row['Assignement_H'], k_i + row['Charge2'])
        if k in nucleic_acid.keys():
            k = nucleic_acid[k]
        elem_loss = clean_up(k)
        elem_loss_c = Counter(elem_loss)
        new_comp = row['Elemental_composition_dict'] - elem_loss_c
        new_comp_mass_d = get_mass(new_comp.copy(), mass_dict)
        new_comp_mass = sum(new_comp_mass_d.values())
        mass = round((new_comp_mass - row['Charge2'] * mass_dict['Hp']) / row['Charge2'], 6)
        elem_comp = combine_elemental_form(new_comp)

        d[i] = {'Assignement': Assignment, 'Elemental composition': elem_comp,
                                 'Mass': mass, 'Charge': row['Charge']}
        i = i + 1
    df = pd.DataFrame.from_dict(d, "index")
    tmp_Array_prec.append(df)

def other_ions_loss_frag(row):
    # Charge2 is the int of Charge
    i = 0
    d = {}
    for k in l_l:
        Assignment = '%s - %s' % (row['Assignement'], k)
        if k == 'H' or k == '2H' or k == '4H':
            k, k_i = swap(k)
            Assignment = '%s - %sH' % (row['Assignement'], k_i)
        if k in nucleic_acid.keys():
            k = nucleic_acid[k]
        elem_loss = clean_up(k)
        elem_loss_c = Counter(elem_loss)
        new_comp = row['Elemental_composition_dict'] - elem_loss_c
        new_comp_mass_d = get_mass(new_comp.copy(), mass_dict)
        new_comp_mass = sum(new_comp_mass_d.values())
        mass = round((new_comp_mass - row['Charge2'] * mass_dict['Hp']) / row['Charge2'], 6)
        elem_comp = combine_elemental_form(new_comp)

        d[i] = {'Assignement': Assignment, 'Elemental composition': elem_comp,
                                 'Mass': mass, 'Charge': row['Charge']}
        i = i + 1
    df = pd.DataFrame.from_dict(d, "index")
    tmp_Array_frag.append(df)


def add_protons_loss_frag(row):
    # Charge2 is the int of Charge
    i = 0
    d = {}
    l_l = ['H', '2H']
    for k in l_l:
        k, k_i = swap(k)
        Assignment = '%s + %sH' % (row['Assignement'], k_i)
        elem_loss = clean_up(k)
        elem_loss_c = Counter(elem_loss)
        new_comp = row['Elemental_composition_dict'] + elem_loss_c
        new_comp_mass_d = get_mass(new_comp.copy(), mass_dict)
        new_comp_mass = sum(new_comp_mass_d.values())
        mass = round((new_comp_mass - row['Charge2'] * mass_dict['Hp']) / row['Charge2'], 6)
        elem_comp = combine_elemental_form(new_comp)
        d[i] = {'Assignement': Assignment, 'Elemental composition': elem_comp,
                                 'Mass': mass, 'Charge': row['Charge']}
        i = i + 1
    df = pd.DataFrame.from_dict(d, "index")
    tmp_Array_frag.append(df)

def get_mass(d, mass_dict):
    '''
    Insert dictionary  of element and number, returns mass of compound
    '''
    for i in d:
        if (i not in mass_dict):
            print('Atom %s not found' % i)
        d[i] = float(d[i] * mass_dict[i])
    return d

def get_mass_df(d):
    '''
    Apply get mass to df and sum the values
    '''
    for i in d:
        if (i not in mass_dict):
            print('Atom %s not found' % i)
        d[i] = float(d[i] * mass_dict[i])
    d = sum(d.values())
    return d


def update_dict(d, n):
    for key in d:
        d[key] *= n
    return d

def get_base_formula(d1,d2):
    '''
    d1 = {A: 1...}
    d2 = {A : CH3PS...}
    '''
    d4 = Counter({})
    for key, value in d1.items():
        elem_comp = d2[key]
        d3 = clean_up(elem_comp)
        d3 = Counter(update_dict(d3, value))
        d4 = d3 + d4
    return d4

def swap(k):
    '''
    Swaps the element and their numbers to be used in mass list  from assignement
    '''
    try:
        a = re.match('\d+', k).group(0)
    except:
        a = 1
    b = re.split('\d+', k)
    b = ''.join(b)
    ba = b, str(a)
    return ''.join(ba), int(a)

def match_H(row):
    '''
    Finds the number of H to do subtraction
    '''
    return re.match('[A-Z]*\d*[a-z]?\d*', row).group(0)

def get_precusror_mass(sequence, nucleic_acid, ring, sample, modification, phos, end_3, end_5, mass_dict):
    '''
    Returns df of precursor mass for each charge state
    Only return the loss of proton
    '''
    base = get_base(sequence) # returns bases from sequence
    seq_list = sequence_builder(base, nucleic_acid, ring, sample, modification, phos, end_3, end_5)
    precursor = ''.join(seq_list) # combines the elements from bases
    element_count = clean_up(precursor) # converts elem comp into dict
    mass_element = get_mass(element_count.copy(), mass_dict) # masses for each element
    mass_precursor = sum(mass_element.values()) # returns the neutral mass
    charge_list = np.arange(1, max_charges, 1) # get a list of charge states upto number of bases
    prec_lab = np.char.add(charge_list.astype(str), 'H')  # add H to number e.g. 1H, 2H
    prec_lab = np.char.add('Mp - ',prec_lab) # Add the Mp so.. Mp-1H etc
    elem_comp = combine_elemental_form(element_count) # cleaned up elem comp
    charge_H = np.multiply(mass_dict['Hp'], charge_list) # multiple charge by proton mass
    mass = np.divide(np.subtract(mass_precursor,charge_H), charge_list).round(6) # neutral mass - nH / n, n=charge
    charge_lab = np.char.add(charge_list.astype(str), '-') # convert charge to number with '-' e.g. 1 to 1-
    df = pd.DataFrame(data={'Assignement': prec_lab, 'Elemental composition': elem_comp,
                             'Mass': mass, 'Charge': charge_lab})
    return df


def sort_fragments(d, l, i, mass_dict):
    '''
    Sort the fragments into elemental list, mass list and labels
    '''
    d = clean_up(d)
    md = get_mass(d.copy(), mass_dict)
    m = sum(md.values())
    l = l + i
    return d, m, l


def combine_elemental_form(d):
    '''
    Combines the elemental formula
    '''
    new_d = dict(OrderedDict(sorted(d.items(), key=lambda t: t[0])))
    return ''.join('{}{}'.format(k, v) for k, v in new_d.items())

def combine_elemental_form_r(d):
    '''
    Combines the elemental formula reverse order
    '''
    return ''.join('{}{}'.format(v, k) for k, v in d.items())


def get_fragment_mass(mass, charge):
    '''
    Returns mass from elemental sequence for each charge
    '''
    return round((mass - charge * mass_dict['Hp']) / charge, 6)


def get_fragments(sequence, ring, sample, nucleic_acid, end, mass_dict, label):
    '''
    a/z ions are negative ions not radical, use radical function to make all ions radical
    '''
    each_lab = list(label)
    base = get_base(sequence)
    print(base)
    if 'w' in each_lab:
        base.reverse()
        # print(base)
    seq_a = ''
    frag_a_B = []
    frag_a = []
    for base_i in range(len(base) - 1):
        seq_a_B = ring + sample + seq_a + 'H-1'  # First one contains no base
        seq_a = nucleic_acid[base[base_i]] + ring + sample + seq_a  # Add base each iteration
        # No phosph linker added yet, this will be added sequentially to generate fragments
        frag_a_B.append(seq_a_B)
        if 'dU' in base[base_i] or 'dA' in base[base_i] or 'dG' in base[base_i] or 'dC' in base[base_i]:
            seq_a = modification + seq_a
        frag_a.append(seq_a)
        if len(frag_a) > 0:
            seq_a = frag_a[base_i] + phos  # Add phos linker for multiple bases
    df = pd.DataFrame()
    for i in range(len(frag_a)):
        '''
        First for loop will go through a1,b1,c1 and d1 then a2,b2,c2 and d2 ...
        '''
        label_i = str(i + 1)
        fragaz = frag_a[i] + 'H' + end # add proton
        fragby = frag_a[i] + 'OH' + end # add OH
        fragcx = frag_a[i] + phos + 'O-1' + end # remove an oxygen
        fragdw = frag_a[i] + phos + end + 'H'
        daz, maz, laz = sort_fragments(fragaz, each_lab[0], label_i, mass_dict)
        dby, mby, lby = sort_fragments(fragby, each_lab[1], label_i, mass_dict)
        dcx, mcx, lcx = sort_fragments(fragcx, each_lab[2], label_i, mass_dict)
        ddw, mdw, ldw = sort_fragments(fragdw, each_lab[3], label_i, mass_dict)
        df1 = pd.DataFrame()
        for charge in range(1):
            '''
            This loop splits everything into its own charge states
            '''
            charge += 1

            maz_i = get_fragment_mass(maz, charge)
            if maz_i < min_mass:
                # Ignores everything less than 100 mz
                continue
            mby_i = get_fragment_mass(mby, charge)
            mcx_i = get_fragment_mass(mcx, charge)
            mdw_i = get_fragment_mass(mdw, charge)

            daz_i = combine_elemental_form(daz)
            dby_i = combine_elemental_form(dby)
            dcx_i = combine_elemental_form(dcx)
            ddw_i = combine_elemental_form(ddw)

            charge_lab = str(charge) + '-'
            df2 = pd.DataFrame(
                data={'Assignement': [laz, lby, lcx, ldw], 'Elemental composition': [daz_i, dby_i, dcx_i, ddw_i],
                      'Mass': [maz_i, mby_i, mcx_i, mdw_i], 'Charge': [charge_lab, charge_lab, charge_lab, charge_lab]})
            df1 = pd.concat([df1, df2], ignore_index=True)
        df = pd.concat([df, df1], ignore_index=True)
    # print(df)
    return df

def generate_base_loss_precursor(df, base_comp, base_label, mass_dict):
    '''
    Generates base loss for precursor, need to correct for each fragment e.g. a1 will only have 1 but a2  can have 2 etc
    base_comp = elemental comp for the bases
    base_label = bases e.g. C, U, A etc...
    '''
    df['Charge2'] = df.Charge.str.extract('(\d+)').astype(int)  # new col containing charge numbers
    df = df[df['Mass'] < max_mass]
    df = df[df['Mass'] > min_mass]
    df['Elemental_composition_dict'] = df['Elemental composition'].apply(
        clean_up_counter)  # Converts elem comp into dict
    l = []
    for i, j in zip(base_comp, base_label):
        base_comp_c = Counter(i)
        new_comp = df['Elemental_composition_dict'] - base_comp_c
        apple = new_comp.copy()
        elem_comp = new_comp.apply(combine_elemental_form) # combines dict into string
        new_comp_mass_d = apple.apply(get_mass_df) # returns masses from counter
        mass = round((new_comp_mass_d - df['Charge2'] * mass_dict['Hp']) / df['Charge2'],
                     6)  # minus proton from each mass
        new_lab = df['Assignement'] + ' - ' + j

        df2 = pd.DataFrame(data={'Assignement': new_lab, 'Elemental composition': elem_comp,
                                     'Mass': mass, 'Charge': df['Charge']})
        l.append(df2)
    df1 = pd.concat(l).reset_index(drop=True)
    df1.drop_duplicates(subset='Mass', keep='first', inplace=True)
    return df1


def get_further_ions_precurosr(df, k):
    '''
    Will add additional losses to each subsequent loss, can build something like M- H2O - PO2S - H etc
    '''
    df['Charge2'] = df.Charge.str.extract('(\d+)').astype(int) # new col containing charge numbers
    # df = df[df['Mass'] < max_mass]
    df = df[df['Mass'] > min_mass]
    df.pop('Charge2') # remove the new column
    for index, row in df.iterrows():
        charge = int(re.findall('\d+', row['Charge'])[0])
        charge_lab = str(charge) + '-'
        elem_com = str(row['Elemental composition'])
        assing_i = row['Assignement']
        # print(assing_i)
        elem_com_d = clean_up(elem_com)
        elem_com_d_c = Counter(elem_com_d)
        elem_loss = clean_up(k)
        elem_loss_c = Counter(elem_loss)
        new_comp = elem_com_d_c - elem_loss_c
        new_comp_mass_d = get_mass(new_comp.copy(), mass_dict)
        new_comp_mass = sum(new_comp_mass_d.values())
        mass = round((new_comp_mass + charge * mass_dict['elec']) / charge, 6)
        Assignment = '%s - %s' % (assing_i, k)
        elem_comp = combine_elemental_form(new_comp)
        df1 = pd.DataFrame(data={'Assignement': [Assignment], 'Elemental composition': elem_comp,
                                 'Mass': mass, 'Charge': charge_lab})
        df = df.append(df1, ignore_index=True)
    # print(df)
    return df


def get_other_ions_precursor(df, sequence):
    '''
    Returns other ions such as loss of proton  for precursor
    This is not ion activation i.e. M-nH20-nH but single loss of a species e.g. M-H20
    '''
    # print(df)
    base = get_base(sequence)
    d_base = {i: base.count(i) for i in base}
    base_comp = []
    base_label = []
    for key, value in d_base.items():
        if value > 1:
            for v in range(value):
                v += 1
                form = nucleic_acid[key]
                form_comp = clean_up(form)
                form_comp['H'] = form_comp['H']
                form_comp.update({n: v * form_comp[n] for n in form_comp.keys()})
                base_comp.append(form_comp)
                if v == 1:
                    label_base = key
                else:
                    label_base = str(v) + key
                base_label.append(label_base)
        else:
            form = nucleic_acid[key]
            label_base = key
            form_comp = clean_up(form)
            base_label.append(label_base)
            base_comp.append(form_comp)
    df['Charge2'] = df.Charge.str.extract('(\d+)').astype(int) # new col containing charge numbers
    # df = df[df['Mass'] < max_mass]
    df = df[df['Mass'] > min_mass]
    df['Elemental_composition_dict'] = df['Elemental composition'].apply(clean_up_counter) # Converts elem comp into dict
    df['Assignement_H'] = df['Assignement'].apply(match_H)
    df.apply(other_ions_loss_prec, axis=1)
    df1 = pd.concat(tmp_Array_prec).reset_index(drop=True)
    df1.drop_duplicates(subset='Mass', keep='first', inplace=True)
    return df1, base_comp, base_label



def get_other_ions_frag(df):
    '''
    Returns other ions such as loss of proton  for each fragment
    This is not ion activation i.e. M-nH20-nH but single loss of a species e.g. M-H20
    '''
    df['Charge2'] = df.Charge.str.extract('(\d+)').astype(int) # new col containing charge numbers
    df = df[df['Mass'] < max_mass]
    df['Elemental_composition_dict'] = df['Elemental composition'].apply(clean_up_counter) # Converts elem comp into dict
    df['Assignement_H'] = df['Assignement'].apply(match_H) # gets the label e.g. a from a4
    df_2 = df.copy()


    df_2.apply(add_protons_loss_frag, axis=1)
    df.apply(other_ions_loss_frag, axis=1)
    df1 = pd.concat(tmp_Array_frag).reset_index(drop=True)
    return df1


def add_salts(df, mass_dict):
    l_l = ['Na', '2Na', '3Na', '4Na', 'K', '2K', '3K', '4K']
    H_dict = {}
    for index, row in df.iterrows():
        charge = int(re.findall('\d+', row['Charge'])[0])
        charge_lab = str(charge) + '-'
        elem_com = str(row['Elemental composition'])
        assing_i = row['Assignement']
        assing_lab = re.match('[A-Z]*[a-z]*\d*', assing_i).group(0)  # assign lab is M or a or w ...
        elem_com_d = clean_up(elem_com)
        elem_com_d_c = Counter(elem_com_d)
        df1 = pd.DataFrame()
        for k in l_l:
            Assignment = '%s + %s' % (assing_lab, k)  # k is the element from the list
            k, k_i = swap(k)  # k_i is the number of atoms added
            elem_loss = clean_up(k)
            elem_loss_c = Counter(elem_loss)
            H_dict['H'] = k_i
            new_comp = elem_com_d_c + elem_loss_c - Counter(H_dict)
            new_comp_mass_d = get_mass(new_comp.copy(), mass_dict)
            new_comp_mass = sum(new_comp_mass_d.values())
            # Generate M- ions only. Technically they're all M-H but -H accounted for in elemental comp
            #
            mass = round((new_comp_mass + charge * mass_dict['elec']) / charge, 6)

            elem_comp = combine_elemental_form(new_comp)
            df2 = pd.DataFrame(data={'Assignement': [Assignment], 'Elemental composition': elem_comp,
                                     'Mass': mass, 'Charge': charge_lab})
            df1 = pd.concat([df1, df2], ignore_index=True)
        df = df.append(df1, ignore_index=True)
        # print(df)
        return df


def ion_loss_mp_charged(df, mass_dict, atom):
    '''
    Charge reduced loss of a species, [M-5H]5- to [M-NCO]4-
    Only works for multiply charged ions
    Abstraction of charge therefore -2 e-.
    '''
    for index, row in df.iterrows():
        charge = int(re.findall('\d+', row['Charge'])[0])
        if charge == 1: # If its -1, loss of an ion would neutralise it
            continue
        else:
            charge -= 1 # Add one for each charge i.e. -2 goes to -1
        charge_lab = str(charge) + '-'
        elem_com = str(row['Elemental composition'])
        assing_i = row['Assignement']
        # assing_lab = re.match('[A-Z]*[a-z]*\d*', assing_i).group(0)  # assign lab is M or a or w ...
        # print(assing_i)
        elem_com_d = clean_up(elem_com)
        elem_com_d_c = Counter(elem_com_d)
        prec_lab = '%s - %s%s' % (assing_i,charge, atom)
        elem_loss = clean_up(atom)
        elem_loss_c = Counter(elem_loss)
        new_comp = elem_com_d_c - elem_loss_c
        new_comp_mass_d = get_mass(new_comp.copy(), mass_dict)
        new_comp_mass = sum(new_comp_mass_d.values())
        mass = round((new_comp_mass - charge * mass_dict['elec']) / charge, 6) # abstraction of a charge not proton therefore -2 e- i.e. add one e- then -2e- or just -1e-
        elem_comp = combine_elemental_form(new_comp)
        df1 = pd.DataFrame(data={'Assignement': [prec_lab], 'Elemental composition': elem_comp,
                                 'Mass': mass, 'Charge': charge_lab})
        df = df.append(df1, ignore_index=True)
    # print(df)
    return df

def base_loss_fragments_terminal(df, sequence):
    '''
    Loss of last base only
    '''
    base = get_base(sequence)
    if df['Assignement'].str.contains('w|x|y|z').any():
        base.reverse()
    # part_seq_d = Counter()
    df['Charge2'] = df.Charge.str.extract('(\d+)').astype(int) # new col containing charge numbers
    df = df[df['Charge2'] < max_charges] # filter charge to select all charges less than ...
    # df = df[df['Mass'] < max_mass]
    df = df[df['Mass'] > min_mass]
    df['Elemental_composition_dict'] = df['Elemental composition'].apply(clean_up_counter) # Converts elem comp into dict
    df['Assignement_val'] = df.Assignement.str.extract('(\d+)').astype(int) # get cleavage number e.g. a3 returns 3
    df['part_seq'] = df['Assignement_val'].apply(seq_part) # get number of bases upto a3
    df.apply(base_loss, axis=1) # calculates each base loss
    df.pop('Charge2')
    df = pd.concat(tmp_Array).reset_index(drop=True) # combine base loss
    df.drop_duplicates(subset='Elemental composition', keep='first', inplace=True) #Removes duplicates but keeps first, a7-T7 >>> a7-T4
    # print(df)
    return df

def opp_symp_abcd(row):
    # gets the complementary cleavage
    if row == 'a':
        row = 'w'
    if row == 'b':
        row = 'x'
    if row == 'c':
        row = 'y'
    if row == 'd':
        row = 'z'
    return row

def opp_symp_wxyz(row):
    # gets the complementary cleavage
    if row == 'w':
        row = 'a'
    if row == 'x':
        row = 'b'
    if row == 'y':
        row = 'c'
    if row == 'z':
        row = 'd'
    return row


def shift_df(df, shift):
    return df.reindex(np.roll(df.index, shift))

def diff_df_func(row):
    return row['Elemental_composition_dict_grp'] - row['Elemental_composition_dict']


def internal_seq_part(row):
    '''
    a3w6 - selects the first 3-6 bases
    '''
    base = get_base(sequence)
    a = row['Assignement_val']
    b = row['R_Assignement_val']
    part_seq = base[int(a):int(b)]
    part_seq = list(dict.fromkeys(part_seq))  # remove dup bases as u cant tell which one is lost
    part_seq = list(reversed(part_seq)) # reverse it so it starts from right to left
    return part_seq

def internal_base_loss(row):
    i = 0
    d = {}
    for each_base in row['part_seq']:
        part_seq_d = Counter()
        part_seq_d[each_base] = 1 # select the base and set the num of bases to 1
        assing_i_lab = combine_elemental_form_r(part_seq_d)
        assing_i_lab_no_int = re.findall(r'[a-zA-Z]+', assing_i_lab)[0]
        elem_comp_i = get_base_formula(part_seq_d, nucleic_acid)
        prec_lab = row['Assignement'] + ' - B(%s)'%(assing_i_lab_no_int)
        new_comp = row['Elemental_composition_dict'] - elem_comp_i #- Counter({'H': 2}) # remove 2 H from base loss
        if 'a' in row['Assignement']:
            new_comp = row['Elemental_composition_dict'] - elem_comp_i - Counter({'H':1})
        # else:
        #     new_comp = elem_com_d_c -  elem_comp_i + Counter({'H': 1})
        new_comp_mass_d = get_mass(new_comp.copy(), mass_dict)
        new_comp_mass = sum(new_comp_mass_d.values())
        mass = round((new_comp_mass - row['Charge2'] * mass_dict['Hp']) / row['Charge2'],
                     6)  # abstraction of a charge not proton therefore -2 e- i.e. add one e- then -2e- or just -1e-
        elem_comp = combine_elemental_form(new_comp)

        d[i] = {'Assignement': prec_lab, 'Elemental composition': elem_comp,
                                 'Mass': mass, 'Charge': row['Charge']}
        i = i + 1

    df = pd.DataFrame.from_dict(d, "index")
    tmp_Array_int.append(df)

def collate_int_frag_abcd(row):
    '''
    Collates multiple internal fragments of single bases into 1
    these consist of e + base-loss + v where e/v represent a/w frag grouped together
    number does not matter here as its a single base loss
    '''
    base = get_base(sequence)
    base_name = base[row['Base_val']]
    return base_name

def collate_int_frag_wxyz(row):
    '''
    Collates multiple internal fragments of single bases into 1
    these consist of e + base-loss + v where e/v represent a/w frag grouped together
    number does not matter here as its a single base loss
    '''
    base = get_base(sequence)
    base.reverse()
    base_name = base[row['Base_val']]
    return  base_name

def get_internal_fragments(df, type):
    '''
    for a sequence of len 10, a1, a2 .... a9
    for a a9w1 fragment, you need to do a9w9 = a9 - a1
    cant have a9w1 as they are the same cleavage (complementary)
    a2w15 = b2x15 = c2y15 = d2z15 because they contain same elements from phos group
    16 bases means that a2w15 is equivalent to 1 base. there wll be loads like this so need to group them by base
    some bases can be same masses e.g. A + O = G
    '''
    base = get_base(sequence)
    #sequence, ring, sample, nucleic_acid, end, mass_dict, label
    df['Charge2'] = df.Charge.str.extract('(\d+)').astype(int) # new col containing charge numbers
    df = df[df['Charge2'] == 1] # reset charge states to 1
    df['Elemental_composition_dict'] = df['Elemental composition'].apply(clean_up_counter) # Converts elem comp into dict
    df['Assignement_val'] = df.Assignement.str.extract('(\d+)').astype(int) # get cleavage number e.g. a3 returns 3
    assignement_values = list(dict.fromkeys(df['Assignement_val'].to_list())) # list of assigned values
    df['Assignement_char'] = df.Assignement.str.extract('(\w?)').astype(str) # get cleavage number e.g. a3 returns 3
    df['R_Assignement_val'] = len(base) - df['Assignement_val']  # gets the reverse number
    df['R_Assignement_val'] =  df['R_Assignement_val'].astype(str)
    if type == 'a': # check if its abcd or not
        df['R_Assignement_char'] = df['Assignement_char'].apply(opp_symp_abcd)
    else:
        df['R_Assignement_char'] = df['Assignement_char'].apply(opp_symp_wxyz)
        base.reverse()
    grouped = df.groupby("Assignement_val") # groups data by the cleavage number e.g. a3 will return a3,b3,c3,d3 etc

    # print(grouped.get_group(2))
    l2 = []
    for i in assignement_values:
        cleavage_chk = df['Assignement_val'] - i # len(base) - 1,2,3...
        diff_df = df[cleavage_chk < 0] # only selects smallest bases, num of cleavages = num of bases - 1

        if len(diff_df) ==0:
            continue
        df_group = grouped.get_group(i)
        # print(df_group[['Assignement', 'Assignement_val', 'R_Assignement_char', 'R_Assignement_val']])
        l = []
        size = int(len(diff_df)/len(df_group)) # find the size diff between group and orig df
        df_group = pd.concat(size*[df_group]) # make the 2 dfs same size
        # print(df_group[['Assignement','Assignement_val','R_Assignement_char','R_Assignement_val']])
        for j in range(len(diff_df)):
            df_shift = shift_df(diff_df,j) # shift the df so you can subtract each row
            df_shift['Elemental_composition_dict_grp'] = df_group['Elemental_composition_dict'].values
            df_shift['Elemental_composition_dict_diff'] = df_shift.apply(diff_df_func, axis=1)
            df_shift['Assignement'] = df_group['Assignement'].values +  df_shift['R_Assignement_char'] + df_shift['R_Assignement_val']
            df_shift['Assignement_val'] =  df_group['Assignement_val'].values
            df_shift['Assignement_char'] = df_group['Assignement_char'].values
            # print(df_shift['A ssignement'])
            l.append(df_shift)
        df_conc = pd.concat(l,ignore_index=True, sort=True).reset_index(drop=True)

        l2.append(df_conc)
    df_int = pd.concat(l2, ignore_index=True, sort=True).reset_index(drop=True)
    df_int['int_lab'] = df_int['Assignement'].replace('(\d)', '', regex=True) # extracts the cleavage letters e.g. d3w1 = dw
    # select subsets of df
    # remove the H based on chemdraw
    '''
    A3y15 has extra H (-H)

    b3Y15 has extra H (-H)

    C3W15 less 1H (+H)

    C3X15 less 1H (+H)

    C3Z15 less 1H (+H)

    D3y15 extra H (-H)
    '''
    df_int.loc[df_int.int_lab == 'ay', 'Elemental_composition_dict_diff'] = df_int['Elemental_composition_dict_diff'] - Counter({'H':1})
    df_int.loc[df_int.int_lab == 'by', 'Elemental_composition_dict_diff'] = df_int['Elemental_composition_dict_diff'] - Counter({'H': 1})
    df_int.loc[df_int.int_lab == 'cw', 'Elemental_composition_dict_diff'] = df_int['Elemental_composition_dict_diff'] + Counter({'H':1})
    df_int.loc[df_int.int_lab == 'cx', 'Elemental_composition_dict_diff'] = df_int['Elemental_composition_dict_diff'] + Counter({'H': 1})
    df_int.loc[df_int.int_lab == 'cz', 'Elemental_composition_dict_diff'] = df_int['Elemental_composition_dict_diff'] + Counter({'H':1})
    df_int.loc[df_int.int_lab == 'dy', 'Elemental_composition_dict_diff'] = df_int['Elemental_composition_dict_diff'] - Counter({'H': 1})

    # need to make a copy otherwise you cant use both functions
    apple = df_int['Elemental_composition_dict_diff']
    apple2 = apple.copy()
    df_int['Elemental composition'] = apple.apply(combine_elemental_form)
    df_int['Mass'] = apple2.apply(get_mass_df)
    df_int['Mass'] = round(df_int['Mass'] - mass_dict['Hp'],6) # Make it -1 charge
    df_int.drop_duplicates(subset='Assignement', keep='first', inplace=True)

    df_int['tot_ass_val'] = df_int['Assignement_val'].astype(int) + df_int['R_Assignement_val'].astype(int)
    cleavage_length = len(base) + 1  # select enough cleavages to equal a single base
    sub_df  = df_int.loc[df_int.tot_ass_val == cleavage_length] # take subset eqal to a single base
    df_int = df_int.loc[~df_int.tot_ass_val == cleavage_length]

    sub_df['Base_val'] = sub_df['Assignement_val'] -1 #extract the base, starts from 0
    group_mass = sub_df.groupby("Mass")
    uniq_mass = sub_df.Mass.unique()
    uniqu_l = []
    for mass in uniq_mass:
        df_uniq = group_mass.get_group(mass)
        if len(df_uniq) ==1: #only 1 mass, then leave it
            uniqu_l.append(df_uniq)
            continue
        if type == 'a':
            df_uniq['base'] = df_uniq.apply(collate_int_frag_abcd, axis=1) # combine the multiple single internal fragments
        else:
            df_uniq['base'] = df_uniq.apply(collate_int_frag_wxyz, axis=1)

        df_uniq['Assignement'] = df_uniq['Assignement_char'] + df_uniq['Assignement_val'].astype(str) + '(' + df_uniq['base']+')' +  df_uniq['R_Assignement_char'] + df_uniq['R_Assignement_val'].astype(str)
        uniqu_l.append(df_uniq)

    df_uniq = pd.concat(uniqu_l, ignore_index=True, sort=True).reset_index(drop=True)
    df_uniq.drop_duplicates(subset='Assignement', keep='first', inplace=True) # drop the duplicates now

    df_int = pd.concat([df_int,df_uniq], ignore_index=True, sort=True)
    df_int = df_int[['Assignement','Elemental composition', 'Mass', 'Charge']]
    df_int = df_int[df_int['Mass'] > 300]
    return df_int

def add_charges(df):
    df['Charge2'] = df.Charge.str.extract('(\d+)').astype(int) # new col containing charge numbers
    charge_l = []
    charge_l.append(df) # append the original df
    max_charges2 = max_charges + 1 # add 1 to be used in range func
    for i in range(2,max_charges2,1): # start from 2
        df_charge = df.copy()
        df_charge['Charge2'] = i
        df_charge['Charge'] = str(i)+'-'
        charge_l.append(df_charge)
    df = pd.concat(charge_l).reset_index(drop=True)
    df['Mass'] = df['Mass'] +mass_dict['Hp'] # accounted for -1 charge already, therefore make it neutral
    df['Mass'] = round((df['Mass'] - df['Charge2']*mass_dict['Hp'])/df['Charge2'],6)
    df = df[df['Mass']>min_mass]
    return df


def generate_peak_list():
    global l_l, tmp_Array,tmp_Array_int

    if calc_neutral == 'y':
        l_l = neutral_losses
    else:
        l_l = []
    if calc_fragments !='y':
        df_prec = get_precusror_mass(sequence, nucleic_acid, ring, sample, modification, phos, end_3, end_5, mass_dict)
        df_prec.to_csv('Precursor_mass.csv',index=None)
        print(df_prec)
    else:
        # calculate precursor masses for each charge state
        df_prec = get_precusror_mass(sequence, nucleic_acid, ring, sample, modification, phos, end_3, end_5, mass_dict)
        print(df_prec)

        #calculate mcluckey fragments and internal fragments
        df_frag_abcd = get_fragments(sequence, ring, sample, nucleic_acid, end_5, mass_dict, 'abcd')  # zyxw
        df_frag_wxyz = get_fragments(sequence, ring, sample, nucleic_acid, end_3, mass_dict, 'zyxw')


        #create empty dfs to be filled with conditions specified above
        #precursors
        df_prec_loss = df_prec.drop(df_prec.index)
        df_prec_base_loss = df_prec.drop(df_prec.index)
        #fragments
        df_frag_aB =df_frag_abcd.drop(df_frag_abcd.index)
        df_frag_wB =df_frag_abcd.drop(df_frag_abcd.index)
        df_int_abcd =df_frag_abcd.drop(df_frag_abcd.index)
        df_int_wxyz = df_frag_abcd.drop(df_frag_abcd.index)

        #base loss for mcluckey fragments and precursors
        if calc_base == 'y':
            df_frag_aB = base_loss_fragments_terminal(df_frag_abcd, sequence)
            df_frag_wB = base_loss_fragments_terminal(df_frag_wxyz, sequence) # a and w fragments are merging
        if calc_internal =='y':
            df_int_abcd = get_internal_fragments(df_frag_abcd,'a')
            df_int_wxyz = get_internal_fragments(df_frag_wxyz,'w')

        #combine all the fragments into one df
        df_frag = pd.concat([df_frag_abcd,df_frag_wxyz,df_frag_aB,df_frag_wB,df_int_abcd,df_int_wxyz], axis=0, ignore_index=True, sort=True).reset_index(drop=True)
        df_frag.drop_duplicates(subset='Assignement', keep='first', inplace=True) # drop the duplicates from int frag of abcd - wxyz
        df_frag = add_charges(df_frag)

        if calc_neutral =='y':
            df_prec_loss, base_comp, base_label = get_other_ions_precursor(df_prec, sequence)
            df_prec_base_loss = generate_base_loss_precursor(df_prec_loss, base_comp, base_label, mass_dict)

            df_frag_Io = get_other_ions_frag(df_frag)
            df_frag_Io = pd.concat([df_frag_Io,df_frag], axis=0, ignore_index=True, sort=True).reset_index(drop=True) #
        else:
            df_frag_Io =df_frag
        # precursor and fragments
        df_all = pd.concat([df_prec,df_prec_loss,df_prec_base_loss,df_frag_Io], sort=True, ignore_index=True).reset_index(drop=True) # df_prec_base_loss,
        df_all.sort_values(by='Mass', inplace=True)
        df_all = df_all[df_all['Mass']>min_mass]
        df_all.pop('Charge2')
        #save the mass list
        df_all.to_csv('Theoretical_masses.csv', index=None)
        print(df_all.sort_values(by='Mass').head(50))
        print('Total number of fragments %d'%(len(df_all)))

generate_peak_list()

end = time.time()
print('Total time taken %.1fs'%(end - start))
