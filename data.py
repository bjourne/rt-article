import numpy as np

data = {}
data['comp1', 'spo6', 'lucy1'] = [
    ['embree', 'n', [20.675]],
    ['embree', 'y', [44.555]],
    ['hh',     'n', [20.314]],
    ['hh',     'y', [50.823]],
    ['hh2',    'n', [19.600]],
    ['hh2',    'y', [52.305]],
    ['sf01',   'n', [19.593]],
    ['sf01',   'y', [48.557]],
    ['mt',     'n', [20.478]],
    ['mt',     'y', [52.540]],
    ['bw12',   'n', [19.855]],
    ['bw12',   'y', [50.741]],
    ['bw9',    'n', []],
    ['bw9',    'y', [50.608]],
    ['shev',   'n', [17.176]],
    ['shev',   'y', [53.112]],
    ['ds',     'n', [18.863]],
    ['ds',     'y', [44.779]]
]

data['comp1', 'spo6', 'lucy2'] = [
    ['embree', 'n', [38.680, 38.519]],
    ['embree', 'y', [62.914, 63.036]],
    ['hh',     'n', [38.098, 37.664]],
    ['hh',     'y', [73.297, 73.640]],
    ['hh2',    'n', [36.661, 36.488]],
    ['hh2',    'y', [75.403, 75.012]],
    ['sf01',   'n', [37.332, 36.622]],
    ['sf01',   'y', [69.444, 69.658]],
    ['mt',     'n', [38.014, 38.101]],
    ['mt',     'y', [75.129, 75.503]],
    ['bw12',   'n', [37.713, 37.203]],
    ['bw12',   'y', [72.561, 72.501]],
    ['bw9',    'n', [25.839]],
    ['bw9',    'y', [72.183, 71.996]],
    ['shev',   'n', [32.386, 32.892]],
    ['shev',   'y', [74.630, 75.271]],
    ['ds',     'n', [35.789, 36.248]],
    ['ds',     'y', [64.187, 65.068]]
    ]

data['comp1', 'spo1', 'sanm1'] = [
    ['embree', 'n', [12.506]],
    ['embree', 'y', [21.516]],
    ['hh',     'n', [12.293]],
    ['hh',     'y', [25.035]],
    ['hh2',    'n', [11.866]],
    ['hh2',    'y', [25.129]],
    ['sf01',   'n', [11.767]],
    ['sf01',   'y', [23.506]],
    ['mt',     'n', [12.139]],
    ['mt',     'y', [25.400]],
    ['bw12',   'n', [12.219]],
    ['bw12',   'y', [24.728]],
    ['bw9',    'n', [9.053]],
    ['bw9',    'y', [24.311]],
    ['shev',   'n', [10.383]],
    ['shev',   'y', [24.932]],
    ['ds',     'n', [11.579]],
    ['ds',     'y', [21.716]]
    ]

data['comp1', 'spe11', 'sanm2'] = [
    ['embree', 'n', [3.423]],
    ['embree', 'y', [7.032]],
    ['hh',     'n', [3.407]],
    ['hh',     'y', [8.000]],
    ['hh2',    'n', [3.326]],
    ['hh2',    'y', [8.002]],
    ['sf01',   'n', [3.284]],
    ['sf01',   'y', [7.534]],
    ['mt',     'n', [3.373]],
    ['mt',     'y', [8.174]],
    ['bw12',   'n', [3.376]],
    ['bw12',   'y', [7.927]],
    ['bw9',    'n', [2.569]],
    ['bw9',    'y', [7.756]],
    ['shev',   'n', [2.901]],
    ['shev',   'y', [8.027]],
    ['ds',     'n', [3.176]],
    ['ds',     'y', [7.021]]
    ]

data['comp1', 'spe5', 'crown'] = [
    ['embree', 'n', [7.388]],
    ['embree', 'y', [13.011]],
    ['hh',     'n', [7.273]],
    ['hh',     'y', [15.306]],
    ['hh2',    'n', [7.169]],
    ['hh2',    'y', [15.432]],
    ['sf01',   'n', [7.092]],
    ['sf01',   'y', [14.099]],
    ['mt',     'n', [7.374]],
    ['mt',     'y', [15.539]],
    ['bw12',   'n', [7.138]],
    ['bw12',   'y', [14.992]],
    ['bw9',    'n', [5.228]],
    ['bw9',    'y', [14.998]],
    ['shev',   'n', [6.195]],
    ['shev',   'y', [15.648]],
    ['ds',     'n', [6.875]],
    ['ds',     'y', [13.108]]
    ]

def print_line(algo, rates, div):
    if not rates:
        fps = '?'
    else:
        fps = np.mean(rates)
        fps = fps / div
        # if packets == 'y':
        #     fps = fps / embreeK
        # else:
        #     fps = fps / embree1
        fps = int(round(fps*100))
    print('%6s %7s' % (algo, fps))

def compute_pcts(comp, id, algo):
    seq = data[comp, id, algo]
    header = 'single rays: %s, %s, %s' % (comp, id, algo)
    print(header)
    print('=' * len(header))
    assert seq[0][0] == 'embree'
    assert seq[1][0] == 'embree'
    embree1 = np.mean(seq[0][2])
    print('%6s %7s' % (seq[0][0], round(embree1, 2)))
    for algo, packets, rates in seq[2:]:
        if packets == 'n':
            print_line(algo, rates, embree1)
    header = 'ray packets: %s, %s, %s' % (comp, id, algo)
    print(header)
    print('=' * len(header))
    embreeK = np.mean(seq[1][2])
    print('%6s %7s' % (seq[1][0], round(embreeK, 1)))
    for algo, packets, rates in seq[2:]:
        if packets == 'y':
            print_line(algo, rates, embreeK)

if __name__ == '__main__':
    compute_pcts('comp1', 'spo1', 'sanm1')
