import numpy as np

data = {}

data['comp1', 'spe3', 'lucy1'] = [
    ['embree',  'n', [19.249]],
    ['embree',  'y', [41.448]],
    ['embree2', 'n', [19.024]],
    ['embree2', 'y', [40.806]],
    ['hh',      'n', [18.852]],
    ['hh',      'y', [47.837]],
    ['hh2',     'n', [18.472]],
    ['hh2',     'y', [48.255]],
    ['sf01',    'n', [18.174]],
    ['sf01',    'y', [44.847]],
    ['mt',      'n', [19.009]],
    ['mt',      'y', [48.908]],
    ['bw12',    'n', [18.467]],
    ['bw12',    'y', [46.584]],
    ['bw9',     'n', [13.470]],
    ['bw9',     'y', [46.698]],
    ['shev',    'n', [15.595]],
    ['shev',    'y', [49.222]],
    ['ds',      'n', [17.577]],
    ['ds',      'y', [41.826]]
    ]

data['comp1', 'spo8', 'lucy2'] = [
    ['embree',  'n', [35.805]],
    ['embree',  'y', [58.251]],
    ['embree2', 'n', [35.480]],
    ['embree2', 'y', [57.936]],
    ['bw9',     'n', [23.942]],
    ['bw9',     'y', [66.728]],
    ['bw12',    'n', [34.706]],
    ['bw12',    'y', [67.331]],
    ['ds',      'n', [33.582]],
    ['ds',      'y', [60.375]],
    ['hh',      'n', [35.208]],
    ['hh',      'y', [58.809]],
    ['hh2',     'n', [34.168]],
    ['hh2',     'y', [70.037]],
    ['mt',      'n', [35.548]],
    ['mt',      'y', [70.003]],
    ['sf01',    'n', [34.649]],
    ['sf01',    'y', [64.638]],
    ['shev',    'n', [30.248]],
    ['shev',    'y', [69.932]],
    ]

data['comp1', 'spo9', 'sanm1'] = [
    ['embree',  'n', [13.012]],
    ['embree',  'y', [22.598]],
    ['embree2', 'n', [12.839]],
    ['embree2', 'y', [22.388]],
    ['bw9',     'n', [9.526]],
    ['bw9',     'y', [25.112]],
    ['bw12',    'n', [12.808]],
    ['bw12',    'y', [25.672]],
    ['ds',      'n', [12.179]],
    ['ds',      'y', [22.985]],
    ['hh',      'n', [12.803]],
    ['hh',      'y', [26.163]],
    ['hh2',     'n', [12.638]],
    ['hh2',     'y', [26.022]],
    ['mt',      'n', [12.942]],
    ['mt',      'y', [26.663]],
    ['sf01',    'n', [12.474]],
    ['sf01',    'y', [24.598]],
    ['shev',    'n', [10.889]],
    ['shev',    'y', [25.928]],

data['comp2', 'x', 'crown'] = [
    ['embree', 'n', [3.945]],
    ['embree', 'y', [6.367]]
    ]

data['comp1', 'spo8', 'sanm2'] = [
    ['embree',  'n', [3.343]],
    ['embree',  'y', [6.847]],
    ['embree2', 'n', [3.348]],
    ['embree2', 'y', [6.805]],
    ['bw9',     'n', [2.508]],
    ['bw9',     'y', [7.576]],
    ['bw12',    'n', [3.300]],
    ['bw12',    'y', [7.731]],
    ['ds',      'n', [3.086]],
    ['ds',      'y', [6.860]],
    ['hh',      'n', [3.320]],
    ['hh',      'y', [7.864]],
    ['hh2',     'n', [3.265]],
    ['hh2',     'y', [7.808]],
    ['mt',      'n', [3.275]],
    ['mt',      'y', [7.984]],
    ['sf01',    'n', [3.220]],
    ['sf01',    'y', [7.348]],
    ['shev',    'n', [2.832]],
    ['shev',    'y', [7.844]],
    ]

data['comp1', 'spo8', 'crown'] = [
    ['embree',  'n', [7.643]],
    ['embree',  'y', [13.601]],
    ['embree2', 'n', [7.490]],
    ['embree2', 'y', [13.410]],
    ['bw9',     'n', [5.415]],
    ['bw9',     'y', [15.318]],
    ['bw12',    'n', [7.526]],
    ['bw12',    'y', [15.660]],
    ['ds',      'n', [7.086]],
    ['ds',      'y', [13.808]],
    ['hh',      'n', [7.600]],
    ['hh',      'y', [15.926]],
    ['hh2',     'n', [7.449]],
    ['hh2',     'y', [16.092]],
    ['mt',      'n', [7.642]],
    ['mt',      'y', [16.249]],
    ['sf01',    'n', [7.319]],
    ['sf01',    'y', [14.700]],
    ['shev',    'n', [6.279]],
    ['shev',    'y', [16.196]],
    ]

# 500 frames
data['comp1', 'spo6', 'dragon'] = [
    ['embree',  'n', [0.4120]],
    ['embree',  'y', [1.1004]],
    ['embree2', 'n', [0.4168]],
    ['embree2', 'y', [1.0952]],
    ['bw9',     'n', [0.3137]],
    ['bw9',     'y', [1.1402]],
    ['bw12',    'n', [0.4060]],
    ['bw12',    'y', [1.1466]],
    ['ds',      'n', [0.3889]],
    ['ds',      'y', [1.0640]],
    ['hh',      'n', [0.4133]],
    ['hh',      'y', [1.2012]],
    ['hh2',     'n', [0.4015]],
    ['hh2',     'y', [1.2203]],
    ['mt',      'n', [0.4160]],
    ['mt',      'y', [1.2264]],
    ['sf01',    'n', [0.3891]],
    ['sf01',    'y', [1.0300]],
    ['shev',    'n', [0.3573]],
    ['shev',    'y', [1.2394]],
    ]


def print_line(algo, rates, div):
    if not rates:
        fps = '?'
    else:
        fps = np.mean(rates)
        fps = fps / div
        fps = int(round(fps*100))
    print('%8s %7s' % (algo, fps))

def compute_pcts(comp, id, view):
    seq = data[comp, id, view]
    header = 'single rays: %s, %s, %s' % (comp, id, view)
    print(header)
    print('=' * len(header))
    assert seq[0][0] == 'embree'
    assert seq[1][0] == 'embree'
    embree1 = np.mean(seq[0][2])
    print('%8s %7s' % (seq[0][0], round(embree1, 2)))
    for algo, packets, rates in seq[2:]:
        if packets == 'n':
            print_line(algo, rates, embree1)
    header = 'ray packets: %s, %s, %s' % (comp, id, view)
    print(header)
    print('=' * len(header))
    embreeK = np.mean(seq[1][2])
    print('%8s %7s' % (seq[1][0], round(embreeK, 2)))
    for algo, packets, rates in seq[2:]:
        if packets == 'y':
            print_line(algo, rates, embreeK)

if __name__ == '__main__':
    compute_pcts('comp1', 'spo6', 'dragon')
