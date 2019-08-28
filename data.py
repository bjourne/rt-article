from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

data = {}

data['comp1', 'spe3', 'lucy1'] = [
    ['embree',  'n', [19.249]],
    ['embree',  'y', [41.448]],
    ['embree2', 'n', [19.024]],
    ['embree2', 'y', [40.806]],
    ['bw9',     'n', [13.470]],
    ['bw9',     'y', [46.698]],
    ['bw12',    'n', [18.467]],
    ['bw12',    'y', [46.584]],
    ['hh',      'n', [18.852]],
    ['hh',      'y', [47.837]],
    ['hh2',     'n', [18.472]],
    ['hh2',     'y', [48.255]],
    ['sf01',    'n', [18.174]],
    ['sf01',    'y', [44.847]],
    ['mt',      'n', [19.009]],
    ['mt',      'y', [48.908]],
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
    ['shev',    'y', [25.928]]
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


########################################################################
data['comp2', 'x', 'crown'] = [
    # 10k frames
    ['embree', 'n', [3.945]],
    ['embree', 'y', [6.367]],
    # 10k frames
    ['embree2',     'n', [3.896]],
    ['embree2',     'y', [6.444]],
    # 10k frames
    ['bw9',     'n', [2.820]],
    ['bw9',     'y', [7.248]],
    # 10k frames
    ['bw12',     'n', [3.837]],
    ['bw12',     'y', [7.266]],
    # 10k frames
    ['ds',     'n', [3.585]],
    ['ds',     'y', [6.036]],
    # 5k frames
    ['hh',     'n', [3.888]],
    ['hh',     'y', [7.420]],
    # 5k frames
    ['hh2',    'n', [3.826]],
    ['hh2',    'y', [7.539]],
    # 10k frames
    ['mt',     'n', [3.903]],
    ['mt',     'y', [7.440]],
    # 5k frames
    ['sf01',   'n', [3.695]],
    ['sf01',   'y', [6.388]],
    # 10k frames
    ['shev',   'n', [3.332]],
    ['shev',   'y', [7.639]],
    ]

# 5k frames
data['comp2', 'x', 'lucy1'] = [
    ['embree', 'n', [10.369]],
    ['embree', 'y', [20.838]],
    ['hh',     'n', [10.102]],
    ['hh',     'y', [23.460]],
    ['hh2',    'n', [9.943]],
    ['hh2',    'y', [23.834]],
    ['sf01',   'n', [9.603]],
    ['sf01',   'y', [20.477]],
    ['mt',     'n', [9.977]],
    ['mt',     'y', [23.713]],
    ['bw12',   'n', [9.906]],
    ['bw12',   'y', [23.385]],
    ['bw9',   'n', [7.423]],
    ['bw9',   'y', [23.111]],
    ['shev',   'n', [8.758]],
    ['shev',   'y', [24.990]],
    ['ds',   'n', [9.245]],
    ['ds',   'y', [20.033]],
    ['embree2', 'n', [10.375]],
    ['embree2', 'y', [20.989]]
]

# 5k frames
data['comp2', 'x', 'lucy2'] = [
    ['embree',  'n', [19.969]],
    ['embree',  'y', [29.512]],
    ['hh',      'n', [19.543]],
    ['hh',      'y', [34.014]],
    ['hh2',     'n', [18.959]],
    ['hh2',     'y', [34.535]],
    ['sf01',    'n', [18.892]],
    ['sf01',    'y', [29.575]],
    ['mt',      'n', [19.620]],
    ['mt',      'y', [34.240]],
    ['bw12',    'n', [19.307]],
    ['bw12',    'y', [34.079]],
    ['bw9',     'n', [13.415]],
    ['bw9',     'y', [34.643]],
    ['shev',    'n', [17.123]],
    ['shev',    'y', [36.595]],
    ['ds',      'n', [18.276]],
    ['ds',      'y', [29.214]],
    ['embree2', 'n', [19.940]],
    ['embree2', 'y', [29.875]]
]

# 5k frames
data['comp2', 'x', 'sanm1'] = [
    ['embree',  'n', [6.591]],
    ['embree',  'y', [10.737]],
    ['hh',      'n', [6.4874]],
    ['hh',      'y', [12.5437]],
    ['hh2',     'n', [6.322]],
    ['hh2',     'y', [12.362]],
    ['sf01',    'n', [6.239]],
    ['sf01',    'y', [10.920]],
    ['mt',      'n', [6.368]],
    ['mt',      'y', [12.532]],
    ['bw12',    'n', [6.405]],
    ['bw12',    'y', [12.402]],
    ['bw9',     'n', [4.853]],
    ['bw9',     'y', [12.349]],
    ['shev',    'n', [5.660]],
    ['shev',    'y', [12.688]],
    ['ds',      'n', [5.954]],
    ['ds',      'y', [10.540]],
    ['embree2', 'n', [6.525]],
    ['embree2', 'y', [10.873]]
]

# 1k frames
data['comp2', 'x', 'sanm2'] = [
    ['embree',  'n', [1.778]],
    ['embree',  'y', [3.432]],
    ['hh',      'n', [1.752]],
    ['hh',      'y', [3.872]],
    ['hh2',     'n', [1.721]],
    ['hh2',     'y', [3.873]],
    ['sf01',    'n', [1.664]],
    ['sf01',    'y', [3.385]],
    ['mt',      'n', [1.748]],
    ['mt',      'y', [3.855]],
    ['bw12',    'n', [1.728]],
    ['bw12',    'y', [3.824]],
    ['bw9',     'n', [1.350]],
    ['bw9',     'y', [3.752]],
    ['shev',    'n', [1.531]],
    ['shev',    'y', [3.907]],
    ['ds',      'n', [1.615]],
    ['ds',      'y', [3.266]],
    ['embree2', 'n', [1.761]],
    ['embree2', 'y', [3.422]]
]

# 500 frames
data['comp2', 'x', 'dragon'] = [
    ['embree',  'n', [0.2078]],
    ['embree',  'y', [0.5255]],
    ['hh',      'n', [0.2047]],
    ['hh',      'y', [0.5553]],
    ['hh2',     'n', [0.2030]],
    ['hh2',     'y', [0.5693]],
    ['sf01',    'n', [0.1923]],
    ['sf01',    'y', [0.4424]],
    ['mt',      'n', [0.2052]],
    ['mt',      'y', [0.5515]],
    ['bw12',    'n', [0.1995]],
    ['bw12',    'y', [0.5357]],
    ['bw9',     'n', [0.1571]],
    ['bw9',     'y', [0.5366]],
    ['shev',    'n', [0.1802]],
    ['shev',    'y', [0.5859]],
    ['ds',      'n', [0.1892]],
    ['ds',      'y', [0.4581]],
    ['embree2', 'n', [0.2082]],
    ['embree2', 'y', [0.5188]]
]


def compute_avg_pct():
    series = defaultdict(list)
    for key, value in data.items():
        for this_algo, this_packets, fps in value:
            # Fix this to support multiple data series
            fps = fps[0]
            if this_algo == 'embree':
                if this_packets == 'n':
                    pct_1ray = fps
                else:
                    pct_Krays = fps
            if this_packets == 'n':
                pct = pct_1ray
            else:
                pct = pct_Krays
            series[this_algo, this_packets].append(fps / pct)
    series = {k : round(np.mean(v)*100) for k, v in series.items()}

    algo_order = ('embree2',
                  'bw9', 'bw12', 'ds', 'hh', 'hh2',
                  'mt', 'sf01', 'shev')

    means_1rays = [series[algo, 'n'] for algo in algo_order]
    means_Krays = [series[algo, 'y'] for algo in algo_order]

    fix, ax = plt.subplots(figsize = (8, 4))
    ax.axhline(y = 100, linestyle = 'dotted', color = 'k', linewidth = 0.8)
    index = np.arange(len(algo_order))
    bar_width = 0.35
    opacity = 0.8
    print(means_1rays, means_Krays)

    rects1 = plt.bar(index, means_1rays, bar_width, opacity,
                     color = 'b',
                     label = 'Single rays')

    rects2 = plt.bar(index + bar_width + 0.04, means_Krays, bar_width, opacity,
                     color = 'g',
                     label = 'Ray packets')

    plt.ylabel('Relative framerate')
    #plt.title("Algorithm performance relative to Embree's default")
    plt.legend(framealpha = 0.4)
    plt.xticks(index + bar_width / 2 + 0.02, algo_order, rotation = 45)

    plt.tight_layout()
    plt.savefig('algos.png')


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
    compute_pcts('comp2', 'x', 'crown')
    compute_avg_pct()
