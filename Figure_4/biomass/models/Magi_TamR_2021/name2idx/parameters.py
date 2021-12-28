NAMES = [
    'v1_k',
    'v1_i',
    'v2_c',
    'v2_a',
    'v2_b',
    'v3_c',
    'v4_k',
    'v5_k',
    'v6_k',
    'v7_c',
    'v7_a',
    'v7_b',
    'v8_c',
    'v9_k',
    'v10_c',
    'v10_a',
    'v10_b',
    'v11_c',
    'v12_k',
]

for idx, name in enumerate(NAMES):
    exec(
        '{} = {:d}'.format(
            name, idx
        )
    )

NUM = len(NAMES)
