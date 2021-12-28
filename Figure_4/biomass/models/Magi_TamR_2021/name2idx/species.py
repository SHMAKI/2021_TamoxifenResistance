NAMES = [
"S",
"P",
"R1",
"R2",
"Tam",
"Tam_integ",
]

for idx, name in enumerate(NAMES):
    exec(
        '{} = {:d}'.format(
            name, idx
        )
    )

NUM = len(NAMES)
