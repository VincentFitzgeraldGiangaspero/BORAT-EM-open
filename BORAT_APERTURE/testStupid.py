import itertools
from tqdm import tqdm

ip=100
it=50

for i1, i2 in tqdm(itertools.product(range(ip), range(it)), total=ip * it, desc="Processing", leave=True):
    print(i1)
    print(i2)