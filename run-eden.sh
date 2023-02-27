set -ex
rm -fr eden
mkdir eden
cd eden
eden nml ../NeuroML2/LEMS_SimMaexDeSchutter1998.xml

python3 << EOF
import pandas as pd
import matplotlib.pyplot as plt

a = pd.read_csv('SimMaexDeSchutter1998.Gols.v.dat', delim_whitespace=True, header=None)
b = pd.read_csv('SimMaexDeSchutter1998.GrCs.v.dat', delim_whitespace=True, header=None)

a, b = a * 1000, b * 1000

a = a.set_index(0).rename(columns='Gol{}'.format)
b = b.set_index(0).rename(columns='GrCs{}'.format)

df = pd.merge(a, b, left_index=True, right_index=True)
df.plot()

plt.show()
EOF
