set -ex
rm -fr eden
mkdir eden
cd eden
eden nml ../NeuroML2/LEMS_GranCellLayer.xml

python3 << EOF
import pandas as pd
import matplotlib.pyplot as plt

a = pd.read_csv('./Golgis_0.0.dat', delim_whitespace=True, header=None)
b = pd.read_csv('./Golgis_1.0.dat', delim_whitespace=True, header=None)
c = pd.read_csv('./Golgis_2.0.dat', delim_whitespace=True, header=None)
d = pd.read_csv('./Golgis_3.0.dat', delim_whitespace=True, header=None)
df = pd.merge(pd.merge(a, b, on=0), pd.merge(c, d, on=0), on=0)
df.set_index(0, inplace=True)
df.plot()
plt.show()
EOF
