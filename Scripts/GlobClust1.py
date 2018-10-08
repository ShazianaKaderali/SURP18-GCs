from gaia_tools import query
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt

globquery="""select ra, dec, l, b, pmra, pmdec from gaiadr2.gaia_source where abs(l-305.35370)<0.1 and a\
bs(b-(-40.47324))<0.1 limit 100"""

out=query.query(globquery,local=True,timeit=True)

plt.figure()
plt.scatter(out['l'], out['b'])
plt.show()
