from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

h1pdm = AlignIO.read("h1pdm.max2ambichar.trimmed.realigned", "fasta")
with open("h1pdm.gappedcols", "w") as f:
	AlignIO.write(h1pdm[:, 33:34]+h1pdm[:, 1262:1263]+h1pdm[:, 1519:1520]+h1pdm[:, 1540:1542]+h1pdm[:, 1681:1687]+h1pdm[:, 1689:1692]+h1pdm[:, 1699:1702]+h1pdm[:, 1703:1704]+h1pdm[:, 1714:1715]+h1pdm[:, 1716:1717], f, "fasta")


