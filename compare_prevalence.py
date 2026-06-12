import sys; sys.path.insert(0, '/home/welcome3/Viscacha_pipeline')
import numpy as np
from layer0.utils.sparse_utils import donor_detection_vector
from layer0.utils.qc_log import QCLogger
import layer0.step01_metadata as step01
import layer0.step02_barcode_merge as step02
import layer0.step03_barcode_filter as step03
from layer0.config import COND_AD, COND_CTRL

qc = QCLogger()
meta     = step01.run(qc)
adata_tx = step02.run(meta, qc)
adata_f  = step03.run(adata_tx, qc)

print("\n{:<22} {:>7} {:>7} {:>7}  {}".format(
    'Cell type', '30%', '40%', 'extra', 'zero-ctrl / zero-AD among extras'))
print('-' * 75)

for ct in sorted(adata_f.obs['cell_type'].dropna().unique()):
    sub    = adata_f[adata_f.obs['cell_type'] == ct]
    donors = sub.obs['donor'].unique()
    n      = len(donors)

    detection = np.zeros(adata_f.n_vars, dtype=np.int32)
    for donor in donors:
        rows = np.where(sub.obs['donor'] == donor)[0]
        detection += donor_detection_vector(sub.X[rows, :]).astype(np.int32)

    prev      = detection / n
    n30       = int((prev >= 0.30).sum())
    n40       = int((prev >= 0.40).sum())
    extra_idx = np.where((prev >= 0.30) & (prev < 0.40))[0]

    ctrl_donors = [d for d in donors
                   if sub.obs[sub.obs['donor'] == d]['condition'].iloc[0] == COND_CTRL]
    ad_donors   = [d for d in donors
                   if sub.obs[sub.obs['donor'] == d]['condition'].iloc[0] == COND_AD]

    zero_ctrl = zero_ad = 0
    for idx in extra_idx:
        c_det = sum(1 for d in ctrl_donors
                    if sub[sub.obs['donor'] == d].X[:, idx].sum() > 0)
        a_det = sum(1 for d in ad_donors
                    if sub[sub.obs['donor'] == d].X[:, idx].sum() > 0)
        if c_det == 0: zero_ctrl += 1
        if a_det == 0: zero_ad   += 1

    if len(extra_idx):
        pct  = (zero_ctrl + zero_ad) / len(extra_idx) * 100
        note = "{} / {}  ({:.0f}% one-sided)".format(zero_ctrl, zero_ad, pct)
    else:
        note = '-'

    print("{:<22} {:>7,} {:>7,} {:>7,}  {}".format(ct, n30, n40, n30-n40, note))
