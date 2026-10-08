"""build4d/slim.pkl for a newly built patient, from build/stage2.pkl.

ofj4d/inputs.py recomputes the runtime inputs (mesiobuccal cusps, crown centres, 16 axes,
crest points) from a full model when slim.pkl is absent; stage 2 carries every field it
reads, so it stands in for build/model.pkl here.
"""
import sys, os, pickle
sys.path.insert(0, os.path.join(os.getcwd(), 'ofj4d'))
import inputs

assert not os.path.exists(inputs.SLIM), f'{inputs.SLIM} exists; remove it to rebuild'
IN, _ = inputs.load('build/stage2.pkl')
os.makedirs('build4d', exist_ok=True)
pickle.dump(IN, open(inputs.SLIM, 'wb'))
print('->', inputs.SLIM, {k: (sorted(v) if isinstance(v, dict) else None) for k, v in IN.items()})
