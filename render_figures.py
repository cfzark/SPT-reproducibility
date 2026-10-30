"""Re-render the current 3 main + 4 supplement figures from aggregate coordinates."""
from pathlib import Path
import subprocess,sys,os
import tempfile
cache=Path(tempfile.gettempdir())/'spt_paper_plot_cache';cache.mkdir(exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(cache/'matplotlib'))
os.environ.setdefault('XDG_CACHE_HOME',str(cache))
root=Path(__file__).resolve().parent
for f in ['ablation.py','large_n.py','clinical.py','large_n_diagnostics.py']:
 subprocess.run([sys.executable,str(root/'plots'/f)],cwd=root,check=True)
subprocess.run([os.environ.get('RSCRIPT','Rscript'),'plots/clinical_curves.R'],cwd=root,check=True)
