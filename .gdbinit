add-auto-load-safe-path /home/markg/Dropbox/blma/.gdbinit
python
import sys
sys.path.insert(0, '/home/markg/eigen')
from printers import register_eigen_printers
register_eigen_printers (None)
end
