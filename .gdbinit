add-auto-load-safe-path /Users/markg/github.com/blma/.gdbinit
python
import sys
sys.path.insert(0, '/Users/markg/github.com/blma/eigen')
from printers import register_eigen_printers
register_eigen_printers (None)
end
