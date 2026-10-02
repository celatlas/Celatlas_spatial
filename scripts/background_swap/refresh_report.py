"""Compatibility module name for the transactional report refresh CLI.

The old A1-only --sample-dir/--sample/--yes interface is intentionally removed.
Use prepare --swap-run ..., then apply --run ... --yes. See README.md.
"""
from .report_refresh import main

if __name__ == '__main__':
    main()
