# -*- coding: utf-8 -*-
"""Module providing subcommands for pdp.

Each file within this module exposes a function with the same
name as the file; this is the code that runs the corresponding
subcommand when called via `pdp <SUBCOMMAND>`

For example, the subcmd_dedupe.py file exposes the function
subcmd_dedupe(), which is called (with appropriate arguments)
to execute the `pdp dedupe` script subcommand.
"""

from .subcmd_config import subcmd_config
from .subcmd_filter import subcmd_filter
from .subcmd_eprimer3 import subcmd_eprimer3
from .subcmd_primer3 import subcmd_primer3
from .subcmd_primersearch import subcmd_primersearch
from .subcmd_dedupe import subcmd_dedupe
from .subcmd_blastscreen import subcmd_blastscreen
from .subcmd_classify import subcmd_classify
from .subcmd_extract import subcmd_extract
from .subcmd_plot import subcmd_plot
