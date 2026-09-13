#!/usr/bin/env python3
"""Telomere-independent estimate of the TR gene's core Template sequence:
start 2bp after the G-rich-5'-end -> Template boundary (HMM column 128,
found by pure cross-species conservation, no repeat data), take 12bp.
Validated at mean IoU=0.836 against 1073 telomere-verified cores - see
notes/tr_template_boundary_prediction.md.

Usage: predict_tr_template.py <padded_template_window>
  (the ~47bp string from outputs/tr_template_loci or outputs/tr_domains)
"""
import sys

START_OFFSET = 2
LENGTH = 12


def predict_core(template_window):
    return template_window[START_OFFSET:START_OFFSET + LENGTH]


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    print(predict_core(sys.argv[1]))
