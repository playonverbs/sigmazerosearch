"""
Provides default definitions for parameters and configuration options.

Currently the only additions here are the equivalent frozen 'tunes' created in the
Lambda analysis. These are prepended with `cthorpe_` to label as them being
from that analysis.
"""

from sigmazerosearch.selection import ParameterSet

cthorpe_FHC_Tune_325 = ParameterSet(
    min_length=10,
    max_separation=1,
    separation_cut=3,
    pid_cut=0.6,
    proton_pid_cut=0.1,
    pion_pid_cut=-0.1,
    w_lambda_min=1,
    w_lambda_max=2,
)
"""
Corresponds to [Frozen FHC Tune 325](https://github.com/cthorpe123/HyperonAnalysis/blob/master/ParameterSets/Parameters.h#L29)
"""

cthorpe_RHC_Tune_397 = ParameterSet(
    min_length=10,
    max_separation=1,
    separation_cut=3,
    pid_cut=0.6,
    proton_pid_cut=0.1,
    pion_pid_cut=-0.1,
    w_lambda_min=1,
    w_lambda_max=2,
)
"""
Corresponds to [Frozen RHC Tune 397](https://github.com/cthorpe123/HyperonAnalysis/blob/master/ParameterSets/Parameters.h#L53)
"""
