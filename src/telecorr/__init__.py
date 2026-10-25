"""telecorr: two-qubit quantum teleportation through a noisy channel, repaired by a correction learned from counts.

The package follows the nine steps of the study (one notebook each):

1. teleportation ............ ``protocol`` (exact six-qubit simulation, qiskit circuits), ``states``, ``paulis``
2. channel noise ............ ``channels`` (Kraus / Choi / transfer matrices), ``twirl`` (what reaches Bob)
3. correction from counts ... ``probes`` (the measurable cost), ``calibrate`` (COBYLA)
4. automatic recalibration .. ``calibrate`` (``AutoCalibrator``)
5. incoherent noise ......... ``ceiling`` (Procrustes), ``calibrate`` (Rotosolve, ``RelativeAutoCalibrator``)
6. correlated noise ......... ``correlated`` (two-qubit twirl, ceilings, TOMO-75), ``detector``
7. efficient learning ....... ``tomography`` (TOMO-9 / TOMO-12), ``tracking`` (Kalman tracker, ABC rule)
8. real hardware ............ ``hardware`` (the IBM runs of 4 October, from the raw counts)
9. replication .............. ``hardware.replication`` (a new day and a second device, 25 October)

``rotations``, ``reference``, ``compare``, ``plotting`` and ``report`` are shared helpers.
"""
from ._version import __version__

__all__ = ["__version__"]
