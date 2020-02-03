# PyVCF-P
Parallel multiprocess version of PyVCF.

## Current Status
Interprocess communication (serialization/deserialization) is bottlenecking the main process and slowing down worker processes.  I am looking into using shared memory reference passing to fix this.
