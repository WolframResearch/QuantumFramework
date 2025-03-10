[The Wolfram Quantum Framework](https://resources.wolframcloud.com/PacletRepository/resources/Wolfram/QuantumFramework/) brings a broad, coherent design for quantum computation, together with a host of leading-edge capabilities, and full integration into Mathematica and Wolfram Language. Starting from discrete quantum mechanics, the Framework provides a high-level symbolic representation of quantum basis, states, and operations. The Framework can perform measurements and is equipped with various well-known states and operations, such as Bell states, multiplexers and quantum channels. Using such simulation capabilities, one can use the Framework to model and simulate the time evolution of isolated or open quantum systems, quantum circuits and algorithms.

Install and use the development version:
```
PacletInstall["https://www.wolfr.am/DevWQCF", ForceVersionInstall -> True]
<< Wolfram`QuantumFramework`
```

Load the local installation:
```
PacletDirectoryLoad["path_to_the_paclet_source/QuantumFramework"]
<< Wolfram`QuantumFramework`
```

To contact [Wolfram quantum team](https://www.wolfram.com/quantum-computation-framework/), please contact us at quantum [at] wolfram [dot] com.
