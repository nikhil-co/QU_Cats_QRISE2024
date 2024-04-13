# Submission for Microsoft Quantum research challenge at QRISE 2024

Team: Qu-Cats

We have implemented Quantum Approximate Optimization Algorithm using **Q#** and **Python**. We also have produced resource estimates for QAOA using the Azure Quantum Resource Estimator.

## File and Folder guide

### Folders

* `src`
  * *Main.qs* - **Q#** code for running the QAOA circuit. It can't be run on its own but requires a **Python** host to run.
* `resource estimation`
  * *qaoa_resource_estimator_graphs.ipynb* - Notebook for resource estimation graphs.
  * *qaoa_resource_estimator.ipynb* - Notebook for data using the **Resource Estimator**.
* `graphs` - has all the graphs we generated in the *PNG* format.

### Files

* *qaoa_host.py* - is a **Python** script for running QAOA on the NPP problem. This script uses the *Main.qs* file to run the QAOA circuit.
* *qaoa_notebook.ipynb* - is a standalone notebook to run QAOA on the NPP problem.

Even though we have implemented QAOA for the Number Partitioning Problem, we have designed the code so that it can be used to solve any other *QUBO* problem by changing the `build_qubo` function.

* *qaoa_steps.pdf* - is a *PDF* we have made to out line the steps involved in implementing and run the QAOA algorithm. This file serves as a reference for our implementation.

### Further Instructions

    To run the Python script, you must set the project root to the **qaoa** folder.

    ``` python
    qsharp.init(project_root = './qaoa')

    ```

    > Note: ***Use*** the forward slash (/), the import will fail with backward slash (\ ).

* Notebook should run without any issues if you have ```qsharp``` installed.
* For resource estimation notebooks you need to have ```qsharp_widgets``` installed.

> We have provided requirements.txt file consisting of all packages you require to run the code we have provided.
